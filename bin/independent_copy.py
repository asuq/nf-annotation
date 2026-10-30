"""Independent copy-on-write staging shared by publication and native inputs."""

import ctypes
import errno
import functools
import os
import shutil
import sys
from pathlib import Path


@functools.cache
def darwin_clonefile():
    """Resolve the optional platform primitive once per allocation."""
    clonefile = ctypes.CDLL(None, use_errno=True).clonefile
    clonefile.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int]
    clonefile.restype = ctypes.c_int
    return clonefile


def copy_independent(source: Path, destination: Path) -> str:
    """Use filesystem copy-on-write when supported; never share writable inodes."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    unsupported = {
        errno.EXDEV,
        errno.EOPNOTSUPP,
        errno.ENOTTY,
        errno.EINVAL,
        errno.ENOSYS,
    }
    try:
        if sys.platform == "linux":
            import fcntl

            with source.open("rb") as src, destination.open("xb") as dst:
                fcntl.ioctl(dst.fileno(), 0x40049409, src.fileno())  # FICLONE
        elif sys.platform == "darwin":
            if darwin_clonefile()(os.fsencode(source), os.fsencode(destination), 0):
                code = ctypes.get_errno()
                raise OSError(code, os.strerror(code))
        else:
            raise OSError(errno.ENOSYS, "Copy-on-write is unavailable")
    except OSError as error:
        if error.errno not in unsupported:
            raise
        destination.unlink(missing_ok=True)
        shutil.copy2(source, destination)
        return "copied"
    shutil.copystat(source, destination)
    return "cloned"
