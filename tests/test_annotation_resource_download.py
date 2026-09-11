"""Acquisition integrity and restart controls using a local HTTP server."""

from __future__ import annotations

import hashlib
import json
import sys
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))
import download_annotation_resources as acquisition


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.destination = Path(self.temporary.name) / "resource"
        self.payload = b"public resource data\n" * 100
        self.etag = '"version-one"'
        self.headers = []
        self.ignore_range = False
        self.short_response = False
        self.omit_get_length = False
        owner = self

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_HEAD(self):
                self.send_response(200)
                self.send_header("Content-Length", str(len(owner.payload)))
                if owner.etag:
                    self.send_header("ETag", owner.etag)
                self.end_headers()

            def do_GET(self):
                owner.headers.append(dict(self.headers))
                if self.headers.get("If-Match") not in (None, owner.etag):
                    self.send_error(412)
                    return
                start = int(self.headers.get("Range", "bytes=0-")[6:-1])
                self.send_response(206 if start and not owner.ignore_range else 200)
                if not owner.omit_get_length:
                    self.send_header("Content-Length", str(len(owner.payload) - start))
                if start:
                    self.send_header(
                        "Content-Range",
                        f"bytes {start}-{len(owner.payload) - 1}/{len(owner.payload)}",
                    )
                if owner.etag:
                    self.send_header("ETag", owner.etag)
                self.end_headers()
                self.wfile.write(
                    owner.payload[start:-3]
                    if owner.short_response
                    else owner.payload[start:]
                )
                self.close_connection = True

        self.server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self.thread = threading.Thread(target=self.server.serve_forever, daemon=True)
        self.thread.start()
        self.url = f"http://127.0.0.1:{self.server.server_port}/resource"

    def tearDown(self):
        self.server.shutdown()
        self.server.server_close()
        self.thread.join()
        self.temporary.cleanup()

    def make_partial(self, checksum=None, etag='"version-one"'):
        self.destination.with_name("resource.part").write_bytes(self.payload[:20])
        self.destination.with_name("resource.part.json").write_text(
            json.dumps(
                {
                    "url": self.url,
                    "bytes": len(self.payload),
                    "etag": etag,
                    "checksum": checksum,
                }
            )
        )

    def test_fresh_md5_and_receipt(self):
        checksum = {"type": "md5", "value": hashlib.md5(self.payload).hexdigest()}
        record = acquisition.acquire_file(self.url, self.destination, checksum)
        self.assertEqual(self.destination.read_bytes(), self.payload)
        self.assertEqual(record["sha256"], hashlib.sha256(self.payload).hexdigest())
        self.assertEqual(record["checksum"], ["md5", checksum["value"]])
        self.assertEqual(self.headers[0]["If-Match"], self.etag)
        self.assertEqual(self.headers[0]["User-Agent"], "nf-annotation/0.4")
        self.assertEqual(
            acquisition.acquire_file(self.url, self.destination, checksum), record
        )
        self.assertEqual(len(self.headers), 1)

    def test_changed_existing_file_fails(self):
        acquisition.acquire_file(self.url, self.destination)
        self.destination.write_bytes(b"tampered")
        with self.assertRaisesRegex(acquisition.ResourceDownloadError, "changed"):
            acquisition.acquire_file(self.url, self.destination)

    def test_response_without_content_length_is_verified_by_total_size(self):
        self.omit_get_length = True
        acquisition.acquire_file(self.url, self.destination)
        self.assertEqual(self.destination.read_bytes(), self.payload)

    def test_resume_requires_matching_strong_validator(self):
        self.make_partial()
        acquisition.acquire_file(self.url, self.destination)
        self.assertEqual(self.headers[0]["Range"], "bytes=20-")
        self.assertEqual(self.headers[0]["If-Match"], self.etag)
        self.assertEqual(self.destination.read_bytes(), self.payload)

    def test_changed_source_and_unrecorded_partial_fail(self):
        self.make_partial()
        self.etag = '"version-two"'
        with self.assertRaisesRegex(
            acquisition.ResourceDownloadError, "Source changed"
        ):
            acquisition.acquire_file(self.url, self.destination)
        self.destination.with_name("resource.part.json").unlink()
        with self.assertRaisesRegex(
            acquisition.ResourceDownloadError, "Unrecorded partial"
        ):
            acquisition.acquire_file(self.url, self.destination)
        self.assertFalse(self.destination.exists())
        self.assertEqual(self.headers, [])

    def test_weak_validator_cannot_authorize_resume(self):
        self.etag = 'W/"weak"'
        self.make_partial(etag=None)
        with self.assertRaisesRegex(
            acquisition.ResourceDownloadError, "Cannot verify resumption"
        ):
            acquisition.acquire_file(self.url, self.destination)

    def test_pinned_checksum_authorizes_resume_without_etag(self):
        self.etag = None
        checksum = {"type": "md5", "value": hashlib.md5(self.payload).hexdigest()}
        self.make_partial(["md5", checksum["value"]], etag=None)
        acquisition.acquire_file(self.url, self.destination, checksum)
        self.assertEqual(self.destination.read_bytes(), self.payload)

    def test_checksum_mismatch_never_promotes_input(self):
        with self.assertRaisesRegex(
            acquisition.ResourceDownloadError, "checksum mismatch"
        ):
            acquisition.acquire_file(
                self.url, self.destination, {"type": "md5", "value": "0" * 32}
            )
        self.assertFalse(self.destination.exists())
        self.assertFalse(
            self.destination.with_name("resource.acquisition.json").exists()
        )

    def test_ignored_range_cannot_append_a_full_response(self):
        self.make_partial()
        self.ignore_range = True
        with self.assertRaisesRegex(acquisition.ResourceDownloadError, "HTTP status"):
            acquisition.acquire_file(self.url, self.destination)
        self.assertEqual(
            self.destination.with_name("resource.part").read_bytes(), self.payload[:20]
        )
        self.assertFalse(self.destination.exists())

    def test_short_response_is_retained_as_partial(self):
        self.short_response = True
        with self.assertRaisesRegex(
            acquisition.ResourceDownloadError, "Incomplete download"
        ):
            acquisition.acquire_file(self.url, self.destination)
        self.assertFalse(self.destination.exists())
        self.short_response = False
        acquisition.acquire_file(self.url, self.destination)
        self.assertEqual(self.destination.read_bytes(), self.payload)

    def test_invalid_checksum_config_fails_before_transfer(self):
        for config in (
            {},
            {"type": "md5"},
            {"type": "md5", "value": "bad"},
            {"type": "sha1", "value": "0" * 40},
        ):
            with (
                self.subTest(config=config),
                self.assertRaises(acquisition.ResourceDownloadError),
            ):
                acquisition.acquire_file(self.url, self.destination, config)
        self.assertEqual(self.headers, [])

    def test_bundle_schema_and_lock_fail_before_transfer(self):
        for files in (
            [{"name": "../escape", "url": "https://example.org/x"}],
            [{"name": "a", "url": "http://example.org/x"}],
            [],
        ):
            with (
                self.subTest(files=files),
                self.assertRaises(acquisition.ResourceDownloadError),
            ):
                acquisition.acquire_component(
                    "pfam",
                    "test",
                    {"kind": "annotation_bundle", "files": files},
                    self.destination,
                )
        self.destination.mkdir()
        (self.destination / ".acquisition.lock").write_text("busy\n")
        with self.assertRaisesRegex(
            acquisition.ResourceDownloadError, "lock already exists"
        ):
            acquisition.acquire_component(
                "pfam",
                "test",
                {
                    "kind": "annotation_bundle",
                    "files": [{"name": "a", "url": "https://example.org/x"}],
                },
                self.destination,
            )


if __name__ == "__main__":
    unittest.main()
