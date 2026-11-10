# Code quality

The repository-owned Python code, tests and container helpers use Python 3.12,
Ruff and ty. The pinned `external/nf-helper` submodule remains an upstream
dependency; run its integration checks without reformatting its checkout.

Install the locked runtime and quality environments, then check the whole project:

```bash
pixi install --locked
pixi install --locked -e quality
pixi run --locked -e quality quality
pixi run --locked test
```

`quality` runs Ruff linting, Ruff's formatting check and ty. Type checking resolves
the runtime's installed scientific libraries. Native-only package imports have
local, documented exclusions because those packages live in their pinned tool
images; native controls verify those interfaces. Negative tests mark deliberately
invalid arguments explicitly rather than loosening production types.

Update the relevant behaviour documentation and tests before implementation.
Keep missing, failed, filtered and valid-zero outcomes distinct. Use precise
types for stable state; heterogeneous JSON records retain explicit dictionary
types at their boundaries. Formatting-only commits should preserve the parsed
Python syntax tree and remain separate from functional changes.

For changed Bash scripts, use `shfmt -w path/to/script.sh`, then
`shfmt -d path/to/script.sh` and `bash -n path/to/script.sh`. Run the relevant
behaviour tests on Linux as well as macOS when portability matters. The GC
backfill helper preserves the caller's `PATH` when locating and running `seqtk`.

Nextflow command, workflow integration and Viper checks use the locked Nextflow
version. A release also requires a small native pilot with recorded source,
runtime and resource identities. Report skipped and unavailable checks separately
from passed checks; production reruns require their own authorisation.
