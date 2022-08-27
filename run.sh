#!/bin/bash
set -euo pipefail
make bin/main
chmod u+x bin/main
bin/main