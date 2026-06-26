#!/usr/bin/env bash
set -euo pipefail

cmake -S . -B build
ln -sf build/compile_commands.json .
cmake --build build

