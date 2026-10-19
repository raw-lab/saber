#!/usr/bin/env bash
# limited.sh KB cmd...  -- run cmd with an address-space limit (ulimit -v, KB)
ulimit -v "$1"; shift; exec "$@"
