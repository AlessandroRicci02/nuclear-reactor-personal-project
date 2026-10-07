#!/usr/bin/env bash
set -u

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="/home/aleri/miniforge3/envs/openmc-env/bin/python"
cd "$HERE" || exit 1

while [[ ! -f solution.json ]]; do
    if "$PYTHON" -u scripts/auto_search.py >> campaign/search_controller.log 2>&1; then
        [[ -f solution.json ]] && exit 0
        printf '%s Search returned without a qualifying solution; restarting.\n' "$(date -Is)" >> campaign/launcher.log
    else
        printf '%s Search process failed; retrying after 60 seconds.\n' "$(date -Is)" >> campaign/launcher.log
    fi
    sleep 60
done
