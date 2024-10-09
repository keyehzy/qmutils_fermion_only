#/usr/bin/env bash

set -e
set -o pipefail
set -u

if ! command -v clang-format &> /dev/null
then
    echo "clang-format could not be found"
    exit 1
fi

ROOT="$(git rev-parse --show-toplevel)"

find_repository_files() {
    git ls-files --cached --exclude-standard -- "${ROOT}"
}

FILES="$(find_repository_files | grep -v '^vendor' | grep -E '\.(c|cpp|h)$')"

echo "${FILES}" | xargs clang-format -i

echo "Formatted the following files:"
echo "${FILES}"
