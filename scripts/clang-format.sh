#!/usr/bin/env bash
set -euo pipefail

CLANG_FORMAT_MAJOR=22

find_clang_format() {
  if command -v "clang-format-${CLANG_FORMAT_MAJOR}" &>/dev/null; then
    echo "clang-format-${CLANG_FORMAT_MAJOR}"
    return
  fi
  if command -v clang-format &>/dev/null; then
    local version
    version=$(clang-format --version | grep -oE '[0-9]+' | head -1)
    if [[ "${version}" == "${CLANG_FORMAT_MAJOR}" ]]; then
      echo "clang-format"
      return
    fi
    echo "error: found clang-format ${version}, need ${CLANG_FORMAT_MAJOR}" >&2
    echo "  Linux: sudo apt-get install clang-format-${CLANG_FORMAT_MAJOR}" >&2
    echo "  macOS: brew install llvm@${CLANG_FORMAT_MAJOR}" >&2
    exit 1
  fi
  echo "error: clang-format-${CLANG_FORMAT_MAJOR} not found." >&2
  echo "  Linux: sudo apt-get install clang-format-${CLANG_FORMAT_MAJOR}" >&2
  echo "  macOS: brew install llvm@${CLANG_FORMAT_MAJOR}" >&2
  exit 1
}

CLANG_FORMAT=$(find_clang_format)
ROOT="$(git rev-parse --show-toplevel)"
FILES=$(find "${ROOT}/src" "${ROOT}/test" -name '*.cpp' -o -name '*.hpp')

if [[ "${1:-}" == "--fix" ]]; then
  echo "${FILES}" | xargs "${CLANG_FORMAT}" -i
  echo "Done."
else
  if ! echo "${FILES}" | xargs "${CLANG_FORMAT}" --dry-run --Werror 2>&1; then
    echo ""
    echo "Run './scripts/clang-format.sh --fix' to apply formatting."
    exit 1
  fi
fi
