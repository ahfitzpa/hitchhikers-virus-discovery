#!/usr/bin/env bash
# ---------------------------------------------------------------------
# stage_reads.sh  MODE  OUTPUT.fastq.gz  INPUT [INPUT ...]
#
#   MODE = symlink | copy
#
# One input, already gzipped  -> symlink (or copy)
# One input, plain fastq      -> gzip into place
# Several inputs              -> concatenate in the given order
#                                (gzip members concatenate legally, so a
#                                 cat of .gz files is a valid .gz file)
#
# Used for split lanes (L001/L002...) and for ONT barcode folders that
# contain many fastq_runid_*.fastq.gz chunks.
# ---------------------------------------------------------------------
set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "usage: stage_reads.sh MODE OUTPUT INPUT [INPUT ...]" >&2
    exit 2
fi

mode="$1"
out="$2"
shift 2

mkdir -p "$(dirname "$out")"
rm -f "$out"

# Prefer pigz when available; harmless if not.
if command -v pigz >/dev/null 2>&1; then
    GZ="pigz -c"
else
    GZ="gzip -c"
fi

if [ "$#" -eq 1 ]; then
    f="$1"
    case "$f" in
        *.gz)
            if [ "$mode" = "copy" ]; then
                cp -f -- "$f" "$out"
            else
                ln -sfn -- "$(readlink -f -- "$f")" "$out"
            fi
            ;;
        *)
            $GZ -- "$f" > "$out"
            ;;
    esac
else
    : > "$out"
    for f in "$@"; do
        case "$f" in
            *.gz) cat -- "$f" >> "$out" ;;
            *)    $GZ -- "$f" >> "$out" ;;
        esac
    done
fi

# Cheap integrity check: a truncated or non-gzip file fails here rather
# than 3 rules downstream inside an aligner.
if command -v gzip >/dev/null 2>&1; then
    gzip -t "$out" 2>/dev/null || {
        echo "[stage_reads] ERROR: $out is not a valid gzip file" >&2
        exit 1
    }
fi
