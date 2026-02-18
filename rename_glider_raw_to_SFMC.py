# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 10:05:21 2026

@author: OllieTwigge
"""

#!/usr/bin/env python3
"""
Rename Slocum .dbd/.ebd files to the value in the header tag `full_filename:`.

Example:
  01320000.dbd  ->  selkie-2025-250-0-0.dbd

Usage:
  python rename_glider_raw_to_SFMC.py "G:/path/to/folder" --dry-run
  python rename_glider_raw_to_SFMC.py "G:/path/to/folder"
"""

#from __future__ import annotations

import argparse
from pathlib import Path
import re


FULL_FILENAME_RE = re.compile(r"^\s*full_filename\s*:\s*(.+?)\s*$", re.IGNORECASE)
FILENAME_EXT_RE = re.compile(r"^\s*filename_extension\s*:\s*(.+?)\s*$", re.IGNORECASE)


def read_header_lines(path: Path, max_bytes: int = 64 * 1024) -> list[str]:
    """
    Read the first `max_bytes` of the file and decode as text (best-effort),
    then split into lines. This works for these files because the header is ASCII.
    """
    with path.open("rb") as f:
        chunk = f.read(max_bytes)
    text = chunk.decode("utf-8", errors="ignore")
    return text.splitlines()


def extract_full_filename(lines: list[str]) -> str | None:
    for line in lines:
        m = FULL_FILENAME_RE.match(line)
        if m:
            return m.group(1).strip()
    return None


def extract_filename_extension(lines: list[str]) -> str | None:
    for line in lines:
        m = FILENAME_EXT_RE.match(line)
        if m:
            return m.group(1).strip().lstrip(".")
    return None


def pick_nonconflicting_path(target: Path) -> Path:
    """
    If target exists, append _1, _2, ... before extension.
    """
    if not target.exists():
        return target

    stem = target.stem
    suffix = target.suffix
    parent = target.parent
    i = 1
    while True:
        candidate = parent / f"{stem}_{i}{suffix}"
        if not candidate.exists():
            return candidate
        i += 1


def rename_one(path: Path, dry_run: bool = False, prefer_header_ext: bool = False) -> tuple[bool, str]:
    """
    Returns (renamed?, message)
    """
    try:
        lines = read_header_lines(path)
    except Exception as e:
        return False, f"SKIP (read error): {path.name} -> {e}"

    full_name = extract_full_filename(lines)
    if not full_name:
        return False, f"SKIP (no full_filename tag): {path.name}"

    # Extension choice:
    # - default: keep the original file extension (safer)
    # - optional: use filename_extension from header
    if prefer_header_ext:
        header_ext = extract_filename_extension(lines)
        ext = header_ext if header_ext else path.suffix.lstrip(".")
    else:
        ext = path.suffix.lstrip(".")

    if not ext:
        return False, f"SKIP (no extension): {path.name}"

    new_name = f"{full_name}.{ext}"
    if path.name == new_name:
        return False, f"OK (already named): {path.name}"

    target = pick_nonconflicting_path(path.with_name(new_name))

    if dry_run:
        return True, f"DRY: {path.name} -> {target.name}"

    try:
        path.rename(target)
        return True, f"RENAMED: {path.name} -> {target.name}"
    except Exception as e:
        return False, f"FAIL (rename error): {path.name} -> {e}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("folder", help="Folder containing .dbd/.ebd files")
    ap.add_argument("--dry-run", action="store_true", help="Print changes without renaming")
    ap.add_argument("--recursive", action="store_true", help="Also search subfolders")
    ap.add_argument(
        "--prefer-header-ext",
        action="store_true",
        help="Use filename_extension from header if present (otherwise keep original extension)",
    )
    args = ap.parse_args()

    folder = Path(args.folder).expanduser().resolve()
    if not folder.exists() or not folder.is_dir():
        raise SystemExit(f"Not a folder: {folder}")

    patterns = ["*.dbd", "*.ebd", "*.DBD", "*.EBD"]
    files: list[Path] = []
    for pat in patterns:
        files.extend(folder.rglob(pat) if args.recursive else folder.glob(pat))

    files = sorted(set(files))

    renamed = 0
    skipped = 0
    for p in files:
        did, msg = rename_one(p, dry_run=args.dry_run, prefer_header_ext=args.prefer_header_ext)
        print(msg)
        if did:
            renamed += 1
        else:
            skipped += 1

    print(f"\nDone. Candidate changes: {renamed}, skipped/errors: {skipped}")
    if args.dry_run:
        print("Note: --dry-run was set; no files were actually renamed.")
#%%
if __name__ == "__main__":
    main()

