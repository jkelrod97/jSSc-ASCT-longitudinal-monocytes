#!/usr/bin/env python3
"""
Move files/folders that are not needed by the pipeline (README order) into an
archive folder, keeping their relative paths. Used before the H1-H4 <-> 594
label correction (2026-11-03): archived items are NOT relabelled.

    python3 archive_unused.py ROOT LIST [--dry-run] [--archive-dir DIR] [--log-dir DIR]

LIST is a text file with one path per line, relative to ROOT. Blank lines and
lines starting with '#' are ignored.

Safety:
  * moves only (os.rename); refuses if the archive is on a different
    filesystem, so nothing is ever copied or deleted
  * refuses if any destination already exists
  * refuses if a symlink outside the moved items points into them, or a
    relative symlink inside them points outside (those would break)
  * counts files and bytes before and after each move
  * running it again with --restore moves everything back
"""

import argparse
import datetime as dt
import os
import sys

ARCHIVE_NAME = "_archive_pre_label_fix"
README = """\
Archived {date} by archive_unused.py.

These files are not used by the pipeline described in the repository README
(0_A ... 6_*). They were moved here before the sample-label correction and
were NOT relabelled: anything named H1-H4 / H here may actually be patient 594
(baseline/six/twelve/twen4), and vice versa. See sample_swap.csv in the
repository for the mapping. Do not use these files without applying it.
"""


def tree_stats(path):
    """(number of regular files, total bytes); symlinks are not followed."""
    if os.path.islink(path) or os.path.isfile(path):
        return 1, os.lstat(path).st_size
    n = size = 0
    for dirpath, _, files in os.walk(path, followlinks=False):
        for f in files:
            n += 1
            size += os.lstat(os.path.join(dirpath, f)).st_size
    return n, size


def inside(path, parents):
    return any(path == p or path.startswith(p + os.sep) for p in parents)


def breaking_symlinks(root, moved):
    """Symlinks that would stop resolving if MOVED (absolute paths) move."""
    bad = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        for name in dirnames + filenames:
            link = os.path.join(dirpath, name)
            if not os.path.islink(link) or not os.path.exists(link):
                continue
            target = os.readlink(link)
            resolved = os.path.normpath(target if os.path.isabs(target)
                                        else os.path.join(dirpath, target))
            link_moves = inside(link, moved)
            target_moves = inside(resolved, moved)
            if os.path.isabs(target):
                broken = target_moves            # absolute: breaks if target moves
            else:
                # relative: fine only if link and target move together
                same = any(inside(link, [m]) and inside(resolved, [m]) for m in moved)
                broken = (link_moves or target_moves) and not same
            if broken:
                bad.append(f"{os.path.relpath(link, root)} -> {target}")
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("list")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--restore", action="store_true",
                    help="move the listed items back out of the archive")
    ap.add_argument("--archive-dir", help=f"default: ROOT/{ARCHIVE_NAME}")
    ap.add_argument("--log-dir", help="default: ARCHIVE_DIR")
    args = ap.parse_args()

    root = os.path.abspath(args.root)
    archive = os.path.abspath(args.archive_dir or os.path.join(root, ARCHIVE_NAME))
    rels = []
    with open(args.list) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                rels.append(os.path.normpath(line.rstrip("/")))

    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    mode = ("RESTORE" if args.restore else "ARCHIVE") + ("_DRYRUN" if args.dry_run else "")
    log_dir = args.log_dir or archive
    os.makedirs(log_dir, exist_ok=True)
    log_path = os.path.join(log_dir, f"archive_{mode}_{stamp}.log")
    log_fh = open(log_path, "w")

    def log(msg=""):
        print(msg)
        log_fh.write(msg + "\n")

    log(f"archive_unused.py  mode={mode}  {dt.datetime.now().isoformat(timespec='seconds')}")
    log(f"root:    {root}\narchive: {archive}\nlist:    {os.path.abspath(args.list)}\n")

    pairs, problems = [], []
    for rel in rels:
        if rel.startswith("..") or os.path.isabs(rel) or rel == ".":
            problems.append(f"not a path inside ROOT: {rel}")
            continue
        src, dst = os.path.join(root, rel), os.path.join(archive, rel)
        if args.restore:
            src, dst = dst, src
        if inside(src, [archive]) and not args.restore:
            problems.append(f"already inside the archive: {rel}")
        elif not os.path.lexists(src):
            problems.append(f"does not exist: {src}")
        elif os.path.lexists(dst):
            problems.append(f"destination already exists: {dst}")
        else:
            pairs.append((rel, src, dst))
    # No item listed twice or nested inside another listed item
    srcs = [s for _, s, _ in pairs]
    for s in srcs:
        if any(s != o and s.startswith(o + os.sep) for o in srcs):
            problems.append(f"listed inside another listed item: {os.path.relpath(s, root)}")
    # Same filesystem only (rename, never copy)
    probe = archive if os.path.exists(archive) else os.path.dirname(archive)
    dev = os.stat(probe).st_dev
    for _, s, _ in pairs:
        if os.lstat(s).st_dev != dev:
            problems.append(f"different filesystem from archive (would need copying): {s}")
    if not args.restore:
        for b in breaking_symlinks(root, srcs):
            problems.append(f"symlink would break: {b}")

    stats = {}
    for rel, s, d in pairs:
        stats[rel] = tree_stats(s)
        n, b = stats[rel]
        log(f"  {rel}   ({n} files, {b / 1e9:.2f} GB)")
    log()
    if problems:
        log("ABORTING, problems:")
        for p in problems:
            log(f"  - {p}")
        log_fh.close()
        return 1
    log(f"Checks passed: {len(pairs)} item(s) to {'restore' if args.restore else 'archive'}.")
    if args.dry_run:
        log(f"\nDRY RUN: nothing moved.\nLog: {log_path}")
        log_fh.close()
        return 0

    if not args.restore:
        os.makedirs(archive, exist_ok=True)
        readme = os.path.join(archive, "README.txt")
        if not os.path.exists(readme):
            with open(readme, "w") as fh:
                fh.write(README.format(date=dt.date.today().isoformat()))
    ok = True
    for rel, s, d in pairs:
        os.makedirs(os.path.dirname(d), exist_ok=True)
        os.rename(s, d)
        after = tree_stats(d)
        good = after == stats[rel] and not os.path.lexists(s)
        ok &= good
        log(f"  {'moved' if good else 'MISMATCH'}: {os.path.relpath(s, root)} -> {os.path.relpath(d, root)}")
    log("\nAll moves verified (same file count and bytes)." if ok else "\nMISMATCH, see above.")
    log(f"Log: {log_path}")
    log_fh.close()
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
