#!/usr/bin/env python3
"""
Correct the H1-H4 <-> 594 sample-label swap (discovered 2026-11-03).

All rules live in sample_swap.csv. Each rule rewrites one path *segment*
(a file/folder name, or one '/'- or ','-delimited piece of a config line):

  match = exact   segment == old                      -> new
  match = prefix  segment == old, or starts with old + '_' / '.' / '-'
                                                      -> new + rest of segment
  scope           regex that the *parent* path must match ('' = anywhere)

Rules are checked in file order; the first matching rule wins.
Every swap is done in two passes (old -> unique placeholder -> new), so
A<->B swaps never collide.

Subcommands
-----------
  rename ROOT [--dry-run] [--skip-md5] [--log-dir DIR]
      Rename files/folders under ROOT (e.g. data/ locally, or the FASTQ /
      cellranger folders on the server). Records md5 of every affected
      file before and after, and writes a log + TSV manifest.

  text FILE [FILE ...] [--dry-run]
      Rewrite label segments inside text files (cellranger config CSVs).
"""

import argparse
import csv
import datetime as dt
import hashlib
import os
import re
import sys
import uuid
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_MAPPING = HERE / "sample_swap.csv"
SKIP_NAMES = {".DS_Store"}
PREFIX_SEPS = ("_", ".", "-")


# ----------------------------------------------------------------------------
# Rules
# ----------------------------------------------------------------------------
def load_rules(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    rules = []
    for i, r in enumerate(rows, start=2):
        if r["match"] not in ("exact", "prefix"):
            sys.exit(f"{path}:{i}: match must be 'exact' or 'prefix'")
        rules.append({
            "line": i,
            "level": r["level"],
            "old": r["old"],
            "new": r["new"],
            "match": r["match"],
            "scope": re.compile(r["scope"]) if r["scope"] else None,
        })
    # Every rule must have an inverse (it is a swap, not a one-way rename)
    pairs = {(r["old"], r["new"], r["scope"] and r["scope"].pattern) for r in rules}
    for old, new, scope in pairs:
        if (new, old, scope) not in pairs:
            sys.exit(f"Mapping is not symmetric: {old} -> {new} has no inverse")
    return rules


def new_segment(segment, parent, rules):
    """Return (new_segment, rule) or (None, None) if no rule applies."""
    for r in rules:
        if r["scope"] is not None and not r["scope"].search(parent):
            continue
        old = r["old"]
        if segment == old:
            return r["new"], r
        if r["match"] == "prefix" and segment.startswith(old) \
                and segment[len(old):len(old) + 1] in PREFIX_SEPS:
            return r["new"] + segment[len(old):], r
    return None, None


# ----------------------------------------------------------------------------
# md5
# ----------------------------------------------------------------------------
def md5_file(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def md5_tree(path, skip_md5):
    """{relative file path: md5} for a file or every file under a directory."""
    path = Path(path)
    files = [path] if path.is_file() else sorted(
        p for p in path.rglob("*") if p.is_file() and p.name not in SKIP_NAMES)
    out = {}
    for f in files:
        rel = "." if f == path else str(f.relative_to(path))
        out[rel] = "SKIPPED" if skip_md5 else md5_file(f)
    return out


# ----------------------------------------------------------------------------
# rename
# ----------------------------------------------------------------------------
def plan_renames(root, rules):
    """Walk ROOT (without following symlinks) and collect matching paths."""
    plan = []
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        rel_parent = os.path.relpath(dirpath, root)
        rel_parent = "" if rel_parent == "." else rel_parent
        # Scope regexes are matched against the absolute parent path so that
        # e.g. "(^|/)FASTQs$" works no matter where ROOT is.
        abs_parent = str(Path(dirpath).resolve())
        for name in sorted(dirnames + filenames):
            if name in SKIP_NAMES:
                continue
            new, rule = new_segment(name, abs_parent, rules)
            if new is not None and new != name:
                plan.append({
                    "rel": os.path.join(rel_parent, name),
                    "parent": rel_parent,
                    "old": name,
                    "new": new,
                    "rule": rule,
                    "is_dir": os.path.isdir(os.path.join(dirpath, name)),
                })
    return plan


def final_rel(rel, renamed):
    """Final relative path of REL after all renames (including ancestors)."""
    parts = Path(rel).parts
    return os.path.join(*[renamed.get(os.path.join(*parts[: i + 1]), parts[i])
                          for i in range(len(parts))])


def check_plan(root, plan):
    problems = []
    renamed = {p["rel"]: p["new"] for p in plan}
    sources = set(renamed)
    for p in plan:
        if os.path.islink(os.path.join(root, p["rel"])):
            problems.append(f"refusing to rename symlink: {p['rel']}")
        # Top-level swaps (parent not renamed) must have their partner
        # present, otherwise this would be a one-way rename (e.g. H exists
        # but p594 is missing). Nested items move with their parent instead.
        parent_renamed = any(p["rel"].startswith(s + os.sep) for s in sources)
        tgt = os.path.join(p["parent"], p["new"])
        if not parent_renamed and tgt not in sources:
            problems.append(f"{p['rel']} -> {tgt}: swap partner not found (one-way rename)")
    # Every path in the tree must still be unique after renaming
    seen = {}
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        for name in dirnames + filenames:
            rel = os.path.relpath(os.path.join(dirpath, name), root)
            fin = final_rel(rel, renamed)
            if fin in seen:
                problems.append(f"collision: {seen[fin]} and {rel} would both become {fin}")
            seen[fin] = rel
    return problems


def current_path(root, rel, state):
    """Resolve REL (original relative path) given the names currently on disk."""
    parts = Path(rel).parts
    cur = []
    for i in range(len(parts)):
        orig = os.path.join(*parts[: i + 1])
        cur.append(state.get(orig, parts[i]))
    return os.path.join(root, *cur)


def do_rename(args):
    root = os.path.abspath(args.root)
    if not os.path.isdir(root):
        sys.exit(f"Not a directory: {root}")
    rules = load_rules(args.mapping)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    mode = "DRYRUN" if args.dry_run else "APPLY"
    log_dir = Path(args.log_dir or os.path.join(root, "sample_swap_logs"))
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"rename_{mode}_{stamp}.log"
    manifest_path = log_dir / f"rename_{mode}_{stamp}_md5.tsv"
    log_fh = open(log_path, "w")

    def log(msg=""):
        print(msg)
        log_fh.write(msg + "\n")

    log(f"sample_swap.py rename  mode={mode}  {dt.datetime.now().isoformat(timespec='seconds')}")
    log(f"root:    {root}")
    log(f"mapping: {os.path.abspath(args.mapping)}  (md5 {md5_file(args.mapping)})")
    log()

    plan = plan_renames(root, rules)
    if not plan:
        log("Nothing matched. No changes.")
        log_fh.close()
        return 0

    log(f"{len(plan)} item(s) to rename:")
    for p in plan:
        kind = "dir " if p["is_dir"] else "file"
        log(f"  [{p['rule']['level']:9s}] {kind}  {p['rel']}  ->  {os.path.join(p['parent'], p['new'])}")
    log()

    problems = check_plan(root, plan)
    if problems:
        log("ABORTING, plan problems:")
        for pr in problems:
            log(f"  - {pr}")
        log_fh.close()
        return 1
    log("Plan check passed (no collisions, every top-level swap has its partner).")
    log()

    # md5 before. Only hash top-level plan items; a nested item's files are
    # covered by its ancestor's hash, so skip anything inside another item.
    tops = [p for p in plan
            if not any(p["rel"].startswith(q["rel"] + os.sep) for q in plan)]
    log("Computing md5 (before)..." if not args.skip_md5 else "md5 skipped (--skip-md5)")
    before = {}
    for p in tops:
        before[p["rel"]] = md5_tree(os.path.join(root, p["rel"]), args.skip_md5)

    if args.dry_run:
        write_manifest(manifest_path, root, tops, plan, before, None)
        log(f"\nDRY RUN: nothing renamed.\nLog:      {log_path}\nManifest: {manifest_path}")
        log_fh.close()
        return 0

    # Two-pass rename, deepest paths first so parents stay valid.
    token = uuid.uuid4().hex[:8]
    order = sorted(plan, key=lambda p: p["rel"].count(os.sep), reverse=True)
    state = {}  # original rel path -> current name on disk
    log(f"Pass 1: old -> temporary names (__swaptmp_{token}__)")
    for p in order:
        src = current_path(root, p["rel"], state)
        tmp = f"{p['old']}.__swaptmp_{token}__"
        dst = os.path.join(os.path.dirname(src), tmp)
        os.rename(src, dst)
        state[p["rel"]] = tmp
        log(f"  {p['rel']}  ->  {tmp}")
    log("Pass 2: temporary names -> new names")
    for p in order:
        src = current_path(root, p["rel"], state)
        dst = os.path.join(os.path.dirname(src), p["new"])
        if os.path.lexists(dst):
            sys.exit(f"Unexpected: {dst} exists during pass 2; temp files left in place, see log")
        os.rename(src, dst)
        state[p["rel"]] = p["new"]
        log(f"  {p['old']}.__swaptmp_{token}__  ->  {p['new']}")
    log()

    # md5 after: compare each moved tree to its hash before the move, using
    # the renamed relative paths for files inside renamed subfolders.
    log("Verifying md5 (after)...")
    after, ok = {}, True
    for p in tops:
        new_top = current_path(root, p["rel"], state)
        got = md5_tree(new_top, args.skip_md5)
        expected = {rename_rel(k, p["rel"], plan): v for k, v in before[p["rel"]].items()}
        after[p["rel"]] = got
        if got != expected:
            ok = False
            log(f"  MISMATCH: {p['rel']} -> {os.path.relpath(new_top, root)}")
    write_manifest(manifest_path, root, tops, plan, before, after, state)
    log("All md5 checksums match." if ok else "md5 MISMATCH, see above.")
    log(f"\nLog:      {log_path}\nManifest: {manifest_path}")
    log_fh.close()
    return 0 if ok else 1


def rename_rel(inner_rel, top_rel, plan):
    """Map a path inside TOP (original names) to the names after renaming."""
    if inner_rel == ".":
        return "."
    renamed = {p["rel"]: p["new"] for p in plan}
    parts = Path(inner_rel).parts
    out = []
    for i in range(len(parts)):
        orig = os.path.join(top_rel, *parts[: i + 1])
        out.append(renamed.get(orig, parts[i]))
    return os.path.join(*out)


def write_manifest(path, root, tops, plan, before, after, state=None):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["old_path", "new_path", "md5_before", "md5_after", "match"])
        for p in tops:
            new_top = os.path.join(p["parent"], p["new"])
            for inner, md5b in before[p["rel"]].items():
                old_f = p["rel"] if inner == "." else os.path.join(p["rel"], inner)
                new_inner = rename_rel(inner, p["rel"], plan)
                new_f = new_top if inner == "." else os.path.join(new_top, new_inner)
                md5a = after[p["rel"]].get(new_inner, "MISSING") if after else ""
                w.writerow([old_f, new_f, md5b, md5a,
                            "" if after is None else str(md5a == md5b)])


# ----------------------------------------------------------------------------
# text
# ----------------------------------------------------------------------------
SPLIT = re.compile(r"([/,\n])")


def swap_text(text, rules):
    """Two-pass segment swap on text. Returns (new_text, list of changes)."""
    token = uuid.uuid4().hex[:8]
    if "__SWAP" in text:
        sys.exit("Text already contains a __SWAP placeholder; aborting.")
    pieces = SPLIT.split(text)
    placeholders, changes = {}, []
    # Pass 1: matched segments -> placeholders. The parent path is everything
    # since the last ',' or newline (one config field), matching how paths
    # are scoped in the rename subcommand.
    field_start = 0
    for i, piece in enumerate(pieces):
        if piece in (",", "\n"):
            field_start = i + 1
            continue
        if piece == "/" or piece == "":
            continue
        parent = "".join(pieces[field_start:i]).rstrip("/")
        new, rule = new_segment(piece, parent, rules)
        if new is not None and new != piece:
            ph = f"__SWAP{len(placeholders)}_{token}__"
            placeholders[ph] = new
            changes.append((piece, new, rule["level"]))
            pieces[i] = ph
    # Pass 2: placeholders -> new values
    out = "".join(pieces)
    for ph, new in placeholders.items():
        out = out.replace(ph, new)
    assert "__SWAP" not in out
    return out, changes


def do_text(args):
    rules = load_rules(args.mapping)
    for f in args.files:
        text = Path(f).read_text()
        new, changes = swap_text(text, rules)
        print(f"{f}: {len(changes)} change(s)")
        for old, nw, level in changes:
            print(f"    [{level}] {old} -> {nw}")
        if not args.dry_run and new != text:
            Path(f).write_text(new)
    if args.dry_run:
        print("DRY RUN: no files written.")
    return 0


# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mapping", default=str(DEFAULT_MAPPING))
    sub = ap.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("rename", help="rename files/folders under ROOT")
    r.add_argument("root")
    r.add_argument("--dry-run", action="store_true")
    r.add_argument("--skip-md5", action="store_true",
                   help="don't hash files (e.g. very large FASTQs); not recommended")
    r.add_argument("--log-dir", help="default: ROOT/sample_swap_logs")

    t = sub.add_parser("text", help="swap labels inside text files")
    t.add_argument("files", nargs="+")
    t.add_argument("--dry-run", action="store_true")

    args = ap.parse_args()
    sys.exit(do_rename(args) if args.cmd == "rename" else do_text(args))


if __name__ == "__main__":
    main()
