"""Sum token usage from Claude Code session transcript(s).

Reads the JSONL transcripts Claude Code writes per session and totals the usage
across all assistant turns: output, fresh input, cache writes and cache reads.

Stdlib only (no dependencies) so it runs anywhere:  python3 analysis/session_tokens.py

Usage:
    python3 analysis/session_tokens.py                # this project's transcripts
    python3 analysis/session_tokens.py <file.jsonl>   # one session
    python3 analysis/session_tokens.py <dir>          # all *.jsonl in a directory

With no argument it locates this project's transcript directory under
~/.claude/projects/<cwd-with-slashes-as-dashes>/ and aggregates every session.
Note cache-read tokens dominate the raw total but are billed at a steep discount;
`/cost` in Claude Code remains the source of truth for dollar cost.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def project_transcript_dir() -> Path:
    """~/.claude/projects/<abs-cwd-with-/-and-.-as-dashes>/ for the current repo."""
    slug = str(Path.cwd().resolve()).replace("/", "-").replace(".", "-")
    return Path.home() / ".claude" / "projects" / slug


def _cache_creation(usage: dict) -> int:
    if "cache_creation_input_tokens" in usage:
        return usage.get("cache_creation_input_tokens") or 0
    cc = usage.get("cache_creation")  # newer nested form
    return sum(v for v in cc.values() if isinstance(v, int)) if isinstance(cc, dict) else 0


def tally(path: Path) -> dict:
    totals = {"turns": 0, "input": 0, "output": 0, "cache_write": 0, "cache_read": 0, "models": set()}
    for line in path.open():
        line = line.strip()
        if not line:
            continue
        try:
            msg = (json.loads(line).get("message") or {})
        except json.JSONDecodeError:
            continue
        usage = msg.get("usage")
        if not usage:
            continue
        totals["turns"] += 1
        totals["input"] += usage.get("input_tokens", 0) or 0
        totals["output"] += usage.get("output_tokens", 0) or 0
        totals["cache_write"] += _cache_creation(usage)
        totals["cache_read"] += usage.get("cache_read_input_tokens", 0) or 0
        if msg.get("model"):
            totals["models"].add(msg["model"])
    return totals


def _report(name: str, t: dict) -> None:
    grand = t["input"] + t["output"] + t["cache_write"] + t["cache_read"]
    non_cached = t["input"] + t["output"] + t["cache_write"]
    print(f"\n{name}")
    print(f"  model(s):                  {', '.join(sorted(t['models'])) or 'n/a'}")
    print(f"  assistant turns:           {t['turns']:>14,}")
    print(f"  output:                    {t['output']:>14,}")
    print(f"  input (fresh, non-cached): {t['input']:>14,}")
    print(f"  input (cache write):       {t['cache_write']:>14,}")
    print(f"  input (cache read):        {t['cache_read']:>14,}")
    print(f"  non-cached (in+out+write): {non_cached:>14,}")
    print(f"  GRAND TOTAL (all tokens):  {grand:>14,}")


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("target", nargs="?", default=None, help="a .jsonl file, a directory, or omit for this project")
    args = ap.parse_args(argv)

    if args.target:
        target = Path(args.target)
        files = [target] if target.is_file() else sorted(target.glob("*.jsonl"))
    else:
        d = project_transcript_dir()
        files = sorted(d.glob("*.jsonl"))
        if not files:
            sys.exit(f"No transcripts found in {d}")

    grand = {"turns": 0, "input": 0, "output": 0, "cache_write": 0, "cache_read": 0, "models": set()}
    for f in files:
        t = tally(f)
        _report(f.name, t)
        for k in ("turns", "input", "output", "cache_write", "cache_read"):
            grand[k] += t[k]
        grand["models"] |= t["models"]

    if len(files) > 1:
        _report(f"ALL {len(files)} SESSIONS", grand)


if __name__ == "__main__":
    main()
