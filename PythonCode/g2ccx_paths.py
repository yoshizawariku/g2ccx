"""Machine-specific lookups kept out of the code: the CalculiX executable and timestamped output folders."""
from __future__ import annotations

import json
import os
from datetime import datetime
from pathlib import Path
import re
import shutil

ENV_SOLVER = "G2CCX_SOLVER"
LOCAL_FILE = ".g2ccx.local.json"  # git-ignored; {"solver": "C:/path/to/ccx_dynamic.exe"}
CANDIDATES = ("ccx_dynamic", "ccx", "ccx_static", "ccx_2.22")
STAMP = re.compile(r"\d{12}")


def _local_file_solver():
    for base in (Path.cwd(), Path(__file__).resolve().parent.parent):
        p = base / LOCAL_FILE
        if p.is_file():
            value = json.loads(p.read_text(encoding="utf-8")).get("solver")
            if value:
                return value, p
    return None, None


def resolve_solver(configured=None):
    """CalculiX executable: config `solver` > $G2CCX_SOLVER > .g2ccx.local.json > PATH."""
    tried = []
    for source, value in (("config 'solver'", configured), (f"${ENV_SOLVER}", os.environ.get(ENV_SOLVER)),
                          (LOCAL_FILE, _local_file_solver()[0])):
        if value:
            p = Path(value).expanduser()
            if p.is_file():
                return str(p.resolve())
            tried.append(f"{source} = {value} (not a file)")
    for name in CANDIDATES:
        found = shutil.which(name)
        if found:
            return found
    raise FileNotFoundError(
        "CalculiX executable not found. Set one of: the 'solver' key in the JSON config, the environment variable "
        f"{ENV_SOLVER}, a {LOCAL_FILE} file ({{\"solver\": \"/path/to/ccx\"}}), or put ccx on PATH."
        + ("".join(f"\n  tried {t}" for t in tried)))


def resolve_output(pattern, create):
    """Expand `{timestamp}` in an output path.

    create=True  -> current time (yyyymmddhhmm), a fresh folder name.
    create=False -> the newest existing folder matching the pattern (so build / solve / view find each other)."""
    pattern = str(pattern)
    if "{timestamp}" not in pattern:
        return pattern
    if create:
        return pattern.replace("{timestamp}", datetime.now().strftime("%Y%m%d%H%M"))
    path = Path(pattern)
    parent = path.parent
    suffix = path.name.replace("{timestamp}", "")
    matches = sorted(p for p in (parent.iterdir() if parent.is_dir() else [])
                     if p.is_dir() and STAMP.fullmatch(p.name[:12]) and p.name[12:] == suffix)
    if not matches:
        raise FileNotFoundError(f"No existing output folder matches {pattern}; run the build step first")
    return str(matches[-1])
