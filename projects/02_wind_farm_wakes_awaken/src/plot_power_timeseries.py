#!/usr/bin/env python3
"""
Question 5: Time series of domain-total POWER from WFP wrfout (2023-08-24, 00–23 UTC).

Sums the POWER field over all horizontal grid cells each hour (WRF: cell-integrated watts).
"""
from __future__ import annotations

import argparse
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from wrf_utils import list_wrfout_d01, read_power_2d, project_root


def parse_time_from_wrfout_name(path: Path) -> datetime | None:
    """
    wrfout_d01_2023-08-24_06_00_00 -> UTC datetime
    """
    m = re.search(r"wrfout_d\d+_(\d{4}-\d{2}-\d{2})_(\d{2})_00_00$", path.name)
    if not m:
        return None
    d, h = m.group(1), int(m.group(2))
    y, mo, da = map(int, d.split("-"))
    return datetime(y, mo, da, h, 0, 0, tzinfo=timezone.utc)


def main() -> None:
    root = project_root()
    ap = argparse.ArgumentParser(description="Total POWER time series")
    ap.add_argument("--wfp-dir", type=str, default=str(root / "data" / "wfp"))
    ap.add_argument("--date", type=str, default="2023-08-24")
    ap.add_argument("--output", type=str, default=str(root / "figures" / "farm_power_timeseries.png"))
    args = ap.parse_args()

    wfp_dir = Path(args.wfp_dir).expanduser()
    if not wfp_dir.is_absolute():
        wfp_dir = root / wfp_dir
    out = Path(args.output).expanduser()
    if not out.is_absolute():
        out = root / out
    out.parent.mkdir(parents=True, exist_ok=True)

    pat = args.date.replace("-", "-")  # keep YYYY-MM-DD
    files = [f for f in list_wrfout_d01(wfp_dir, date_glob=pat + "_*") if parse_time_from_wrfout_name(f)]
    files.sort(key=lambda p: parse_time_from_wrfout_name(p) or datetime.min.replace(tzinfo=timezone.utc))

    if not files:
        print(f"No hourly wrfout files found for date {args.date} in {wfp_dir}", file=sys.stderr)
        sys.exit(1)

    times = []
    total_w = []
    for f in files:
        t = parse_time_from_wrfout_name(f)
        if t is None:
            continue
        if t.date().isoformat() != args.date:
            continue
        p = read_power_2d(f)
        times.append(t)
        total_w.append(float(np.sum(p)))

    if len(times) != 24 and len(times) < 2:
        print(f"Warning: expected 24 hours, got {len(times)}", file=sys.stderr)

    times_np = mdates.date2num(times)
    mw = np.array(total_w) * 1e-6

    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.plot(times, mw, color="darkgreen", linewidth=1.8, marker="o", markersize=4)
    ax.set_xlabel("Time (UTC)")
    ax.set_ylabel("Total power production (MW)")
    ax.set_title(f"WFP simulation — domain-integrated POWER — {args.date}")
    ax.grid(True, alpha=0.35)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=2))
    fig.autofmt_xdate()
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
