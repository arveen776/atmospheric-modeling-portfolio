#!/usr/bin/env python3
"""
Question 3: Terrain elevation (AWAKEN region) + wind turbine locations.

Uses HGT from WFP wrfout. Map extent defaults to a zoom around turbine locations (bounding
box + padding). Turbine positions from optional windturbines.txt, else grid cells where
POWER was ever > 0 across all wrfout files in the WFP directory.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from wrf_utils import (
    list_wrfout_d01,
    open_wrf,
    parse_windturbines_txt,
    project_root,
    turbine_mask_from_power,
)

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    _HAS_CARTOPY = True
except Exception:  # pragma: no cover
    _HAS_CARTOPY = False


def main() -> None:
    root = project_root()
    ap = argparse.ArgumentParser(description="Terrain + turbines map")
    ap.add_argument(
        "--wfp-dir",
        type=str,
        default=str(root / "data" / "wfp"),
        help="Directory with WFP wrfout_d01_* files",
    )
    ap.add_argument(
        "--windturbines",
        type=str,
        default="",
        help="Path to windturbines.txt (lat lon type, or index lat lon). "
        "Default: inputs/windturbines.txt if it exists, else infer turbine cells from POWER.",
    )
    ap.add_argument(
        "--output",
        type=str,
        default=str(root / "figures" / "terrain_and_turbines.png"),
        help="Output PNG path",
    )
    ap.add_argument(
        "--pad-lon",
        type=float,
        default=0.07,
        help="Degrees longitude padding around turbine bounding box (zoom)",
    )
    ap.add_argument(
        "--pad-lat",
        type=float,
        default=0.05,
        help="Degrees latitude padding around turbine bounding box (zoom)",
    )
    ap.add_argument(
        "--turbine-size",
        type=float,
        default=32.0,
        help="Scatter marker area scale (matplotlib ``s``)",
    )
    ap.add_argument(
        "--turbine-color",
        type=str,
        default="#0066cc",
        help="Face color for turbine markers (hex or matplotlib name)",
    )
    args = ap.parse_args()
    wfp_dir = Path(args.wfp_dir).expanduser()
    if not wfp_dir.is_absolute():
        wfp_dir = root / wfp_dir
    out = Path(args.output).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)

    files = list_wrfout_d01(wfp_dir)
    if not files:
        print(f"No wrfout_d01_* files in {wfp_dir}", file=sys.stderr)
        sys.exit(1)

    ds = open_wrf(files[0])
    hgt = ds["HGT"].isel(Time=0).values
    lats = ds["XLAT"].isel(Time=0).values
    lons = ds["XLONG"].isel(Time=0).values
    ds.close()

    wt_arg = args.windturbines.strip()
    default_wt = root / "inputs" / "windturbines.txt"
    if wt_arg:
        wt_path = Path(wt_arg).expanduser()
        if not wt_path.is_absolute():
            wt_path = root / wt_path
    elif default_wt.is_file():
        wt_path = default_wt
    else:
        wt_path = None

    if wt_path is not None:
        if not wt_path.is_file():
            print(f"windturbines file not found: {wt_path}", file=sys.stderr)
            sys.exit(1)
        tlat, tlon = parse_windturbines_txt(wt_path)
        if len(tlat) == 0:
            print("No turbine rows parsed from windturbines.txt", file=sys.stderr)
            sys.exit(1)
    else:
        # One day is enough to flag turbine cells; avoids scanning the full multi-day archive.
        mask = turbine_mask_from_power(wfp_dir, date_glob="2023-08-24_*")
        tlat = lats[mask]
        tlon = lons[mask]

    # Zoom to turbine footprint (+ padding) so terrain and farm fit the frame
    lon_min = float(tlon.min()) - args.pad_lon
    lon_max = float(tlon.max()) + args.pad_lon
    lat_min = float(tlat.min()) - args.pad_lat
    lat_max = float(tlat.max()) + args.pad_lat

    fig = plt.figure(figsize=(10, 8))
    if _HAS_CARTOPY:
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent(
            [lon_min, lon_max, lat_min, lat_max],
            crs=ccrs.PlateCarree(),
        )
        ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.6)
        ax.add_feature(cfeature.BORDERS.with_scale("50m"), linewidth=0.4, linestyle=":")
        ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.35, edgecolor="0.45")
        gl = ax.gridlines(draw_labels=True, linewidth=0.5, color="gray", alpha=0.55, linestyle="--")
        gl.top_labels = False
        gl.right_labels = False
        cf = ax.contourf(
            lons,
            lats,
            hgt,
            levels=25,
            transform=ccrs.PlateCarree(),
            cmap="terrain",
            extend="both",
        )
        ax.scatter(
            tlon,
            tlat,
            s=args.turbine_size,
            c=args.turbine_color,
            edgecolors="white",
            linewidths=0.65,
            transform=ccrs.PlateCarree(),
            label="Wind turbines",
            zorder=6,
        )
    else:
        ax = fig.add_subplot(111)
        cf = ax.contourf(lons, lats, hgt, levels=25, cmap="terrain", extend="both")
        ax.scatter(
            tlon,
            tlat,
            s=args.turbine_size,
            c=args.turbine_color,
            edgecolors="white",
            linewidths=0.65,
            label="Wind turbines",
        )
        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)
        ax.set_xlabel("Longitude (°E)")
        ax.set_ylabel("Latitude (°N)")
        ax.set_aspect("equal")

    cb = fig.colorbar(cf, ax=ax, orientation="vertical", fraction=0.035, pad=0.02)
    cb.set_label("Terrain elevation (m)")

    ax.set_title("AWAKEN region: terrain elevation and wind turbine locations")
    ax.legend(loc="lower left", framealpha=0.92)

    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
