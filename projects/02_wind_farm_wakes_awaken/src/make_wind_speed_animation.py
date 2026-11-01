#!/usr/bin/env python3
"""
Create an animated GIF of WFP wind speed evolution across all available hours.

This script reads all hourly ``wrfout_d01_{date}_HH_00_00`` files in the WFP directory,
plots wind speed at a fixed AGL height (default 113 m), overlays turbine sites (red),
and saves a GIF with consistent color scaling across frames.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from wrf_utils import (
    destagger_wind_speed_at_agl,
    destagger_wind_speed_at_level,
    hub_mass_level_index,
    open_wrf,
    parse_windturbines_txt,
    project_root,
)

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    _HAS_CARTOPY = True
except Exception:  # pragma: no cover
    _HAS_CARTOPY = False


def _turbine_map_extent(tlon: np.ndarray, tlat: np.ndarray, pad_lon: float, pad_lat: float) -> tuple[float, float, float, float]:
    return (
        float(tlon.min()) - pad_lon,
        float(tlon.max()) + pad_lon,
        float(tlat.min()) - pad_lat,
        float(tlat.max()) + pad_lat,
    )


def _available_hours(wfp_dir: Path, date: str) -> list[int]:
    pat = re.compile(rf"^wrfout_d01_{re.escape(date)}_(\d\d)_00_00$")
    hrs: list[int] = []
    for p in sorted(wfp_dir.glob(f"wrfout_d01_{date}_*_00_00")):
        m = pat.match(p.name)
        if not m:
            continue
        try:
            hrs.append(int(m.group(1)))
        except ValueError:
            continue
    return sorted(set(hrs))


def _plot_frame(
    lons: np.ndarray,
    lats: np.ndarray,
    wspd: np.ndarray,
    *,
    title: str,
    out_png: Path,
    levels: np.ndarray,
    extent: tuple[float, float, float, float] | None,
    tlon: np.ndarray | None,
    tlat: np.ndarray | None,
    cmap: str,
) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    if extent is None:
        extent = (float(lons.min()), float(lons.max()), float(lats.min()), float(lats.max()))
    lon_min, lon_max, lat_min, lat_max = extent
    show_t = tlon is not None and tlat is not None and len(tlon) > 0

    fig = plt.figure(figsize=(8.6, 6.2))
    if _HAS_CARTOPY:
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.5)
        ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.35, edgecolor="0.45")
        gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="gray", alpha=0.5, linestyle="--")
        gl.top_labels = False
        gl.right_labels = False
        cf = ax.contourf(
            lons,
            lats,
            wspd,
            levels=levels,
            transform=ccrs.PlateCarree(),
            cmap=cmap,
            extend="neither",
        )
        if show_t:
            ax.scatter(
                tlon,
                tlat,
                s=28,
                c="red",
                edgecolors="white",
                linewidths=0.7,
                transform=ccrs.PlateCarree(),
                zorder=8,
            )
    else:
        ax = fig.add_subplot(111)
        cf = ax.contourf(lons, lats, wspd, levels=levels, cmap=cmap, extend="neither")
        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)
        ax.set_xlabel("Longitude (°E)")
        ax.set_ylabel("Latitude (°N)")
        ax.set_aspect("equal")
        if show_t:
            ax.scatter(tlon, tlat, s=28, c="red", edgecolors="white", linewidths=0.7, zorder=8)

    # Straight-ended colorbar (no extend caps); keep it compact for GIF frames.
    cb = fig.colorbar(cf, ax=ax, orientation="vertical", fraction=0.04, pad=0.02)
    cb.set_label("Wind speed (m s⁻¹)")
    ax.set_title(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    root = project_root()
    ap = argparse.ArgumentParser(description="WFP wind speed evolution GIF (all available hours)")
    ap.add_argument("--wfp-dir", type=str, default=str(root / "data" / "wfp"))
    ap.add_argument("--date", type=str, default="2023-08-24")
    ap.add_argument("--hub-agl-m", type=float, default=113.0, help="Target hub height AGL (m); 100–126 typical")
    ap.add_argument(
        "--mass-level",
        type=int,
        default=None,
        help="If set, use this bottom_top index instead of AGL interpolation",
    )
    ap.add_argument("--fig-dir", type=str, default=str(root / "figures"))
    ap.add_argument(
        "--full-domain",
        action="store_true",
        help="Plot full d01 extent instead of zoom on turbine cluster",
    )
    ap.add_argument(
        "--zoom-pad-deg",
        type=float,
        nargs=2,
        default=(0.12, 0.08),
        metavar=("LON", "LAT"),
        help="Extra degrees lon/lat around turbine bbox when zooming",
    )
    ap.add_argument("--windturbines", type=str, default=str(root / "inputs" / "windturbines.txt"))
    ap.add_argument("--cmap", type=str, default="viridis")
    ap.add_argument("--duration-ms", type=int, default=220, help="Frame duration in ms")
    args = ap.parse_args()

    wfp_dir = Path(args.wfp_dir).expanduser()
    if not wfp_dir.is_absolute():
        wfp_dir = root / wfp_dir
    fig_dir = Path(args.fig_dir).expanduser()
    if not fig_dir.is_absolute():
        fig_dir = root / fig_dir

    if not wfp_dir.is_dir():
        print(f"WFP directory not found: {wfp_dir}", file=sys.stderr)
        sys.exit(1)

    hours = _available_hours(wfp_dir, args.date)
    if not hours:
        print(f"No hourly WFP files found for {args.date} under {wfp_dir}", file=sys.stderr)
        print("Expected: wrfout_d01_YYYY-MM-DD_HH_00_00", file=sys.stderr)
        sys.exit(1)

    wt_path = Path(args.windturbines).expanduser()
    if not wt_path.is_absolute():
        wt_path = root / wt_path
    if not wt_path.is_file():
        print(f"Missing windturbines file: {wt_path}", file=sys.stderr)
        sys.exit(1)
    tlat, tlon = parse_windturbines_txt(wt_path)
    pad_lon, pad_lat = float(args.zoom_pad_deg[0]), float(args.zoom_pad_deg[1])
    zoom_extent = _turbine_map_extent(tlon, tlat, pad_lon, pad_lat)
    map_extent: tuple[float, float, float, float] | None = None if args.full_domain else zoom_extent

    # Load grid from first hour
    ds0 = open_wrf(wfp_dir / f"wrfout_d01_{args.date}_{hours[0]:02d}_00_00")
    lats = ds0["XLAT"].isel(Time=0).values
    lons = ds0["XLONG"].isel(Time=0).values

    use_level = args.mass_level is not None
    if use_level:
        k = int(args.mass_level)
        nz = int(ds0.sizes["bottom_top"])
        if not (0 <= k < nz):
            print(f"--mass-level {k} out of range [0, {nz - 1}]", file=sys.stderr)
            sys.exit(1)
    else:
        k = hub_mass_level_index(ds0, target_agl_m=args.hub_agl_m)
    ds0.close()

    wspd_by_hour: dict[int, np.ndarray] = {}
    for h in hours:
        p = wfp_dir / f"wrfout_d01_{args.date}_{h:02d}_00_00"
        ds = open_wrf(p)
        if use_level:
            wspd_by_hour[h] = destagger_wind_speed_at_level(ds, k)
        else:
            wspd_by_hour[h] = destagger_wind_speed_at_agl(ds, args.hub_agl_m)
        ds.close()

    all_stack = np.stack([wspd_by_hour[h] for h in hours], axis=0)
    mn = float(np.nanmin(all_stack))
    mx = float(np.nanmax(all_stack))
    if not np.isfinite(mn) or not np.isfinite(mx) or mx <= mn:
        mn, mx = 0.0, 1.0
    levels = np.linspace(mn, mx, 22)

    height_note = f"level {k}" if use_level else f"{args.hub_agl_m:g} m AGL"
    frame_dir = fig_dir / "_animation_frames"
    frame_dir.mkdir(parents=True, exist_ok=True)
    frame_paths: list[Path] = []
    for h in hours:
        title = f"WFP wind speed • {height_note} • {args.date} {h:02d} UTC"
        out_png = frame_dir / f"frame_{h:02d}utc.png"
        _plot_frame(
            lons,
            lats,
            wspd_by_hour[h],
            title=title,
            out_png=out_png,
            levels=levels,
            extent=map_extent,
            tlon=tlon,
            tlat=tlat,
            cmap=args.cmap,
        )
        frame_paths.append(out_png)

    try:
        from PIL import Image
    except Exception as e:  # pragma: no cover
        print("Missing Pillow. Install it (e.g., pip install pillow) to write GIFs.", file=sys.stderr)
        raise e

    imgs = [Image.open(p).convert("P", palette=Image.Palette.ADAPTIVE) for p in frame_paths]
    first, rest = imgs[0], imgs[1:]
    out_gif = fig_dir / f"wind_speed_with_farm_{hours[0]:02d}to{hours[-1]:02d}utc.gif"
    first.save(
        out_gif,
        save_all=True,
        append_images=rest,
        duration=int(args.duration_ms),
        loop=0,
        optimize=False,
        disposal=2,
    )
    print(f"Saved GIF: {out_gif}")


if __name__ == "__main__":
    main()

