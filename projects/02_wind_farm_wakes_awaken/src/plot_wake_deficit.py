#!/usr/bin/env python3
"""
Question 4: Destagger U/V, wind speed near hub height; contour WFP & NWF at 0600 and 1800 UTC
as side-by-side pairs (one shared colorbar per figure); difference maps (WFP − NWF) likewise paired.

By default, wind speed is interpolated to a **fixed height AGL** at each horizontal grid point
(using each run's own PH/PHB + HGT). That matches hub-height sampling over terrain and avoids
terrain-parallel artifacts from differencing a single eta level (which is not constant AGL).

Maps default to a **zoom around the turbine cluster** (from windturbines.txt) so the wake is visible.
Use --full-domain to plot the entire d01 domain.
"""
from __future__ import annotations

import argparse
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
    wrfout_path_for_utc,
)

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature

    _HAS_CARTOPY = True
except Exception:  # pragma: no cover
    _HAS_CARTOPY = False


def _style_colorbar_flat_ends(cb) -> None:
    """
    Thin, rounded colorbar frame.
    Slightly rounded corners on the colorbar frame where Matplotlib allows it.
    """
    try:
        outline = getattr(cb, "outline", None)
        if outline is not None:
            outline.set_linewidth(1.0)
            outline.set_edgecolor("0.35")
            for name in ("set_capstyle", "set_joinstyle"):
                fn = getattr(outline, name, None)
                if callable(fn):
                    fn("round")
    except Exception:
        pass


def _turbine_map_extent(tlon: np.ndarray, tlat: np.ndarray, pad_lon: float, pad_lat: float) -> tuple[float, float, float, float]:
    return (
        float(tlon.min()) - pad_lon,
        float(tlon.max()) + pad_lon,
        float(tlat.min()) - pad_lat,
        float(tlat.max()) + pad_lat,
    )


def _scatter_turbines_cartopy(ax, tlon: np.ndarray, tlat: np.ndarray, *, legend: bool) -> None:
    kw: dict = dict(
        s=22,
        c="white",
        edgecolors="black",
        linewidths=0.9,
        transform=ccrs.PlateCarree(),
        zorder=8,
    )
    if legend:
        kw["label"] = "Turbine sites"
    ax.scatter(tlon, tlat, **kw)


def _scatter_turbines_plain(ax, tlon: np.ndarray, tlat: np.ndarray, *, legend: bool) -> None:
    kw: dict = dict(s=22, c="white", edgecolors="black", linewidths=0.9, zorder=8)
    if legend:
        kw["label"] = "Turbine sites"
    ax.scatter(tlon, tlat, **kw)


def _plot_wspd_pair(
    lons: np.ndarray,
    lats: np.ndarray,
    wspd_06: np.ndarray,
    wspd_18: np.ndarray,
    suptitle: str,
    panel_titles: tuple[str, str],
    out: Path,
    levels: np.ndarray,
    extent: tuple[float, float, float, float] | None = None,
    tlon: np.ndarray | None = None,
    tlat: np.ndarray | None = None,
    cmap: str = "viridis",
    show_turbines: bool = True,
) -> None:
    """Side-by-side 06 and 18 UTC wind speed maps; one shared colorbar."""
    out.parent.mkdir(parents=True, exist_ok=True)
    if extent is None:
        extent = (float(lons.min()), float(lons.max()), float(lats.min()), float(lats.max()))
    lon_min, lon_max, lat_min, lat_max = extent
    show_t = show_turbines and tlon is not None and tlat is not None and len(tlon) > 0
    fig_w, fig_h = (14.0, 5.8) if _HAS_CARTOPY else (12.0, 5.2)

    if _HAS_CARTOPY:
        fig, axes = plt.subplots(
            1,
            2,
            figsize=(fig_w, fig_h),
            subplot_kw={"projection": ccrs.PlateCarree()},
            constrained_layout=True,
        )
        ax0, ax1 = axes[0], axes[1]
        for ax, left_lbl in ((ax0, True), (ax1, False)):
            ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.5)
            ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.35, edgecolor="0.45")
            gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="gray", alpha=0.5, linestyle="--")
            gl.top_labels = False
            gl.right_labels = False
            gl.left_labels = left_lbl
            gl.bottom_labels = True
        cf0 = ax0.contourf(
            lons,
            lats,
            wspd_06,
            levels=levels,
            transform=ccrs.PlateCarree(),
            cmap=cmap,
            extend="neither",
        )
        ax1.contourf(
            lons,
            lats,
            wspd_18,
            levels=levels,
            transform=ccrs.PlateCarree(),
            cmap=cmap,
            extend="neither",
        )
        if show_t:
            _scatter_turbines_cartopy(ax0, tlon, tlat, legend=True)
            _scatter_turbines_cartopy(ax1, tlon, tlat, legend=False)
    else:
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(fig_w, fig_h), constrained_layout=True)
        for ax in (ax0, ax1):
            ax.set_xlim(lon_min, lon_max)
            ax.set_ylim(lat_min, lat_max)
            ax.set_aspect("equal")
        cf0 = ax0.contourf(lons, lats, wspd_06, levels=levels, cmap=cmap, extend="neither")
        ax1.contourf(lons, lats, wspd_18, levels=levels, cmap=cmap, extend="neither")
        ax0.set_ylabel("Latitude (°N)")
        for ax in (ax0, ax1):
            ax.set_xlabel("Longitude (°E)")
        if show_t:
            _scatter_turbines_plain(ax0, tlon, tlat, legend=True)
            _scatter_turbines_plain(ax1, tlon, tlat, legend=False)

    ax0.set_title(panel_titles[0])
    ax1.set_title(panel_titles[1])

    cb = fig.colorbar(cf0, ax=[ax0, ax1], orientation="vertical", fraction=0.046, pad=0.02)
    _style_colorbar_flat_ends(cb)
    cb.set_label("Wind speed (m s⁻¹)")
    if show_t:
        ax0.legend(loc="lower left", fontsize=8, framealpha=0.92)
    fig.suptitle(suptitle, fontsize=12)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _deficit_vmax_for_hour(
    diff: np.ndarray,
    lons: np.ndarray,
    lats: np.ndarray,
    map_extent: tuple[float, float, float, float] | None,
) -> float:
    if map_extent is not None:
        lon0, lon1, la0, la1 = map_extent
        in_zoom = (lons >= lon0) & (lons <= lon1) & (lats >= la0) & (lats <= la1)
        local = np.abs(diff[in_zoom])
        local = local[np.isfinite(local)]
        vmax = float(np.percentile(local, 98)) if local.size > 0 else 1.0
        vmax = max(vmax, 0.15)
        vmax = min(vmax, 8.0)
    else:
        vmax = float(np.percentile(np.abs(diff), 99))
        if vmax < 1e-6:
            vmax = max(float(np.abs(diff).max()), 0.1)
        vmax = min(vmax, 10.0)
    return vmax


def _plot_deficit_pair(
    lons: np.ndarray,
    lats: np.ndarray,
    diff_06: np.ndarray,
    diff_18: np.ndarray,
    suptitle: str,
    panel_titles: tuple[str, str],
    out: Path,
    vmax: float,
    extent: tuple[float, float, float, float] | None = None,
    tlon: np.ndarray | None = None,
    tlat: np.ndarray | None = None,
    n_levels: int = 45,
) -> None:
    """Side-by-side 06 and 18 UTC WFP−NWF difference maps; one shared diverging colorbar."""
    out.parent.mkdir(parents=True, exist_ok=True)
    if extent is None:
        extent = (float(lons.min()), float(lons.max()), float(lats.min()), float(lats.max()))
    lon_min, lon_max, lat_min, lat_max = extent
    levels = np.linspace(-vmax, vmax, n_levels)
    show_t = tlon is not None and len(tlon) > 0
    fig_w, fig_h = (14.0, 5.8) if _HAS_CARTOPY else (12.0, 5.2)

    if _HAS_CARTOPY:
        fig, axes = plt.subplots(
            1,
            2,
            figsize=(fig_w, fig_h),
            subplot_kw={"projection": ccrs.PlateCarree()},
            constrained_layout=True,
        )
        ax0, ax1 = axes[0], axes[1]
        for ax, left_lbl in ((ax0, True), (ax1, False)):
            ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.COASTLINE.with_scale("50m"), linewidth=0.5)
            ax.add_feature(cfeature.STATES.with_scale("50m"), linewidth=0.35, edgecolor="0.45")
            gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="gray", alpha=0.5, linestyle="--")
            gl.top_labels = False
            gl.right_labels = False
            gl.left_labels = left_lbl
            gl.bottom_labels = True
        cf0 = ax0.contourf(
            lons,
            lats,
            diff_06,
            levels=levels,
            transform=ccrs.PlateCarree(),
            cmap="RdBu_r",
            extend="both",
        )
        ax1.contourf(
            lons,
            lats,
            diff_18,
            levels=levels,
            transform=ccrs.PlateCarree(),
            cmap="RdBu_r",
            extend="both",
        )
        if show_t:
            _scatter_turbines_cartopy(ax0, tlon, tlat, legend=True)
            _scatter_turbines_cartopy(ax1, tlon, tlat, legend=False)
    else:
        fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(fig_w, fig_h), constrained_layout=True)
        for ax in (ax0, ax1):
            ax.set_xlim(lon_min, lon_max)
            ax.set_ylim(lat_min, lat_max)
            ax.set_aspect("equal")
        cf0 = ax0.contourf(lons, lats, diff_06, levels=levels, cmap="RdBu_r", extend="both")
        ax1.contourf(lons, lats, diff_18, levels=levels, cmap="RdBu_r", extend="both")
        ax0.set_ylabel("Latitude (°N)")
        for ax in (ax0, ax1):
            ax.set_xlabel("Longitude (°E)")
        if show_t:
            _scatter_turbines_plain(ax0, tlon, tlat, legend=True)
            _scatter_turbines_plain(ax1, tlon, tlat, legend=False)

    ax0.set_title(panel_titles[0])
    ax1.set_title(panel_titles[1])

    cb = fig.colorbar(cf0, ax=[ax0, ax1], orientation="vertical", fraction=0.046, pad=0.02)
    _style_colorbar_flat_ends(cb)
    cb.set_label("Wind speed difference (m s⁻¹), WFP − NWF")
    if show_t:
        ax0.legend(loc="lower left", fontsize=8, framealpha=0.92)
    fig.suptitle(suptitle, fontsize=12)
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    root = project_root()
    ap = argparse.ArgumentParser(description="Wind speed maps and WFP−NWF deficit")
    ap.add_argument("--wfp-dir", type=str, default=str(root / "data" / "wfp"))
    ap.add_argument("--nwf-dir", type=str, default=str(root / "data" / "nwf"))
    ap.add_argument("--date", type=str, default="2023-08-24")
    ap.add_argument("--hub-agl-m", type=float, default=113.0, help="Target hub height AGL (m); 100–126 typical")
    ap.add_argument(
        "--mass-level",
        type=int,
        default=None,
        help="If set, use this bottom_top index for both runs instead of AGL interpolation",
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
    ap.add_argument(
        "--windturbines",
        type=str,
        default=str(root / "inputs" / "windturbines.txt"),
        help="Turbine lat/lon for map zoom and markers on WFP and deficit maps (NWF maps have no turbines)",
    )
    args = ap.parse_args()

    wfp_dir = Path(args.wfp_dir).expanduser()
    nwf_dir = Path(args.nwf_dir).expanduser()
    if not wfp_dir.is_absolute():
        wfp_dir = root / wfp_dir
    if not nwf_dir.is_absolute():
        nwf_dir = root / nwf_dir
    fig_dir = Path(args.fig_dir).expanduser()
    if not fig_dir.is_absolute():
        fig_dir = root / fig_dir

    if not nwf_dir.is_dir() or not list(nwf_dir.glob("wrfout_d01_*")):
        print(
            "NWF (no wind farm) wrfout files not found.\n"
            "Provide the no-turbine MYNN run (same grid and physics as the WFP run) in:\n"
            f"  {nwf_dir}\n"
            "Expected files: wrfout_d01_*.\n"
            "Then re-run this script (or pass --nwf-dir PATH).",
            file=sys.stderr,
        )
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

    sample_wfp = wrfout_path_for_utc(wfp_dir, args.date, 6)
    ds0 = open_wrf(sample_wfp)
    lats = ds0["XLAT"].isel(Time=0).values
    lons = ds0["XLONG"].isel(Time=0).values

    use_level = args.mass_level is not None
    if use_level:
        k = int(args.mass_level)
        nz = int(ds0.sizes["bottom_top"])
        if not (0 <= k < nz):
            print(f"--mass-level {k} out of range [0, {nz - 1}]", file=sys.stderr)
            sys.exit(1)
        print(f"Using fixed mass level k={k} (same index for WFP and NWF)")
    else:
        k = hub_mass_level_index(ds0, target_agl_m=args.hub_agl_m)
        print(
            f"Using wind speed at {args.hub_agl_m:g} m AGL (interpolated per column; "
            f"center-column nearest mass index ~ {k})"
        )
    ds0.close()

    hours = [6, 18]
    wfp_wspd: dict[int, np.ndarray] = {}
    nwf_wspd: dict[int, np.ndarray] = {}

    for h in hours:
        pw = wrfout_path_for_utc(wfp_dir, args.date, h)
        pn = wrfout_path_for_utc(nwf_dir, args.date, h)
        d_w = open_wrf(pw)
        d_n = open_wrf(pn)
        if use_level:
            wfp_wspd[h] = destagger_wind_speed_at_level(d_w, k)
            nwf_wspd[h] = destagger_wind_speed_at_level(d_n, k)
        else:
            wfp_wspd[h] = destagger_wind_speed_at_agl(d_w, args.hub_agl_m)
            nwf_wspd[h] = destagger_wind_speed_at_agl(d_n, args.hub_agl_m)
        d_w.close()
        d_n.close()

    height_note = f"level {k}" if use_level else f"{args.hub_agl_m:g} m AGL"

    # One color scale per run type across 06 and 18 UTC
    wfp_mn = min(float(wfp_wspd[6].min()), float(wfp_wspd[18].min()))
    wfp_mx = max(float(wfp_wspd[6].max()), float(wfp_wspd[18].max()))
    levels_wfp = np.linspace(wfp_mn, wfp_mx, 22)

    nwf_mn = min(float(nwf_wspd[6].min()), float(nwf_wspd[18].min()))
    nwf_mx = max(float(nwf_wspd[6].max()), float(nwf_wspd[18].max()))
    levels_nwf = np.linspace(nwf_mn, nwf_mx, 22)

    _plot_wspd_pair(
        lons,
        lats,
        wfp_wspd[6],
        wfp_wspd[18],
        suptitle=f"Wind speed (WFP), {height_note} — {args.date}",
        panel_titles=("06:00 UTC", "18:00 UTC"),
        out=fig_dir / "wind_speed_with_farm_06_18utc.png",
        levels=levels_wfp,
        extent=map_extent,
        tlon=tlon,
        tlat=tlat,
    )
    _plot_wspd_pair(
        lons,
        lats,
        nwf_wspd[6],
        nwf_wspd[18],
        suptitle=f"Wind speed (NWF), {height_note} — {args.date}",
        panel_titles=("06:00 UTC", "18:00 UTC"),
        out=fig_dir / "wind_speed_no_farm_06_18utc.png",
        levels=levels_nwf,
        extent=map_extent,
        show_turbines=False,
    )

    diff6 = wfp_wspd[6] - nwf_wspd[6]
    diff18 = wfp_wspd[18] - nwf_wspd[18]
    vmax6 = _deficit_vmax_for_hour(diff6, lons, lats, map_extent)
    vmax18 = _deficit_vmax_for_hour(diff18, lons, lats, map_extent)
    vmax_pair = max(vmax6, vmax18)

    _plot_deficit_pair(
        lons,
        lats,
        diff6,
        diff18,
        suptitle=f"Wind speed difference (WFP − NWF), {height_note} — {args.date}",
        panel_titles=("06:00 UTC", "18:00 UTC"),
        out=fig_dir / "wake_deficit_06_18utc.png",
        vmax=vmax_pair,
        extent=map_extent,
        tlon=tlon,
        tlat=tlat,
    )
    print(
        f"Saved deficit pair: shared color scale +/-{vmax_pair:.2f} m/s "
        f"(per-hour 98th-pct suggestions: 06 UTC {vmax6:.2f}, 18 UTC {vmax18:.2f}); "
        f"max|diff| ~ {float(np.abs(diff6).max()):.3f} (06), {float(np.abs(diff18).max()):.3f} (18) m/s"
    )

    print(f"Figures written under {fig_dir}")


if __name__ == "__main__":
    main()
