#!/usr/bin/env python3
"""
Katrina nesting comparison plots (d01 vs d02)
Two variables:
  1) RAINNC (accumulated precip, mm)
  2) WS10 (10 m wind speed = sqrt(U10^2 + V10^2), m/s)

Outputs (PNG):
  - MAP_RAINNC_*.png
  - MAP_WS10_*.png
  - TS_RAINNC_*.png
  - TS_WS10_*.png
"""

import os
import warnings
warnings.filterwarnings("ignore")

from pathlib import Path
import re
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

import cartopy.crs as ccrs
import cartopy.feature as cfeature


# =========================
# USER SETTINGS (EDIT THESE)
# =========================
# Folder containing wrfout_d01_* and wrfout_d02_* (set KATRINA_DATA, or use ./data)
PROJECT_DIR = Path(__file__).resolve().parent.parent
CASE_DIR = Path(os.environ.get("KATRINA_DATA", PROJECT_DIR / "data"))
OUT_DIR  = PROJECT_DIR / "figures"

# time-series point: near the Mississippi River delta, along the storm's approach
LAT0 = 29.350
LON0 = -89.520

# If you want incremental rain instead of cumulative RAINNC:
PLOT_INCREMENTAL_RAIN = False


# =========================
# HELPERS
# =========================
def list_wrf_files(case_dir: Path, domain: str):
    files = sorted(case_dir.glob(f"wrfout_{domain}_*"))
    if not files:
        raise FileNotFoundError(f"No files found like wrfout_{domain}_* in {case_dir}")
    return files

def timestamp_from_name(p: Path, domain: str):
    # expects: wrfout_d01_2005-08-29_00_00_00
    m = re.match(rf"wrfout_{domain}_(\d{{4}}-\d{{2}}-\d{{2}}_\d{{2}}_\d{{2}}_\d{{2}})$", p.name)
    if not m:
        return None
    return m.group(1)

def choose_common_time(files_d01, files_d02):
    t1 = {timestamp_from_name(p, "d01") for p in files_d01}
    t2 = {timestamp_from_name(p, "d02") for p in files_d02}
    t1.discard(None); t2.discard(None)
    common = sorted(t1.intersection(t2))
    if not common:
        raise RuntimeError(
            "No common timestamps found between d01 and d02 filenames. "
            "Double-check you downloaded overlapping times for both domains."
        )
    # choose a “middle-ish” time for a representative map
    return common[len(common)//2]

def open_single(path: Path) -> xr.Dataset:
    return xr.open_dataset(path, decode_times=False)

def open_all(files):
    # Build combined dataset without dask so it works even when dask is missing or fails
    parts = []
    for f in files:
        ds = xr.open_dataset(str(f), decode_times=False, engine="netcdf4")
        parts.append(ds)
    out = xr.concat(parts, dim="Time")
    for ds in parts:
        ds.close()
    return out

def parse_times_from_ds(ds: xr.Dataset) -> pd.DatetimeIndex:
    """
    Parse WRF Times char array into pandas times.
    Handles both:
      2005-08-29_00:00:00
      2005-08-29_00_00_00
    """
    if "Times" not in ds.variables:
        raise ValueError("Dataset has no 'Times' variable (needed for time axis).")

    arr = ds["Times"].values  # (Time, DateStrLen)
    times_str = []
    for row in arr:
        s = "".join([x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in row]).strip()
        times_str.append(s)

    # Try both common formats (underscores in time first to match typical WRF filenames)
    t = pd.to_datetime(times_str, format="%Y-%m-%d_%H_%M_%S", errors="coerce")
    if t.isna().all():
        t = pd.to_datetime(times_str, format="%Y-%m-%d_%H:%M:%S", errors="coerce")

    if t.isna().any():
        t = pd.to_datetime(times_str, errors="coerce")

    return t


def times_from_file_list(files, domain: str) -> pd.DatetimeIndex:
    """
    Build time index from filenames so it exactly matches the concatenated data order.
    Files are in sorted order; concat preserves that order, so these times match the data.
    """
    stamps = [timestamp_from_name(p, domain) for p in files]
    if any(s is None for s in stamps):
        raise ValueError(f"Some {domain} filenames did not match expected wrfout_{domain}_YYYY-MM-DD_HH_MM_SS pattern.")
    return pd.to_datetime(stamps, format="%Y-%m-%d_%H_%M_%S")

def nearest_ij(lat2d, lon2d, lat0, lon0):
    d2 = (lat2d - lat0)**2 + (lon2d - lon0)**2
    j, i = np.unravel_index(np.argmin(d2), d2.shape)
    return int(i), int(j)

def add_map_features(ax):
    ax.coastlines(resolution="50m", linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5)
    ax.add_feature(cfeature.STATES, linewidth=0.3, edgecolor="gray")
    gl = ax.gridlines(draw_labels=True, linewidth=0.2, alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False

def plot_maps(lat1, lon1, v1, lat2, lon2, v2, extent, title, cbar_label, outpath):
    proj = ccrs.PlateCarree()
    fig, axs = plt.subplots(1, 2, figsize=(14, 6), subplot_kw={"projection": proj},
                            constrained_layout=True)

    vmin = np.nanmin([np.nanmin(v1), np.nanmin(v2)])
    vmax = np.nanmax([np.nanmax(v1), np.nanmax(v2)])

    ax = axs[0]
    ax.set_extent(extent, crs=proj)
    add_map_features(ax)
    m = ax.contourf(lon1, lat1, v1, levels=20, transform=proj, vmin=vmin, vmax=vmax)
    ax.set_title("d01 (coarse)")

    ax = axs[1]
    ax.set_extent(extent, crs=proj)
    add_map_features(ax)
    ax.contourf(lon2, lat2, v2, levels=20, transform=proj, vmin=vmin, vmax=vmax)
    ax.set_title("d02 (fine)")

    cbar = fig.colorbar(m, ax=axs, orientation="vertical", shrink=0.9, pad=0.02)
    cbar.set_label(cbar_label)

    fig.suptitle(title, fontsize=14)
    fig.savefig(outpath, dpi=200)
    plt.close(fig)

def plot_timeseries(time1, y1, time2, y2, title, ylabel, outpath):
    plt.figure(figsize=(12, 5))
    plt.plot(time1, y1, label="d01 (coarse)")
    plt.plot(time2, y2, label="d02 (fine)")
    plt.xlabel("Time")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()


# =========================
# MAIN
# =========================
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    files_d01 = list_wrf_files(CASE_DIR, "d01")
    files_d02 = list_wrf_files(CASE_DIR, "d02")
    print(f"Found {len(files_d01)} d01 files and {len(files_d02)} d02 files.")

    # Pick a common time automatically for the MAP plots
    t_common = choose_common_time(files_d01, files_d02)
    f1 = CASE_DIR / f"wrfout_d01_{t_common}"
    f2 = CASE_DIR / f"wrfout_d02_{t_common}"
    print(f"Using common timestamp for MAP plots: {t_common}")
    print(f"  d01 file: {f1.name}")
    print(f"  d02 file: {f2.name}")

    # ---------- MAPS ----------
    ds1 = open_single(f1)
    ds2 = open_single(f2)

    lat1 = ds1["XLAT"].isel(Time=0).values
    lon1 = ds1["XLONG"].isel(Time=0).values
    lat2 = ds2["XLAT"].isel(Time=0).values
    lon2 = ds2["XLONG"].isel(Time=0).values

    # d02 bounds for extent + d01 masking
    latmin, latmax = np.nanmin(lat2), np.nanmax(lat2)
    lonmin, lonmax = np.nanmin(lon2), np.nanmax(lon2)
    extent = [lonmin, lonmax, latmin, latmax]
    mask_overlap = (lat1 >= latmin) & (lat1 <= latmax) & (lon1 >= lonmin) & (lon1 <= lonmax)

    # Variable 1 (MAP): RAINNC
    rain1 = ds1["RAINNC"].isel(Time=0).values.astype(float)
    rain2 = ds2["RAINNC"].isel(Time=0).values.astype(float)
    rain1 = np.where(mask_overlap, rain1, np.nan)

    out_map_rain = OUT_DIR / f"MAP_RAINNC_{t_common}.png"
    plot_maps(
        lat1, lon1, rain1,
        lat2, lon2, rain2,
        extent=extent,
        title=f"RAINNC (Accumulated Precip) at {t_common}",
        cbar_label="mm",
        outpath=out_map_rain
    )
    print(f"Saved: {out_map_rain.name}")

    # Variable 2 (MAP): WS10 = sqrt(U10^2 + V10^2)
    u10_1 = ds1["U10"].isel(Time=0).values.astype(float)
    v10_1 = ds1["V10"].isel(Time=0).values.astype(float)
    ws1 = np.sqrt(u10_1**2 + v10_1**2)

    u10_2 = ds2["U10"].isel(Time=0).values.astype(float)
    v10_2 = ds2["V10"].isel(Time=0).values.astype(float)
    ws2 = np.sqrt(u10_2**2 + v10_2**2)

    ws1 = np.where(mask_overlap, ws1, np.nan)

    out_map_ws = OUT_DIR / f"MAP_WS10_{t_common}.png"
    plot_maps(
        lat1, lon1, ws1,
        lat2, lon2, ws2,
        extent=extent,
        title=f"10 m Wind Speed (WS10) at {t_common}",
        cbar_label="m/s",
        outpath=out_map_ws
    )
    print(f"Saved: {out_map_ws.name}")

    # ---------- TIME SERIES ----------
    ds1_all = open_all(files_d01)
    ds2_all = open_all(files_d02)

    # Use times from filenames so the time axis exactly matches the concatenated data order
    time1 = times_from_file_list(files_d01, "d01")
    time2 = times_from_file_list(files_d02, "d02")

    # nearest indices
    lat1_0 = ds1_all["XLAT"].isel(Time=0).values
    lon1_0 = ds1_all["XLONG"].isel(Time=0).values
    lat2_0 = ds2_all["XLAT"].isel(Time=0).values
    lon2_0 = ds2_all["XLONG"].isel(Time=0).values

    i1, j1 = nearest_ij(lat1_0, lon1_0, LAT0, LON0)
    i2, j2 = nearest_ij(lat2_0, lon2_0, LAT0, LON0)
    print(f"Nearest grid point to ({LAT0:.3f},{LON0:.3f}): d01(i,j)=({i1},{j1}), d02(i,j)=({i2},{j2})")

    # TS 1: RAINNC
    rain_ts1 = ds1_all["RAINNC"][:, j1, i1].values.astype(float)
    rain_ts2 = ds2_all["RAINNC"][:, j2, i2].values.astype(float)

    if PLOT_INCREMENTAL_RAIN:
        rain_ts1 = np.diff(rain_ts1, prepend=rain_ts1[0])
        rain_ts2 = np.diff(rain_ts2, prepend=rain_ts2[0])
        ylabel = "Incremental rain per output interval (mm)"
        tag = "RAININC"
        title = f"Incremental Rain (from RAINNC) at ({LAT0:.3f},{LON0:.3f})"
    else:
        ylabel = "Accumulated rain (mm)"
        tag = "RAINNC"
        title = f"RAINNC (Accumulated) at ({LAT0:.3f},{LON0:.3f})"

    out_ts_rain = OUT_DIR / f"TS_{tag}_{LAT0:.3f}_{LON0:.3f}.png"
    plot_timeseries(time1, rain_ts1, time2, rain_ts2, title, ylabel, out_ts_rain)
    print(f"Saved: {out_ts_rain.name}")

    # TS 2: WS10
    u1 = ds1_all["U10"][:, j1, i1].values.astype(float)
    v1 = ds1_all["V10"][:, j1, i1].values.astype(float)
    ws_ts1 = np.sqrt(u1**2 + v1**2)

    u2 = ds2_all["U10"][:, j2, i2].values.astype(float)
    v2 = ds2_all["V10"][:, j2, i2].values.astype(float)
    ws_ts2 = np.sqrt(u2**2 + v2**2)

    out_ts_ws = OUT_DIR / f"TS_WS10_{LAT0:.3f}_{LON0:.3f}.png"
    plot_timeseries(
        time1, ws_ts1, time2, ws_ts2,
        title=f"10 m Wind Speed (WS10) at ({LAT0:.3f},{LON0:.3f})",
        ylabel="Wind speed (m/s)",
        outpath=out_ts_ws
    )
    print(f"Saved: {out_ts_ws.name}")

    print("\nDONE ✅  All plots saved to:")
    print(OUT_DIR)

if __name__ == "__main__":
    main()