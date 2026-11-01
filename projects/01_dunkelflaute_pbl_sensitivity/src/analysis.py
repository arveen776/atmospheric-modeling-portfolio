# =============================================================================
# Dunkelflaute PBL Scheme Sensitivity — EN.530.462 Final Project
# YSU vs MYNN2.5 | Maryland | Jan 17–19 2024
#
# Approach:
#   - xarray + numpy only (NO wrf-python)
#   - Height from (PH+PHB)/9.81 - HGT
#   - Winds from destaggered U, V
#   - Times parsed from WRF Times byte array
#   - ERA5 100m winds as the reference (loaded via cfgrib)
#
# Required packages: see requirements.txt at the repository root.
#
# Run from this project's folder:  python src/analysis.py
# Point DUNKELFLAUTE_DATA at the data folder if it is not ./data:
# data/
# ├── geo_em.d01.nc
# ├── geo_em.d02.nc
# ├── era5_100m_winds.grib       ← ERA5 u100/v100, Jan 16-20 2024
# ├── YSU/
# │   └── wrfout_d02_2024-01-16_00_00_00
# └── MYNN/
#     ├── wrfout_d02_2024-01-16_00_00_00
#     └── wrfout_d02_2024-01-19_01_00_00
# =============================================================================

# %%
# =============================================================================
# CELL 1 — Imports
# =============================================================================
import os

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import TwoSlopeNorm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pathlib import Path
from scipy.interpolate import interp1d
import warnings
warnings.filterwarnings("ignore")

print("Imports OK")


# %%
# =============================================================================
# CELL 2 — File paths and settings  (only edit this cell)
# =============================================================================

PROJECT_DIR    = Path(__file__).resolve().parent.parent
DATA_DIR       = Path(os.environ.get("DUNKELFLAUTE_DATA", PROJECT_DIR / "data"))
FIG_DIR        = PROJECT_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)
YSU_DIR        = DATA_DIR / "YSU"
MYNN_DIR       = DATA_DIR / "MYNN"
GEO_D01        = DATA_DIR / "geo_em.d01.nc"
GEO_D02        = DATA_DIR / "geo_em.d02.nc"
ERA5_GRIB_FILE = DATA_DIR / "era5_100m_winds.grib"   # ERA5 u100/v100

# Analysis window (after spin-up, trimmed to MYNN end)
T_START = pd.Timestamp("2024-01-17 00:00:00")
T_END   = pd.Timestamp("2024-01-19 17:00:00")

# Evaluation point — Chesapeake Bay / MD coastal plain
LAT0 = 38.80
LON0 = -76.00

# Hub height of interest (m AGL)
HUB_HT = 100.0

# Colors — consistent across all plots
COL_YSU  = "#2471A3"   # blue
COL_MYNN = "#C0392B"   # red
COL_ERA5 = "#1A1A1A"   # near-black

print(f"Period : {T_START} → {T_END}")
print(f"Point  : {LAT0}N, {abs(LON0)}W")
print(f"Hub Ht : {HUB_HT} m AGL")


# %%
# =============================================================================
# CELL 3 — Open WRF files with xarray
# =============================================================================

def open_wrf_ds(directory: Path, domain="d02"):
    """Open all wrfout files for a domain, concatenated along Time."""
    files = sorted(directory.glob(f"wrfout_{domain}_*"))
    if not files:
        raise FileNotFoundError(f"No wrfout_{domain}_* in {directory}")
    ds = xr.open_mfdataset(
        files,
        combine="nested",
        concat_dim="Time",
        decode_times=False,
        parallel=False,
        engine="netcdf4",
    )
    return ds

def add_time_coord(ds):
    """Parse WRF Times byte array → pandas DatetimeIndex, add as 'time' coord."""
    tr = np.asarray(ds["Times"].values)

    def _row_to_str(r):
        r = np.asarray(r).ravel()
        if r.size == 0:
            return ""
        if r.dtype.kind in "iu":
            return r.astype(np.uint8, copy=False).tobytes().decode(
                "ascii", errors="replace"
            ).strip("\x00").strip()
        if r.dtype.kind == "S" and r.dtype.itemsize > 1 and r.size == 1:
            return np.bytes_(r.flat[0]).decode("utf-8", errors="replace").strip()
        if r.dtype.kind in "SU":
            return b"".join(r.tolist()).decode("utf-8", errors="replace").strip("\x00").strip()
        return str(r.flat[0]).strip()

    if tr.ndim == 1 and tr.dtype.kind == "S" and tr.dtype.itemsize > 1:
        times_str = [np.bytes_(x).decode("utf-8", errors="replace").strip() for x in tr]
    elif tr.ndim == 1 and tr.dtype.kind == "U":
        times_str = [str(x).strip() for x in tr]
    elif tr.ndim == 1 and tr.dtype.kind == "O":
        times_str = [str(x).strip() for x in tr]
    elif tr.ndim >= 2:
        times_str = [_row_to_str(tr[i]) for i in range(tr.shape[0])]
    else:
        times_str = [_row_to_str(tr)]

    times_dt = pd.to_datetime(times_str, format="%Y-%m-%d_%H:%M:%S")
    ds = ds.assign_coords(time=("Time", times_dt))

    # A restart can write the same valid times to a second wrfout file. Keep the
    # first copy of each so no hour is counted twice in means or error metrics.
    keep = ~times_dt.duplicated(keep="first")
    if not keep.all():
        print(f"  Dropping {int((~keep).sum())} duplicate timesteps")
        ds = ds.isel(Time=np.where(keep)[0])
    if not ds["time"].to_index().is_monotonic_increasing:
        raise ValueError("WRF times are not monotonically increasing")
    return ds

print("Opening YSU dataset...")
ds_ysu = add_time_coord(open_wrf_ds(YSU_DIR))
print("Opening MYNN dataset...")
ds_mynn = add_time_coord(open_wrf_ds(MYNN_DIR))

print(f"\nYSU  timesteps: {ds_ysu.sizes['Time']}  "
      f"({ds_ysu.time.values[0]} → {ds_ysu.time.values[-1]})")
print(f"MYNN timesteps: {ds_mynn.sizes['Time']}  "
      f"({ds_mynn.time.values[0]} → {ds_mynn.time.values[-1]})")


# %%
# =============================================================================
# CELL 4 — Helper functions
# =============================================================================

def time_mask(ds, t0, t1):
    """Boolean mask for the Time dimension."""
    t = ds["time"].values
    return (t >= np.datetime64(t0)) & (t <= np.datetime64(t1))

def nearest_ij(ds, lat0, lon0):
    """Return (j, i) nearest grid point."""
    lats = ds["XLAT"].isel(Time=0).values
    lons = ds["XLONG"].isel(Time=0).values
    dist = np.sqrt((lats - lat0)**2 + (lons - lon0)**2)
    jj, ii = np.unravel_index(np.argmin(dist), dist.shape)
    print(f"  Nearest grid: j={jj}, i={ii}  "
          f"→ {float(lats[jj,ii]):.3f}N, {float(lons[jj,ii]):.3f}E")
    return jj, ii

def geopotential_height_agl(ds, j, i):
    """
    Height AGL at mass levels for a single grid point.
    z = (PH+PHB)/g  destaggered  minus terrain.
    Returns DataArray(Time, bottom_top).
    """
    g   = 9.81
    ph  = ds["PH"].isel(south_north=j, west_east=i)
    phb = ds["PHB"].isel(south_north=j, west_east=i)
    z_stag = (ph + phb) / g
    z_mass = 0.5 * (z_stag.isel(bottom_top_stag=slice(0, -1)) +
                    z_stag.isel(bottom_top_stag=slice(1,  None)))
    hgt = ds["HGT"].isel(Time=0, south_north=j, west_east=i)
    return z_mass - hgt

def geopotential_height_agl_3d(ds):
    """Full-domain height AGL at mass levels. Returns (Time, nz, ny, nx)."""
    g   = 9.81
    ph  = ds["PH"].values
    phb = ds["PHB"].values
    z_stag = (ph + phb) / g
    z_mass = 0.5 * (z_stag[:, :-1, :, :] + z_stag[:, 1:, :, :])
    hgt    = ds["HGT"].values[0, :, :]
    return z_mass - hgt[np.newaxis, np.newaxis, :, :]

def destagger_uv(ds, j, i):
    """
    Destagger U and V winds to mass point at (j, i).
    Returns (wspd, wdir) DataArrays (Time, bottom_top).
    """
    U = ds["U"].isel(south_north=j, west_east_stag=slice(i, i+2))
    U = 0.5 * (U.isel(west_east_stag=0) + U.isel(west_east_stag=1))
    V = ds["V"].isel(south_north_stag=slice(j, j+2), west_east=i)
    V = 0.5 * (V.isel(south_north_stag=0) + V.isel(south_north_stag=1))
    wspd = np.hypot(U, V)
    wdir = (270.0 - np.degrees(np.arctan2(V, U))) % 360.0
    return wspd, wdir

def interp_to_height(z_agl, field, z_target):
    """
    Interpolate (Time, nz) field to z_target m AGL at each timestep.
    Returns (Time,) array.
    """
    out = []
    for tt in range(field.shape[0]):
        z = z_agl.isel(Time=tt).values
        f = field.isel(Time=tt).values
        out.append(float(np.interp(z_target, z, f)))
    return np.array(out)

def interp_profile_mean(z_agl, field, mask, z_common):
    """
    Time-average profile interpolated to z_common grid.
    mask selects timesteps. Returns (nz,) array.
    """
    profiles = []
    t_indices = np.where(mask)[0]
    for tt in t_indices:
        z = z_agl.isel(Time=int(tt)).values
        f = field.isel(Time=int(tt)).values
        profiles.append(np.interp(z_common, z, f))
    return np.mean(profiles, axis=0)

def add_map_features(ax):
    ax.add_feature(cfeature.STATES.with_scale("10m"),
                   linewidth=0.8, edgecolor="black", zorder=3)
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=1.0, zorder=3)
    ax.add_feature(cfeature.LAKES.with_scale("10m"),
                   facecolor="#d6eaf8", edgecolor="steelblue",
                   linewidth=0.3, zorder=2)
    ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#d6eaf8", zorder=1)
    ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#f5f5f0", zorder=1)

print("Helper functions ready.")


# %%
# =============================================================================
# CELL 5 — Find nearest WRF grid point and extract time info
# =============================================================================
print("YSU nearest point:")
jj, ii = nearest_ij(ds_ysu, LAT0, LON0)

ysu_times  = ds_ysu["time"].values
mynn_times = ds_mynn["time"].values

mask_ysu  = time_mask(ds_ysu,  T_START, T_END)
mask_mynn = time_mask(ds_mynn, T_START, T_END)

ysu_times_trim  = ysu_times[mask_ysu]
mynn_times_trim = mynn_times[mask_mynn]

print(f"\nYSU  in window: {mask_ysu.sum()} timesteps")
print(f"MYNN in window: {mask_mynn.sum()} timesteps")


# %%
# =============================================================================
# CELL 5b — Load ERA5 100-m winds as the reference
# Uses xarray + cfgrib backend — no wrf-python
# =============================================================================
print("\nLoading ERA5 100-m wind reference...")

ds_era5 = xr.open_dataset(ERA5_GRIB_FILE, engine="cfgrib",
                          backend_kwargs={"indexpath": ""})

# ERA5 lat/lon arrays — lons already negative in this file (-90 to -65)
era5_lats = ds_era5.latitude.values    # 1-D
era5_lons = ds_era5.longitude.values   # 1-D, negative convention

# Find nearest ERA5 grid point to evaluation location
era5_dist = np.sqrt(
    (era5_lats[:, None] - LAT0)**2 + (era5_lons[None, :] - LON0)**2
)
jj_e, ii_e = np.unravel_index(np.argmin(era5_dist), era5_dist.shape)
print(f"  Nearest ERA5 point: {float(era5_lats[jj_e]):.3f}N, "
      f"{float(era5_lons[ii_e]):.3f}E")

# Extract u100, v100 at that point and compute wind speed
u100 = ds_era5["u100"].values[:, jj_e, ii_e]   # (time,)
v100 = ds_era5["v100"].values[:, jj_e, ii_e]
era5_wspd = np.hypot(u100, v100)

# ERA5 time coordinate (valid_time is the actual datetime for forecast data)
era5_times_raw = (ds_era5.valid_time.values
                  if "valid_time" in ds_era5.coords
                  else ds_era5.time.values)
era5_times = pd.to_datetime(era5_times_raw)

# Trim to analysis window
era5_mask       = (era5_times >= T_START) & (era5_times <= T_END)
era5_times_trim = era5_times[era5_mask]
era5_ws_trim    = era5_wspd[era5_mask]

print(f"  ERA5 100m mean (in window) = {era5_ws_trim.mean():.2f} m/s")
print(f"  ERA5 timesteps in window   = {era5_mask.sum()}")


# %%
# =============================================================================
# CELL 6 — Extract WRF hub-height time series at evaluation point
# =============================================================================
print("\nExtracting WRF hub-height time series...")

# YSU
z_agl_ysu   = geopotential_height_agl(ds_ysu,  jj, ii)
wspd_ysu, _ = destagger_uv(ds_ysu,  jj, ii)
ysu_hub_ts  = interp_to_height(
    z_agl_ysu.isel(Time=mask_ysu),
    wspd_ysu.isel(Time=mask_ysu),
    HUB_HT
)

# MYNN
z_agl_mynn   = geopotential_height_agl(ds_mynn, jj, ii)
wspd_mynn, _ = destagger_uv(ds_mynn, jj, ii)
mynn_hub_ts  = interp_to_height(
    z_agl_mynn.isel(Time=mask_mynn),
    wspd_mynn.isel(Time=mask_mynn),
    HUB_HT
)

print(f"YSU  100m mean = {ysu_hub_ts.mean():.2f} m/s")
print(f"MYNN 100m mean = {mynn_hub_ts.mean():.2f} m/s")
print(f"ERA5 100m mean = {era5_ws_trim.mean():.2f} m/s")


# %%
# =============================================================================
# CELL 7 — PLOT 0: Domain map using geo_em files
# =============================================================================
print("\nPlot 0: Domain map...")

geo_d01 = xr.open_dataset(GEO_D01, engine="netcdf4")
geo_d02 = xr.open_dataset(GEO_D02, engine="netcdf4")

lat_d01 = geo_d01["XLAT_M"].values[0]
lon_d01 = geo_d01["XLONG_M"].values[0]
hgt_d01 = geo_d01["HGT_M"].values[0]
lat_d02 = geo_d02["XLAT_M"].values[0]
lon_d02 = geo_d02["XLONG_M"].values[0]

def domain_box(lat2d, lon2d):
    lx = [lon2d[0,0],  lon2d[0,-1],  lon2d[-1,-1], lon2d[-1,0],  lon2d[0,0]]
    ly = [lat2d[0,0],  lat2d[0,-1],  lat2d[-1,-1], lat2d[-1,0],  lat2d[0,0]]
    return lx, ly

# Build Lambert projection from WRF global attributes
truelat1  = float(ds_ysu.attrs.get("TRUELAT1",  38.8))
truelat2  = float(ds_ysu.attrs.get("TRUELAT2",  38.8))
stand_lon = float(ds_ysu.attrs.get("STAND_LON", -76.0))
proj = ccrs.LambertConformal(
    central_longitude=stand_lon,
    central_latitude=truelat1,
    standard_parallels=(truelat1, truelat2)
)

fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={"projection": proj})
cf = ax.pcolormesh(lon_d01, lat_d01, hgt_d01,
                   transform=ccrs.PlateCarree(),
                   cmap="terrain", vmin=0, vmax=900, alpha=0.65, zorder=2)
plt.colorbar(cf, ax=ax, label="Terrain Elevation (m)", shrink=0.6, pad=0.02)
add_map_features(ax)

lxd1, lyd1 = domain_box(lat_d01, lon_d01)
lxd2, lyd2 = domain_box(lat_d02, lon_d02)
ax.plot(lxd1, lyd1, "b-", linewidth=2.5, transform=ccrs.PlateCarree(),
        label="d01 — 12 km (100×100)", zorder=5)
ax.plot(lxd2, lyd2, "r-", linewidth=2.5, transform=ccrs.PlateCarree(),
        label="d02 — 4 km (100×100)", zorder=5)

ax.text(lon_d01[0,0]+0.4, lat_d01[-1,0]-1.2, "d01  12 km",
        color="blue", fontsize=11, fontweight="bold",
        transform=ccrs.PlateCarree(),
        bbox=dict(facecolor="white", alpha=0.75, edgecolor="none"), zorder=6)
ax.text(lon_d02[0,0]+0.3, lat_d02[-1,0]-0.9, "d02  4 km",
        color="red", fontsize=11, fontweight="bold",
        transform=ccrs.PlateCarree(),
        bbox=dict(facecolor="white", alpha=0.75, edgecolor="none"), zorder=6)

ax.plot(LON0, LAT0, "k*", markersize=16, transform=ccrs.PlateCarree(),
        zorder=7, label=f"Eval. point ({LAT0}N, {abs(LON0)}W)")

ax.set_extent([lon_d01[0,0]-0.5, lon_d01[0,-1]+0.5,
               lat_d01[0,0]-0.5, lat_d01[-1,0]+0.5],
              crs=ccrs.PlateCarree())
ax.legend(loc="lower right", fontsize=10, framealpha=0.9)
ax.set_title("WRF Model Domains — Mid-Atlantic Dunkelflaute Study\n"
             "ERA5 IC/BC | Lambert Conformal | Jan 2024",
             fontsize=12, fontweight="bold")
gl = ax.gridlines(draw_labels=True, linewidth=0.4, color="gray",
                  alpha=0.5, linestyle="--")
gl.top_labels = False
gl.right_labels = False

plt.tight_layout()
plt.savefig(FIG_DIR / "plot0_domain_map.png", dpi=150, bbox_inches="tight")
plt.show()
print("  Saved: plot0_domain_map.png")


# %%
# =============================================================================
# CELL 8 — PLOT 1: Time series with ERA5 reference
# Three lines: YSU (blue), MYNN2.5 (red), ERA5 (black dashed)
# =============================================================================
print("Plot 1: Time series with ERA5 reference...")

fig, ax = plt.subplots(figsize=(13, 5.5))

# ERA5 first so it sits behind WRF lines but above shading
ax.plot(era5_times_trim, era5_ws_trim, color=COL_ERA5,
        linewidth=2.0, linestyle="-.", label="ERA5 (reference)", zorder=3)

ax.plot(ysu_times_trim,  ysu_hub_ts,  color=COL_YSU,
        linewidth=2.2, label="YSU", zorder=4)
ax.plot(mynn_times_trim, mynn_hub_ts, color=COL_MYNN,
        linewidth=2.2, label="MYNN2.5", zorder=4)

# Dunkelflaute threshold
ax.axhline(4.0, color="gray", linestyle="--", linewidth=1.2, alpha=0.7,
           label="Low-wind threshold (4 m s$^{-1}$)")

# Shade nighttime hours (~23Z–11Z for January Maryland)
for d in [17, 18, 19]:
    ns = pd.Timestamp(f"2024-01-{d-1:02d} 23:00") if d > 17 else T_START
    ne = pd.Timestamp(f"2024-01-{d:02d} 11:00")
    if ns >= T_START and ne <= T_END:
        ax.axvspan(ns, ne, alpha=0.06, color="navy",
                   label="Nighttime" if d == 18 else "")

# Annotate max difference between the two WRF schemes
min_len  = min(len(ysu_hub_ts), len(mynn_hub_ts))
diff_arr = np.abs(ysu_hub_ts[:min_len] - mynn_hub_ts[:min_len])
idx_max  = np.argmax(diff_arr)
max_time = pd.Timestamp(ysu_times_trim[idx_max]).strftime("%b %d %H:00Z")
ax.text(0.02, 0.05,
        f"Max |Δ YSU–MYNN| = {diff_arr[idx_max]:.2f} m s$^{{-1}}$ at {max_time}",
        transform=ax.transAxes, fontsize=10,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.85))

ax.set_xlabel("Date / Time (UTC)", fontsize=12)
ax.set_ylabel("100-m Wind Speed (m s$^{-1}$)", fontsize=12)
ax.set_title(f"100-m Hub-Height Wind Speed During Mid-Atlantic Dunkelflaute\n"
             f"YSU vs MYNN2.5 vs ERA5 reference | Jan 17–19 2024 | "
             f"{LAT0}°N, {abs(LON0)}°W",
             fontsize=11)
ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%H:00Z"))
ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
ax.set_xlim(T_START, T_END)
ax.set_ylim(bottom=0)
ax.legend(fontsize=10, framealpha=0.9, loc="upper right")
ax.grid(True, alpha=0.3, linestyle=":")

plt.tight_layout()
plt.savefig(FIG_DIR / "plot1_timeseries_100m_wind.png", dpi=150, bbox_inches="tight")
plt.show()
print("  Saved: plot1_timeseries_100m_wind.png")


# %%
# =============================================================================
# CELL 9 — Extract time-averaged vertical profiles
# =============================================================================
print("Extracting vertical profiles...")

z_common = np.arange(0, 801, 10)   # 0–800 m every 10 m

ws_prof_ysu = interp_profile_mean(
    z_agl_ysu,
    wspd_ysu.rename({"bottom_top": "bottom_top"}),
    mask_ysu,
    z_common
)
ws_prof_mynn = interp_profile_mean(
    z_agl_mynn,
    wspd_mynn.rename({"bottom_top": "bottom_top"}),
    mask_mynn,
    z_common
)

# ERA5 reference: time-averaged 100m wind speed
# ERA5 only provides one level (100m) so we add a single-point marker
era5_mean_100m = float(era5_ws_trim.mean())

print(f"  YSU  at {HUB_HT}m : {np.interp(HUB_HT, z_common, ws_prof_ysu):.2f} m/s")
print(f"  MYNN at {HUB_HT}m : {np.interp(HUB_HT, z_common, ws_prof_mynn):.2f} m/s")
print(f"  ERA5 at {HUB_HT}m : {era5_mean_100m:.2f} m/s")


# %%
# =============================================================================
# CELL 10 — PLOT 2: Vertical profile with ERA5 marker at hub height
# ERA5 has only 100m — plotted as a horizontal diamond showing "target" value
# =============================================================================
print("Plot 2: Vertical profile...")

fig, ax = plt.subplots(figsize=(6, 9))

ax.plot(ws_prof_ysu,  z_common, color=COL_YSU,  linewidth=2.5, label="YSU")
ax.plot(ws_prof_mynn, z_common, color=COL_MYNN, linewidth=2.5, label="MYNN2.5")

# ERA5 reference marker — diamond at 100 m
ax.plot(era5_mean_100m, HUB_HT, marker="D", markersize=11,
        color=COL_ERA5, zorder=6, label=f"ERA5 at 100 m ({era5_mean_100m:.2f} m s$^{{-1}}$)")
# Vertical dotted line from 0 to 100m at ERA5 value to guide the eye
ax.vlines(era5_mean_100m, 0, HUB_HT, color=COL_ERA5,
          linewidth=0.8, linestyle=":", alpha=0.5)

# Hub height line
ax.axhline(HUB_HT, color="black", linestyle="--", linewidth=1.3, alpha=0.8)
ax.text(0.02, HUB_HT + 12, f"Hub height ({int(HUB_HT)} m)",
        fontsize=9, transform=ax.get_yaxis_transform())

# Rotor layer shading 40–160 m
ax.axhspan(40, 160, alpha=0.07, color="green")
ax.text(0.97, 100, "Rotor\nlayer", transform=ax.get_yaxis_transform(),
        fontsize=8, color="darkgreen", ha="right", va="center")

# Annotate difference at hub height
ws_y  = float(np.interp(HUB_HT, z_common, ws_prof_ysu))
ws_m  = float(np.interp(HUB_HT, z_common, ws_prof_mynn))
diff_h = ws_y - ws_m
higher = "YSU higher" if diff_h > 0 else "MYNN2.5 higher"
ax.annotate(f"Δ = {abs(diff_h):.2f} m s$^{{-1}}$\n({higher})",
            xy=(max(ws_y, ws_m), HUB_HT),
            xytext=(25, -40), textcoords="offset points",
            arrowprops=dict(arrowstyle="->", lw=1.2),
            fontsize=9,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="wheat", alpha=0.85))

ax.set_xlabel("Wind Speed (m s$^{-1}$)", fontsize=12)
ax.set_ylabel("Height AGL (m)", fontsize=12)
ax.set_title(f"Time-Averaged Vertical Wind Speed Profile\n"
             f"YSU vs MYNN2.5 vs ERA5 | Jan 17–19 2024\n"
             f"{LAT0}°N, {abs(LON0)}°W", fontsize=11)
ax.set_ylim(0, 800)
ax.set_xlim(left=0)
ax.yaxis.set_major_locator(MultipleLocator(100))
ax.yaxis.set_minor_locator(MultipleLocator(50))
ax.legend(fontsize=11, framealpha=0.9)
ax.grid(True, alpha=0.3, linestyle=":")

plt.tight_layout()
plt.savefig(FIG_DIR / "plot2_vertical_profile.png", dpi=150, bbox_inches="tight")
plt.show()
print("  Saved: plot2_vertical_profile.png")


# %%
# =============================================================================
# CELL 11 — Extract spatial hub-height maps (vectorized with numpy)
# =============================================================================
print("\nExtracting spatial hub-height maps (vectorized)...")

def spatial_hub_map_fast(ds, mask, hub_ht=100.0):
    print("  Computing 3-D height field...")
    z3d = geopotential_height_agl_3d(ds)
    z3d = z3d[mask]

    print("  Destaggering U...")
    U = ds["U"].values
    U = 0.5 * (U[:, :, :, :-1] + U[:, :, :, 1:])
    U = U[mask]

    print("  Destaggering V...")
    V = ds["V"].values
    V = 0.5 * (V[:, :, :-1, :] + V[:, :, 1:, :])
    V = V[mask]

    wspd = np.hypot(U, V)

    print(f"  Interpolating to {hub_ht} m AGL across domain...")
    nT, nz, ny, nx = wspd.shape
    ws_hub = np.zeros((nT, ny, nx))

    for t in range(nT):
        for j in range(ny):
            for i in range(nx):
                ws_hub[t, j, i] = float(np.interp(
                    hub_ht, z3d[t, :, j, i], wspd[t, :, j, i]
                ))
        if (t+1) % 10 == 0:
            print(f"    timestep {t+1}/{nT}")

    return ws_hub.mean(axis=0)

lats2d = ds_ysu["XLAT"].isel(Time=0).values
lons2d = ds_ysu["XLONG"].isel(Time=0).values

print("YSU:")
ws_map_ysu  = spatial_hub_map_fast(ds_ysu,  mask_ysu,  HUB_HT)
print("MYNN:")
ws_map_mynn = spatial_hub_map_fast(ds_mynn, mask_mynn, HUB_HT)

ws_diff = ws_map_mynn - ws_map_ysu
print(f"\nDomain-mean diff (MYNN−YSU): {ws_diff.mean():.3f} m/s")


# %%
# =============================================================================
# CELL 12 — PLOT 3: Spatial maps — YSU, MYNN, difference (diverging colormap)
# =============================================================================
print("Plot 3: Spatial maps...")

fig, axes = plt.subplots(1, 3, figsize=(19, 6),
                          subplot_kw={"projection": proj})

vmin, vmax = 0, 10
max_abs = max(abs(np.percentile(ws_diff, 2)),
              abs(np.percentile(ws_diff, 98)))

for ax in axes:
    add_map_features(ax)

# Panel A — YSU
cf = axes[0].pcolormesh(lons2d, lats2d, ws_map_ysu,
                         transform=ccrs.PlateCarree(),
                         cmap="viridis", vmin=vmin, vmax=vmax, zorder=2)
plt.colorbar(cf, ax=axes[0], label="Wind Speed (m s$^{-1}$)", shrink=0.82, pad=0.02)
axes[0].plot(LON0, LAT0, "r*", markersize=12,
             transform=ccrs.PlateCarree(), zorder=5, label="Eval. point")
axes[0].set_title("(a) YSU\nTime-averaged 100-m wind speed", fontsize=10)
axes[0].legend(loc="lower right", fontsize=8)

# Panel B — MYNN
cf = axes[1].pcolormesh(lons2d, lats2d, ws_map_mynn,
                         transform=ccrs.PlateCarree(),
                         cmap="viridis", vmin=vmin, vmax=vmax, zorder=2)
plt.colorbar(cf, ax=axes[1], label="Wind Speed (m s$^{-1}$)", shrink=0.82, pad=0.02)
axes[1].plot(LON0, LAT0, "r*", markersize=12,
             transform=ccrs.PlateCarree(), zorder=5)
axes[1].set_title("(b) MYNN2.5\nTime-averaged 100-m wind speed", fontsize=10)

# Panel C — Difference with diverging RdBu_r colormap
norm_diff = TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs)
cf = axes[2].pcolormesh(lons2d, lats2d, ws_diff,
                         transform=ccrs.PlateCarree(),
                         cmap="RdBu_r", norm=norm_diff, zorder=2)
plt.colorbar(cf, ax=axes[2],
             label="Δ Wind Speed (m s$^{-1}$)\n(MYNN2.5 − YSU)",
             shrink=0.82, pad=0.02)
axes[2].contour(lons2d, lats2d, ws_diff, levels=[0],
                colors="black", linewidths=0.8,
                transform=ccrs.PlateCarree(), linestyles="--", zorder=4)
axes[2].plot(LON0, LAT0, "k*", markersize=12,
             transform=ccrs.PlateCarree(), zorder=5)
axes[2].set_title("(c) MYNN2.5 − YSU\nDifference in 100-m wind speed",
                  fontsize=10)

for ax in axes:
    ax.set_extent([lons2d.min()-0.2, lons2d.max()+0.2,
                   lats2d.min()-0.2, lats2d.max()+0.2],
                  crs=ccrs.PlateCarree())

fig.suptitle("Time-Averaged 100-m Hub-Height Wind Speed | "
             "Jan 17–19 2024 | Maryland d02 (4 km)",
             fontsize=12, fontweight="bold", y=1.01)
plt.tight_layout()
plt.savefig(FIG_DIR / "plot3_spatial_maps.png", dpi=150, bbox_inches="tight")
plt.show()
print("  Saved: plot3_spatial_maps.png")


# %%
# =============================================================================
# CELL 13 — PLOT 4: Time-height cross section (as in ICA07)
# =============================================================================
print("Plot 4: Time-height cross section...")

def time_height_wspd(z_agl, wspd, mask, z_max=600):
    z_vals = z_agl.values
    w_vals = wspd.values
    z_mean = z_vals[mask].mean(axis=0)
    zmask  = z_mean <= z_max
    return z_mean[zmask], w_vals[mask][:, zmask]

z_ysu_th,  ws_th_ysu  = time_height_wspd(z_agl_ysu,  wspd_ysu,  mask_ysu)
z_mynn_th, ws_th_mynn = time_height_wspd(z_agl_mynn, wspd_mynn, mask_mynn)

levels_th = np.linspace(0, 12, 25)
fig, axes = plt.subplots(2, 1, figsize=(14, 9))

for ax, times, z_th, ws_th, lbl, col in [
    (axes[0], ysu_times_trim,  z_ysu_th,  ws_th_ysu,  "(a) YSU",     COL_YSU),
    (axes[1], mynn_times_trim, z_mynn_th, ws_th_mynn, "(b) MYNN2.5", COL_MYNN),
]:
    cf = ax.contourf(times, z_th, ws_th.T,
                     levels=levels_th, cmap="plasma", extend="max")
    ax.axhline(HUB_HT, color="white", linestyle="--", linewidth=1.8,
               label=f"{int(HUB_HT)}-m hub height")
    plt.colorbar(cf, ax=ax, label="Wind Speed (m s$^{-1}$)",
                 pad=0.01, shrink=0.95)
    ax.set_ylabel("Height AGL (m)", fontsize=11)
    ax.set_title(f"{lbl} — Wind Speed (m s$^{{-1}}$)", fontsize=11, color=col)
    ax.set_ylim(0, 600)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%H:00Z"))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=12))
    ax.legend(loc="upper right", fontsize=9, framealpha=0.85)
    ax.grid(True, alpha=0.15, color="white", linestyle=":")

fig.suptitle(f"Time-Height Wind Speed | Jan 17–19 2024 | "
             f"{LAT0}°N, {abs(LON0)}°W | d02 4-km domain",
             fontsize=12, fontweight="bold")
plt.tight_layout()
plt.savefig(FIG_DIR / "plot4_timeheight.png", dpi=150, bbox_inches="tight")
plt.show()
print("  Saved: plot4_timeheight.png")


# %%
# =============================================================================
# CELL 13b — Interpolate ERA5 to WRF timestamps for bias calculation
# ERA5 is hourly; WRF is hourly — but timestamps may be offset slightly.
# We interpolate ERA5 onto WRF time grids using numeric timestamps.
# =============================================================================

# Convert to numeric (nanoseconds since epoch) for interpolation
era5_t_ns  = era5_times_trim.astype(np.int64)
ysu_t_ns   = pd.to_datetime(ysu_times_trim).astype(np.int64)
mynn_t_ns  = pd.to_datetime(mynn_times_trim).astype(np.int64)

# Linear interpolation — only interpolate within ERA5 range
era5_interp_fn = interp1d(era5_t_ns, era5_ws_trim,
                           kind="linear", bounds_error=False,
                           fill_value=np.nan)

era5_at_ysu  = era5_interp_fn(ysu_t_ns)    # ERA5 sampled at YSU timestamps
era5_at_mynn = era5_interp_fn(mynn_t_ns)   # ERA5 sampled at MYNN timestamps

# Bias = WRF - ERA5 (positive = WRF overestimates)
bias_ysu  = ysu_hub_ts  - era5_at_ysu
bias_mynn = mynn_hub_ts - era5_at_mynn

# Remove NaNs for stats
valid_y = ~np.isnan(era5_at_ysu)
valid_m = ~np.isnan(era5_at_mynn)

mae_ysu  = float(np.abs(bias_ysu[valid_y]).mean())
mae_mynn = float(np.abs(bias_mynn[valid_m]).mean())
mbe_ysu  = float(bias_ysu[valid_y].mean())
mbe_mynn = float(bias_mynn[valid_m].mean())
rmse_ysu  = float(np.sqrt((bias_ysu[valid_y]**2).mean()))
rmse_mynn = float(np.sqrt((bias_mynn[valid_m]**2).mean()))

winner = "MYNN2.5" if mae_mynn < mae_ysu else "YSU"
print(f"\nComparison vs ERA5 100m:")
print(f"  YSU  MAE={mae_ysu:.2f}  MBE={mbe_ysu:+.2f}  RMSE={rmse_ysu:.2f} m/s")
print(f"  MYNN MAE={mae_mynn:.2f}  MBE={mbe_mynn:+.2f}  RMSE={rmse_mynn:.2f} m/s")
print(f"  Closer to ERA5: {winner}")


# %%
# =============================================================================
# CELL 13c — PLOT 5: Bias vs ERA5 reference
# Two panels:
#   Top: WRF bias time series (YSU - ERA5 and MYNN - ERA5)
#   Bottom: Error metric summary bar chart
# =============================================================================
print("Plot 5: Bias vs ERA5 reference...")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 8),
                                gridspec_kw={"height_ratios": [2.2, 1]})

# ── Top: bias time series ────────────────────────────────────────────────
ax1.axhline(0, color="black", linewidth=1.0, alpha=0.5, linestyle="-")
ax1.axhline(+1, color="gray", linewidth=0.7, linestyle=":", alpha=0.5)
ax1.axhline(-1, color="gray", linewidth=0.7, linestyle=":", alpha=0.5)

ax1.plot(pd.to_datetime(ysu_times_trim),  bias_ysu,
         color=COL_YSU,  linewidth=1.8, label="YSU − ERA5", zorder=4)
ax1.plot(pd.to_datetime(mynn_times_trim), bias_mynn,
         color=COL_MYNN, linewidth=1.8, label="MYNN2.5 − ERA5", zorder=4)

# Shade nighttime
for d in [17, 18, 19]:
    ns = pd.Timestamp(f"2024-01-{d-1:02d} 23:00") if d > 17 else T_START
    ne = pd.Timestamp(f"2024-01-{d:02d} 11:00")
    if ns >= T_START and ne <= T_END:
        ax1.axvspan(ns, ne, alpha=0.06, color="navy",
                    label="Nighttime" if d == 18 else "")

ax1.set_ylabel("Bias vs ERA5 (m s$^{-1}$)", fontsize=12)
ax1.set_title(f"WRF Hub-Height Wind Speed Bias vs ERA5 Reference\n"
              f"Jan 17–19 2024 | {LAT0}°N, {abs(LON0)}°W | positive = WRF overestimates",
              fontsize=11)
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%b %d\n%H:00Z"))
ax1.xaxis.set_major_locator(mdates.HourLocator(interval=12))
ax1.set_xlim(T_START, T_END)
ax1.legend(fontsize=10, framealpha=0.9, loc="upper right")
ax1.grid(True, alpha=0.25, linestyle=":")

# Annotate mean bias values on the plot
ax1.text(0.02, 0.92, f"YSU  mean bias = {mbe_ysu:+.2f} m s$^{{-1}}$",
         transform=ax1.transAxes, fontsize=10, color=COL_YSU, fontweight="bold")
ax1.text(0.02, 0.82, f"MYNN mean bias = {mbe_mynn:+.2f} m s$^{{-1}}$",
         transform=ax1.transAxes, fontsize=10, color=COL_MYNN, fontweight="bold")

# ── Bottom: error metric bar chart ───────────────────────────────────────
metrics   = ["MAE\n(m s⁻¹)", "MBE\n(m s⁻¹)", "RMSE\n(m s⁻¹)"]
vals_ysu  = [mae_ysu,  mbe_ysu,  rmse_ysu]
vals_mynn = [mae_mynn, mbe_mynn, rmse_mynn]

x     = np.arange(len(metrics))
width = 0.32

bars_y = ax2.bar(x - width/2, vals_ysu,  width, label="YSU",
                  color=COL_YSU,  alpha=0.88, edgecolor="white", linewidth=0.5)
bars_m = ax2.bar(x + width/2, vals_mynn, width, label="MYNN2.5",
                  color=COL_MYNN, alpha=0.88, edgecolor="white", linewidth=0.5)

# Label each bar
for bar in bars_y:
    h = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2, h + 0.02,
             f"{h:+.2f}", ha="center", va="bottom", fontsize=9,
             color=COL_YSU, fontweight="bold")
for bar in bars_m:
    h = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2, h + 0.02,
             f"{h:+.2f}", ha="center", va="bottom", fontsize=9,
             color=COL_MYNN, fontweight="bold")

# Winner annotation — upper center, above MBE (avoids covering RMSE bar labels on the right)
better_col  = COL_MYNN if winner == "MYNN2.5" else COL_YSU
ax2.text(0.5, 0.97,
         f"Closer to ERA5: {winner}  (MAE Δ = {abs(mae_ysu - mae_mynn):.2f} m s$^{{-1}}$)",
         transform=ax2.transAxes, ha="center", va="top", fontsize=11,
         color=better_col, fontweight="bold",
         bbox=dict(boxstyle="round,pad=0.35", facecolor="white",
                   edgecolor=better_col, alpha=0.9))

ax2.axhline(0, color="black", linewidth=0.8, alpha=0.5)
ax2.set_ylim(top=max(vals_ysu + vals_mynn) * 1.35)
ax2.set_xticks(x)
ax2.set_xticklabels(metrics, fontsize=11)
ax2.set_ylabel("Error (m s$^{-1}$)", fontsize=11)
ax2.set_title("Error Metrics vs ERA5 Reference", fontsize=11)
ax2.legend(fontsize=10, loc="upper left", bbox_to_anchor=(1.0, 1.0))
ax2.grid(True, axis="y", alpha=0.25, linestyle=":")

plt.tight_layout()
plt.savefig(FIG_DIR / "plot5_bias_vs_era5.png", dpi=150, bbox_inches="tight")
plt.show()
print("  Saved: plot5_bias_vs_era5.png")


# %%
# =============================================================================
# CELL 14 — Summary statistics (conclusions slide)
# =============================================================================
min_len   = min(len(ysu_hub_ts), len(mynn_hub_ts))
diff_ts   = ysu_hub_ts[:min_len] - mynn_hub_ts[:min_len]
mad       = float(np.abs(diff_ts).mean())
mean_bias = float(diff_ts.mean())

print("\n" + "=" * 60)
print("SUMMARY STATISTICS")
print("=" * 60)
print(f"\n--- Inter-scheme comparison ---")
print(f"YSU  mean 100m wind : {ysu_hub_ts.mean():.2f} m/s")
print(f"MYNN mean 100m wind : {mynn_hub_ts.mean():.2f} m/s")
print(f"ERA5 mean 100m wind : {era5_ws_trim.mean():.2f} m/s")
print(f"Mean bias (YSU−MYNN): {mean_bias:+.2f} m/s")
print(f"MAD (|YSU−MYNN|)    : {mad:.2f} m/s")
print(f"Max |diff|          : {np.abs(diff_ts).max():.2f} m/s")

print(f"\n--- Comparison vs ERA5 100m ---")
print(f"YSU  MAE  : {mae_ysu:.2f} m/s")
print(f"MYNN MAE  : {mae_mynn:.2f} m/s")
print(f"YSU  MBE  : {mbe_ysu:+.2f} m/s")
print(f"MYNN MBE  : {mbe_mynn:+.2f} m/s")
print(f"YSU  RMSE : {rmse_ysu:.2f} m/s")
print(f"MYNN RMSE : {rmse_mynn:.2f} m/s")
print(f"Closer to ERA5: {winner}")

print(f"\n--- Hours below the 4 m/s low-wind threshold ({len(era5_ws_trim)} h window) ---")
print(f"ERA5 : {int((era5_ws_trim < 4).sum())}")
print(f"YSU  : {int((ysu_hub_ts < 4).sum())}")
print(f"MYNN : {int((mynn_hub_ts < 4).sum())}")

print(f"\n--- Profile at hub height ---")
print(f"YSU  at {int(HUB_HT)}m : {np.interp(HUB_HT, z_common, ws_prof_ysu):.2f} m/s")
print(f"MYNN at {int(HUB_HT)}m : {np.interp(HUB_HT, z_common, ws_prof_mynn):.2f} m/s")
print(f"ERA5 at {int(HUB_HT)}m : {era5_mean_100m:.2f} m/s")

print(f"\n--- Spatial domain ---")
print(f"Domain-mean diff (MYNN−YSU): {ws_diff.mean():.3f} m/s")

print(f"\n--- Hypothesis assessment ---")
hyp = "SUPPORTED" if mad >= 1.0 else "NOT MET (MAD < 1 m/s)"
print(f"Inter-scheme MAD ≥ 1 m/s  : {hyp}  [MAD={mad:.2f}]")
print(f"Comparison with ERA5       : {winner} is closer")
print("=" * 60)
print("\nAll plots saved:")
for p in ["plot0_domain_map.png", "plot1_timeseries_100m_wind.png",
          "plot2_vertical_profile.png", "plot3_spatial_maps.png",
          "plot4_timeheight.png", "plot5_bias_vs_era5.png"]:
    print(f"  {p}")