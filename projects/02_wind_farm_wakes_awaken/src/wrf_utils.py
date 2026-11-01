"""
Shared helpers for the WRF AWAKEN wind-farm plots.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

G = 9.80665


def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def list_wrfout_d01(data_dir: Path, date_glob: str = "*") -> list[Path]:
    return sorted(data_dir.glob(f"wrfout_d01_{date_glob}"))


def open_wrf(path: Path) -> xr.Dataset:
    return xr.open_dataset(path, decode_times=False)


def mass_level_heights_agl_m(ds: xr.Dataset, time_idx: int = 0) -> np.ndarray:
    """Height above ground (m) at WRF mass levels, shape (nz, ny, nx)."""
    ph = (ds["PH"] + ds["PHB"]).isel(Time=time_idx).values / G
    hgt = ds["HGT"].isel(Time=time_idx).values
    z_mass = 0.5 * (ph[:-1, ...] + ph[1:, ...])
    return z_mass - hgt[np.newaxis, ...]


def hub_mass_level_index(ds: xr.Dataset, target_agl_m: float = 113.0, time_idx: int = 0) -> int:
    """Pick mass level whose AGL height at domain center is closest to target_agl_m (hub ~100–126 m)."""
    agl = mass_level_heights_agl_m(ds, time_idx=time_idx)
    ny, nx = agl.shape[1], agl.shape[2]
    j, i = ny // 2, nx // 2
    z_col = agl[:, j, i]
    k = int(np.argmin(np.abs(z_col - target_agl_m)))
    return k


def destagger_uv_at_level(ds: xr.Dataset, k: int, time_idx: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Destagger U,V to mass points; shapes (ny, nx). u eastward, v northward (m/s)."""
    u_stag = ds["U"].isel(Time=time_idx, bottom_top=k).values
    v_stag = ds["V"].isel(Time=time_idx, bottom_top=k).values
    u = 0.5 * (u_stag[:, :-1] + u_stag[:, 1:])
    v = 0.5 * (v_stag[:-1, :] + v_stag[1:, :])
    return u, v


def destagger_wind_speed_at_level(ds: xr.Dataset, k: int, time_idx: int = 0) -> np.ndarray:
    """U,V destaggered to mass points; return wind speed (m/s) shape (ny, nx)."""
    u, v = destagger_uv_at_level(ds, k, time_idx=time_idx)
    return np.sqrt(u**2 + v**2)


def destagger_wind_speed_at_agl(ds: xr.Dataset, target_agl_m: float, time_idx: int = 0) -> np.ndarray:
    """
    Destagger U,V on each mass level, then linearly interpolate to ``target_agl_m`` AGL at each column.

    Uses that column's own (PH+PHB)/g and HGT so hub height follows local terrain. This avoids
    comparing WFP and NWF on a fixed eta surface that sits at different physical heights over
    ridges/valleys (a common source of terrain-parallel artifacts in difference maps).
    """
    agl = mass_level_heights_agl_m(ds, time_idx=time_idx)
    nz, ny, nx = agl.shape
    u_stack = np.empty((nz, ny, nx), dtype=np.float64)
    v_stack = np.empty((nz, ny, nx), dtype=np.float64)
    for k in range(nz):
        u, v = destagger_uv_at_level(ds, k, time_idx=time_idx)
        u_stack[k] = u
        v_stack[k] = v
    ui = np.empty((ny, nx), dtype=np.float64)
    vi = np.empty((ny, nx), dtype=np.float64)
    for j in range(ny):
        for i in range(nx):
            z = agl[:, j, i]
            ui[j, i] = np.interp(
                target_agl_m,
                z,
                u_stack[:, j, i],
                left=u_stack[0, j, i],
                right=u_stack[-1, j, i],
            )
            vi[j, i] = np.interp(
                target_agl_m,
                z,
                v_stack[:, j, i],
                left=v_stack[0, j, i],
                right=v_stack[-1, j, i],
            )
    return np.sqrt(ui**2 + vi**2)


def read_power_2d(path: Path) -> np.ndarray:
    """First time slice of POWER only (netCDF4 avoids xarray metadata overhead)."""
    import netCDF4 as nc

    with nc.Dataset(str(path), "r") as d:
        return np.array(d.variables["POWER"][0, ...], dtype=np.float64)


def turbine_mask_from_power(wfp_dir: Path, date_glob: str = "*") -> np.ndarray:
    """True where any matching hourly file has POWER > 0 (turbine-hosting cells)."""
    files = list_wrfout_d01(wfp_dir, date_glob=date_glob)
    mx = None
    for f in files:
        p = read_power_2d(f)
        mx = p if mx is None else np.maximum(mx, p)
    return mx > 0


def parse_windturbines_txt(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """
    Parse windturbines.txt in either format:
    - AWAKEN: ``lat lon turbine_type`` (first two columns are decimal degrees)
    - Alternate WRF layout: ``index lat lon ...`` (small integer index in column 0)
    """
    lats, lons = [], []
    # AWAKEN benchmark nest (d01) — lats ~35–37°N, lons ~97–98°W
    def _in_awaken(lat: float, lon: float) -> bool:
        return 34.0 <= lat <= 38.5 and -100.0 <= lon <= -95.5

    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.replace(",", " ").split()
        if len(parts) < 2:
            continue
        try:
            a = float(parts[0])
            b = float(parts[1])
        except ValueError:
            continue
        if len(parts) >= 3 and _in_awaken(a, b):
            lat, lon = a, b
        elif len(parts) >= 3:
            try:
                lat = float(parts[1])
                lon = float(parts[2])
            except ValueError:
                continue
        else:
            lat, lon = a, b
        if not (-90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0):
            continue
        lats.append(lat)
        lons.append(lon)
    return np.array(lats), np.array(lons)


def wrfout_path_for_utc(data_dir: Path, date: str, hour: int) -> Path:
    """date='2023-08-24', hour=6 -> wrfout_d01_2023-08-24_06_00_00"""
    name = f"wrfout_d01_{date}_{hour:02d}_00_00"
    p = data_dir / name
    if not p.exists():
        raise FileNotFoundError(f"Missing WRF file: {p}")
    return p
