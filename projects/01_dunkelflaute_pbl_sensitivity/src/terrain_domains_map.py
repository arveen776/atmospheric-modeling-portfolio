#!/usr/bin/env python3
"""
Terrain elevation and WRF nest boundaries from geo_em NetCDF files.

Looks for geo_em.d01.nc, geo_em.d02.nc, ... in the data directory (same layout
as final.py). Plots terrain (HGT_M) on the outer domain and draws rectangular
outlines for every nest found.

Dependencies: numpy, xarray, netCDF4, matplotlib, cartopy
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.colors import LinearSegmentedColormap


def find_geo_em_files(directory: Path) -> list[Path]:
    pat = re.compile(r"geo_em\.d(\d+)\.nc$", re.IGNORECASE)
    found: list[tuple[int, Path]] = []
    for p in directory.glob("geo_em.d*.nc"):
        m = pat.match(p.name)
        if m:
            found.append((int(m.group(1)), p))
    found.sort(key=lambda t: t[0])
    return [p for _, p in found]


def domain_box_corners(lat2d: np.ndarray, lon2d: np.ndarray) -> tuple[list, list]:
    lx = [
        lon2d[0, 0],
        lon2d[0, -1],
        lon2d[-1, -1],
        lon2d[-1, 0],
        lon2d[0, 0],
    ]
    ly = [
        lat2d[0, 0],
        lat2d[0, -1],
        lat2d[-1, -1],
        lat2d[-1, 0],
        lat2d[0, 0],
    ]
    return lx, ly


def projection_from_geo_em(ds: xr.Dataset) -> ccrs.Projection:
    """Build a cartopy CRS from WRF geo_em global attributes."""
    map_proj = int(ds.attrs.get("MAP_PROJ", 1))
    stand_lon = float(ds.attrs.get("STAND_LON", -97.0))
    truelat1 = float(ds.attrs.get("TRUELAT1", 30.0))
    truelat2 = float(ds.attrs.get("TRUELAT2", truelat1))

    if map_proj == 1:
        # Lambert conformal (WRF MAP_PROJ = 1)
        return ccrs.LambertConformal(
            central_longitude=stand_lon,
            central_latitude=truelat1,
            standard_parallels=(truelat1, truelat2),
        )
    if map_proj == 2:
        # Polar stereographic
        return ccrs.NorthPolarStereo(central_longitude=stand_lon)
    if map_proj == 3:
        return ccrs.Mercator(central_longitude=stand_lon)
    if map_proj in (0, 6):
        return ccrs.PlateCarree(central_longitude=stand_lon)
    # Unknown: assume Lambert (common for continental US nests)
    return ccrs.LambertConformal(
        central_longitude=stand_lon,
        central_latitude=truelat1,
        standard_parallels=(truelat1, truelat2),
    )


# Basemap colors: open water reads as blue; dry land backdrop is neutral so
# elevation colors read clearly on top.
OCEAN_BLUE = "#4a90c8"
LAKE_BLUE = "#5c9fd4"
LAND_BACKDROP = "#e6e8e0"


def add_ocean_and_land(ax) -> None:
    ax.add_feature(
        cfeature.OCEAN.with_scale("10m"),
        facecolor=OCEAN_BLUE,
        zorder=1,
    )
    ax.add_feature(
        cfeature.LAND.with_scale("10m"),
        facecolor=LAND_BACKDROP,
        zorder=1,
    )


def add_lakes(ax) -> None:
    """Draw above terrain so inland water matches ocean-like blue."""
    ax.add_feature(
        cfeature.LAKES.with_scale("10m"),
        facecolor=LAKE_BLUE,
        edgecolor="#3a7aaa",
        linewidth=0.3,
        zorder=2.5,
    )


def add_states_and_coast(ax) -> None:
    ax.add_feature(
        cfeature.STATES.with_scale("10m"),
        linewidth=0.8,
        edgecolor="black",
        zorder=3,
    )
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), linewidth=1.0, zorder=3)


# Land-only elevation ramp (low → high): greens and tans, no ocean blues.
LAND_ELEV_CMAP = LinearSegmentedColormap.from_list(
    "land_elev",
    [
        "#c8e6c9",
        "#81c784",
        "#558b2f",
        "#c9a961",
        "#8d6e63",
        "#bcaaa4",
        "#eceff1",
    ],
    N=256,
)


def land_elevation_for_plot(
    hgt: np.ndarray, geo_ds: xr.Dataset
) -> np.ndarray:
    """Return HGT with water cells as NaN so ocean/lake basemap shows through."""
    h = np.asarray(hgt, dtype=float).copy()
    if "LANDMASK" in geo_ds:
        lm = np.asarray(geo_ds["LANDMASK"].values[0])
        h = np.where(lm > 0.5, h, np.nan)
    else:
        h = np.where(h > 0.5, h, np.nan)
    return h


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Terrain map + WRF domain boxes from geo_em files."
    )
    ap.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory containing geo_em.d*.nc (default: script directory)",
    )
    ap.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("terrain_domains_map.png"),
        help="Output PNG path",
    )
    ap.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Figure DPI",
    )
    ap.add_argument(
        "--vmin",
        type=float,
        default=None,
        help="Terrain color scale minimum (m); default: data min",
    )
    ap.add_argument(
        "--vmax",
        type=float,
        default=None,
        help="Terrain color scale maximum (m); default: data max",
    )
    ap.add_argument(
        "--marker",
        nargs=2,
        type=float,
        metavar=("LAT", "LON"),
        default=None,
        help="Optional lat lon for a star marker (e.g. --marker 38.8 -76.0)",
    )
    args = ap.parse_args()

    data_dir = args.data_dir.resolve()
    geo_files = find_geo_em_files(data_dir)
    if not geo_files:
        raise SystemExit(
            f"No geo_em.d*.nc files found in {data_dir}\n"
            "Place geo_em.d01.nc, geo_em.d02.nc, ... next to this script or pass --data-dir."
        )

    # Outer domain = lowest dXX number (first file)
    outer_path = geo_files[0]
    geo_outer = xr.open_dataset(outer_path, engine="netcdf4")
    proj = projection_from_geo_em(geo_outer)

    lat = geo_outer["XLAT_M"].values[0]
    lon = geo_outer["XLONG_M"].values[0]
    hgt = geo_outer["HGT_M"].values[0]

    hgt_land = land_elevation_for_plot(hgt, geo_outer)
    land_vals = hgt_land[np.isfinite(hgt_land)]
    if land_vals.size == 0:
        raise SystemExit("No land cells in HGT/LANDMASK — check geo_em file.")
    vmin = args.vmin if args.vmin is not None else float(np.nanmin(land_vals))
    vmax = args.vmax if args.vmax is not None else float(np.nanmax(land_vals))

    # Distinct colors for nest outlines
    box_colors = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e", "#17becf"]

    fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={"projection": proj})
    add_ocean_and_land(ax)
    h_plot = np.ma.masked_invalid(hgt_land)
    cf = ax.pcolormesh(
        lon,
        lat,
        h_plot,
        transform=ccrs.PlateCarree(),
        cmap=LAND_ELEV_CMAP,
        vmin=vmin,
        vmax=vmax,
        shading="auto",
        zorder=2,
    )
    plt.colorbar(cf, ax=ax, label="Terrain elevation (m MSL)", shrink=0.6, pad=0.02)
    add_lakes(ax)
    add_states_and_coast(ax)

    for i, path in enumerate(geo_files):
        ds = xr.open_dataset(path, engine="netcdf4") if path != outer_path else geo_outer
        lat_i = ds["XLAT_M"].values[0]
        lon_i = ds["XLONG_M"].values[0]
        dom_id = re.search(r"d(\d+)", path.stem, re.I)
        label = f"d{dom_id.group(1)}" if dom_id else path.stem
        color = box_colors[i % len(box_colors)]
        lx, ly = domain_box_corners(lat_i, lon_i)
        ax.plot(
            lx,
            ly,
            color=color,
            linewidth=2.5,
            transform=ccrs.PlateCarree(),
            label=f"{label} ({lat_i.shape[0]}×{lat_i.shape[1]})",
            zorder=5,
        )
        if path != outer_path:
            ds.close()

    if args.marker is not None:
        mlat, mlon = args.marker
        ax.plot(
            mlon,
            mlat,
            "k*",
            markersize=16,
            transform=ccrs.PlateCarree(),
            zorder=7,
            label=f"Point ({mlat}°N, {abs(mlon)}°W)" if mlon < 0 else f"Point ({mlat}, {mlon})",
        )

    pad = 0.5
    ax.set_extent(
        [
            lon[0, 0] - pad,
            lon[0, -1] + pad,
            lat[0, 0] - pad,
            lat[-1, 0] + pad,
        ],
        crs=ccrs.PlateCarree(),
    )
    ax.legend(loc="lower right", fontsize=10, framealpha=0.9)
    ax.set_title(
        "Terrain and WRF domains (geo_em)\n"
        f"Outer grid: {outer_path.name}",
        fontsize=12,
        fontweight="bold",
    )
    gl = ax.gridlines(
        draw_labels=True, linewidth=0.4, color="gray", alpha=0.5, linestyle="--"
    )
    gl.top_labels = False
    gl.right_labels = False

    plt.tight_layout()
    out = args.output.resolve()
    plt.savefig(out, dpi=args.dpi, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out}")
    geo_outer.close()


if __name__ == "__main__":
    main()
