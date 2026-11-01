#!/usr/bin/env python3
"""
Plot manufacturer-style power curves (kW vs m/s) from WRF Fitch ``wind-turbine-*.tbl`` files.

The AWAKEN case uses four GE open-source-style curves (types 1–4 in ``windturbines.txt``).
If ``wind_turbine_tables/`` is missing the tables, this script writes GE 1.72 / 1.79 / 2.3 / 2.82 MW
tables, then saves ``figures/turbine_power_curves.png``.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# -----------------------------------------------------------------------------
# WRF .tbl format (see WRF doc/README.windturbine)
# -----------------------------------------------------------------------------


@dataclass
class TurbineTbl:
    n: int
    hub_m: float
    rotor_d_m: float
    stand_ct: float
    nominal_mw: float
    wind_speed: np.ndarray
    thrust_ct: np.ndarray
    power_kw: np.ndarray


def parse_wrf_turbine_tbl(path: Path) -> TurbineTbl:
    lines = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    n = int(lines[0])
    hub, diam, ct0, pnom = map(float, lines[1].split())
    ws, ct, pkw = [], [], []
    for ln in lines[2 : 2 + n]:
        a, b, c = ln.split()
        ws.append(float(a))
        ct.append(float(b))
        pkw.append(float(c))
    return TurbineTbl(
        n=n,
        hub_m=hub,
        rotor_d_m=diam,
        stand_ct=ct0,
        nominal_mw=pnom,
        wind_speed=np.array(ws),
        thrust_ct=np.array(ct),
        power_kw=np.array(pkw),
    )


def write_wrf_turbine_tbl(path: Path, tbl: TurbineTbl) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(tbl.n), f"{tbl.hub_m} {tbl.rotor_d_m} {tbl.stand_ct} {tbl.nominal_mw}"]
    for i in range(tbl.n):
        lines.append(f"{tbl.wind_speed[i]} {tbl.thrust_ct[i]} {tbl.power_kw[i]}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def build_ge_awaken_tables() -> dict[int, TurbineTbl]:
    """
    Four curves (GE open-source nominal powers).
    Thrust vs wind follows the shape of WRF's example ``wind-turbine-1.tbl`` for ws 4–25 m/s;
    power columns match rated speeds from the figure (10 m/s except 2.82 MW at ~11 m/s).
    """
    # Thrust coefficient template (WRF run/wind-turbine-1.tbl), ws = 4..25
    ct_template = np.array(
        [
            0.805,
            0.805,
            0.805,
            0.805,
            0.805,
            0.805,
            0.790,
            0.740,
            0.700,
            0.400,
            0.300,
            0.250,
            0.200,
            0.160,
            0.140,
            0.120,
            0.100,
            0.080,
            0.070,
            0.060,
            0.055,
            0.050,
        ]
    )
    ws_main = np.arange(4.0, 26.0, 1.0)
    assert len(ws_main) == len(ct_template)

    # Include a cut-in row at 3 m/s (zero power)
    ws = np.concatenate([[3.0], ws_main])
    ct = np.concatenate([[0.820], ct_template])

    def pack(p_main: list[float]) -> TurbineTbl:
        assert len(p_main) == len(ws_main)
        pkw = np.concatenate([[0.0], np.array(p_main)])
        return TurbineTbl(
            n=len(ws),
            hub_m=0.0,
            rotor_d_m=0.0,
            stand_ct=0.0,
            nominal_mw=0.0,
            wind_speed=ws,
            thrust_ct=ct,
            power_kw=pkw,
        )

    # Power (kW) at ws 4..25 — piecewise to match slide plateaus
    p172 = [95, 240, 430, 690, 1000, 1360] + [1720] * 16
    p179 = [100, 250, 450, 720, 1050, 1400] + [1790] * 16
    p230 = [150, 350, 600, 950, 1350, 1850] + [2300] * 16
    p282 = [200, 400, 700, 1100, 1600, 2200, 2650] + [2820] * 15

    specs = {
        1: dict(nom_mw=1.72, hub=80.0, d=87.0, stand=0.130, power=p172),
        2: dict(nom_mw=1.79, hub=90.0, d=100.0, stand=0.130, power=p179),
        3: dict(nom_mw=2.30, hub=100.0, d=116.0, stand=0.130, power=p230),
        4: dict(nom_mw=2.82, hub=110.0, d=127.0, stand=0.130, power=p282),
    }
    out: dict[int, TurbineTbl] = {}
    for tid, s in specs.items():
        t = pack(s["power"])
        t.hub_m = s["hub"]
        t.rotor_d_m = s["d"]
        t.stand_ct = s["stand"]
        t.nominal_mw = s["nom_mw"]
        out[tid] = t
    return out


def main() -> None:
    root = project_root()
    ap = argparse.ArgumentParser(description="Plot GE power curves from WRF wind-turbine-*.tbl")
    ap.add_argument(
        "--tbl-dir",
        type=str,
        default=str(root / "inputs" / "wind_turbine_tables"),
        help="Directory with wind-turbine-1.tbl … wind-turbine-4.tbl",
    )
    ap.add_argument(
        "--write-defaults",
        action="store_true",
        help="Write GE AWAKEN tables into --tbl-dir if missing",
    )
    ap.add_argument(
        "--output",
        type=str,
        default=str(root / "figures" / "turbine_power_curves.png"),
        help="Output figure path",
    )
    args = ap.parse_args()
    tbl_dir = Path(args.tbl_dir).expanduser()
    if not tbl_dir.is_absolute():
        tbl_dir = root / tbl_dir

    need_write = args.write_defaults
    for i in range(1, 5):
        if not (tbl_dir / f"wind-turbine-{i}.tbl").is_file():
            need_write = True
            break
    if need_write:
        for tid, tbl in build_ge_awaken_tables().items():
            write_wrf_turbine_tbl(tbl_dir / f"wind-turbine-{tid}.tbl", tbl)
        print(f"Wrote wind-turbine-{{1..4}}.tbl under {tbl_dir}")

    meta = {
        1: ("GE 1.72 MW open-source", "red"),
        2: ("GE 1.79 MW open-source", "green"),
        3: ("GE 2.3 MW open-source", "tab:blue"),
        4: ("GE 2.82 MW open-source", "darkorange"),
    }

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for tid in range(1, 5):
        path = tbl_dir / f"wind-turbine-{tid}.tbl"
        tbl = parse_wrf_turbine_tbl(path)
        label, color = meta[tid]
        ax.plot(
            tbl.wind_speed,
            tbl.power_kw,
            "-",
            lw=2.2,
            color=color,
            label=f"{label} ({tbl.nominal_mw:.2f} MW nominal)",
        )

    ax.set_xlabel(r"Wind speed (m s$^{-1}$)")
    ax.set_ylabel("Power (kW)")
    ax.set_title("Wind turbine power curves used for WRF wind-farm types 1–4 (AWAKEN)")
    ax.set_xlim(0, 25)
    ax.set_ylim(0, 3000)
    ax.grid(True, alpha=0.35)
    ax.legend(loc="lower right", framealpha=0.95, fontsize=9)

    out = Path(args.output).expanduser()
    if not out.is_absolute():
        out = root / out
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()
