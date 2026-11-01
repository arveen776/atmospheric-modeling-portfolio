#!/usr/bin/env bash
# Optional: extract hourly domain-total POWER (W) from WFP wrfouts using NCO (ncwa).
# Usage:  bash extract_power_nco.sh /path/to/wrf_files 2023-08-24 > power_totals.txt
set -euo pipefail
WFP_DIR=${1:?wrf_files directory}
DATE=${2:-2023-08-24}
TMP=$(mktemp)
trap 'rm -f "$TMP"' EXIT
echo "# hour_utc total_power_W"
for f in "${WFP_DIR}/wrfout_d01_${DATE}"_*; do
  [[ -e "$f" ]] || continue
  hr=$(basename "$f" | sed -n "s/.*_${DATE}_\\([0-9][0-9]\\)_00_00/\\1/p")
  ncwa -O -y sum -v POWER -a south_north,west_east "$f" "$TMP"
  tot=$(ncks -H -v POWER "$TMP" | awk '{print $NF}')
  echo "$hr $tot"
done | sort -n
