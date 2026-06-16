#!/usr/bin/env bash
#
# Mirrors GFS-Wave hourly data from the NOAA S3 bucket to a local directory,
# iterating over every valid calendar date for the requested year/month
# combinations. Invalid dates (e.g. February 30) are no-ops, not errors.
#
# Usage:
#   bash scripts/sync_s3.sh --years 2024 2025 --months 1 2 3 --data-dir ./data
#
# Requires AWS credentials to already be configured in the environment
# (e.g. via `aws configure` or AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY).

set -euo pipefail

BUCKET_ROOT="s3://noaa-weather-service/process"
HOURS=("00" "03" "06" "09" "12" "15" "18" "21")

years=()
months=()
data_dir=""

usage() {
  cat >&2 <<'EOF'
usage: sync_s3.sh --years <y1> [y2 ...] --months <m1> [m2 ...] --data-dir <path>
EOF
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --years)
      shift
      while [[ $# -gt 0 && ! "$1" == --* ]]; do years+=("$1"); shift; done
      ;;
    --months)
      shift
      while [[ $# -gt 0 && ! "$1" == --* ]]; do months+=("$1"); shift; done
      ;;
    --data-dir)
      shift
      [[ $# -gt 0 ]] || usage
      data_dir="$1"
      shift
      ;;
    *)
      usage
      ;;
  esac
done

[[ ${#years[@]} -gt 0 && ${#months[@]} -gt 0 && -n "$data_dir" ]] || usage

days_in_month() {
  local year="$1" month="$2"
  case "$((10#$month))" in
    1|3|5|7|8|10|12) echo 31 ;;
    4|6|9|11) echo 30 ;;
    2)
      if (( (10#$year % 4 == 0 && 10#$year % 100 != 0) || 10#$year % 400 == 0 )); then
        echo 29
      else
        echo 28
      fi
      ;;
    *) echo 0 ;;
  esac
}

for year in "${years[@]}"; do
  for month in "${months[@]}"; do
    month_padded=$(printf "%02d" "$((10#$month))")
    echo "=== Starting ${year}-${month_padded} ==="
    nd=$(days_in_month "$year" "$month_padded")
    for ((day = 1; day <= nd; day++)); do
      day_padded=$(printf "%02d" "$day")
      for hour in "${HOURS[@]}"; do
        aws s3 sync \
          "${BUCKET_ROOT}/${year}/${month_padded}/${day_padded}/${hour}/" \
          "${data_dir}/${year}/${month_padded}/${day_padded}/${hour}/"
      done
    done
    echo "=== Finished ${year}-${month_padded} ==="
  done
done
