#!/usr/bin/env bash
# DrEICH channel-head resolution-sensitivity sweep on the Ohio test case
# (Mid Bailey Run, USGS 1 m LiDAR -> averaged to coarser grids).
#
# Run INSIDE a GRASS session in the MidBaileyRun location, e.g.:
#   REPO=$(pwd) PROJ_DATA=/usr/share/proj PYTHONPATH="$REPO" \
#     grass /home/awickert/grassdata/MidBaileyRun/dreich_restest \
#       --exec bash dev/dreich_resolution_sweep.sh std
#
# Arg 1 = scheme:  std  (window=7*res cells, tan_curv_threshold=0.1; the standard
#                        1 m-tuned LSDTT/DrEICH parameters at every resolution)
#                  adapt (window=7*res cells, tan_curv_threshold = a fixed high
#                        percentile of each resolution's own curvature, anchored
#                        to what 0.1 is at 1 m -- see CURV_PCTL below)
# Routing is the internal FastScape-style D8 (the canonical DrEICH path).
# DEM coarsening is mean aggregation (r.resamp.stats method=average).
#
# Per resolution N it writes (in the current mapset):
#   dem_${N}m, heads_${SCHEME}_${N}m, net_${SCHEME}_${N}m, rnet_${SCHEME}_${N}m
# and appends "N  heads  status" to $SUMMARY.
set -u
SCHEME="${1:-std}"
REPO="${REPO:-$(pwd)}"
CH="$REPO/r.fluvial.channelheads/r.fluvial.channelheads.py"
EXTENT="n=4369656.1523925 s=4362549.1523925 w=398767.32685636 e=405033.32685636"
RESOLUTIONS="${RESOLUTIONS:-1 2 3 4 5 8 10 12 15 20 30}"
SUMMARY="${SUMMARY:-/tmp/dreich_restest_${SCHEME}.txt}"
CURV_PCTL="${CURV_PCTL:-/tmp/dreich_curv_pctl.txt}"   # written by the adapt calibrator

printf "# scheme=%s  res_m  heads  tan_curv_threshold  status\n" "$SCHEME" > "$SUMMARY"
for N in $RESOLUTIONS; do
  g.region $EXTENT res=$N -a 2>/dev/null
  r.resamp.stats input=dem output=dem_${N}m method=average --o --q 2>/dev/null
  # window mode: grow7 = 7*res cells (default); const = hillslope-scale W_h held
  # constant where valid, else the smallest valid window (~1.5*res, kr>=2).
  case "${WINMODE:-grow7}" in
    const) WR=$(python3 -c "import math;print(max(${WH:-7}, math.ceil(1.5*$N)))");;
    *)     WR=$((7 * N));;
  esac
  if [ "$SCHEME" != "std" ]; then          # any adaptive scheme uses the calibration
    TC=$(grep -E "^${N} " "$CURV_PCTL" 2>/dev/null | awk '{print $2}')
    [ -z "$TC" ] && TC=0.1
  else
    TC=0.1
  fi
  log=$(python3 "$CH" method=dreich elevation=dem_${N}m \
        points=heads_${SCHEME}_${N}m network=net_${SCHEME}_${N}m \
        raster_network=rnet_${SCHEME}_${N}m \
        window_radius=$WR tan_curv_threshold=$TC threshold=${SRCTHR:-100} a_0=1000 \
        m_over_n=0.525 n_connecting_nodes=${NCON:-10} \
        min_segment_length=${MINSEG:-10} --o 2>&1)
  n=$(echo "$log" | grep -oE "Found [0-9]+" | grep -oE "[0-9]+")
  if [ -n "$n" ]; then st="ok"; else n=0; st="no_heads"; fi
  printf "%s  %s  %s  %s\n" "$N" "$n" "$TC" "$st" | tee -a "$SUMMARY"
done
echo "DONE ($SCHEME). Summary: $SUMMARY"
