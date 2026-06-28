#!/usr/bin/env bash
# Mean distance from a hillslope cell to the channel network, vs resolution, on
# the Ohio (Mid Bailey Run) DrEICH networks. Quantifies the network retreat as
# heads are lost: as resolution coarsens, hillslopes effectively lengthen.
#
# ALONG-FLOW-PATH distance via r.stream.distance (overland flow length to the
# stream), with the flow direction from r.watershed -s (SFD, ~matches the DrEICH
# routing). Truth network = 1 m; coarse networks = adapt2 (adaptive threshold).
#
# Run in the MidBaileyRun/dreich_restest GRASS session:
#   PROJ_DATA=/usr/share/proj grass <...>/dreich_restest \
#     --exec bash dev/dreich_distance_to_channel.sh
set -u
EXT="n=4369656.1523925 s=4362549.1523925 w=398767.32685636 e=405033.32685636"
echo "res  network          mean_flowdist(m)  median(m)"
for spec in "1:rnet_std_1m" "2:rnet_adapt2_2m" "3:rnet_adapt2_3m" "4:rnet_adapt2_4m" \
            "5:rnet_adapt2_5m" "8:rnet_adapt2_8m" "10:rnet_adapt2_10m" \
            "12:rnet_adapt2_12m" "15:rnet_adapt2_15m" "20:rnet_adapt2_20m" "30:rnet_adapt2_30m"; do
  N=${spec%%:*}; R=${spec##*:}
  g.region $EXT res=$N -a 2>/dev/null
  if ! g.findfile element=cell file=$R mapset=dreich_restest >/dev/null 2>&1; then
    echo "$N  $R  (absent)"; continue; fi
  r.watershed -s elevation=dem_${N}m drainage=_dir --o --q 2>/dev/null
  r.stream.distance stream_rast=$R direction=_dir method=downstream distance=_fd --o --q 2>/dev/null
  r.mapcalc "_hd = if(isnull($R), _fd, null())" --o --q 2>/dev/null      # hillslope cells only
  read mn md < <(r.univar -g -e map=_hd 2>/dev/null | awk -F= '/^mean=/{m=$2}/^median=/{d=$2}END{print m, d}')
  printf "%-3s %-16s %16.1f %11.1f\n" "$N" "$R" "$mn" "$md"
done
