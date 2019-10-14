DEM=DEM
slope=slope
accum=accum
draindir=draindir # not needed
Athresh_km2=10
cellsize_km2=cellsize_km2
streams=streams
window_meters=200
outfile_basename=Toro
nullmask=nullmask # if DEM is not perfectly rectangular

g.region -p rast=$DEM --o
r.cell.area output=$cellsize_km2 units=km2 --o
r.slope.aspect elevation=$DEM slope=$slope --o
r.watershed elevation=$DEM accumulation=$accum flow=$cellsize_km2 --o
r.mapcalc "$nullmask = $DEM * 0 + 1" --o
r.mapcalc "$accum = $accum * $nullmask" --o
r.stream.extract elevation=$DEM accumulation=$accum threshold=$Athresh_km2 \
                 stream_vector=$streams --o # Note: MFD
v.stream.network map=$streams

# This gives different slopes and drainage areas than expected -- perhaps SFD needed

start_segment_cat=101
v.stream.profiler cat=$start_segment_cat streams=$streams direction=downstream \
                  elevation=$DEM accumulation=$accum slope=$slope units=km2 \
                  window=$window_meters \
                  plots=LongProfile,SlopeAccum,SlopeDistance,AccumDistance \
                  outfile_original=${outfile_basename}_fullres.out \
                  outfile_smoothed=${outfile_basename}_${window_meters}m_average.out
                  
# Update this to save first.
