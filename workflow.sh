# Required add-ons
#   r.cell.area
#   r.stream.extract (now mainstream?)
#   v.stream.network
#   v.stream.profiler


######################
# COMPUTER VARIABLES #
######################

memory_Mb=8000 # Mb for analysis


############
# SETTINGS #
############

Athresh_km2=10 # Threshold drainage area to define a channel
window_meters=1000 # Smoothing window (moving average)
dx_target=100 # Downsample DEM to this dx when constructing profiles (speed)


###############
# INPUT NAMES #
###############

DEM=DEM


################
# OUTPUT NAMES #
################

slope=slope
accum=accum
cellsize_km2=cellsize_km2
streams=streams
outfile_basename=Cottonwood
nullmask=nullmask # if DEM is not perfectly rectangular
outstream=LittleCottonwood

###########
# PROCESS #
###########

g.region -p rast=$DEM --o
r.cell.area output=$cellsize_km2 units=km2 --o
r.slope.aspect elevation=$DEM slope=$slope format=percent --o
r.mapcalc "slope = slope / 100." --o
# Have to run r.watershed separately because r.stream.extract does not support 
# creation of its own accumulation grids or weighting of flow accumulation
r.watershed elevation=$DEM accumulation=$accum flow=$cellsize_km2 -s \
            memory=$memory_Mb --o
r.mapcalc "$nullmask = $DEM * 0 + 1" --o
r.mapcalc "$accum = $accum * $nullmask" --o
r.stream.extract elevation=$DEM accumulation=$accum threshold=$Athresh_km2 \
                 stream_vector=$streams d8cut=0 memory=$memory_Mb --o
v.stream.network map=$streams

# This gives different slopes and drainage areas than expected -- perhaps SFD needed

#!!!!!!!!!!!!!!!!!!!!
start_segment_cat=135 # !!! USE GSFLOW--GRASS TOOLS TO HELP SELECT THIS!
#!!!!!!!!!!!!!!!!!!!!
v.stream.profiler cat=$start_segment_cat streams=$streams direction=upstream \
                  elevation=$DEM accumulation=$accum slope=$slope units=km2 \
                  plots=LongProfile,SlopeAccum,SlopeDistance,AccumDistance \
                  window=$window_meters dx_target=$dx_target \
                  outstream=$outstream \
                  outfile=${outfile_basename}.out \
                  --o
                  
# Update this to save first.
