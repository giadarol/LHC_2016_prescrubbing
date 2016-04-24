Some example commands:

/afs/cern.ch/work/l/lhcscrub/LHC_2016_prescrubbing

python 016a_heat_load_period.py tstart=22-04-2016,12:00 tend=22-04-2016,19:00 plotaverage=n plotall=y timein=hourtime tagfills=y varlists=AVG_ARC mode=norm_to_intensity normlength=none hourtickspac=1

run 016_heat_load_arcs_vs_SAMs.py 4615  --noaverage --zeroat=1.2
