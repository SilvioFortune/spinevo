#!/bin/bash
julia   <<EOD
include("/home/moon/sfortune/spinevo/pkg/meta.jl")
flush(stdout)
write_halo_stories(start=519, iterstep=-1,stop=350, central_switch=false, root_mass_thr=1e10, subtype = "central", root_list="/home/moon/sfortune/spinevo/data/root_list_20220306.jld", outdir="/home/moon/sfortune/spinevo/data/silvio_stories_spinmap_noswitch")

EOD

#1 start=1, stop=49, 
#2 start=120, stop=519, 
#3 start=520, stop=919, 
#4 start=920, stop=1319, 

#11 start=50, stop=119, 
#21 start=519, iterstep=-1,stop=350, 



#1 start=1, stop=300, 
#2 start=600, iterstep=-1,stop=301, 
#3 start=601, stop=900, 
#4 start=1319, iterstep=-1,stop=901, 