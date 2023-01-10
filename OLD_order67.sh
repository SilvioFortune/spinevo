#!/bin/bash
julia   <<EOD
include("/home/moon/sfortune/spinevo/pkg/meta.jl")
flush(stdout)
tree2story(start=1160, stop=1220, iterstep=1, skipexisting=false, outdir="/home/moon/sfortune/spinevo/data/centralstories_nobreak_physical_v20220127_min0.0Gyr")

EOD