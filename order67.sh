#!/bin/bash
julia   <<EOD
include("/home/moon/sfortune/spinevo/pkg/meta.jl")
flush(stdout)
assemble_halostories(mass_ST_thr=1e10, indir="/home/moon/sfortune/spinevo/data/silvio_stories_spinmap", central_only=true, outdir="/home/moon/sfortune/phd/bathtub/data", outfile="assembly_spinmap_centrals_switch_03Gyr.jld")

EOD

# endstate_quickndirty()
# update_mergercollections(story_dir="/home/moon/sfortune/spinevo/data/silvio_stories_spinmap_noswitch", min_time=0.0)