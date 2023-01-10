#!/bin/bash
julia   <<EOD
include("/home/moon/sfortune/spinevo/pkg/meta.jl")
flush(stdout)
endstate_quickndirty()

EOD


#endstate_quickndirty()
#assemble_halostories(indir="/home/moon/sfortune/spinevo/data/silvio_stories_spinmap_09Gyr", central_only=true, outdir="/home/moon/sfortune/spinevo/data", outfile="assembly_spinmap_09Gyr.jld")