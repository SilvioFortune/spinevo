#!/bin/bash
julia   <<EOD
include("/home/moon/sfortune/spinevo/pkg/meta.jl")
flush(stdout)
felix2jld(start=500, stop=100, iterstep=-1, outdir="/home/moon/sfortune/spinevo/halostories_v20211204_min0.0Gyr")

EOD