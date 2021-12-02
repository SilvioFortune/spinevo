#!/bin/bash
julia   <<EOD
include("/home/moon/sfortune/spinevo/pkg/meta.jl")
flush(stdout)
felix2jld(start=1501, outdir="/home/moon/sfortune/spinevo/halostories_v20211127_min0.0Gyr")

EOD