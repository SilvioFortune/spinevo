#!/bin/bash
julia   <<EOD
include("/home/moon/sfortune/spinevo/pkg/meta.jl")
flush(stdout)
assemble_halostories(central_only=true, outdir="/home/moon/sfortune/spinevo/mergerimpact")

EOD