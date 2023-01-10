#!/bin/bash
julia   <<EOD
include("/home/moon/sfortune/spinevo/pkg/meta.jl")
flush(stdout)
quickadd_feature("hs"; target="/home/moon/sfortune/spinevo/data/silvio_stories_spinmap", min_time=0.0)

EOD


#