"""
Read features from snapshot and store them in file via JLD package
"""


# Packages

using LinearAlgebra
using JLD
using QuadGK
using Statistics
using PyCall
using PyPlot
using LaTeXStrings
using GadgetIO
using GadgetUnits
using GadgetGalaxies
using UnitfulAstro
using Missings



# Initial Parameters

snap_id     = 136
box4        = "/HydroSims/Magneticum/Box4/uhr_test/"
filepath    = string(box4, "groups_$snap_id/sub_$snap_id")
isub        = 166478



# Prerequisites

snapshot    = Snapshot(box4, snap_id)
snapshot.snapbase
snapshot.subbase

isub        = Int.(read_subfind(filepath, "FSUB"))
nsub        = Int.(read_subfind(filepath, "NSUB"))
nsub_sum    = sum(nsub)
#neighbors   = missings(UInt8, nsub_sum)
neighbors   = zeros(nsub_sum)
# k = 0
Threads.@threads for i in 1:length(nsub)
    for ii in 1:nsub[i]
        temp = i-1
        #current = nsub[1:temp]
        neighbors[ii+sum(nsub[1:temp])] = nsub[i] - 1
    end
    #k += nsub[i]
end

file_success    = open(joinpath(@__DIR__, "success.txt"), "a+")
file_error      = open(joinpath(@__DIR__, "errors.txt"), "a+")



### The actual thing

# Reading and saving the subhalos
last_id = sum(nsub) - 1
println("Reading from subid 0 to $last_id")
#Threads.@threads for i in 0:last_id
@time for i in 0:last_id
    g = Galaxy(snapshot, i)
    try
        #read_halo!(g, props=((:stars,["MASS", "POS", "VEL"]),), units=:physical)
        #save("/home/moon/sfortune/nsub_vs_b-value/subs_snap_$snap_id/subid_$ii.jld",    "subid", g.subid, 
                                                    #"mass", g.stars.mass, 
                                                    #"pos", g.stars.pos, 
                                                    #"vel", g.stars.vel)
        read_halo!(g, units=:physical)
        save(joinpath(@__DIR__, "subs_snap$snap_id/subid_$i.jld"), 
                "subhalo",      g, 
                "neighbors",    neighbors[i+1])
        write(file_success, "Saved SubID $i using Thread $(Threads.threadid()).\n")
        #println("Saved SubID $i using Thread $(Threads.threadid()).")
    catch e
        write(file_error, "SubID $i Thread $(Threads.threadid())\n$e\n\n\n")
        #println("Error caught on SubID $i using Thread $(Threads.threadid()).")
    end
end

close(file_success)
close(file_error)


# Reading and saving the subhalos
#last_id = 10#sum(nsub) - 1
#println(last_id)
#for i in 0:last_id
#    g = Galaxy(snapshot, i)
#    read_halo!(g, units=:physical)
#    save("/home/moon/sfortune/nsub_vs_b-value/subs_snap_$snap_id/subid_$i.jld", "g", g)
#end


# Output

println("\n---------------------------\n
Now witness the firepower of this fully armed and operational battle station!")