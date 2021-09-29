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


### Read in

bvalues     = load(joinpath(@__DIR__, "results.jld"), "bvalues")
neighbors   = load(joinpath(@__DIR__, "results.jld"), "neighbors")



### Figure

fig, ax = subplots()

ax.plot(neighbors, bvalues, "b.", label="> 1e9 M⊙", alpha= 0.2)

ax.set_xlabel("N neighbors")
ax.set_ylabel("b-value")
ax.grid()
ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)

scale=0.5

fig.set_size_inches(16scale, 9scale)
fig.tight_layout()


fig.savefig(joinpath(@__DIR__, "plot.pdf"), bbox_inches="tight", pad_inches=.1)



# Output

println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!")
