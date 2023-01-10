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



# Functions

function bval(m,j)
    return log10(norm(j.*1u"s/kpc/km")) - (2/3)log10(m/(1u"Msun"))
end



# Initial Parameters

snap_id     = 136
box4        = "/HydroSims/Magneticum/Box4/uhr_test/"
filepath    = string(box4, "groups_$snap_id/sub_$snap_id")
isub        = 5001
igroup      = 10000



# Calculations

nsub        = read_subfind(filepath, "NSUB")
spinvecs    = read_subfind(filepath, "SPIN")
spins       = zeros(length(spinvecs[1,:]))
for i in 1:length(spins)
    spins[i] = norm(spinvecs[:,i])
end
bvalues     = zeros(length(spins))


snapshot    = Snapshot(box4, snap_id)
snapshot.snapbase
snapshot.subbase

mass_limit  = 1e10u"Msun"
fraction    = 12/length(bvalues)

id_limit    = Int(floor(length(bvalues) * fraction))
below_mass_limit_counter    = 0
Threads.@threads for i in 1:id_limit
    #println("Progress: ", 100*i/id_limit, " %   ---   $i / ", id_limit)
    g = Galaxy(snapshot, i-1)
    try
        read_halo!(g, props=((:stars,["MASS", "POS", "VEL"]),), units=:full)
        # sfc: DONE IN save2save.jl: implement mass cut to leave out halos that are too small
        # sfc: radius für masse auswählen
        m_star_tot  = sum(g.stars.mass)
        if m_star_tot > mass_limit
            #zeros(typeof(1u"km*kpc*Msun/s"), 3) # alternative initialization
            j_stars     = (g.stars.mass[1]/m_star_tot) .* (g.stars.pos[:,1] × g.stars.vel[:,1])
            for ii in 2:length(g.stars.mass)
                j_stars .+= (g.stars.mass[ii]/m_star_tot) .* (g.stars.pos[:,ii] × g.stars.vel[:,ii])
            end
            bvalues[i]  = bval(m_star_tot,j_stars)
        else
            # sfc: done in append mass limit id array
        end
    catch error
        println("Skipping Error in ID: $(i-1)")
        # sfc: append skipped id
    end
end

neighbors       = zeros(length(spins))
b_stats         = zeros(2, length(nsub))

k = 0
for i in 1:length(nsub)
    for j in 1:nsub[i]
        neighbors[j+k] = nsub[i] - 1
    end
    b_stats[1, i]   = Statistics.mean(bvalues[k+1:k+nsub[i]])
    b_stats[2, i]   = Statistics.std(bvalues[k+1:k+nsub[i]])
    k += nsub[i]
end



# Figure

fig, ax = subplots()

ax.plot(neighbors[1:id_limit], bvalues[1:id_limit], "b.", label="Subhalos", alpha= 0.1)
ax.plot(nsub .- 1, b_stats[1,:], "r.", label="Mean value", alpha= 1)

ax.set_xlabel("N neighbors")
ax.set_ylabel("b-value")
ax.grid()
ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)

scale=0.5

fig.set_size_inches(16scale, 9scale)
fig.tight_layout()
fig.savefig(joinpath(@__DIR__, "test.pdf"), bbox_inches="tight", pad_inches=.1)


# Save Calculations
save("~/nsub_vs_b-value/save.jld", "neighbors", neighbors, "bvalues", bvalues)

# Output
println(size(bvalues))
println(maximum(bvalues))
println(minimum(bvalues))
println("$below_mass_limit_counter subhalos below mass limit of $mass_limit")