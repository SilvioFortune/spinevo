# Packages

using LinearAlgebra
using Printf
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



### Functions

function bval(m,j)
    return log10(norm(j.*1u"s/kpc/km")) - (2/3)log10(m/(1u"Msun"))
end

"""
function jₖ() # Specific Angular momentum of component k (Teklu,Rhemus)
    out = zeros(3)

    return ()
end
#"""



### Initial Parameters

snap_id     = 136
box4        = "/HydroSims/Magneticum/Box4/uhr_test/"
boxdir      = string(box4, "groups_$snap_id/sub_$snap_id")
savedir     = "/home/moon/sfortune/nsub_vs_b-value/files1/"

mass_limit  = 1e9u"Msun"



### Calculations

nsub        = read_subfind(boxdir, "NSUB")
nsub_sum = sum(nsub)

neighbors   = zeros(nsub_sum)
bvalues     = zeros(nsub_sum)

k = 0
for i in 1:length(nsub)
    for ii in 1:nsub[i]
        neighbors[ii+k] = nsub[i] - 1
    end
    k += nsub[i]
end

# sfc: radius für masse auswählen
# Loading Data
for i in 1:nsub_sum # Loop over expected number of subhalos
    try
        masses      = load("/home/moon/sfortune/nsub_vs_b-value/files1/subid_$(i-1).jld", "mass") * 1u"Msun"
        m_star_tot  = sum(masses)
        lnth        = length(masses)
        if m_star_tot > mass_limit
            print("$(@sprintf("%.2f", 100*i/(nsub_sum))) %   ---   $m_star_tot\r")
            pos         = load("/home/moon/sfortune/nsub_vs_b-value/files1/subid_$(i-1).jld", "pos") .* 1u"kpc"
            vel         = load("/home/moon/sfortune/nsub_vs_b-value/files1/subid_$(i-1).jld", "vel") .* 1u"km/s"
            
            j_stars = (masses[1]/m_star_tot) .* (pos[:,1] × vel[:,1])
            for ii in 2:lnth
                j_stars .+= (masses[ii]/m_star_tot) .* (pos[:,ii] × vel[:,ii])
            end
            bvalues[i]  = bval(m_star_tot, j_stars)
        else
            println("Halo with ID $(i-1) has stellar mass of $m_star_tot < $mass_limit")
            # store ID in array via cat
        end
    catch e
        println("Halo with ID $(i-1) has been destroyed!")
        #sfc: store ID in array via cat
    end
end

save("/home/moon/sfortune/nsub_vs_b-value/results.jld", "bvalues", bvalues, 
                                                    "neighbors", neighbors)

# Output

println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!")
