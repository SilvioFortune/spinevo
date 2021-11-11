# Collection of Functions and Packages 

##################################################################
### Packages, sfc: needed?

using Printf
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
using Unitful
using UnitfulAstro
using Missings
using Distributions
using CSV
using DataFrames
using Cosmology
using Suppressor
using StatsBase
using ColorSchemes

##################################################################



##################################################################
### Functions
println("Loading Functions:")

function align_plane(vone, vtwo, axis=[0,0,1])
    """
    DESCRIPTION:
        - Rotate two vectors by aligning their plane with a given axis in a way that the third dimension is 0 for both
        - sfc: check consistencies: Lengths, Angles
    INPUT:
        - vone: first vector
        - vtwo: second vector
        - axis: axis to align their plane; DEFAULT: [0,0,1]
    OUTPUT:
        - rotated first vector
        - rotated second vector
    """
    start       = normalize(vone × vtwo)
    axis_norm   = normalize(axis)
    vx      =  start × axis_norm
    cx      = transpose(start) * axis_norm
    Vx      = [0 -vx[3] vx[2]; vx[3] 0 -vx[1]; -vx[2] vx[1] 0]
    R = I + Vx + (Vx * Vx ./ (1+cx))
    return R*vone, R*vtwo
end
print("'align_plane'   ")

function orbit_j(subID, centID, snapNR; boxNR=4) # Find central sub maybe using SUBFIND -> centID = FSUB[GRNR[subID+1]+1]
    """
    DESCRIPTION:
        - Calculate the orbital angular momentum with respect to the central subhalo
        - sfc: WARNING! This does not feel right. instead, it should be necessary to calculate everything with respect to the common center of mass
    INPUT:
        - subID:    sub halo ID
        - centID:   central halo ID
        - snapNR:   snap number
        - boxNR:    box number
    OUTPUT:
        - mass-independent orbital angular momentum vector
    """
    if boxNR == 4
        pathtofile  = "/HydroSims/Magneticum/Box4/uhr_test/groups_$(@sprintf("%03i", snapNR))/sub_$(@sprintf("%03i", snapNR))"
    else
        error("boxNR not known.")
    end
    # Read subhalo features
    head        = read_header(pathtofile)
    spos        = convert_units_physical(read_subfind(pathtofile, "SPOS"), :pos, head)
    svel        = read_subfind(pathtofile, "SVEL") # sfc: no conversion since already physical units??
    # Return with masses added
    return ((spos[:,subID+1] .- spos[:,centID+1]) × (svel[:,subID+1] .- svel[:,centID+1]))
end
print("'orbit_j'   ")


# Calculate Total Halo velocity
function avg_velocity(m_ar, vel_ar)
    """
    DESCRIPTION:
        - mean velocity by mass-weighted velocities
    INPUT:
        - m_ar:     mass array dimensions: n
        - vel_ar:   velocity array dimensions: (3,n)
    OUTPUT:
        - mean velocity vector
    """
    tot_vel = zeros(3)
    for i in 1:length(m_ar)
        tot_vel += m_ar[i] .* vel_ar[:,i]
    end
    tot_vel ./= sum(m_ar)
    return tot_vel
end
print("'avg_velocity'   ")


function borders_tan(coord, steplength, flim)
    """
    DESCRIPTION:
        - construct bin limits for an array with a tangential function for higher resolution in the median central part
    INPUT:
        - coord:        1D-array that is to be divided into a grid
        - steplength:   scale length for bin size
        - flim:         upper (and lower) limit of the tangens x-range, 0 < flim < pi/2, the larger the more center-focused
    OUTPUT:
        - array of grid (bin) borders
    """
    a           = median(coord)
    coord_ext   = maximum(coord)-minimum(coord)
    b           = coord_ext / ( 2*tan(flim) )
    return a .+ b .* tan.( LinRange(-flim, flim, Int(round(coord_ext/steplength))+1) )
end
print("'borders_tan'   ")

function borders_arctan(coord, steplength, flim)
    """
    DESCRIPTION:
        - construct bin limits for an array with an arctan function for higher resolution in the outer rim
    INPUT:
        - coord:        1D-array that is to be divided into a grid
        - steplength:   scale length for bin size
        - flim:         upper (and lower) limit of the arctan x-range, the larger the more outer-rim-focused
    OUTPUT:
        - array of grid (bin) borders
    """
    a           = median(coord)
    coord_ext   = maximum(coord)-minimum(coord)
    b           = coord_ext / ( 2*atan(flim) )
    return a .+ b .* atan.( LinRange(-flim, flim, Int(round(coord_ext/steplength))+1) )
end
print("'borders_arctan'   ")

function reduce_n_particles(n_cells, positions, masses; focus_type="center", focus_strength=0.7, m_th=1e-6, verbose=true)
    """
    DESCRIPTION:
        - reduce particle number by distribution to a 3D-grid
    INPUT:
        - n_cells:          approximate number of cells for the final grid
        - positions:        (3 x n)D array for the positions of n reduce_n_particles
        - masses:           nD array for the particle masses
        - focus_type:       "center", "outerrim" or "linear"; determines where the cell resolution should be highest
        - focus_strength:   [0. , 1.0] determines degree of focus
        - m_th:             minimum cell mass fraction of total mass to accept the cell in the final output, used to get rid of empty cells
        - verbose:          true or false; progress output
    OUTPUT:
        - pos_new:  (3 x m)D array of mass-weighted cell positions
        - mnew:     mD array of cell masses
        - n_part:   mD array of number of particles in each cell
    """
    
    pos_new = Array{Float64}(undef, 3, 0)
    mnew    = Array{Float64}(undef, 0)
    n_part      = Array{Int64}(undef, 0)
    mtot    = sum(masses)


    x_ext   = maximum(positions[1,:])-minimum(positions[1,:])
    y_ext   = maximum(positions[2,:])-minimum(positions[2,:])
    z_ext   = maximum(positions[3,:])-minimum(positions[3,:])
    cell_length = maximum([z_ext,y_ext,x_ext]) / cbrt(n_cells)
    if focus_type == "center"
        limits  = 0.7 + focus_strength*0.8    # x-range for tan function: 0.7 to 1.5
        xborders    = borders_tan(positions[1,:], cell_length, limits)
        yborders    = borders_tan(positions[2,:], cell_length, limits)
        zborders    = borders_tan(positions[3,:], cell_length, limits)
    elseif focus_type == "outerrim"
        limits  = 0.5 + focus_strength*4.5    # x-range for arctan function: 0.5 to 5
        xborders    = borders_arctan(positions[1,:], cell_length, limits)
        yborders    = borders_arctan(positions[2,:], cell_length, limits)
        zborders    = borders_arctan(positions[3,:], cell_length, limits)
    elseif focus_type == "linear"
        xborders    = LinRange(minimum(positions[1,:]), maximum(positions[1,:]), Int(round(x_ext/cell_length))+1)
        yborders    = LinRange(minimum(positions[2,:]), maximum(positions[2,:]), Int(round(y_ext/cell_length))+1)
        zborders    = LinRange(minimum(positions[3,:]), maximum(positions[3,:]), Int(round(z_ext/cell_length))+1)
    else
        error("Unknown Focus type")
    end
    #cell_IDs    = Array{Int64}(undef, 3, length(masses)) # assign each particle to a cell

    # looping over upper cell limits to assign to cells
    if verbose == true
        println("$(length(xborders))   ---   $(length(yborders))   ---   $(length(zborders))")
    end
    #Threads.@threads for i in 2:length(xborders)
    for i in 2:length(xborders)
        if verbose == true
            println("$i ")
        end
        for ii in 2:length(yborders)
            for iii in 2:length(zborders)
                m_temp      = Array{Float64}(undef, 0)
                pos_temp    = Array{Float64}(undef, 3, 0)
                n_temp      = 0
                for iiii in 1:length(masses)
                    if xborders[i-1] <= positions[1,iiii] <= xborders[i] && yborders[ii-1] <= positions[2,iiii] <= yborders[ii] && zborders[iii-1] <= positions[3,iiii] <= zborders[iii]
                        #cell_IDs[:,iiii] = [i,ii,iii]
                        m_temp           = vcat(m_temp, masses[iiii])
                        pos_temp         = hcat(pos_temp, positions[:,iiii])
                        n_temp          += 1
                    end
                end
                if sum(m_temp) > mtot*m_th
                    n_part  = vcat(n_part, n_temp)
                    mnew    = vcat(mnew, sum(m_temp))
                    #pos_new = hcat(pos_new, [ mean(pos_temp[1,:], weights(m_temp)), mean(pos_temp[2,:], weights(m_temp)), mean(pos_temp[3,:], weights(m_temp)) ])
                    pos_new = hcat(pos_new, mean(pos_temp, weights(m_temp); dims=2))
                end
            end
        end
    end
    println("$(size(pos_new))   $(size(mnew))   $(size(n_part))")
    return pos_new, mnew, n_part
end
print("'reduce_n_particles'   ")

function find_felixID(isub; dir="/home/moon/sfortune/spinevo/halostories_v20211007_min0.0Gyr")
    """
    DESCRIPTION:
        - find halo ID according to Felix' notation based on the ID in the final snap shot
    INPUT:
        - isub: subfind ID of the halo in the final snap shot
        - dir:  directory with all the processed Silvio-processed halo files
    OUTPUT:
        - isub:             Int; halo ID according to Felix' notation
        - halo_filestring:  String; jld filename of that halo
    """
    storyfilelist   = readdir(dir)
    halo_filestring = " "
    ifelix          = Int64
    for i in 1:length(storyfilelist)
        if occursin("_$(isub).jld", storyfilelist[i])
            halo_filestring = storyfilelist[i]
            ifelix          = parse( Int64, chop(replace(halo_filestring, "$(isub)" => "" ), head=5,tail=5) )
        end
    end
    return ifelix, halo_filestring
end
print("'find_felixID'   ")

function find_subID(ifelix; dir="/home/moon/sfortune/spinevo/halostories_v20211007_min0.0Gyr")
    """
    DESCRIPTION:
        - find sub ID in the final snap shot from to Felix' notation ID
    INPUT:
        - ifelix:   halo ID from Felix
        - dir:      directory with all the processed Silvio-processed halo files
    OUTPUT:
        - isub:             Int64; sub ID in the final snap shot
        - halo_filestring:  String; jld filename of that halo
    """
    storyfilelist   = readdir(dir)
    halo_filestring = " "
    isub            = Int64
    for i in 1:length(storyfilelist)
        if occursin("halo_$(ifelix)_", storyfilelist[i])
            halo_filestring = storyfilelist[i]
            isub            = parse( Int64, chop( replace( halo_filestring, "$(ifelix)" => "" ), head=6,tail=4 ) )
        end
    end
    return isub, halo_filestring
end
print("'find_ISUB'   ")

function ssFPfinder(mID, halsto)
    """
    DESCRIPTION:
        - Find same-snap First Progenitor index depending on merger index in Felix' list
    INPUT:
        - mID:      merger index
        - halsto:   specific halostory dictionary
    OUTPUT:
        - result_ID:    Int64; same-snap First Progenitor index
    """
    result_ID = Int64
    for finder_ID in mID:-1:2    # walking the array towards the first progenitor
        if halsto["SNAP"][finder_ID] < halsto["SNAP"][finder_ID-1]
            result_ID   = finder_ID
            break
        end
    end
    return result_ID
end
print("'ssFPfinder'   ")

function cartesian_to_spherical(x)
    """
    DESCRIPTION:
        - Tranform Cartesian vector into [r,θ,ϕ] Spherical vector
    INPUT:
        - x:    cartesian vector
    OUTPUT:
        - s:    spherical vector in [r,θ,ϕ]
    """
    s = zeros(3)
    s[1]    = norm(x) # Radius
    s[2]    = acosd(x[3] / s[1])# * 180 / π # θ[°]
    s[3]    = acosd(x[1] / sqrt(x[1]*x[1] + x[2]*x[2]))# * 180 / π # ϕ[°]
    return s
end
print("'cartesian_to_spherical'   ")

function aligner(x; ref=[0,0,1])
    """
    DESCRIPTION:
        - Rotation Matrix to align x with reference vector
    INPUT:
        - x:    cartesian vector
        - ref:  reference vector
    OUTPUT:
        - rotmat:   rotation matrix to align x with ref
    """
    x       = x ./ norm(x)
    ref     = ref ./ norm(ref)
    vx      = x × ref
    cx      = transpose(x) * ref
    Vx      = [0 -vx[3] vx[2]; vx[3] 0 -vx[1]; -vx[2] vx[1] 0]
    rotmat  = I + Vx + (Vx * Vx ./ (1+cx))
    return rotmat
end
print("'aligner'   ")


println()
##################################################################