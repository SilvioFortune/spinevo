println("   Processing...")


@doc """
DESCRIPTION:\n
    - reduce particle number by distribution to a 3D-grid
INPUT:\n
    - n_cells:          approximate number of cells for the final grid
    - positions:        (3 x n)D array for the positions of n reduce_n_particles
    - masses:           nD array for the particle masses
    - focus_type:       "center", "outerrim" or "linear"; determines where the cell resolution should be highest
    - focus_strength:   [0. , 1.0] determines degree of focus
    - m_th:             minimum cell mass fraction of total mass to accept the cell in the final output, used to get rid of empty cells
    - verbose:          true or false; progress output
OUTPUT:\n
    - pos_new:  (3 x m)D array of mass-weighted cell positions
    - mnew:     mD array of cell masses
    - n_part:   mD array of number of particles in each cell
""" ->
function reduce_n_particles(n_cells, positions, masses; focus_type="center", focus_strength=0.7, m_th=1e-6, verbose=true)
    
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



########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################



@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function fake_flip_finder(indict, idx; verbose=false)
    flip1 = 3
    if 1 < idx < length(indict["M_felix"]) && 
    count(ismissing, indict["j_main"][:,idx-1]) + count(ismissing, indict["j_main"][:,idx]) + 
    count(ismissing, indict["j_main"][:,idx+1]) == 0 # can be checked
        flip1 = ifelse( norm(indict["j_main"][:,idx+1]-indict["j_main"][:,idx-1]) < norm(indict["j_main"][:,idx]-indict["j_main"][:,idx-1]) && # && indict["M_felix"][idx+1]-indict["M_felix"][idx-1] < indict["M_felix"][idx]-indict["M_felix"][idx-1], 
                        norm(indict["j_main"][:,idx+1]-indict["j_main"][:,idx-1]) < norm(indict["j_main"][:,idx+1]-indict["j_main"][:,idx]),
                        1, 0 
                        )
    end

    flip2 = 3
    if 2 < idx && 
    count(ismissing, indict["j_main"][:,idx-2]) + count(ismissing, indict["j_main"][:,idx-1]) + 
    count(ismissing, indict["j_main"][:,idx]) == 0 # can be checked
        flip2 = ifelse( norm(indict["j_main"][:,idx]-indict["j_main"][:,idx-2]) < norm(indict["j_main"][:,idx-1]-indict["j_main"][:,idx-2]) && # && indict["M_felix"][idx]-indict["M_felix"][idx-2] < indict["M_felix"][idx-1]-indict["M_felix"][idx-2], 
                        norm(indict["j_main"][:,idx]-indict["j_main"][:,idx-2]) < norm(indict["j_main"][:,idx]-indict["j_main"][:,idx-1]),
                        1, 0 
                        )
    end
    if verbose
        println("Snap $(indict["SNAP"][idx])   ---   in = $flip1, out = $flip2")
    end

    if flip1 + flip2 == 0 # both not fake
        return 0
    elseif flip1 + flip2 < 3 # any is fake
        return 1
    else
        return missing
    end
end

print("'fake_flip_finder'   ")



########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################



# Functions that are too long and have own file
include("/home/moon/sfortune/spinevo/pkg/tree2story.jl")
include("/home/moon/sfortune/spinevo/pkg/assemble_halostories.jl")
include("/home/moon/sfortune/spinevo/pkg/treewalker.jl")



########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################



println()