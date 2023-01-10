println("   Auxiliary...")


@doc """
DESCRIPTION:\n
    - construct bin limits for an array with a tangential function for higher resolution in the median central part
INPUT:\n
    - coord:        1D-array that is to be divided into a grid
    - steplength:   scale length for bin size
    - flim:         upper (and lower) limit of the tangens x-range, 0 < flim < pi/2, the larger the more center-focused
OUTPUT:\n
    - array of grid (bin) borders
""" ->
function borders_tan(coord, steplength, flim)
    a           = median(coord)
    coord_ext   = maximum(coord)-minimum(coord)
    b           = coord_ext / ( 2*tan(flim) )
    return a .+ b .* tan.( LinRange(-flim, flim, Int(round(coord_ext/steplength))+1) )
end

print("'borders_tan'   ")


@doc """
DESCRIPTION:
    - construct bin limits for an array with an arctan function for higher resolution in the outer rim
INPUT:
    - coord:        1D-array that is to be divided into a grid
    - steplength:   scale length for bin size
    - flim:         upper (and lower) limit of the arctan x-range, the larger the more outer-rim-focused
OUTPUT:
    - array of grid (bin) borders
""" ->
function borders_arctan(coord, steplength, flim)
    a           = median(coord)
    coord_ext   = maximum(coord)-minimum(coord)
    b           = coord_ext / ( 2*atan(flim) )
    return a .+ b .* atan.( LinRange(-flim, flim, Int(round(coord_ext/steplength))+1) )
end

print("'borders_arctan'   ")


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


@doc """
DESCRIPTION:\n
    - find halo ID according to Felix' notation based on the ID in the final snap shot
INPUT:\n
    - isub: subfind ID of the halo in the final snap shot
    - dir:  directory with all the processed Silvio-processed halo files
OUTPUT:\n
    - isub:             Int; halo ID according to Felix' notation
    - halo_filestring:  String; jld filename of that halo
""" ->
function find_felixID(rootID; dir=current_dir_stories)
    storyfilelist   = readdir(dir)
    halo_filestring = " "
    ifelix          = Int64
    for i in 1:length(storyfilelist)
        if occursin("_$(rootID).jld", storyfilelist[i])
            halo_filestring = storyfilelist[i]
            ifelix          = parse( Int64, chop(replace(halo_filestring, "_$(rootID)." => "_." ), head=5,tail=5) )
        end
    end
    return ifelix, halo_filestring
end

print("'find_felixID'   ")


@doc """
DESCRIPTION:\n
    - find sub ID in the final snap shot from to Felix' notation ID
INPUT:\n
    - ifelix:   halo ID from Felix
    - dir:      directory with all the processed Silvio-processed halo files
OUTPUT:\n
    - isub:             Int64; sub ID in the final snap shot
    - halo_filestring:  String; jld filename of that halo
""" ->
function find_rootID(felixID; dir=current_dir_stories)
    storyfilelist   = readdir(dir)
    halo_filestring = " "
    isub            = Int64
    for i in 1:length(storyfilelist)
        if occursin("halo_$(felixID)_", storyfilelist[i])
            halo_filestring = storyfilelist[i]
            isub            = parse( Int64, chop( replace( halo_filestring, "_$(felixID)_" => "__" ), head=6,tail=4 ) )
        end
    end
    return isub, halo_filestring
end

print("'find_ISUB'   ")


@doc """
DESCRIPTION:\n
    - Find same-snap First Progenitor index depending on merger index in Felix' list
INPUT:\n
    - mID:      merger index
    - halsto:   specific halostory dictionary
OUTPUT:\n
    - result_ID:    Int64; same-snap First Progenitor index
""" ->
function ssFPfinder(mID, halsto)
    result_ID = Int64
    for finder_ID in mID:-1:2    # walking the array towards the first progenitor
        if halsto["snapFP"][finder_ID] < halsto["snapFP"][finder_ID-1]
            result_ID   = finder_ID
            break
        end
    end
    return result_ID
end

print("'ssFPfinder'   ")


@doc """
DESCRIPTION:\n
    - Tranform Cartesian vector into [r,θ,ϕ] Spherical vector
INPUT:\n
    - x:    cartesian vector
OUTPUT:\n
    - s:    spherical vector in [r,θ,ϕ]
""" ->
function cartesian_to_spherical(x)
    s = zeros(3)
    s[1]    = norm(x) # Radius
    s[2]    = acosd(x[3] / s[1])# * 180 / π # θ[°]
    s[3]    = acosd(x[1] / sqrt(x[1]*x[1] + x[2]*x[2]))# * 180 / π # ϕ[°]
    return s
end

print("'cartesian_to_spherical'   ")


@doc """
DESCRIPTION:\n
    - Rotation Matrix to align x with reference vector
INPUT:\n
    - x:    cartesian vector
    - ref:  reference vector
OUTPUT:\n
    - rotmat:   rotation matrix to align x with ref
""" ->
function aligner(x; ref=[0,0,1])
    nx=normalize(x)
    nref=normalize(nref)
    vx      = nx × nref
    cx      = transpose(nx) * nref
    Vx      = [0 -vx[3] vx[2]; vx[3] 0 -vx[1]; -vx[2] vx[1] 0]
    return I + Vx + (Vx * Vx ./ (1+cx))
end

print("'aligner'   ")


function find_merging_progenitors(snapFP; felixID=" ", rootID=" ", subID=" ", 
    path_to_halostories=current_dir_stories, check_ambiguity=false)
    #subIDlist = Array{Int64}(undef, 0)
    if typeof(felixID)==Int  # simple case since tree is provided
        halo_story  = load(joinpath(path_to_halostories, "halo_$(felixID)_$(find_rootID(felixID; dir=path_to_halostories)[1]).jld"), "halo_story")
        #subIDlist   = halo_story["subID"][findall(x->x.==snapFP, halo_story["snapFP"])]
        return halo_story["subID"][findall(x->x.==snapFP, halo_story["snapFP"])]
    elseif typeof(rootID)==Int  # simple case since tree is provided
        halo_story  = load(joinpath(path_to_halostories, "halo_$(rootID).jld"), "halo_story")
        #subIDlist   = halo_story["subID"][findall(x->x.==snapFP, halo_story["snapFP"])]
        return halo_story["subID"][findall(x->x.==snapFP, halo_story["snapFP"])]
    elseif typeof(subID)==Int   # here we need to find the right tree first
        storyfilelist   = readdir(path_to_halostories)
        if check_ambiguity
            subIDlist = Array{Int64}(undef, 0)
            for i in 1:length(storyfilelist)
                print("$i ")
                halo_story  = load(joinpath(path_to_halostories, storyfilelist[i]), "halo_story")
                if length(findall(x->x.==subID, halo_story["subID"][findall(x->x.==snapFP, halo_story["snapFP"])])) > 0
                    if length(subIDlist) == 0
                        subIDlist = halo_story["subID"][findall(x->x.==snapFP, halo_story["snapFP"])]
                    else
                        error("Combination of subID and snapFP is ambiguous.")
                    end
                end
            end
            return subIDlist
        else
            for i in 1:length(storyfilelist)
                print("$i ")
                halo_story  = load(joinpath(path_to_halostories, storyfilelist[i]), "halo_story")
                if length(findall(x->x.==subID, halo_story["subID"][findall(x->x.==snapFP, halo_story["snapFP"])])) > 0
                    return halo_story["subID"][findall(x->x.==snapFP, halo_story["snapFP"])]
                    break
                end
            end
        end
    else
        error("Either felixID, subID or rootID has to be provided.")
    end
end

print("'find_merging_progenitors'   ")


@doc """
DESCRIPTION:\n
    - Scan merger map for indices with condition
INPUT:\n
    - asmbly:       dictionary
    -- condtition:  
OUTPUT:\n
    - dictionary
""" ->
function mergermap_indices(asmly; 
    condition="mergers"
    , min=-1, max=1e16
    , mtype=convert_parttype_to_idx(str="stars"), snap=136,
    )
    mm_i    = Dict{String, Vector{Union{Missing, Int64}}}()
    mm_i["main"]            = missings(Int64, 0)
    if condition == "mergers"
        mm_i["most_massive"]    = missings(Int64, 0)
        mm_i["last"]            = missings(Int64, 0)
        mm_i["first"]           = missings(Int64, 0)
        for i in 1:length(asmly["snapNR"])
            #if sum( asmly["merger map"][2, 1+sum(asmly["N_MERGERS"][1:i-1]) : sum(asmly["N_MERGERS"][1:i])] ) > 0. # is there a merger?
            if asmly["N_MERGERS"][i] > 0 # is there a merger?
                mm_i["most_massive"]    = vcat( mm_i["most_massive"], sum(asmly["N_MERGERS"][1:i-1]) + argmax( asmly["merger map"][mtype,1+sum(asmly["N_MERGERS"][1:i-1]):sum(asmly["N_MERGERS"][1:i])] ) )
                mm_i["main"] = vcat( mm_i["main"], i )
                mm_i["last"] = vcat( mm_i["last"], sum(asmly["N_MERGERS"][1:i]) )
                mm_i["first"] = vcat( mm_i["first"], sum(asmly["N_MERGERS"][1:i-1]) + 1 )
            end
        end
    elseif condition == "mass diff" 
        for i in 1:length(asmly["snapNR"])
            if !ismissing(asmly[string("Δ",ifelse(mtype==1,"M","M2"))][i]) && min < abs(asmly[string("Δ",ifelse(mtype==1,"M","M2"))][i])/asmly[ifelse(mtype==1,"M_felix","M2_felix")][i] < max # mass change
                mm_i["main"] = vcat( mm_i["main"], i )
            end
        end
    elseif condition == "snap" 
        for i in 1:length(asmly["snapNR"])
            if asmly["snapNR"][i] == snap
                mm_i["main"] = vcat( mm_i["main"], i )
            end
        end
    else
        error("Unknown contition")
    end
    return mm_i
end

print("'mergermap_indices'   ")


@doc """
DESCRIPTION:\n
    - Collect merger incidences with given ratio
    INPUT:\n
        - asmbly:       dictionary
    OUTPUT:\n
        - dictionary
""" ->
function merger_rates(asmly; 
    min=-1, max=1e10, 
    indices=" ", type="most massive",
    verbose=true
    )
    if typeof(indices) == String
        indices = mergermap_indices(asmly; condition="mergers")["main"]
    end

    incidences  = Dict{String, Vector{Union{Missing, Any}}}()
    incidences["subID"]     = missings(Int64, 0)
    incidences["rootID"]    = missings(Int64, 0)
    incidences["snapNR"]    = missings(Int64, 0)
    incidences["ratio"]     = missings(Float64, 0)
    if type == "most massive"
        for i in indices
            # select merger ratios
            temp_ratio  = asmly["merger map"][4,sum(asmly["N_MERGERS"][1:i-1])+argmax(asmly["merger map"][2,1+sum(asmly["N_MERGERS"][1:i-1]):sum(asmly["N_MERGERS"][1:i])])] / asmly["merger map"][2,sum(asmly["N_MERGERS"][1:i-1])+argmax(asmly["merger map"][2,1+sum(asmly["N_MERGERS"][1:i-1]):sum(asmly["N_MERGERS"][1:i])])]
            if min < temp_ratio < max
                incidences["subID"]     = vcat( incidences["subID"]     , asmly["subID"][i] )
                incidences["snapNR"]    = vcat( incidences["snapNR"]    , asmly["snapNR"][i] )
                incidences["rootID"]    = vcat( incidences["rootID"]    , asmly["ID_ISUB"][i] )
                incidences["ratio"]     = vcat( incidences["ratio"]     , temp_ratio )
                if verbose
                    println("Halo $(asmly["subID"][i])  =  $(asmly["ID_ISUB"][i]) @ Snap $(asmly["snapNR"][i])   ---   Ratio $temp_ratio")
                end
            end
        end
    else
        error("Unknown Type")
    end
    return incidences
end

print("'merger_rates'   ")


@doc """
DESCRIPTION:\n
    INPUT:\n
    OUTPUT:\n
""" ->
function find_halo_file(; indir=current_dir_stories
    , felixID=" ", rootID=" "
    )

    storyfilelist   = readdir(indir)
    halo_filestring = " "
    if typeof(felixID)==Int  # simple case since tree is provided
        for i in 1:length(storyfilelist)
            if occursin("halo_$(felixID)_", storyfilelist[i])
                println(storyfilelist[i])
                return storyfilelist[i]
            end
        end
    elseif typeof(rootID)==Int  # simple case since tree is provided
        for i in 1:length(storyfilelist)
            if occursin("_$(rootID).jld", storyfilelist[i])
                println(storyfilelist[i])
                return storyfilelist[i]
            end
        end
    else
        error("Either felixID or rootID has to be provided.")
    end
end

print("'find_halo_file'   ")







@doc """
DESCRIPTION:\n
    INPUT:\n
    OUTPUT:\n
""" ->
function pyplottable(arr; 
                    calc_norm=true, 
                    )
    if sum(size(arr)) == sum(size(arr[:]))+1 # 2d array thats actually 1d
        return replace(replace(replace(arr[:], missing => NaN), Inf => NaN), -Inf => NaN)
    elseif length(arr) > sum(size(arr)) && calc_norm # 2d array to be normalized
        return norm.( eachcol( replace(replace(replace(arr, missing => NaN), Inf => NaN), -Inf => NaN) ) )
    else
        return replace(replace(replace(arr, missing => NaN), Inf => NaN), -Inf => NaN)
    end
end

print("'pyplottable'   ")







@doc """
DESCRIPTION:\n
    INPUT:\n
    OUTPUT:\n
""" ->
function findcs(arr;
                    eq=" ", gt=" ", lt=" ", geq=" ", leq=" ",
                    calc_norm=true,
                    comparewith=" "
                    )
    if gt != " " && lt != " "
        found =  findall(x-> gt .< x .< lt, pyplottable(arr,calc_norm=calc_norm))
    elseif eq === missing
        found =  findall(x-> x .=== eq, arr)
    elseif eq === NaN
        found =  findall(x-> x .=== eq, arr)
    elseif eq != " "
        found =  findall(x-> x .== eq, pyplottable(arr,calc_norm=calc_norm))
    elseif geq != " " && leq != " "
        found =  findall(x-> geq .<= x .<= leq, pyplottable(arr,calc_norm=calc_norm))
    elseif geq != " " && lt != " "
        found =  findall(x-> geq .<= x .< lt, pyplottable(arr,calc_norm=calc_norm))
    elseif gt != " " && leq != " "
        found =  findall(x-> gt .< x .<= leq, pyplottable(arr,calc_norm=calc_norm))
    elseif gt != " "
        found =  findall(x-> x .> gt, pyplottable(arr,calc_norm=calc_norm))
    elseif lt != " "
        found =  findall(x-> x .< lt, pyplottable(arr,calc_norm=calc_norm))
    elseif geq != " "
        found =  findall(x-> x .>= geq, pyplottable(arr,calc_norm=calc_norm))
    elseif leq != " "
        found =  findall(x-> x .<= leq, pyplottable(arr,calc_norm=calc_norm))
    else
        found = findall(x-> x .== 1, pyplottable(arr,calc_norm=calc_norm))
    end

    if comparewith == " "
        return found
    else
        return found[found .∈ Ref(Set(comparewith))]
    end
end

print("'findcs'   ")







@doc """
DESCRIPTION:\n
    INPUT:\n
    OUTPUT:\n
""" ->
function idxmatch(arr1, arr2)
    return arr1[ arr1 .∈ Ref(Set(arr2 ))]
end

print("'idxmatch'   ")

@doc """
DESCRIPTION:\n
    INPUT:\n
    OUTPUT:\n
""" ->
function idxclude(arr1, arr2)
    return arr1[ arr1 .∉ Ref(Set(arr2 ))]
end

print("'idxclude'   ")







@doc """
DESCRIPTION:\n
    INPUT:\n
    OUTPUT:\n
""" ->
function pctl68upper(x)
    return StatsBase.percentile(x,168.27/2)
end

print("'pctl68upper'   ")

@doc """
DESCRIPTION:\n
    INPUT:\n
    OUTPUT:\n
""" ->
function pctl68lower(x)
    return StatsBase.percentile(x,(100-68.27)/2)
end

print("'pctl68lower'   ")







@doc """
DESCRIPTION:\n 
        - Use as @Name(arg)
    INPUT:\n
    OUTPUT:\n
""" ->
macro Name(arg)
    string(arg)
 end

print("'Name'   ")






@doc """
DESCRIPTION:\n
    INPUT:\n
    OUTPUT:\n
""" ->
function convert_parttype_to_idx(;str=" ", sym=:default)
    if lowercase(str) == "stars" || sym == :stars
        return 5
    elseif lowercase(str) == "gas" || sym == :gas
        return 1
    elseif lowercase(str) == "dm" || sym == :dm
        return 2
    else
        error("Unknown particle type: $str or $sym ")
    end
end

print("'convert_parttype_to_idx'   ")







@doc """
DESCRIPTION:\n sfc: WIP
    INPUT:\n
    OUTPUT:\n
""" ->
function running_stats(x,y;
                        resolution=50, width=(maximum(x)-minimum(x))/resolution
                        )
                        
end

print("'running_stats'   ")







@doc """
DESCRIPTION:\n sfc: WIP
    INPUT:\n
    OUTPUT:\n
""" ->
function export_halostories_to_csv(;
    story_dir = current_dir_stories, outdir="./TEMP_csv",
    rootID_list = " ", 
    )

    # Select Halos to export
    storyfilelist   = readdir(story_dir)
    if rootID_list == " "
        list_idx    = Vector{Int64}(undef, 0)
        for i in 1:length(storyfilelist)
            if occursin(".jld", storyfilelist[i])
                list_idx    = vcat( list_idx,   i   )
            end
        end
        storyfilelist   = storyfilelist[list_idx]
    else
        list_idx    = Vector{Int64}(undef, 0)
        for i in 1:length(storyfilelist)
            for ii in 1:length(rootID_list)
                #if occursin("$(rootID_list[ii]).jld", storyfilelist[i])
                if "halo_$(rootID_list[ii]).jld" == storyfilelist[i]
                    list_idx    = vcat( list_idx,   i   )
                end
            end
        end
        storyfilelist   = storyfilelist[list_idx]
    end

    # Exort loop
    for i in storyfilelist
        println("   Writing $i")
        flush(stdout)
        halo_story = load(joinpath(story_dir, i), "halo_story")
        #CSV.write( joinpath(outdir, replace(i, ".jld" => ".csv" )), CSV.RowWriter(halo_story, delim=" ") )
        #fid    = XLSX.openxlsx( joinpath(outdir, replace(i, ".jld" => ".xlsx" )), mode="w" )
        #header = ["keys";"values"]
        #data   = [[collect(keys(halo_story))];[collect(values(halo_story))]]
        #XLSX.writetable(joinpath(outdir, replace(i, ".jld" => ".xlsx" )), data,header)

        # Create DataFrame
        df_halo_story
        for ii in keys(halo_story)
            
        end
    end

    println("\n\n\n---------------------------\n\nNow witness the firepower of this fully armed and operational battle station!\n---------------------------\n\n")
    return nothing
end

print("'export_halostories_to_csv'   ")






@doc """
DESCRIPTION:\n
    INPUT:\n
    OUTPUT:\n
""" ->
function find_mintime_idx(  lbt_array   , idx_array , idx   , mintime   ,
    )
    for i in idx:-1:1
        if lbt_array[idx_array[i]] - lbt_array[idx_array[idx]] > mintime
            #println("found $(idx_array[i]) with $(lbt_array[idx_array[i]] - lbt_array[idx_array[idx]]) Gyr")
            return idx_array[i]
            break
        end
        if i == 1
            #println("nothing found")
            return idx_array[idx]
        end
    end
end


print("'find_mintime_idx'   ")






@doc """
DESCRIPTION:\n
- only defined for abs(x) ≥ 1
    INPUT:\n
    OUTPUT:\n
""" ->
function symlog10(x)
    if !ismissing(x) && !isnan(x)
        if x ≥ 1
            return log10(x)
        elseif x ≤ -1
            return -log10(abs(x))
        else
            return missing
        end
    else
        return missing
    end
end


print("'symlog10'   ")


println()