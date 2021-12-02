println("   Plotting...")
using PyPlot
using ColorSchemes
using JSServe
using WGLMakie
WGLMakie.activate!()


@doc """
DESCRIPTION:\n
    - Use WGLMakie to plot a galaxy, its mergers and its surroundings
INPUT:\n
    - snapNR:       snapshot number
    -- subID:       subhalo ID, dont ever use this one!
    -- lastID:      subID in the last snapshot
    -- felixID:     felix ID to quickly find the right group
    -- simbox:         path to simulation box
    -- res:         plot resolution (m, n)
    -- rad:         radius for surrounding medium from GadgetIO in kpc
    -- ptsize:      size of scatterplot particles
    -- arsize:      size of arrow heads
    -- property:    "pos" or "vel"
    -- port:        ports used on local and remote
    -- stepsize:    plot only every nth particle when there are too many
OUTPUT:\n
    - scene:    plot
""" ->
function plot_group( snapNR; 
    simbox="/HydroSims/Magneticum/Box4/uhr_test", res=(1600,900), felixID="", subID="", lastID="", 
    property = "pos",
    rad=200, ptsize=3000, arsize=4., port=1688, stepsize=1)
    JSServe.configure_server!(listen_port=port, forwarded_port=port)

    if typeof(felixID)==Int  # simple case since tree is provided
        galaxyID_list   = find_merging_progenitors(snapNR; felixID=felixID)
    elseif typeof(lastID)==Int  # simple case since tree is provided
        galaxyID_list   = find_merging_progenitors(snapNR; lastID=lastID)
    elseif typeof(subID)==Int   # here we need to find the right tree first
        galaxyID_list   = find_merging_progenitors(snapNR; subID=subID)
    else
        error("Either felixID, subID or lastID has to be provided.")
    end


    ### Data processing
    snapshot    = Snapshot(simbox, snapNR)
    galaxies        = Dict{Int64, Any}()
    particleID_list = Array{UInt64}(undef,0)
    for i in 1:length(galaxyID_list)
        galaxies[i] = Galaxy(snapshot, galaxyID_list[i])
        rvir_group  = read_galaxy_prop(get_group(galaxies[1]), "RVIR", :physical)
        #sph_small   = GadgetGalaxies.Sphere(0.1*rvir_group)
        #sph_large   = GadgetGalaxies.Sphere(rvir_group)
        read_halo!(galaxies[i], units=:physical, props=((:stars, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=0.1*rvir_group)
        galaxies[i].stars.pos .+= read_galaxy_pos(galaxies[i], :physical)
        particleID_list         = vcat(particleID_list, galaxies[i].stars.id)
        #read_halo!(galaxies[i], units=:physical, props=((:gas, ["POS", "VEL", "MASS"]),), radius_units=:physical, radius=0.1*rvir_group)
        #read_halo!(galaxies[i], units=:physical, props=((:dm, ["POS", "VEL"]),), radius_units=:physical, radius=rvir_group)
    end
    # all particles
    head        = read_header("$simbox/groups_$(@sprintf("%03i", snapNR))/sub_$(@sprintf("%03i", snapNR))")
    filepath    = "$simbox/snapdir_$(@sprintf("%03i", snapNR))/snap_$(@sprintf("%03i", snapNR))"
    blocks      = ["MASS", "POS", "ID", "VEL"]
    radius      = rad * ( read_galaxy_pos(galaxies[1], :sim) ./ read_galaxy_pos(galaxies[1], :physical) )[1] # conversion into simulation units, mean is for 
    position    = read_galaxy_pos(galaxies[1], :sim)
    ALL_STARS   = read_particles_in_volume(filepath, blocks, position, radius; parttype=4, verbose=true, use_keys=true)
    #ALL_STARS["POS"]    .-= position
    ALL_STARS["MASS"]   = convert_units_physical(ALL_STARS["MASS"], :mass, head)
    ALL_STARS["POS"]    = convert_units_physical(ALL_STARS["POS"], :pos, head)
    ALL_STARS["VEL"]    = convert_units_physical(ALL_STARS["VEL"], :vel, head)
    # Crop duplicates with target
    all_notin   = ALL_STARS["ID"] .∉ Ref(Set(particleID_list))
    size_factor = ptsize / maximum(galaxies[1].stars.mass[1:stepsize:end])

    #Plotting
    if property == "pos"
        set_theme!(resolution=res, backgroundcolor = :black)
        #println(typeof(ALL_STARS["POS"][:,all_notin][:,1:stepsize:end]))
        scene = Scene()
        # Halos
        for i in 1:length(galaxyID_list)
            if i == 1
                scatter!( scene,
                    galaxies[i].stars.pos[:,1:stepsize:end],
                    color = galaxies[i].stars.mass[1:stepsize:end] .* size_factor,
                    markersize = galaxies[i].stars.mass[1:stepsize:end] .* size_factor,#markersize = 0.1,
                    colormap = :winter,
                    opacity = true,
                    transparency = true
                    )
            else
                scatter!( scene,
                    galaxies[i].stars.pos[:,1:stepsize:end],
                    color = galaxies[i].stars.mass[1:stepsize:end] .* size_factor,
                    markersize = galaxies[i].stars.mass[1:stepsize:end] .* size_factor,#markersize = 0.1,
                    colormap = :summer,
                    opacity = true,
                    transparency = true
                    )
            end
        end
        scatter!( scene,
            ALL_STARS["POS"][:,all_notin][:,1:stepsize:end],
            color = ALL_STARS["MASS"][all_notin][1:stepsize:end] .* size_factor,
            markersize = ALL_STARS["MASS"][all_notin][1:stepsize:end] .* size_factor,
            colormap = :autumn1,
            opacity = true,
            transparency = true
            )
    elseif property == "vel"
        set_theme!(resolution=res, backgroundcolor = :white)
        scene = Scene()
        for i in 1:length(galaxyID_list)
            if i == 1
                arrows!(scene, 
                    galaxies[i].stars.pos[1,1:stepsize:end], 
                    galaxies[i].stars.pos[2,1:stepsize:end], 
                    galaxies[i].stars.pos[3,1:stepsize:end],
                    galaxi499es[i].stars.vel[1,1:stepsize:end], 
                    galaxies[i].stars.vel[2,1:stepsize:end], 
                    galaxies[i].stars.vel[3,1:stepsize:end],
                    # convert kms/s to kpc/10Myr by 0.010227047347441423
                    lengthscale = 0.010227047347441423,#arsize*1e-2,#arsize * (maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]),
                    arrowsize   = ( arsize / maximum(galaxies[1].stars.mass[1:stepsize:end]) ) .* galaxies[i].stars.mass[1:stepsize:end],
                    #arrowsize   = arsize,# (arsize*(maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["MASS"][all_notin][1:stepsize:end])) .* galaxies[i].stars.mass[1:stepsize:end],
                    linewidth   = 0.1,
                    linecolor   = galaxies[i].stars.mass[1:stepsize:end] .* norm.(eachcol(galaxies[i].stars.vel[:,1:stepsize:end])) ./ maximum(galaxies[1].stars.mass[1:stepsize:end] .* norm.(eachcol(galaxies[1].stars.vel[:,1:stepsize:end]))),
                    arrowcolor  = galaxies[i].stars.mass[1:stepsize:end] .* norm.(eachcol(galaxies[i].stars.vel[:,1:stepsize:end])) ./ maximum(galaxies[1].stars.mass[1:stepsize:end] .* norm.(eachcol(galaxies[1].stars.vel[:,1:stepsize:end]))),
                    colormap    = :winter,
                    align       = :head,
                    quality     = 3,
                    shading     = false,
                    ssao        = false,
                    normalize   = false,
                    transparency = true
                )
            else
                arrows!(scene, 
                    galaxies[i].stars.pos[1,1:stepsize:end], 
                    galaxies[i].stars.pos[2,1:stepsize:end], 
                    galaxies[i].stars.pos[3,1:stepsize:end],
                    galaxies[i].stars.vel[1,1:stepsize:end], 
                    galaxies[i].stars.vel[2,1:stepsize:end], 
                    galaxies[i].stars.vel[3,1:stepsize:end],
                    # convert kms/s to kpc/10Myr by 0.010227047347441423
                    lengthscale = 0.010227047347441423,#arsize*1e-2,#arsize * (maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]),
                    arrowsize   = ( arsize / maximum(galaxies[1].stars.mass[1:stepsize:end]) ) .* galaxies[i].stars.mass[1:stepsize:end],
                    #arrowsize   = arsize,# (arsize*(maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["MASS"][all_notin][1:stepsize:end])) .* galaxies[i].stars.mass[1:stepsize:end],
                    linewidth   = 0.1,
                    linecolor   = galaxies[i].stars.mass[1:stepsize:end] .* norm.(eachcol(galaxies[i].stars.vel[:,1:stepsize:end])) ./ maximum(galaxies[1].stars.mass[1:stepsize:end] .* norm.(eachcol(galaxies[1].stars.vel[:,1:stepsize:end]))),
                    arrowcolor  = galaxies[i].stars.mass[1:stepsize:end] .* norm.(eachcol(galaxies[i].stars.vel[:,1:stepsize:end])) ./ maximum(galaxies[1].stars.mass[1:stepsize:end] .* norm.(eachcol(galaxies[1].stars.vel[:,1:stepsize:end]))),
                    colormap    = :summer,
                    align       = :head,
                    quality     = 3,
                    shading     = false,
                    ssao        = false,
                    normalize   = false,
                    transparency = true
                )
            end
        end
        arrows!(scene, 
            ALL_STARS["POS"][1,all_notin][1:stepsize:end], 
            ALL_STARS["POS"][2,all_notin][1:stepsize:end], 
            ALL_STARS["POS"][3,all_notin][1:stepsize:end],
            ALL_STARS["VEL"][1,all_notin][1:stepsize:end], 
            ALL_STARS["VEL"][2,all_notin][1:stepsize:end], 
            ALL_STARS["VEL"][3,all_notin][1:stepsize:end],
            # convert kms/s to kpc/10Myr by 0.010227047347441423
            lengthscale = 0.010227047347441423,#arsize*1e-2,#arsize * (maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]),
            arrowsize   = (arsize/maximum(galaxies[1].stars.mass[1:stepsize:end])) .* ALL_STARS["MASS"][all_notin][1:stepsize:end],
            linewidth   = 0.1,
            linecolor   = ALL_STARS["MASS"][all_notin][1:stepsize:end] .* norm.(eachcol(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end])) ./ maximum(galaxies[1].stars.mass[1:stepsize:end] .* norm.(eachcol(galaxies[1].stars.vel[:,1:stepsize:end]))),
            arrowcolor  = ALL_STARS["MASS"][all_notin][1:stepsize:end] .* norm.(eachcol(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end])) ./ maximum(galaxies[1].stars.mass[1:stepsize:end] .* norm.(eachcol(galaxies[1].stars.vel[:,1:stepsize:end]))),
            colormap    = :autumn1,
            align       = :head,
            quality     = 3,
            shading     = false,
            ssao        = false,
            normalize   = false,
            transparency = true
        )
    end
    
    return scene
end

print("'plot_group'   ")



include("/home/moon/sfortune/spinevo/pkg/plot_group_direct.jl")



println()