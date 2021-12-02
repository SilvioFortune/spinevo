
@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_group_subID( snapNR, subID; 
    simbox="/HydroSims/Magneticum/Box4/uhr_test", res=(1600,900), 
    property = "pos",
    rad=200, ptsize=3000, arsize=4., port=1688, stepsize=1)
    JSServe.configure_server!(listen_port=port, forwarded_port=port)


    ### Data processing
    snapshot    = Snapshot(simbox, snapNR)
    g = Galaxy(snapshot, subID)
    rvir_group  = read_galaxy_prop(get_group(g), "RVIR", :physical)
    #sph_small   = GadgetGalaxies.Sphere(0.1*rvir_group)
    #sph_large   = GadgetGalaxies.Sphere(rvir_group)
    read_halo!(g, units=:physical, props=((:stars, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=0.1*rvir_group)
    g.stars.pos .+= read_galaxy_pos(g, :physical)
    particleID_list         = vcat(particleID_list, g.stars.id)
    #read_halo!(g, units=:physical, props=((:gas, ["POS", "VEL", "MASS"]),), radius_units=:physical, radius=0.1*rvir_group)
    #read_halo!(g, units=:physical, props=((:dm, ["POS", "VEL"]),), radius_units=:physical, radius=rvir_group)
    
    # all particles
    head        = read_header("$simbox/groups_$(@sprintf("%03i", snapNR))/sub_$(@sprintf("%03i", snapNR))")
    filepath    = "$simbox/snapdir_$(@sprintf("%03i", snapNR))/snap_$(@sprintf("%03i", snapNR))"
    blocks      = ["MASS", "POS", "ID", "VEL"]
    radius      = rad * ( read_galaxy_pos(g, :sim) ./ read_galaxy_pos(g, :physical) )[1] # conversion into simulation units, mean is for 
    position    = read_galaxy_pos(g, :sim)
    ALL_STARS   = read_particles_in_volume(filepath, blocks, position, radius; parttype=4, verbose=true, use_keys=true)
    #ALL_STARS["POS"]    .-= position
    ALL_STARS["MASS"]   = convert_units_physical(ALL_STARS["MASS"], :mass, head)
    ALL_STARS["POS"]    = convert_units_physical(ALL_STARS["POS"], :pos, head)
    ALL_STARS["VEL"]    = convert_units_physical(ALL_STARS["VEL"], :vel, head)
    # Crop duplicates with target
    all_notin   = ALL_STARS["ID"] .∉ Ref(Set(particleID_list))
    size_factor = ptsize / maximum(g.stars.mass[1:stepsize:end])

    #Plotting
    if property == "pos"
        set_theme!(resolution=res, backgroundcolor = :black)
        #println(typeof(ALL_STARS["POS"][:,all_notin][:,1:stepsize:end]))
        scene = Scene()
        scatter!( scene,
            g.stars.pos[:,1:stepsize:end],
            color = g.stars.mass[1:stepsize:end] .* size_factor,
            markersize = g.stars.mass[1:stepsize:end] .* size_factor,#markersize = 0.1,
            colormap = :winter,
            opacity = true,
            transparency = true
            )
            
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
        arrows!(scene, 
            g.stars.pos[1,1:stepsize:end], 
            g.stars.pos[2,1:stepsize:end], 
            g.stars.pos[3,1:stepsize:end],
            galaxi499es[i].stars.vel[1,1:stepsize:end], 
            g.stars.vel[2,1:stepsize:end], 
            g.stars.vel[3,1:stepsize:end],
            # convert kms/s to kpc/10Myr by 0.010227047347441423
            lengthscale = 0.010227047347441423,#arsize*1e-2,#arsize * (maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]),
            arrowsize   = ( arsize / maximum(g.stars.mass[1:stepsize:end]) ) .* g.stars.mass[1:stepsize:end],
            #arrowsize   = arsize,# (arsize*(maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["MASS"][all_notin][1:stepsize:end])) .* g.stars.mass[1:stepsize:end],
            linewidth   = 0.1,
            linecolor   = g.stars.mass[1:stepsize:end] .* norm.(eachcol(g.stars.vel[:,1:stepsize:end])) ./ maximum(g.stars.mass[1:stepsize:end] .* norm.(eachcol(g.stars.vel[:,1:stepsize:end]))),
            arrowcolor  = g.stars.mass[1:stepsize:end] .* norm.(eachcol(g.stars.vel[:,1:stepsize:end])) ./ maximum(g.stars.mass[1:stepsize:end] .* norm.(eachcol(g.stars.vel[:,1:stepsize:end]))),
            colormap    = :winter,
            align       = :head,
            quality     = 3,
            shading     = false,
            ssao        = false,
            normalize   = false,
            transparency = true
        )
        arrows!(scene, 
            ALL_STARS["POS"][1,all_notin][1:stepsize:end], 
            ALL_STARS["POS"][2,all_notin][1:stepsize:end], 
            ALL_STARS["POS"][3,all_notin][1:stepsize:end],
            ALL_STARS["VEL"][1,all_notin][1:stepsize:end], 
            ALL_STARS["VEL"][2,all_notin][1:stepsize:end], 
            ALL_STARS["VEL"][3,all_notin][1:stepsize:end],
            # convert kms/s to kpc/10Myr by 0.010227047347441423
            lengthscale = 0.010227047347441423,#arsize*1e-2,#arsize * (maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]),
            arrowsize   = (arsize/maximum(g.stars.mass[1:stepsize:end])) .* ALL_STARS["MASS"][all_notin][1:stepsize:end],
            linewidth   = 0.1,
            linecolor   = ALL_STARS["MASS"][all_notin][1:stepsize:end] .* norm.(eachcol(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end])) ./ maximum(g.stars.mass[1:stepsize:end] .* norm.(eachcol(g.stars.vel[:,1:stepsize:end]))),
            arrowcolor  = ALL_STARS["MASS"][all_notin][1:stepsize:end] .* norm.(eachcol(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end])) ./ maximum(g.stars.mass[1:stepsize:end] .* norm.(eachcol(g.stars.vel[:,1:stepsize:end]))),
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

print("'plot_group_subID'   ")


@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_pos( snapNR, position; 
    simbox="/HydroSims/Magneticum/Box4/uhr_test", res=(1600,900), 
    property = "pos",
    rad=200, ptsize=3000, arsize=4., port=1688, stepsize=1)
    JSServe.configure_server!(listen_port=port, forwarded_port=port)


    # all particles
    head        = read_header("$simbox/groups_$(@sprintf("%03i", snapNR))/sub_$(@sprintf("%03i", snapNR))")
    filepath    = "$simbox/snapdir_$(@sprintf("%03i", snapNR))/snap_$(@sprintf("%03i", snapNR))"
    blocks      = ["MASS", "POS", "ID", "VEL"]
    radius      = rad
    ALL_STARS   = read_particles_in_volume(filepath, blocks, position, radius; parttype=4, verbose=true, use_keys=true)
    #ALL_STARS["POS"]    .-= position
    ALL_STARS["MASS"]   = convert_units_physical(ALL_STARS["MASS"], :mass, head)
    ALL_STARS["POS"]    = convert_units_physical(ALL_STARS["POS"], :pos, head)
    ALL_STARS["VEL"]    = convert_units_physical(ALL_STARS["VEL"], :vel, head)


    size_factor = ptsize / maximum(g.stars.mass[1:stepsize:end])

    #Plotting
    if property == "pos"
        set_theme!(resolution=res, backgroundcolor = :black)
        #println(typeof(ALL_STARS["POS"][:,:][:,1:stepsize:end]))
        scene = Scene()
            
        scatter!( scene,
            ALL_STARS["POS"][:,:][:,1:stepsize:end],
            color = ALL_STARS["MASS"][:][1:stepsize:end] .* size_factor,
            markersize = ALL_STARS["MASS"][:][1:stepsize:end] .* size_factor,
            colormap = :autumn1,
            opacity = true,
            transparency = true
            )
    elseif property == "vel"
        set_theme!(resolution=res, backgroundcolor = :white)
        scene = Scene()
        arrows!(scene, 
            ALL_STARS["POS"][1,:][1:stepsize:end], 
            ALL_STARS["POS"][2,:][1:stepsize:end], 
            ALL_STARS["POS"][3,:][1:stepsize:end],
            ALL_STARS["VEL"][1,:][1:stepsize:end], 
            ALL_STARS["VEL"][2,:][1:stepsize:end], 
            ALL_STARS["VEL"][3,:][1:stepsize:end],
            # convert kms/s to kpc/10Myr by 0.010227047347441423
            lengthscale = 0.010227047347441423,#arsize*1e-2,#arsize * (maximum(ALL_STARS["POS"][1,:][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,:][1,1:stepsize:end])) / maximum(ALL_STARS["VEL"][:,:][:,1:stepsize:end]),
            arrowsize   = (arsize/maximum(g.stars.mass[1:stepsize:end])) .* ALL_STARS["MASS"][:][1:stepsize:end],
            linewidth   = 0.1,
            linecolor   = ALL_STARS["MASS"][:][1:stepsize:end] .* norm.(eachcol(ALL_STARS["VEL"][:,:][:,1:stepsize:end])) ./ maximum(g.stars.mass[1:stepsize:end] .* norm.(eachcol(g.stars.vel[:,1:stepsize:end]))),
            arrowcolor  = ALL_STARS["MASS"][:][1:stepsize:end] .* norm.(eachcol(ALL_STARS["VEL"][:,:][:,1:stepsize:end])) ./ maximum(g.stars.mass[1:stepsize:end] .* norm.(eachcol(g.stars.vel[:,1:stepsize:end]))),
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

print("'plot_pos'   ")