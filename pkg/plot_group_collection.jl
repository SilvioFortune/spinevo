


function testmakie(; port=1688, res=(1600,900), bgcolor=:white)
    JSServe.configure_server!(listen_port=port, forwarded_port=port)

    #Plotting
    set_theme!(resolution=res, backgroundcolor = bgcolor)
    pos = [1.1,1.1,1.1]
    vel = [3.1,3.1,3.1]
    ps = [Point3f(x, y, z) for (x,y,z) in zip(-5:2:5,-5:2:5,-5:2:5)]# for y in -5:2:5 for z in -5:2:5]
    ns = map(p -> 0.1 * Vec3f(p[2], p[3], p[1]), ps)
    lengths = norm.(ns)
    scene = WGLMakie.arrows(
        ps, ns, fxaa=true, # turn on anti-aliasing
        color=lengths,
        linewidth = 0.1, arrowsize = Vec3f(0.3, 0.3, 0.4),
        align = :center, axis=(type=Axis3,)
    )
    #scene = WGLMakie.arrows(
    ##arrows!(ax, 
    #    1,
    #    # convert kms/s to 10 kpc/Myr by 0.010227047347441423
    #    #lengthscale = 0.010227047347441423,#arsize*1e-2,#arsize * (maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]),
    #    #arrowsize   = (arsize/maximum(galaxies[1].stars.mass[1:stepsize:end])) .* ALL_STARS["MASS"][all_notin][1:stepsize:end],
    #    #linewidth   = 0.1,
    #    #linecolor   = @view(ALL_STARS["MASS"][all_notin][1:stepsize:end]) .* norm.(eachcol(@view(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]))) ./ maximum(@view(galaxies[1].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[1].stars.vel[:,1:stepsize:end])))),
    #    #arrowcolor  = @view(ALL_STARS["MASS"][all_notin][1:stepsize:end]) .* norm.(eachcol(@view(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]))) ./ maximum(@view(galaxies[1].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[1].stars.vel[:,1:stepsize:end])))),
    #    #colormap    = :autumn1,
    #    align       = :head,
    #    quality     = 3,
    #    shading     = false,
    #    ssao        = false,
    #    normalize   = false,
    #    transparency = true,
    #)
    return scene #fig
end

@doc """
DESCRIPTION:\n
    - Use WGLMakie to plot a galaxy, its mergers and its surroundings
INPUT:\n
    - snapNR:       snapshot number
    -- subID:       subhalo ID, dont ever use this one!
    -- rootID:      subID in the last snapshot
    -- felixID:     felix ID to quickly find the right group
    -- simbox:      path to simulation box
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
    simbox=current_dir_simbox, res=(1600,900), felixID=" ", subID=" ", rootID=" ", 
    property = "pos", phkeys=true, 
    rad=200, ptsize=3000, arsize=4., port=1688, stepsize=1,
    bgcolor=:white
    , indir=current_dir_stories, spin=true, boxfix=false, file = " "
    )
    JSServe.configure_server!(listen_port=port, forwarded_port=port)

    scene = 1
    if typeof(felixID)==Int  # simple case since tree is provided
        galaxyID_list   = find_merging_progenitors(snapNR; felixID=felixID)
        println(galaxyID_list)
    elseif typeof(rootID)==Int  # simple case since tree is provided
        galaxyID_list   = find_merging_progenitors(snapNR; rootID=rootID)
        println(galaxyID_list)
    elseif typeof(subID)==Int   # here we need to find the right tree first
        galaxyID_list   = find_merging_progenitors(snapNR; subID=subID)
        println(galaxyID_list)
    else
        error("Either felixID, subID or rootID has to be provided.")
    end


    #### Data processing
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
    ALL_STARS   = read_particles_in_volume(filepath, blocks, position, radius; parttype=4, verbose=true, use_keys=phkeys)
    #ALL_STARS["POS"]    .-= position
    ALL_STARS["MASS"]   = convert_units_physical(ALL_STARS["MASS"], :mass, head)
    ALL_STARS["POS"]    = convert_units_physical(ALL_STARS["POS"], :pos, head)
    ALL_STARS["VEL"]    = convert_units_physical(ALL_STARS["VEL"], :vel, head)
    # Crop duplicates with target
    all_notin   = ALL_STARS["ID"] .∉ Ref(Set(particleID_list))
    size_factor = ptsize / maximum(@view(galaxies[1].stars.mass[1:stepsize:end]))

    #Plotting
    set_theme!(resolution=res, backgroundcolor = bgcolor)
    if property == "pos"
        #println(typeof(ALL_STARS["POS"][:,all_notin][:,1:stepsize:end]))
        #scene = Scene()
        #fig = WGLMakie.Figure()
        #ax = WGLMakie.Axis3(fig[1,1], viewmode=:fit)


        #scatter( #ax,
        scene = WGLMakie.scatter(
            @view(ALL_STARS["POS"][:,all_notin][:,1:stepsize:end]),
            color = @view(ALL_STARS["MASS"][all_notin][1:stepsize:end]) .* size_factor,
            markersize = @view(ALL_STARS["MASS"][all_notin][1:stepsize:end]) .* size_factor,
            colormap = :autumn1,
            opacity = true,
            transparency = true,
            )
        # Halos
        for i in 1:length(galaxyID_list)    # mergers halo in front
        #for i in length(galaxyID_list):-1:1 # main in front
            if i == 1
                #scatter( #ax,
                scatter!(
                    @view(galaxies[i].stars.pos[:,1:stepsize:end]),
                    color = @view(galaxies[i].stars.mass[1:stepsize:end]) .* size_factor,
                    markersize = @view(galaxies[i].stars.mass[1:stepsize:end]) .* size_factor,#markersize = 0.1,
                    colormap = :winter,
                    opacity = true,
                    transparency = true,
                    )
            elseif !boxfix
                #scatter( #ax,
                scatter!(
                    @view(galaxies[i].stars.pos[:,1:stepsize:end]),
                    color = @view(galaxies[i].stars.mass[1:stepsize:end]) .* size_factor,
                    markersize = @view(galaxies[i].stars.mass[1:stepsize:end]) .* size_factor,#markersize = 0.1,
                    colormap = :summer,
                    opacity = true,
                    transparency = true,
                    )
            end
        end
    elseif property == "vel"
        set_theme!(resolution=res, backgroundcolor = bgcolor)
        #print(size(ALL_STARS["POS"][:,all_notin][:,1:stepsize:end])) 
        pts = [Point3f(p1,p2,p3) for (p1,p2,p3) in zip(ALL_STARS["POS"][1,all_notin][1:stepsize:end], ALL_STARS["POS"][2,all_notin][1:stepsize:end], ALL_STARS["POS"][3,all_notin][1:stepsize:end])]
        vcs = [Vec3f(v1,v2,v3) for (v1,v2,v3) in zip(ALL_STARS["VEL"][1,all_notin][1:stepsize:end], ALL_STARS["VEL"][2,all_notin][1:stepsize:end], ALL_STARS["VEL"][3,all_notin][1:stepsize:end])]
        #masssize = [Vec3f(m1,m2,m3) for (m1,m2,m3) in zip(log10.((arsize/maximum(galaxies[1].stars.mass[1:stepsize:end])) .* ALL_STARS["MASS"][all_notin][1:stepsize:end]), log10.((arsize/maximum(galaxies[1].stars.mass[1:stepsize:end])) .* ALL_STARS["MASS"][all_notin][1:stepsize:end]), log10.((arsize/maximum(galaxies[1].stars.mass[1:stepsize:end])) .* ALL_STARS["MASS"][all_notin][1:stepsize:end]))]
        scene = WGLMakie.arrows(
        #arrows!(ax, 
            pts,#Point3f(1.1,1.2,1.3),#[Point3f(p1,p2,p3) for (p1,p2,p3) in zip(ALL_STARS["POS"][1,all_notin][1:stepsize:end], ALL_STARS["POS"][2,all_notin][1:stepsize:end], ALL_STARS["POS"][3,all_notin][1:stepsize:end])],
            vcs,#Point3f(2.1,2.2,2.3),#[Point3f(v1,v2,v3) for (v1,v2,v3) in zip(ALL_STARS["VEL"][1,all_notin][1:stepsize:end], ALL_STARS["VEL"][2,all_notin][1:stepsize:end], ALL_STARS["VEL"][3,all_notin][1:stepsize:end])],
            # convert kms/s to 10 kpc/Myr by 0.010227047347441423
            lengthscale = 0.010227047347441423,#arsize*1e-2,#arsize * (maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]),
            #arrowsize   = masssize, #Vec3f.((arsize/maximum(galaxies[1].stars.mass[1:stepsize:end])) .* ALL_STARS["MASS"][all_notin][1:stepsize:end],(arsize/maximum(galaxies[1].stars.mass[1:stepsize:end])) .* ALL_STARS["MASS"][all_notin][1:stepsize:end],(arsize/maximum(galaxies[1].stars.mass[1:stepsize:end])) .* ALL_STARS["MASS"][all_notin][1:stepsize:end]),#(arsize/maximum(galaxies[1].stars.mass[1:stepsize:end])) .* ALL_STARS["MASS"][all_notin][1:stepsize:end],
            linewidth   = 0.1,
            linecolor   = @view(ALL_STARS["MASS"][all_notin][1:stepsize:end]) .* norm.(eachcol(@view(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]))) ./ maximum(@view(galaxies[1].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[1].stars.vel[:,1:stepsize:end])))),
            arrowcolor  = @view(ALL_STARS["MASS"][all_notin][1:stepsize:end]) .* norm.(eachcol(@view(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]))) ./ maximum(@view(galaxies[1].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[1].stars.vel[:,1:stepsize:end])))),
            colormap    = :autumn1,
            align       = :head,
            quality     = 3,
            shading     = false,
            ssao        = false,
            normalize   = false,
            transparency = true,
        )
        #scene = Scene()
        #fig = WGLMakie.Figure()
        #ax = WGLMakie.Axis3(fig[1,1], viewmode=:fit)
        for i in 1:length(galaxyID_list)    # mergers halo in front
        #for i in length(galaxyID_list):-1:1 # main in front
            if i == 1
                pts = [Point3f(p1,p2,p3) for (p1,p2,p3) in zip(galaxies[i].stars.pos[1,1:stepsize:end], galaxies[i].stars.pos[2,1:stepsize:end], galaxies[i].stars.pos[3,1:stepsize:end])]
                vcs = [Point3f(v1,v2,v3) for (v1,v2,v3) in zip(galaxies[i].stars.vel[1,1:stepsize:end], galaxies[i].stars.vel[2,1:stepsize:end], galaxies[i].stars.vel[3,1:stepsize:end])]
                print(size(pts)) 
                #arrows!(ax, 
                WGLMakie.arrows!(
                    #@view(galaxies[i].stars.pos[1,1:stepsize:end]), 
                    #@view(galaxies[i].stars.pos[2,1:stepsize:end]), 
                    #@view(galaxies[i].stars.pos[3,1:stepsize:end]),
                    #@view(galaxies[i].stars.vel[1,1:stepsize:end]), 
                    #@view(galaxies[i].stars.vel[2,1:stepsize:end]), 
                    #@view(galaxies[i].stars.vel[3,1:stepsize:end]),
                    pts,#[Point3f(p1,p2,p3) for (p1,p2,p3) in zip(galaxies[i].stars.pos[1,1:stepsize:end], galaxies[i].stars.pos[2,1:stepsize:end], galaxies[i].stars.pos[3,1:stepsize:end])],
                    vcs,#[Vec3f(v1,v2,v3) for (v1,v2,v3) in zip(galaxies[i].stars.vel[1,1:stepsize:end], galaxies[i].stars.vel[2,1:stepsize:end], galaxies[i].stars.vel[3,1:stepsize:end])],
                    # convert kms/s to kpc/10Myr by 0.010227047347441423
                    lengthscale = 0.010227047347441423,#arsize*1e-2,#arsize * (maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]),
                    #arrowsize   = ( arsize / maximum(@view(galaxies[1].stars.mass[1:stepsize:end])) ) .* @view(galaxies[i].stars.mass[1:stepsize:end]),
                    #arrowsize   = arsize,# (arsize*(maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["MASS"][all_notin][1:stepsize:end])) .* galaxies[i].stars.mass[1:stepsize:end],
                    linewidth   = 0.1,
                    linecolor   = @view(galaxies[i].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[i].stars.vel[:,1:stepsize:end]))) ./ maximum(@view(galaxies[1].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[1].stars.vel[:,1:stepsize:end])))),
                    arrowcolor  = @view(galaxies[i].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[i].stars.vel[:,1:stepsize:end]))) ./ maximum(@view(galaxies[1].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[1].stars.vel[:,1:stepsize:end])))),
                    colormap    = :winter,
                    align       = :head,
                    quality     = 3,
                    shading     = false,
                    ssao        = false,
                    normalize   = false,
                    transparency = true
                )
            elseif !boxfix
                #arrows!(ax, 
                pts = [Point3f(p1,p2,p3) for (p1,p2,p3) in zip(galaxies[i].stars.pos[1,1:stepsize:end], galaxies[i].stars.pos[2,1:stepsize:end], galaxies[i].stars.pos[3,1:stepsize:end])]
                vcs = [Point3f(v1,v2,v3) for (v1,v2,v3) in zip(galaxies[i].stars.vel[1,1:stepsize:end], galaxies[i].stars.vel[2,1:stepsize:end], galaxies[i].stars.vel[3,1:stepsize:end])]
                WGLMakie.arrows!(
                    #@view(galaxies[i].stars.pos[1,1:stepsize:end]), 
                    #@view(galaxies[i].stars.pos[2,1:stepsize:end]), 
                    #@view(galaxies[i].stars.pos[3,1:stepsize:end]),
                    #@view(galaxies[i].stars.vel[1,1:stepsize:end]), 
                    #@view(galaxies[i].stars.vel[2,1:stepsize:end]), 
                    #@view(galaxies[i].stars.vel[3,1:stepsize:end]),
                    pts,#[Point3f(p1,p2,p3) for (p1,p2,p3) in zip(galaxies[i].stars.pos[1,1:stepsize:end], galaxies[i].stars.pos[2,1:stepsize:end], galaxies[i].stars.pos[3,1:stepsize:end])],
                    vcs,#[Vec3f(v1,v2,v3) for (v1,v2,v3) in zip(galaxies[i].stars.vel[1,1:stepsize:end], galaxies[i].stars.vel[2,1:stepsize:end], galaxies[i].stars.vel[3,1:stepsize:end])],
                    # convert kms/s to 10 kpc/10Myr by 0.010227047347441423
                    lengthscale = 0.010227047347441423,#arsize*1e-2,#arsize * (maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["VEL"][:,all_notin][:,1:stepsize:end]),
                    #arrowsize   = ( arsize / maximum(galaxies[1].stars.mass[1:stepsize:end]) ) .* galaxies[i].stars.mass[1:stepsize:end],
                    #arrowsize   = arsize,# (arsize*(maximum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])-minimum(ALL_STARS["POS"][1,all_notin][1,1:stepsize:end])) / maximum(ALL_STARS["MASS"][all_notin][1:stepsize:end])) .* galaxies[i].stars.mass[1:stepsize:end],
                    linewidth   = 0.1,
                    linecolor   = @view(galaxies[i].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[i].stars.vel[:,1:stepsize:end]))) ./ maximum(@view(galaxies[1].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[1].stars.vel[:,1:stepsize:end])))),
                    arrowcolor  = @view(galaxies[i].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[i].stars.vel[:,1:stepsize:end]))) ./ maximum(@view(galaxies[1].stars.mass[1:stepsize:end]) .* norm.(eachcol(@view(galaxies[1].stars.vel[:,1:stepsize:end])))),
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
    end
    if spin
        merger_collection   = load( joinpath(indir, find_halo_file(indir=indir, felixID=felixID, rootID=rootID)), "merger_collection_STARS" )
        POSITIONS           = hcat( 
            read_galaxy_pos(galaxies[1], :physical)
            , read_galaxy_pos(galaxies[1], :physical)
            , read_galaxy_pos(galaxies[1], :physical)
             )
        #@show merger_collection["snapNR"]
        #println(maximum([1,findall(x->x.==snapNR, merger_collection["snapNR"])[1]-1]))
        #println(findall(x->x.==snapNR, merger_collection["snapNR"])[1])
        #println(minimum([length(merger_collection["snapNR"]), findall(x->x.==snapNR, merger_collection["snapNR"])[1]+1]))
        sANGMOM             = hcat( 
            @view(merger_collection["j_main"][:,maximum([1,findall(x->x.==snapNR, merger_collection["snapNR"])[1]-1])])
            , @view(merger_collection["j_main"][:,findall(x->x.==snapNR, merger_collection["snapNR"])[1]])
            , @view(merger_collection["j_main"][:,minimum([length(merger_collection["snapNR"]), findall(x->x.==snapNR, merger_collection["snapNR"])[1]+1])])
             )
        #arrows!(ax,
        pts = [Point3f(p1,p2,p3) for (p1,p2,p3) in zip(POSITIONS[1,:], POSITIONS[2,:], POSITIONS[3,:])]
        vcs = [Point3f(v1,v2,v3) for (v1,v2,v3) in zip(sANGMOM[1,:], sANGMOM[2,:], sANGMOM[3,:])]
        WGLMakie.arrows!(
            #@view(POSITIONS[1,:]), 
            #@view(POSITIONS[2,:]), 
            #@view(POSITIONS[3,:]),
            #@view(sANGMOM[1,:]), 
            #@view(sANGMOM[2,:]), 
            #@view(sANGMOM[3,:]),
            pts,#[Point3f(p1,p2,p3) for (p1,p2,p3) in zip(POSITIONS[1,:], POSITIONS[2,:], POSITIONS[3,:])],
            vcs,#[Vec3f(v1,v2,v3) for (v1,v2,v3) in zip(sANGMOM[1,:], sANGMOM[2,:], sANGMOM[3,:])],
            # convert kpc*kms/s to 5 kpc²/Myr by 0.10227047347441423/2
            lengthscale = 0.1227047347441423/5,
            arrowsize   = arsize,#norm(sANGMOM)/500,
            linewidth   = 1,
            linecolor   = 1:3,
            arrowcolor  = 1:3,
            colormap    = :Set1_3,
            quality     = 3,
            shading     = false,
            ssao        = false,
            normalize   = false,
            transparency = true
        )
    end
    return scene #fig
end

print("'plot_group'   ")



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