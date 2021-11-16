### Plot galaxies

include("/home/moon/sfortune/spinevo/spinevo.jl")
JSServe.configure_server!(listen_port=1688, forwarded_port=1688)
set_theme!(resolution=(1920, 1080), backgroundcolor = :black)


# Test plotting window
N = 60
function xy_data(x, y)
    r = sqrt(x^2 + y^2)
    r == 0.0 ? 1f0 : (sin(r)/r)
end
l = range(-10, stop = 10, length = N)
z = Float32[xy_data(x, y) for x in l, y in l]

scene = Scene()
surface!( scene,
        -1..1, -1..1, z,
        colormap = :solar
)
scene

##################################################################
### Settings

snapNR  = 128# 124# 128
subID   = 8127# 8366# 8127
nbID    = 8129
box     = "/HydroSims/Magneticum/Box4/uhr_test"


##################################################################



##################################################################
### Data processing

snapshot    = Snapshot(box, snapNR)
snapshot.snapbase
snapshot.subbase

galaxyID_list   = find_merging_progenitors(snapNR; felixID=18)
galaxies        = Dict{Int64, Any}()
particleID_list = Array{UInt64}(undef,0)
for i in 1:length(galaxyID_list)
    galaxies[i] = Galaxy(snapshot, galaxyID_list[i])
    rvir_group  = read_galaxy_prop(get_group(galaxies[1]), "RVIR", :physical)
    sph_small   = GadgetGalaxies.Sphere(0.1*rvir_group)
    sph_large   = GadgetGalaxies.Sphere(rvir_group)
    read_halo!(galaxies[i], units=:physical, props=((:stars, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=0.1*rvir_group)
    galaxies[i].stars.pos .+= read_galaxy_pos(galaxies[i], :physical)
    particleID_list         = vcat(particleID_list, galaxies[i].stars.id)
    #read_halo!(galaxies[i], units=:physical, props=((:gas, ["POS", "VEL", "MASS"]),), radius_units=:physical, radius=0.1*rvir_group)
    #read_halo!(galaxies[i], units=:physical, props=((:dm, ["POS", "VEL"]),), radius_units=:physical, radius=rvir_group)
end

# all particles
head        = read_header("$box/groups_$(@sprintf("%03i", snapNR))/sub_$(@sprintf("%03i", snapNR))")
filepath    = "$box/snapdir_$(@sprintf("%03i", snapNR))/snap_$(@sprintf("%03i", snapNR))"
blocks      = ["MASS", "POS", "ID", "VEL"]
radius      = 200
position    = read_galaxy_pos(galaxies[1], :sim)
ALL_STARS   = read_particles_in_volume(filepath, blocks, position, radius; parttype=4, verbose=true, use_keys=true)
#ALL_STARS["POS"]    .-= position
ALL_STARS["MASS"]   = convert_units_physical(ALL_STARS["MASS"], :mass, head)
ALL_STARS["POS"]    = convert_units_physical(ALL_STARS["POS"], :pos, head)


# Crop duplicates with target
all_notin   = ALL_STARS["ID"] .∉ Ref(Set(particleID_list))

stepsize = 1
size_factor = 3000 / maximum(galaxies[1].stars.mass)
radius_plot = 500 #kpc
#f = WGLMakie.Figure()
#ax = Axis3(f[1, 1])
#xlims!(ax, (read_galaxy_pos(galaxies[1], :physical)[1]-200, read_galaxy_pos(galaxies[1], :physical)[1]+200))
#ylims!(ax, (read_galaxy_pos(galaxies[1], :physical)[2]-200, read_galaxy_pos(galaxies[1], :physical)[2]+200))
#zlims!(ax, (read_galaxy_pos(galaxies[1], :physical)[3]-200, read_galaxy_pos(galaxies[1], :physical)[3]+200))
#f
scene = Scene()
scatter!( scene,
    ALL_STARS["POS"][:,all_notin],
    color = ALL_STARS["MASS"][all_notin] .* size_factor,
    #markersize = log10.(ALL_STARS["MASS"][1:stepsize:end]) .* size_factor ),
    markersize = ALL_STARS["MASS"][all_notin] .* size_factor,
    colormap = :autumn1,
    transparency = false
    )

# Halos
#scene = Scene()
for i in 1:length(galaxyID_list)
    if i == 1
        scatter!( scene,
            galaxies[i].stars.pos[:,1:stepsize:end],
            color = galaxies[i].stars.mass[1:stepsize:end] .* size_factor,
            markersize = galaxies[i].stars.mass[1:stepsize:end] .* size_factor,#markersize = 0.1,
            colormap = :winter,
            transparency = false
            )
    else
        scatter!( scene,
            galaxies[i].stars.pos[:,1:stepsize:end],
            color = galaxies[i].stars.mass[1:stepsize:end] .* size_factor,
            markersize = galaxies[i].stars.mass[1:stepsize:end] .* size_factor,#markersize = 0.1,
            colormap = :summer,
            transparency = false
            )
    end
end
#cam3d!(scene; lookat = read_galaxy_pos(galaxies[1], :physical))
#Camera3D(scene; rotation_center = read_galaxy_pos(galaxies[1], :physical))
scene

#plotrange   = FRect3D(read_galaxy_pos(galaxies[1], :physical), [100,100,100])
#update_limits!(scene, plotrange)


g           = Galaxy(snapshot, subID)
g_nb        = Galaxy(snapshot, nbID)
rvir_group  = read_galaxy_prop(get_group(g), "RVIR", :physical)
sph_small   = GadgetGalaxies.Sphere(0.1*rvir_group)
sph_large   = GadgetGalaxies.Sphere(rvir_group)
main_pos    = read_galaxy_pos(g, :physical)
nb_pos      = read_galaxy_pos(g_nb,:physical)
read_halo!(g, units=:physical, props=((:gas, ["POS", "VEL", "MASS"]), (:stars, ["POS", "VEL", "MASS"])), radius_units=:physical, radius=0.1*rvir_group) 
read_halo!(g, units=:physical, props=((:dm, ["POS", "VEL"]),), radius_units=:physical, radius=rvir_group)
read_halo!(g_nb, units=:physical, props=((:gas, ["POS", "VEL", "MASS"]), (:stars, ["POS", "VEL", "MASS"])), radius_units=:physical, radius=0.1*rvir_group)
read_halo!(g_nb, units=:physical, props=((:dm, ["POS", "VEL"]),), radius_units=:physical, radius=rvir_group)

# Shift neighbor to FP frame

#g.stars.pos  .+= main_pos
#g.gas.pos    .+= main_pos
#g.dm.pos     .+= main_pos
g_nb.stars.pos  .+= nb_pos - main_pos
g_nb.gas.pos    .+= nb_pos - main_pos
g_nb.dm.pos     .+= nb_pos - main_pos


groupID = GadgetGalaxies.get_group(g)
gr          = GadgetGalaxies.GalaxyGroup(snapshot, 20)

# stack all sub halos
stacked_pos     = hcat(g.stars.pos, g_nb.stars.pos)
stacked_mass    = vcat(g.stars.mass, g_nb.stars.mass)

# Angular Momentum
J_STARS     = convert( AbstractVector{Float64}, GadgetGalaxies.angular_momentum(g.stars, sph_small) )
j_STARS     = convert( AbstractVector{Float64}, GadgetGalaxies.specific_angular_momentum(g.stars, sph_small))
Jj_STARS    = hcat(normalize(J_STARS), normalize(j_STARS))
abs_Jj_STARS= [norm(J_STARS), norm(j_STARS)]

reduction_factor = 2
particle_frac   = length(g.stars.pos[1,:]) / length(g_nb.stars.pos[1,:])
plot_positions2, plot_masses2, plot_n2 = reduce_n_particles(length(g_nb.stars.pos[1,:])/reduction_factor/particle_frac, g_nb.stars.pos, g_nb.stars.mass, "center", 0.1)
plot_positions, plot_masses, plot_n = reduce_n_particles(length(g.stars.pos[1,:])/reduction_factor, g.stars.pos, g.stars.mass, "center", 0.5)
stacked_grid_pos, stacked_grid_mass, stacked_grid_n = reduce_n_particles(length(stacked_pos[1,:])/reduction_factor, stacked_pos, stacked_mass, "center", 0.5)

##################################################################

##################################################################
### Plot

# Original
stepsize = 1
size_factor = 2500
scene = Scene()
scatter!( scene,
    g.stars.pos[:,1:stepsize:end],
    color = g.stars.mass[1:stepsize:end] ./( maximum(g.stars.mass[1:stepsize:end]) ),
    markersize = g.stars.mass[1:stepsize:end] ./( maximum(g.stars.mass[1:stepsize:end])/size_factor ),#markersize = 0.1,
    colormap = :winter,
    transparency = false
    )
scatter!( scene,
    g_nb.stars.pos[:,1:stepsize:end],
    color = g_nb.stars.mass[1:stepsize:end] ./( maximum(g.stars.mass[1:stepsize:end]) ),
    markersize = g_nb.stars.mass[1:stepsize:end] ./( maximum(g.stars.mass[1:stepsize:end])/size_factor ),#markersize = 0.1,
    colormap = :autumn1,
    transparency = false
    )
#arrows!( scene,
#    zeros(3,2),
#    Jj_STARS .* 50,
#    arrowsize = 0.03, linecolor = (:white, 0.6), linewidth = 3
#    )
scene

# Stacked sub halos
stepsize = 1
size_factor = 10000
scene = Scene()
scatter!( scene,
    stacked_pos[:,1:stepsize:end],
    color = stacked_mass[1:stepsize:end] ./( maximum(stacked_mass[1:stepsize:end]) ),
    #markersize = log10.(stacked_mass[1:stepsize:end]) ./( maximum(log10.(stacked_mass[1:stepsize:end]))/size_factor ),
    markersize = stacked_mass[1:stepsize:end] ./( maximum(stacked_mass[1:stepsize:end])/size_factor ),
    colormap = :winter,
    transparency = false
    )
scene

# Processed
size_factor = 1000
scene = Scene()
scatter!( scene,
    plot_positions,
    markersize = log10.(plot_n) ./ (maximum(log10.(plot_n2))/size_factor), #plot_masses ./ (maximum(plot_masses)/1), #
    color = plot_masses,
    colormap = :winter,
    transparency = false
    )
scatter!( scene,
    plot_positions2,
    markersize = log10.(plot_n2) ./ (maximum(log10.(plot_n2))/size_factor), #plot_masses ./ (maximum(plot_masses)/1), #
    color = plot_masses2,
    colormap = :autumn1,
    transparency = false
    )
scene

# Stacked with grid
stepsize = 1
size_factor = 2000
scene = Scene()
scatter!( scene,
    stacked_grid_pos[:,1:stepsize:end],
    color = stacked_grid_mass[1:stepsize:end] ./( maximum(stacked_grid_mass[1:stepsize:end]) ),
    markersize = log10.(stacked_grid_n) ./ (maximum(log10.(stacked_grid_n))/size_factor), 
    colormap = :winter,
    transparency = false
    )
scene


# Testing Earlier snap
stepsize = 1
size_factor = 1000
scene = Scene()
scatter!( scene,
    g.stars.pos[:,1:stepsize:end],
    color = g.stars.mass[1:stepsize:end] ./( maximum(g.stars.mass[1:stepsize:end]) ),
    markersize = g.stars.mass[1:stepsize:end] ./( maximum(g.stars.mass[1:stepsize:end])/size_factor ),#markersize = 0.1,
    colormap = :winter,
    transparency = false
    )
scene


println("Plotting finished")
##################################################################


# Testing Area
head        = read_header("$box/groups_$(@sprintf("%03i", snapNR))/sub_$(@sprintf("%03i", snapNR))")
filepath    = "$box/snapdir_$(@sprintf("%03i", snapNR))/snap_$(@sprintf("%03i", snapNR))"
blocks      = ["MASS", "POS"]
radius      = 500
position    = read_galaxy_pos(g, :sim)
particles   = read_particles_in_volume(filepath, blocks, position, radius; parttype=4, verbose=true, use_keys=true)
particles["POS"]    .-= position
particles["MASS"]   = convert_units_physical(particles["MASS"], :mass, head)
particles["POS"]    = convert_units_physical(particles["POS"], :pos, head)


stepsize = 1
size_factor = 10000
scene = Scene()
scatter!( scene,
    particles["POS"][:,1:stepsize:end],
    color = particles["MASS"][1:stepsize:end] ./( maximum(particles["MASS"][1:stepsize:end]) ),
    #markersize = log10.(particles["MASS"][1:stepsize:end]) ./( maximum(log10.(particles["MASS"][1:stepsize:end]))/size_factor ),
    markersize = particles["MASS"][1:stepsize:end] ./( maximum(particles["MASS"][1:stepsize:end])/size_factor ),
    colormap = :winter,
    transparency = false
    )
#scene


# Stacked sub halos
#stepsize = 1
#size_factor = 10000
#scene = Scene()
#update_limits!(scene, Rect3f(Vec3f(-300, -300, -300), Vec3f(300, 300, 300)))
scatter!( scene,
    g.stars.pos[:,1:stepsize:end],
    color = g.stars.mass[1:stepsize:end] ./( maximum(particles["MASS"][1:stepsize:end]) ),
    #markersize = log10.(g.stars.mass[1:stepsize:end]) ./( maximum(log10.(g.stars.mass[1:stepsize:end]))/size_factor ),
    markersize = g.stars.mass[1:stepsize:end] ./( maximum(particles["MASS"][1:stepsize:end])/size_factor ),
    colormap = :autumn1,
    transparency = false
    )
scene
