println("   Plotting...")
using PyPlot
using ColorSchemes
using JSServe
using WGLMakie
WGLMakie.activate!()

include("/home/moon/sfortune/spinevo/pkg/plot_group_collection.jl")

###########
# Lucas'
###########

function get_bvalue_cmap()
    # cs = [colorant"hsl(0, 100%, 18%)", colorant"hsl(0, 100%, 50%)", colorant"hsl(54, 100%, 81%)", colorant"hsl(212, 86%, 83%)", colorant"hsl(212, 100%, 50%)", colorant"hsl(212, 100%, 18%)"]
    # rs = [529, 170, 85, 85, 131]
    cs = [colorant"hsl(0, 100%, 18%)", colorant"hsl(0, 100%, 50%)", colorant"hsl(54, 100%, 81%)", colorant"hsl(54, 40%, 83%)", colorant"hsl(54, 5%, 83%)", colorant"hsl(212, 5%, 83%)", colorant"hsl(212, 60%, 83%)", colorant"hsl(212, 100%, 50%)", colorant"hsl(212, 100%, 18%)"]
    rs = [529, 170, 65, 12, 2, 15, 76, 131] # 1000 colors sampled

    cmap = ColorScheme(vcat([range(cs[i], cs[i+1]; length=r+1)[1:end-1] for (i, r) in enumerate(rs)]...))
    return cmap
end

get_blackberry() = "#8D1D75"

##########################



@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function testmin(; tc_port=1688
    )
    JSServe.configure_server!(listen_port=tc_port, forwarded_port=tc_port)
    set_theme!(resolution=(1600, 900), backgroundcolor = :white)
    # Test plotting window
    x = rand(10)
    y = rand(10)
    z = rand(10)
    w = rand(3,10)
    sc = WGLMakie.scatter(
        w, markersize=10
        )
    
    #scene = Scene()
    #fig = WGLMakie.Figure()
    #ax = WGLMakie.Axis3(fig[1,1], viewmode=:fit)
    scatter!(
        z,x,y, markersize=50
        )
    
    return sc
end

print("'testmin'   ")


@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function testconnection(; tc_port=1688
    )
    JSServe.configure_server!(listen_port=tc_port, forwarded_port=tc_port)
    set_theme!(resolution=(1600, 900), backgroundcolor = :black)
    # Test plotting window
    N = 60
    function xy_data(x, y)
        r = sqrt(x^2 + y^2)
        r == 0.0 ? 1f0 : (sin(r)/r)
    end
    l = range(-10, stop = 10, length = N)
    z = Float32[xy_data(x, y) for x in l, y in l]
    return WGLMakie.surface(-1..1, -1..1, z, colormap = :winter, shading = false)
end

print("'testconnection'   ")


@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function vec_surface(pos, v; r=1, tc_port=1688
    )
    JSServe.configure_server!(listen_port=tc_port, forwarded_port=tc_port)
    set_theme!(resolution=(1600, 900), backgroundcolor = :black)
    # Test plotting window
    N=10
    function srfc_z(v, pos, x,y)
        return ( pos[3] - v[1]*(x-pos[1]) - v[2]*(y-pos[2]) ) / v[3]
    end
    l = range(-r*norm(v), stop = r*norm(v), length = N)
    z = Float64[srfc_z(v, pos, x,y) for x in l, y in l]
    #idcs    = findcs(z, geq=-r*norm(v)+pos[3], leq=r*norm(v)+pos[3])
    idcs    = findall(dm-> -r*norm(v)+pos[3] .<= dm .<= r*norm(v)+pos[3], z)
    x_idcs  = Array{Int64}(undef, 0)
    y_idcs  = Array{Int64}(undef, 0)
    for i in 1:length(idcs)
        x_idcs  = vcat( x_idcs, idcs[i][1] )
        y_idcs  = vcat( y_idcs, idcs[i][2] )
    end
    return surface(l, l, z, color = :blue, shading = false, transparency = true)
end

print("'vec_surface'   ")




@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_halostory_arrows(; outfile="./OUT_plot_halostory_arrows.png", felixID=" ", rootID=" "
    , parttype="STARS", snap=false
    , scale=1, plotwidth=16, plotheight=3, gridspace=0.6 #Gyr
    , indir=current_dir_stories
    )

    halo_filestring = find_halo_file(indir=indir,felixID=felixID,rootID=rootID)

    merger_collection = load(joinpath(indir, halo_filestring), string("merger_collection_", parttype))
    
    lbt     = Array{Float64}(undef, 1)
    lbt[1]  = merger_collection["lookbacktime"][1]
    rs      = Array{Float64}(undef, 1)
    rs[1]   = merger_collection["redshift"][1]
    snaps    = Array{Int64}(undef, 1)
    snaps[1] = merger_collection["snapNR"][1]
    for i in 2:length(merger_collection["snapNR"])-1
        if lbt[end] - merger_collection["lookbacktime"][i] > gridspace
            snaps    = vcat( snaps, merger_collection["snapNR"        ][i] )
            lbt     = vcat( lbt , merger_collection["lookbacktime"][i] )
            rs      = vcat( rs  , merger_collection["redshift"    ][i] )
        end
    end
    snaps    = vcat( snaps, merger_collection["snapNR"        ][end] )
    lbt     = vcat( lbt , merger_collection["lookbacktime"][end] )
    rs      = vcat( rs  , merger_collection["redshift"    ][end] )
    
    
    J_norm      = vcat( transpose(norm.(eachcol(merger_collection["J_main"]))), transpose(norm.(eachcol(merger_collection["J_main"]))), transpose(norm.(eachcol(merger_collection["J_main"]))) )
    
    proj_vec    = (merger_collection["J_main"]./J_norm)[1,:]
    x_vec       = (merger_collection["J_main"]./J_norm)[2,:]
    y_vec       = (merger_collection["J_main"]./J_norm)[3,:]
    
    pltcm       = pyimport("matplotlib.cm")
    pltcolors   = pyimport("matplotlib.colors")
    colormap    = pltcm.get_cmap(name="coolwarm")
    
    fig, ax = subplots()
    if typeof(rootID) == Int
        ax.set_title("Sub Halo $rootID")
    elseif typeof(felixID) == Int
        ax.set_title("Sub Halo $(find_rootID(felixID)[1])")
    end
    ax.set_ylim(-1, 1)
    ax.set_yticks([-1,0,1])
    ax.set_xlim( minimum(merger_collection["lookbacktime"])-0.2-2*maximum([0, x_vec[end]]), maximum(merger_collection["lookbacktime"])+0.2-2*minimum([0, x_vec[1]]) )
    ax.set_xticks( lbt )
    ax.set_xlabel("Lookback Time [Gyr]")
    
    ax.quiver(merger_collection["lookbacktime"], zeros(length(merger_collection["lookbacktime"])), 
            x_vec, y_vec,
            color=colormap(proj_vec.*0.5.+0.5), scale_units="y", scale=scale, alpha=0.7, width=0.005, headwidth=2, headlength=3, edgecolor="black", lw=0.5)
    
    ax.invert_xaxis()
    ax.grid()
    clb = fig.colorbar(pltcm.ScalarMappable(norm=pltcolors.Normalize(vmin=-1, vmax=1), cmap="coolwarm"), ticks=[-1, 0, 1])
    clb.ax.set_ylabel("Z-Axis Component")
    clb.ax.yaxis.set_label_position("left")
    
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(ax.get_xticks())
    ax2.set_xticklabels(round.(rs, digits=2)[1:1:end])
    ax2.set_xlabel("z")
    
    if snap
        ax3 = ax.twiny()
        ax3.set_xlim(ax.get_xlim())
        ax3.set_xticks(ax.get_xticks())
        ax3.set_xticklabels(snaps[1:1:end])
        ax3.set_xlabel("Snap")
        ax3.xaxis.set_ticks_position("bottom")
        ax3.xaxis.set_label_position("bottom")
        ax3.spines["bottom"].set_position(("outward", 40))
    end
    
    fig.set_size_inches(plotwidth, plotheight)
    
    fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
    
    return nothing
end

print("'plot_halostory_arrows'   ")














@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_halostory_data(; outfile="./OUT_plot_halostory_data.png", felixID=" ", rootID=" "
    , parttype="STARS", snap=false
    , scale=0.7, gridspace=0.6 #Gyr
    , indir=current_dir_stories
    )
    
    mcs = load(joinpath(indir, find_halo_file(indir=indir,felixID=felixID,rootID=rootID)), "merger_collection_$parttype")
    
    lbt     = Array{Float64}(undef, 1)
    lbt[1]  = mcs["lookbacktime"][1]
    rs      = Array{Float64}(undef, 1)
    rs[1]   = mcs["redshift"][1]
    snaps    = Array{Int64}(undef, 1)
    snaps[1] = mcs["snapNR"][1]
    for i in 2:length(mcs["snapNR"])-1
        if lbt[end] - mcs["lookbacktime"][i] > gridspace
            snaps    = vcat( snaps, mcs["snapNR"        ][i] )
            lbt     = vcat( lbt , mcs["lookbacktime"][i] )
            rs      = vcat( rs  , mcs["redshift"    ][i] )
        end
    end
    snaps    = vcat( snaps, mcs["snapNR"        ][end] )
    lbt     = vcat( lbt , mcs["lookbacktime"][end] )
    rs      = vcat( rs  , mcs["redshift"    ][end] )
    
    j_spherical = missings(Float64, 3, 0)
    for i in 1:length(mcs["j_main"][1,:])
            j_spherical = hcat( j_spherical, cartesian_to_spherical(mcs["j_main"][:,i]) )
    end
    
    ### Figure
    N_plots = 5
    fig, ax = subplots(N_plots, sharex=true)
    fig.subplots_adjust(hspace=0)
    if typeof(rootID) == Int
        ax[1].set_title("Sub Halo $rootID")
    elseif typeof(felixID) == Int
        ax[1].set_title("Sub Halo $(find_rootID(felixID)[1])")
    end

    ax[1].set_xticks( lbt )
    
    ax[1].plot(mcs["lookbacktime"], replace(mcs["ϕ_flip"], missing => NaN), "b-", lw=5, label="ϕ_flip", alpha=1, zorder=3)
    ax[1].plot(mcs["lookbacktime"], replace(j_spherical[2,:], missing => NaN), "g-", lw=2, label="θ", alpha=1, zorder=1)
    ax[1].plot(mcs["lookbacktime"], replace(j_spherical[3,:], missing => NaN), "c-", lw=2, label="ϕ", alpha=1, zorder=2)
    ax[1].invert_xaxis()
    ax[1].set_ylabel("Angle [°]")
    ax[1].grid()
    ax[1].axhline(color="black")
    ax[1].legend(loc="upper center", frameon=true, borderpad=1, handlelength=1.8)
    ax2 = ax[1].twiny()
    ax2.set_xlim(ax[1].get_xlim())
    ax2.set_xticks(ax[1].get_xticks())
    ax2.set_xticklabels(round.(rs, digits=2)[1:1:end])
    ax2.set_xlabel("z")
    
    ax[2].plot(mcs["lookbacktime"], mcs["BVAL_0"], "-", lw=2, label="b-value converted to z=0", color="green", alpha=1, zorder=2)
    ax[2].grid()
    ax[2].axhline(b_disk,color="blue", label="Disks") # disks
    ax[2].axhline(b_ell,color="red", label="Ellipticals") # ellipticals
    ax[2].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
    
    ax[3].plot(mcs["lookbacktime"], norm.(eachcol(replace(mcs["j_main"], missing => 0.))), "-", lw=2, label="| j_$parttype |", color="navy", alpha=1, zorder=2)
    ax[3].set_ylabel("[kpc km/s]")
    ax[3].grid()
    ax[3].axhline(color="black")
    ax[3].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
    
    #ax[4].plot(mcs["lookbacktime"], norm.(eachcol(replace(mcs["Δj_main"], missing => 0.))), "-", lw=2, label="| Δj_$parttype |", color="navy", alpha=1, zorder=2)
    #ax[4].set_ylabel("[kpc km/s]")
    #ax[4].grid()
    #ax[4].axhline(color="black")
    #ax[4].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
    
    ax[4].plot(mcs["lookbacktime"], pyplottable(mcs["SFR"]), "-", lw=2, label="SFR", color="deepskyblue", alpha=1, zorder=2)
    ax[4].set_ylabel("[ M⊙ / yr ]")
    ax[4].grid()
    ax[4].axhline(color="black")
    ax[4].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
    
    ax[5].plot(mcs["lookbacktime"], mcs["M2"], "-", lw=2, label="M2_$parttype", color="gold", alpha=1, zorder=4)
    ax[5].plot(mcs["lookbacktime"], mcs["M"], "-", lw=2, label="M_$parttype", color="green", alpha=1, zorder=7)
    ax[5].plot(mcs["lookbacktime"], mcs["M2"], ".", ms=7, color="gold", alpha=0.5, zorder=8)
    ax[5].plot(mcs["lookbacktime"], mcs["M"],  ".", ms=7, color="green", alpha=0.5, zorder=9)
    ax[5].plot(mcs["merger map"][7,:], mcs["merger map"][9,:], ".", ms=7, color="gold", alpha=0.5, zorder=8)
    ax[5].plot(mcs["merger map"][7,:], mcs["merger map"][1,:], ".", ms=7, color="green", alpha=0.5, zorder=9)
    ax[5].set_ylabel("Mass [M⊙]")
    ax[5].set_yscale("log")
    ax[5].grid()
    ax[5].axhline(color="black")
    ax[5].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
    
    ax[N_plots].set_xticks( lbt )

    if snap
        ax3 = ax[N_plots].twiny()
        ax3.set_xlim(ax[N_plots].get_xlim())
        ax3.set_xticks(ax[N_plots].get_xticks())
        ax3.set_xticklabels(snaps)
        ax3.set_xlabel("Snap")
        ax3.xaxis.set_ticks_position("bottom")
        ax3.xaxis.set_label_position("bottom")
        ax3.spines["bottom"].set_position(("outward", 40))
    end
    
    fig.tight_layout()#h_pad=0.0)
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(16scale, 4*N_plots*scale)
    fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)    
    
    return nothing
end

print("'plot_halostory_data'   ")














@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_halostory_full(; outfile="./OUT_plot_halostory_data.png", felixID=" ", rootID=" "
    , parttype="STARS", snap=false, title = false
    , scale=0.7, gridspace=0.6 #Gyr
    , indir=current_dir_stories
    )
    
    mcs = load(joinpath(indir, find_halo_file(indir=indir,felixID=felixID,rootID=rootID)), "merger_collection_$parttype")
    
    lbt     = Array{Float64}(undef, 1)
    lbt[1]  = mcs["lookbacktime"][1]
    rs      = Array{Float64}(undef, 1)
    rs[1]   = mcs["redshift"][1]
    snaps    = Array{Int64}(undef, 1)
    snaps[1] = mcs["snapNR"][1]
    for i in 2:length(mcs["snapNR"])-1
        if lbt[end] - mcs["lookbacktime"][i] > gridspace
            snaps    = vcat( snaps, mcs["snapNR"        ][i] )
            lbt     = vcat( lbt , mcs["lookbacktime"][i] )
            rs      = vcat( rs  , mcs["redshift"    ][i] )
        end
    end
    snaps    = vcat( snaps, mcs["snapNR"        ][end] )
    lbt     = vcat( lbt , mcs["lookbacktime"][end] )
    rs      = vcat( rs  , mcs["redshift"    ][end] )
    
    j_spherical = missings(Float64, 3, 0)
    for i in 1:length(mcs["j_main"][1,:])
            j_spherical = hcat( j_spherical, cartesian_to_spherical(mcs["j_main"][:,i]) )
    end
    
    fake1 = missings(Int64, length(mcs["M"]))
    for i in 1:length(mcs["M"])
        fake1[i] = fake_flip_finder(mcs, i; verbose=false)
    end
    fakeflip = findcs(fake1, eq=1)
    switches = findcs(mcs["switch"], eq=1)
    small_fi = findcs(mcs["M"], lt=2e10)
    small_st = findcs(mcs["M"] .- mcs["ΔM"], lt=2e10)
    excluded = vcat(fakeflip,switches, small_fi, small_st)
    
    # arrows 
    J_norm      = vcat( transpose(norm.(eachcol(mcs["J_main"]))), transpose(norm.(eachcol(mcs["J_main"]))), transpose(norm.(eachcol(mcs["J_main"]))) )
    proj_vec    = (mcs["J_main"]./J_norm)[1,:]
    x_vec       = (mcs["J_main"]./J_norm)[2,:]
    y_vec       = (mcs["J_main"]./J_norm)[3,:]
    pltcm       = pyimport("matplotlib.cm")
    pltcolors   = pyimport("matplotlib.colors")
    colormap    = pltcm.get_cmap(name="coolwarm")

    ### Figure
    N_plots = 6
    fig, ax = subplots(N_plots, sharex=true)
    fig.subplots_adjust(hspace=0)
    if typeof(rootID) == Int && title == true
        ax[1].set_title("Sub Halo $rootID")
    elseif typeof(felixID) == Int && title == true
        ax[1].set_title("Sub Halo $(find_rootID(felixID)[1])")
    end

    ax[1].set_xticks( lbt )
    
    ax[1].set_ylim(0, 180)
    #ax[1].vlines(mcs["lookbacktime"][excluded], ymin=0, ymax=180, color="black")
    ax[1].vlines(mcs["lookbacktime"][small_fi], ymin=0, ymax=180, color="black")
    ax[1].vlines(mcs["lookbacktime"][small_st], ymin=0, ymax=180, color="black")
    ax[1].vlines(mcs["lookbacktime"][fakeflip], ymin=0, ymax=180, color="firebrick")
    ax[1].vlines(mcs["lookbacktime"][switches], ymin=0, ymax=180, color="blueviolet")
    ax[1].plot(mcs["lookbacktime"], replace(mcs["ϕ_flip"], missing => NaN), "b-", lw=5, label="ϕ_flip", alpha=1, zorder=3)
    ax[1].plot(mcs["lookbacktime"], replace(j_spherical[2,:], missing => NaN), "g-", lw=2, label="θ", alpha=1, zorder=1)
    ax[1].plot(mcs["lookbacktime"], replace(j_spherical[3,:], missing => NaN), "c-", lw=2, label="ϕ", alpha=1, zorder=2)
    ax[1].set_ylabel("Angle [°]")
    ax[1].grid()
    ax[1].axhline(color="black")
    ax[1].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
    ax[1].invert_xaxis()
    #ax[1].vlines(mcs["lookbacktime"][switches], ymin=ax[1].get_ylim()[1], ymax=ax[1].get_xlim()[2], color="black")
    ax2 = ax[1].twiny()
    ax2.set_xlim(ax[1].get_xlim())
    ax2.set_xticks(ax[1].get_xticks())
    ax2.set_xticklabels(round.(rs, digits=2)[1:1:end])
    ax2.set_xlabel("z")
    

    ax[2].set_ylim(-1, 1)
    #ax[2].vlines(mcs["lookbacktime"][excluded], ymin=-1, ymax=1, color="black")
    ax[2].vlines(mcs["lookbacktime"][small_fi], ymin=-1, ymax=1, color="black")
    ax[2].vlines(mcs["lookbacktime"][small_st], ymin=-1, ymax=1, color="black")
    ax[2].vlines(mcs["lookbacktime"][fakeflip], ymin=-1, ymax=1, color="firebrick")
    ax[2].vlines(mcs["lookbacktime"][switches], ymin=-1, ymax=1, color="blueviolet")
    ax[2].set_yticks([-1,0,1])
    ax[2].grid()
    ax[2].axhline(color="black")
    #ax[2].set_xlim( minimum(mcs["lookbacktime"])-0.2-2*maximum([0, x_vec[end]]), maximum(mcs["lookbacktime"])+0.2-2*minimum([0, x_vec[1]]) )
    ax[2].quiver(mcs["lookbacktime"], zeros(length(mcs["lookbacktime"])), 
            x_vec, y_vec,
            color=colormap(proj_vec.*0.5.+0.5), scale_units="y", scale=scale, alpha=0.7, width=0.005, headwidth=2, headlength=3, edgecolor="black", lw=0.5)
    #ax[2].invert_xaxis()
    #
    

    ax[3].set_ylim(-6.1, -3.9)
    #ax[3].vlines(mcs["lookbacktime"][excluded], ymin=-6.1, ymax=-3.9, color="black")
    ax[3].vlines(mcs["lookbacktime"][small_fi], ymin=-6.1, ymax=-3.9, color="black")
    ax[3].vlines(mcs["lookbacktime"][small_st], ymin=-6.1, ymax=-3.9, color="black")
    ax[3].vlines(mcs["lookbacktime"][fakeflip], ymin=-6.1, ymax=-3.9, color="firebrick")
    ax[3].vlines(mcs["lookbacktime"][switches], ymin=-6.1, ymax=-3.9, color="blueviolet")
    ax[3].plot(mcs["lookbacktime"], mcs["BVAL_0"], "-", lw=2, label="b-value converted to z=0", color="green", alpha=1, zorder=2)
    ax[3].grid()
    ax[3].axhline(b_disk,color="blue", label="Disks") # disks
    ax[3].axhline(b_ell,color="red", label="Ellipticals") # ellipticals
    ax[3].set_ylim([-6.1,-3.9])
    ax[3].legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
    #ax[3].vlines(mcs["lookbacktime"][switches], ymin=ax[3].get_ylim()[1], ymax=ax[3].get_xlim()[2], color="black")
    #ax[3].invert_xaxis()
    
    ax[4].set_ylim(0, 4)
    #ax[4].vlines(mcs["lookbacktime"][excluded], ymin=0, ymax=4, color="black")
    ax[4].vlines(mcs["lookbacktime"][small_fi], ymin=0, ymax=4, color="black")
    ax[4].vlines(mcs["lookbacktime"][small_st], ymin=0, ymax=4, color="black")
    ax[4].vlines(mcs["lookbacktime"][fakeflip], ymin=0, ymax=4, color="firebrick")
    ax[4].vlines(mcs["lookbacktime"][switches], ymin=0, ymax=4, color="blueviolet")
    ax[4].plot(mcs["lookbacktime"], log10.(pyplottable(mcs["j_main"])), "-", lw=2, label="log(| j_$parttype |)", color="navy", alpha=1, zorder=2)
    ax[4].set_ylabel("log[kpc km/s]")
    #ax[4].set_yscale("log")
    ax[4].grid()
    ax[4].axhline(color="black")
    ax[4].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
    #ax[4].vlines(mcs["lookbacktime"][switches], ymin=ax[4].get_ylim()[1], ymax=ax[4].get_xlim()[2], color="black")
    #ax[4].invert_xaxis()
    
    #ax[4].plot(mcs["lookbacktime"], norm.(eachcol(replace(mcs["Δj_main"], missing => 0.))), "-", lw=2, label="| Δj_$parttype |", color="navy", alpha=1, zorder=2)
    #ax[4].set_ylabel("[kpc km/s]")
    #ax[4].grid()
    #ax[4].axhline(color="black")
    #ax[4].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
    
    ax[5].set_ylim(-5, 1)
    #ax[5].vlines(mcs["lookbacktime"][excluded], ymin=-5, ymax=1, color="black")
    ax[5].vlines(mcs["lookbacktime"][small_fi], ymin=-5, ymax=1, color="black")
    ax[5].vlines(mcs["lookbacktime"][small_st], ymin=-5, ymax=1, color="black")
    ax[5].vlines(mcs["lookbacktime"][fakeflip], ymin=-5, ymax=1, color="firebrick")
    ax[5].vlines(mcs["lookbacktime"][switches], ymin=-5, ymax=1, color="blueviolet")
    ax[5].plot(mcs["lookbacktime"], log10.(pyplottable(1e9 .* mcs["SFR"] ./ mcs["M"])), "-", lw=2, label="sSFR", color="deepskyblue", alpha=1, zorder=2)
    ax[5].set_ylabel("1 / Gyr")
    ax[5].grid()
    #ax[5].axhline(color="black")
    #ax[5].set_yscale("log")
    ax[5].legend(loc="lower left", frameon=true, borderpad=1, handlelength=1.8)
    #ax[5].vlines(mcs["lookbacktime"][switches], ymin=ax[1].get_ylim()[1], ymax=ax[1].get_xlim()[2], color="black")
    #ax[5].invert_xaxis()
    
    ax[6].set_ylim(6, 14)
    #ax[6].vlines(mcs["lookbacktime"][excluded], ymin=6, ymax=14, color="black")
    ax[6].vlines(mcs["lookbacktime"][small_fi], ymin=6, ymax=14, color="black")
    ax[6].vlines(mcs["lookbacktime"][small_st], ymin=6, ymax=14, color="black")
    ax[6].vlines(mcs["lookbacktime"][fakeflip], ymin=6, ymax=14, color="firebrick")
    ax[6].vlines(mcs["lookbacktime"][switches], ymin=6, ymax=14, color="blueviolet")
    ax[6].hlines(log10(2e10), xmin=minimum(lbt), xmax=maximum(lbt), color="black")
    #ax[6].plot(mcs["lookbacktime"], mcs["M2"], "-", lw=2, label="M2_$parttype", color="gold", alpha=1, zorder=4)
    ax[6].plot(mcs["lookbacktime"], log10.(mcs["M"]), "-", lw=2, color="green", alpha=1, zorder=7)
    #ax[6].plot(mcs["lookbacktime"], mcs["M2"], ".", ms=7, color="gold", alpha=0.5, zorder=8)
    ax[6].plot(mcs["lookbacktime"], log10.(mcs["M"]),  ".", ms=7, color="green", alpha=0.5, zorder=9)
    #ax[6].plot(mcs["merger map"][7,:], mcs["merger map"][9,:], ".", ms=7, color="gold", alpha=0.5, zorder=8)
    ax[6].plot(mcs["merger map"][7,:], log10.(mcs["merger map"][1,:]), ".", ms=7, color="green", alpha=0.5, zorder=9)
    ax[6].set_ylabel("log( M_$parttype [M⊙] )")
    #ax[6].set_yscale("log")
    ax[6].grid()
    ax[6].axhline(color="black")
    #ax[6].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
    #ax[6].vlines(mcs["lookbacktime"][switches], ymin=ax[6].get_ylim()[1], ymax=ax[6].get_xlim()[2], color="black")
    #ax[6].invert_xaxis()
    
    ax[N_plots].set_xticks( lbt )

    if snap
        ax3 = ax[N_plots].twiny()
        ax3.set_xlim(ax[N_plots].get_xlim())
        ax3.set_xticks(ax[N_plots].get_xticks())
        ax3.set_xticklabels(snaps)
        ax3.set_xlabel("Snap")
        ax3.xaxis.set_ticks_position("bottom")
        ax3.xaxis.set_label_position("bottom")
        ax3.spines["bottom"].set_position(("outward", 40))
    end
    
    fig.tight_layout()#h_pad=0.0)
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(16scale, 4*N_plots*scale)
    fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)    
    
    return nothing
end

print("'plot_halostory_full'   ")












@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_hexbins(xdata, ydata; selection=:, outfile="./OUT_plot_hexbins.png"
    , xmin=" ", xmax=" ", ymin=" ", ymax=" ", title=" "
    , scale=15, gridres=" ", xscale="linear", yscale="linear", xlabel="", ylabel="", clabel=nothing, grid=false
    , lognorm=true, colorlimits=(0.8,nothing)
    , cmap="BuPu", cfunc=sum
    , binned_median = " ", plot_bval=" ", plot_line=" ", plot_line2=" ", plot_mergers=" ", plot_vline=" "
    , colormod=" "
    , calc_norm=true    # pyplottable
    )

    xd = pyplottable(xdata[selection], calc_norm=calc_norm)
    xsane = findall(x-> !isnan.(x), xd)
    yd = pyplottable(ydata[selection], calc_norm=calc_norm)
    ysane = findall(x-> !isnan.(x), yd)
    allsane = xsane[xsane .∈ Ref(Set(ysane))]

    #xmin = ifelse(xmin == " ", minimum(xd[findall(x->x .!== NaN, xd)]), xmin)
    #xmax = ifelse(xmax == " ", maximum(xd[findall(x->x .!== NaN, xd)]), xmax)
    #ymin = ifelse(xmin == " ", minimum(xd[findall(x->x .!== NaN, yd)]), ymin)
    #ymax = ifelse(xmax == " ", maximum(xd[findall(x->x .!== NaN, yd)]), ymax)
    if xmin == " "
        xmin = minimum(xd[xsane])
    end
    if xmax == " "
        xmax = maximum(xd[xsane])
    end
    if ymin == " "
        ymin = minimum(yd[ysane])
    end
    if ymax == " "
        ymax = maximum(yd[ysane])
    end
    if gridres == " "
        gridres = (50,25)
    elseif typeof(gridres) == Int
        gridres = ( gridres, Int(cld(gridres,2)) )
    end

    ### Figure
    #rc          = pyimport("matplotlib.rc")
    pltcolors   = pyimport("matplotlib.colors")
    PyPlot.matplotlib[:rc]("text", usetex=true)
    PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
    fig, ax = subplots()
    if grid
        ax.grid()
    end
    if colormod == " "
        h = ax.hexbin( xd, yd
            , gridsize=gridres
            , extent=[xmin,xmax, ymin,ymax]
            , cmap=cmap, clim=colorlimits
            , zorder=1
            , norm=ifelse(lognorm, pltcolors.LogNorm(), nothing)
            )
    else
        h = ax.hexbin( xd, yd, C = colormod
            , gridsize=gridres
            , extent=[xmin,xmax, ymin,ymax]
            , cmap=cmap, reduce_C_function=cfunc, clim=colorlimits
            , zorder=1
            , norm=ifelse(lognorm, pltcolors.LogNorm(), nothing)
            )
    end

    #println("Check")
    #a = h[:get_array]
    #println(typeof(a))
    #println(a)
    
    if plot_bval == "y"
        ax.hlines(b_disk, xmin=xmin, xmax=xmax, colors="blue")
        ax.hlines(b_ell, xmin=xmin, xmax=xmax, colors="red")
    elseif plot_bval == "x"
        ax.vlines(b_disk, ymin=ymin, ymax=ymax, colors="blue")
        ax.vlines(b_ell, ymin=ymin, ymax=ymax, colors="red")
    end
        
    if plot_mergers == "y"
        ax.hlines(1, xmin=xmin, xmax=xmax, colors="black")
        ax.hlines(log10(3) , xmin=xmin, xmax=xmax, colors="black")
        ax.hlines(0, xmin=xmin, xmax=xmax, colors="black")
    elseif plot_mergers == "x"
        ax.vlines(0, ymin=ymin, ymax=ymax, colors="black")
        ax.vlines(1, ymin=ymin, ymax=ymax, colors="black")
        ax.vlines(log10(3) , ymin=ymin, ymax=ymax, colors="black")
    end

    inboundaries = idxmatch(findcs(yd[allsane], geq=ymin, leq=ymax), findcs(xd[allsane], geq=xmin, leq=xmax))

    if binned_median == "x"
        if xscale == "log"
            medians, edges, bin_number     = stats.binned_statistic(log10.(yd[allsane][inboundaries]), xd[allsane][inboundaries], statistic="median", bins=gridres[2])
            bin_width   = edges[2] - edges[1]
            centers     = edges[2:end] .- (bin_width/2)
            pctl_68, edges_up, bin_number_up     = stats.binned_statistic(log10.(yd[allsane][inboundaries]), xd[allsane][inboundaries], statistic=pctl68upper, bins=gridres[2])
            pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(log10.(yd[allsane][inboundaries]), xd[allsane][inboundaries], statistic=pctl68lower, bins=gridres[2])
            ax.fill_betweenx(log10.(centers), pctl_32, pctl_68, color="blue", alpha=0.2)#, transform=ax.get_xaxis_transform())
            ax.plot(medians, log10.(centers), "b-", title="Median")
        elseif xscale == "linear" # tested
            medians, edges, bin_number     = stats.binned_statistic(yd[allsane][inboundaries], xd[allsane][inboundaries], statistic="median", bins=gridres[2])
            bin_width   = edges[2] - edges[1]
            centers     = edges[2:end] .- (bin_width/2)
            pctl_68, edges_up, bin_number_up     = stats.binned_statistic(yd[allsane][inboundaries], xd[allsane][inboundaries], statistic=pctl68upper, bins=gridres[2])
            pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(yd[allsane][inboundaries], xd[allsane][inboundaries], statistic=pctl68lower, bins=gridres[2])
            ax.fill_betweenx(centers, pctl_32, pctl_68, color="blue", alpha=0.2)#, transform=ax.get_xaxis_transform())
            ax.plot(medians, centers, "b-", label="Median")
        else
            error("Unknown xscale: $xscale")
        end
    elseif binned_median == "y"
        if yscale == "log"
            medians, edges, bin_number     = stats.binned_statistic(log10.(xd[allsane][inboundaries]), yd[allsane][inboundaries], statistic="median", bins=gridres[1])
            bin_width   = edges[2] - edges[1]
            centers     = edges[2:end] .- (bin_width/2)
            pctl_68, edges_up, bin_number_up     = stats.binned_statistic(log10.(xd[allsane][inboundaries]), yd[allsane][inboundaries], statistic=pctl68upper, bins=gridres[1])
            pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(log10.(xd[allsane][inboundaries]), yd[allsane][inboundaries], statistic=pctl68lower, bins=gridres[1])
            ax.fill_between(log10.(centers), pctl_32, pctl_68, color="blue", alpha=0.2)
            ax.plot(log10.(centers), medians, "b-", title="Median", linewidth=5)
        elseif yscale == "linear" #tested
            medians, edges, bin_number     = stats.binned_statistic(xd[allsane][inboundaries], yd[allsane][inboundaries], statistic="median", bins=gridres[1])
            bin_width   = edges[2] - edges[1]
            centers     = edges[2:end] .- (bin_width/2)
            pctl_68, edges_up, bin_number_up     = stats.binned_statistic(xd[allsane][inboundaries], yd[allsane][inboundaries], statistic=pctl68upper, bins=gridres[1])
            pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(xd[allsane][inboundaries], yd[allsane][inboundaries], statistic=pctl68lower, bins=gridres[1])
            ax.fill_between(centers, pctl_32, pctl_68, color="blue", alpha=0.2)
            ax.plot(centers, medians, "b-", label="Median", linewidth=5)
        else
            error("Unknown yscale: $yscale")
        end
    end

    if plot_line != " "
        ax.plot(plot_line[:,1], plot_line[:,2], "-", color="black")
    end
    if plot_line2 != " "
        ax.plot(plot_line2[:,1], plot_line2[:,2], "-", color="black")
    end
    if plot_vline != " "
        ax.axvline(plot_vline, color="black")
    end


    #ax.set_yticks([1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6])
    #rc("axes", ticksize="40")
    ax.set_xlabel(xlabel, fontsize="$(2*scale)")
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_ylabel(ylabel, fontsize="$(2*scale)")
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.tick_params(axis="both", labelsize="$(2.5*scale)")
    if title == " "
    else
        ax.set_title(title, fontsize="$(2.5*scale)")
    end
    colbar  = fig.colorbar(h, ax=ax)
    colbar.ax.tick_params(axis="y", labelsize="$(2*scale)")
    colbar.ax.set_ylabel(clabel, fontsize="$(2*scale)")
    colbar.ax.yaxis.set_label_position("left")
    
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(scale, 0.7scale)
    fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
    close("all")
    PyPlot.matplotlib[:rc]("text", usetex=false)
    
    return nothing
end

print("'plot_hexbins'   ")












@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_hexbins_nxm(xdata, ydata, selections, labels, rows, cols; 
    outfile="./OUT_plot_hexbins.png", 
    xmin=" ", xmax=" ", ymin=" ", ymax=" ",
    scale=15, gridres=" ", xscale="linear", yscale="linear", xlabel="", ylabel="", grid=false, fontsize=20,  dim_ratio=0.7, tick_size=1,
    lognorm=true, label_pos=" ",
    cmap="BuPu", cbarshift=2.0, colorlimits=(0.8,nothing), cpad=0.05,
    binned_median = " ", plot_bval=" ", plot_mergers=" ",
    calc_norm=true,    # pyplottable
    cfunc=sum, colormod=" ",
    plot_line=" ", plot_line2=" ", plot_vline=" ",
    
    )

    PyPlot.matplotlib[:rc]("text", usetex=true)
    PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
    #   Check for dimensions
    if rows*cols != length(keys(selections))
        error("Mismatch of dimensions $n x $m and datasets $(length(keys(selections))).")
    end
    if label_pos == " "
        label_pos =  Array{Nothing}(undef, rows*cols)
    end

    #   Create sanity selection
    all_sanes   = Vector{Int64}(undef, 0)
    sane        = Dict{Int64, Vector{Int64}}()
    for i in 1:rows*cols
        sane[i]     = findall(x-> !isnan.(x), pyplottable(xdata[selections[i]], calc_norm=calc_norm))[findall(x-> !isnan.(x), pyplottable(xdata[selections[i]], calc_norm=calc_norm)) .∈ Ref(Set(findall(x-> !isnan.(x), pyplottable(ydata[selections[i]], calc_norm=calc_norm))))]
        all_sanes   = vcat( all_sanes, selections[i][sane[i]] )
    end


    if xmin == " "
        xmin = minimum(xdata[all_sanes])
    end
    if xmax == " "
        xmax = maximum(xdata[all_sanes])
    end
    if ymin == " "
        ymin = minimum(ydata[all_sanes])
    end
    if ymax == " "
        ymax = maximum(ydata[all_sanes])
    end
    if gridres == " "
        gridres = (50,25)
    elseif typeof(gridres) == Int
        gridres = ( gridres, Int(cld(gridres,2)) )
    end


    ### Figure
    #rc          = pyimport("matplotlib.rc")
    pltcolors   = pyimport("matplotlib.colors")
    fig, axs    = subplots(rows, cols; sharex=true, sharey=true)

    # workaround for maximum count
    histograms  = Dict{Int64, Any}()
    n   = 0
    m   = 1
    allcounts = Vector{Int64}(undef, 0)
    for i in 1:rows*cols
        n += 1
        if n == rows+1
            n   = 1
            m  += 1
        end
        xd = xdata[selections[i]][sane[i]]
        yd = ydata[selections[i]][sane[i]]
        if colormod == " "
            c = nothing
        else
            c = colormod[selections[i]][sane[i]]
        end
        histograms[i]   = axs[n, m].hexbin( xd, yd, C=c,#colormod[selections[i]][sane[i]]
                                            gridsize=gridres, #label=labels[i],
                                            extent=[xmin,xmax, ymin,ymax], 
                                            cmap=cmap, zorder=0, #clim=(colorlimits[1], ifelse(colorlimits[2]===nothing,maximum(allcounts), colorlimits[2])),
                                            norm=ifelse(lognorm, pltcolors.LogNorm(), nothing), # inside LogNorm(): vmin=colorlimits[1], vmax=ifelse(colorlimits[2]===nothing,maximum(allcounts), colorlimits[2])
                                            reduce_C_function=cfunc, 
                                            )
        #histograms[i]   = axs[n, m].hexbin( xd, yd, 
        #                                    gridsize=gridres, #label=labels[i],
        #                                    extent=[xmin,xmax, ymin,ymax], 
        #                                    cmap=cmap, zorder=1, norm=ifelse(lognorm, pltcolors.LogNorm(), nothing), 
        #                                    )
        allcounts = vcat( allcounts, histograms[i].get_array() )
    end
    close("all")
    fig, axs    = subplots(rows, cols; sharex=true, sharey=true)

    # now the real thing
    histograms  = Dict{Int64, Any}()
    n   = 0
    m   = 1
    for i in 1:rows*cols
        n += 1
        if n == rows+1
            n   = 1
            m  += 1
        end
        xd = xdata[selections[i]][sane[i]]
        yd = ydata[selections[i]][sane[i]]
        if colormod == " "
            c = nothing
        else
            c = colormod[selections[i]][sane[i]]
        end
        histograms[i]   = axs[n, m].hexbin( xd, yd, C=c,
                                            gridsize=gridres, label=labels[i],
                                            extent=[xmin,xmax, ymin,ymax], 
                                            cmap=cmap, zorder=1, clim=(colorlimits[1], ifelse(colorlimits[2]===nothing,maximum(allcounts), colorlimits[2])),
                                            norm=ifelse(lognorm, pltcolors.LogNorm(), nothing), # inside LogNorm(): vmin=colorlimits[1], vmax=ifelse(colorlimits[2]===nothing,maximum(allcounts), colorlimits[2])
                                            reduce_C_function=cfunc, 
                                            )
        axs[n,m].set_xlim([xmin,xmax])
        axs[n,m].set_ylim([ymin,ymax])
        if grid
            axs[n,m].grid()
        end
        if plot_bval == "y"
            axs[n, m].hlines(b_disk, xmin=xmin, xmax=xmax, colors="blue")
            axs[n, m].hlines(b_ell , xmin=xmin, xmax=xmax, colors="red")
        elseif plot_bval == "x"
            axs[n, m].vlines(b_disk, ymin=ymin, ymax=ymax, colors="blue")
            axs[n, m].vlines(b_ell , ymin=ymin, ymax=ymax, colors="red")
        end
        
        if plot_mergers == "y"
            axs[n, m].hlines(1, xmin=xmin, xmax=xmax, colors="black")
            axs[n, m].hlines(log10(3) , xmin=xmin, xmax=xmax, colors="black")
            axs[n, m].hlines(0, xmin=xmin, xmax=xmax, colors="black")
        elseif plot_mergers == "x"
            axs[n, m].vlines(0, ymin=ymin, ymax=ymax, colors="black")
            axs[n, m].vlines(1, ymin=ymin, ymax=ymax, colors="black")
            axs[n, m].vlines(log10(3) , ymin=ymin, ymax=ymax, colors="black")
        end
        axs[n, m].legend(loc=label_pos[i], frameon=true, borderpad=1, handlelength=1.8, fontsize=floor(0.5*fontsize))

        inboundaries = idxmatch(findcs(yd, geq=ymin, leq=ymax), findcs(xd, geq=xmin, leq=xmax))
        if binned_median == "x"
            if xscale == "log"
                medians, edges, bin_number     = stats.binned_statistic(log10.(yd), xd, statistic="median", bins=gridres[2])
                bin_width   = edges[2] - edges[1]
                centers     = edges[2:end] .- (bin_width/2)
                pctl_68, edges_up, bin_number_up     = stats.binned_statistic(log10.(yd), xd, statistic=pctl68upper, bins=gridres[2])
                pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(log10.(yd), xd, statistic=pctl68lower, bins=gridres[2])
                axs[n, m].fill_betweenx(log10.(centers), pctl_32, pctl_68, color="blue", alpha=0.2)#, transform=ax.get_xaxis_transform())
                axs[n, m].plot(medians, log10.(centers), "b-", title="Median")
            elseif xscale == "linear" # tested
                medians, edges, bin_number     = stats.binned_statistic(yd[inboundaries], xd[inboundaries], statistic="median", bins=gridres[2])
                bin_width   = edges[2] - edges[1]
                centers     = edges[2:end] .- (bin_width/2)
                pctl_68, edges_up, bin_number_up     = stats.binned_statistic(yd[inboundaries], xd[inboundaries], statistic=pctl68upper, bins=gridres[2])
                pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(yd[inboundaries], xd[inboundaries], statistic=pctl68lower, bins=gridres[2])
                axs[n, m].fill_betweenx(centers, pctl_32, pctl_68, color="blue", alpha=0.2)#, transform=ax.get_xaxis_transform())
                axs[n, m].plot(medians, centers, "b-", label="Median")
            else
                error("Unknown xscale: $xscale")
            end
        elseif binned_median == "y"
            if yscale == "log"
                medians, edges, bin_number     = stats.binned_statistic(log10.(xd), yd, statistic="median", bins=gridres[1])
                bin_width   = edges[2] - edges[1]
                centers     = edges[2:end] .- (bin_width/2)
                pctl_68, edges_up, bin_number_up     = stats.binned_statistic(log10.(xd), yd, statistic=pctl68upper, bins=gridres[1])
                pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(log10.(xd), yd, statistic=pctl68lower, bins=gridres[1])
                axs[n, m].fill_between(log10.(centers), pctl_32, pctl_68, color="blue", alpha=0.2)
                axs[n, m].plot(log10.(centers), medians, "b-", title="Median", linewidth=5)
            elseif yscale == "linear" #tested
                medians, edges, bin_number     = stats.binned_statistic(xd[inboundaries], yd[inboundaries], statistic="median", bins=gridres[1])
                bin_width   = edges[2] - edges[1]
                centers     = edges[2:end] .- (bin_width/2)
                pctl_68, edges_up, bin_number_up     = stats.binned_statistic(xd[inboundaries], yd[inboundaries], statistic=pctl68upper, bins=gridres[1])
                pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(xd[inboundaries], yd[inboundaries], statistic=pctl68lower, bins=gridres[1])
                axs[n, m].fill_between(centers, pctl_32, pctl_68, color="blue", alpha=0.2)
                axs[n, m].plot(centers, medians, "b-", label="Median", linewidth=5)
            else
                error("Unknown yscale: $yscale")
            end
        end

        if plot_line != " "
            axs[n, m].plot(plot_line[:,1], plot_line[:,2], "-", color="black", zorder = 100)
        end
        if plot_line2 != " "
            axs[n, m].plot(plot_line2[:,1], plot_line2[:,2], "-", color="black", zorder = 100)
        end
        if plot_vline != " "
            axs[n, m].axvline(plot_vline[1], plot_vline[2], plot_vline[3], color="black",zorder=100)
        end
    end
    

    #ax.set_yticks([1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6])
    #rc("axes", ticksize="$(fontsize)")
    fig.supxlabel(xlabel, fontsize="$(fontsize)")
    fig.supylabel(ylabel, fontsize="$(fontsize)")
    #axs[3,1].set_xscale(xscale)
    #axs[3,2].set_xscale(xscale)
    #axs[3,1].set_xlim([xmin, xmax])
    #axs[3,2].set_xlim([xmin, xmax])
    #axs[1,1].set_yscale(yscale)
    #axs[1,1].set_ylabel(ylabel, fontsize="$(2*scale)")
    #axs[1,1].set_ylim([ymin, ymax])

    #ticker   = pyimport("matplotlib.ticker")
    #nbins = length(axs[3,1].get_xticklabels()) # added 
    #axs[1,1].tick_params(axis="y", labelsize="$(scale)")
    #axs[2,1].tick_params(axis="y", labelsize="$(scale)")
    ##axs[2,1].yaxis.set_major_locator(ticker.MaxNLocator(prune="upper")) # added 
    #axs[3,1].tick_params(axis="both", labelsize="$(scale)")
    #axs[3,1].yaxis.set_major_locator(ticker.MaxNLocator(prune="lower")) # added 
    #axs[3,1].xaxis.set_major_locator(ticker.MaxNLocator(prune="upper")) # added 
    for i in 1:cols
        axs[rows,i].tick_params(axis="x", labelsize="$(tick_size*scale)")
    end
    for i in 1:rows
        axs[i,1].tick_params(axis="y", labelsize="$(tick_size*scale)")
    end
    colbar  = fig.colorbar(histograms[rows*cols], ax=axs[:,cols], anchor=(cbarshift, 0.5), pad=cpad)
    colbar.ax.tick_params(axis="y", labelsize="$(tick_size*scale)")
    
    #fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(scale, dim_ratio*scale)
    fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
    close("all")

    PyPlot.matplotlib[:rc]("text", usetex=false)
    
    return histograms
end

print("'plot_hexbins_nxm'   ")















@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_scatter_nxm(xdata, ydata, selections, labels, rows, cols, colordata; 
    outfile="./OUT_plot_hexbins.png", 
    xmin=" ", xmax=" ", ymin=" ", ymax=" ",
    scale=15, gridres=" ", xscale="linear", yscale="linear", xlabel="", ylabel="", grid=false, fontsize=20, 
    lognorm=true, 
    cmap="rainbow_r", cbarshift=2.0, pointsize=5,
    binned_median = " ", plot_bval=" ",
    calc_norm=true,    # pyplottable
    )

    PyPlot.matplotlib[:rc]("text", usetex=true)
    PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
    #   Check for dimensions
    if rows*cols != length(keys(selections))
        error("Mismatch of dimensions $n x $m and datasets $(length(keys(selections))).")
    end

    #   Create sanity selection
    all_sanes   = Vector{Int64}(undef, 0)
    sane        = Dict{Int64, Vector{Int64}}()
    for i in 1:rows*cols
        sane[i]     = findall(x-> !isnan.(x), pyplottable(xdata[selections[i]], calc_norm=calc_norm))[findall(x-> !isnan.(x), pyplottable(xdata[selections[i]], calc_norm=calc_norm)) .∈ Ref(Set(findall(x-> !isnan.(x), pyplottable(ydata[selections[i]], calc_norm=calc_norm))))]
        all_sanes   = vcat( all_sanes, selections[i][sane[i]] )
    end


    if xmin == " "
        xmin = minimum(xdata[all_sanes])
    end
    if xmax == " "
        xmax = maximum(xdata[all_sanes])
    end
    if ymin == " "
        ymin = minimum(ydata[all_sanes])
    end
    if ymax == " "
        ymax = maximum(ydata[all_sanes])
    end
    if gridres == " "
        gridres = (50,25)
    elseif typeof(gridres) == Int
        gridres = ( gridres, Int(cld(gridres,2)) )
    end

    ### Figure
    #rc          = pyimport("matplotlib.rc")
    pltcolors   = pyimport("matplotlib.colors")
    fig, axs    = subplots(rows, cols; sharex=true, sharey=true)

    colr    = pyplottable(colordata)#./as09["Δlookbacktime"] ./-1000)
    maxcolr = maximum(colr[findall(x->!isnan.(x), colr)])
    #slct    = findall(x->!isnan.(x), colr)
    # now the real thing
    scatters  = Dict{Int64, Any}()
    n   = 0
    m   = 1
    for i in 1:rows*cols
        n += 1
        if n == rows+1
            n   = 1
            m  += 1
        end
        xd = xdata[selections[i]][sane[i]]
        yd = ydata[selections[i]][sane[i]]
        scatters[i]   = axs[n, m].scatter( xd, yd, 
                                        s= (pointsize/maxcolr) .* colr[selections[i]][sane[i]]    , alpha=0.6,
                                        cmap="plasma_r", c=colr[selections[i]][sane[i]])
        if plot_bval == "y"
            axs[n, m].hlines(b_disk, xmin=xmin, xmax=xmax, colors="blue")
            axs[n, m].hlines(b_ell, xmin=xmin, xmax=xmax, colors="red")
        elseif plot_bval == "x"
            axs[n, m].vlines(b_disk, ymin=ymin, ymax=ymax, colors="blue")
            axs[n, m].vlines(b_ell, ymin=ymin, ymax=ymax, colors="red")
        end
        axs[n, m].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8, fontsize=floor(0.5*fontsize))

        if binned_median == "y"
            if yscale == "log"
                medians, edges, bin_number     = stats.binned_statistic(log10.(yd), xd, statistic="median", bins=gridres[2])
                bin_width   = edges[2] - edges[1]
                centers     = edges[2:end] .- (bin_width/2)
                pctl_68, edges_up, bin_number_up     = stats.binned_statistic(log10.(yd), xd, statistic=pctl68upper, bins=gridres[2])
                pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(log10.(yd), xd, statistic=pctl68lower, bins=gridres[2])
                axs[n, m].fill_betweenx(log10.(centers), pctl_32, pctl_68, color="blue", alpha=0.2)#, transform=ax.get_xaxis_transform())
                axs[n, m].plot(medians, log10.(centers), "b-", title="Median")
            elseif yscale == "linear" # tested
                medians, edges, bin_number     = stats.binned_statistic(yd, xd, statistic="median", bins=gridres[2])
                bin_width   = edges[2] - edges[1]
                centers     = edges[2:end] .- (bin_width/2)
                pctl_68, edges_up, bin_number_up     = stats.binned_statistic(yd, xd, statistic=pctl68upper, bins=gridres[2])
                pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(yd, xd, statistic=pctl68lower, bins=gridres[2])
                axs[n, m].fill_betweenx(centers, pctl_32, pctl_68, color="blue", alpha=0.2)#, transform=ax.get_xaxis_transform())
                axs[n, m].plot(medians, centers, "b-", label="Median")
            else
                error("Unknown yscale: $yscale")
            end
        elseif binned_median == "x"
            if xscale == "log"
                medians, edges, bin_number     = stats.binned_statistic(log10.(xd), yd, statistic="median", bins=gridres[1])
                bin_width   = edges[2] - edges[1]
                centers     = edges[2:end] .- (bin_width/2)
                pctl_68, edges_up, bin_number_up     = stats.binned_statistic(log10.(xd), yd, statistic=pctl68upper, bins=gridres[1])
                pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(log10.(xd), yd, statistic=pctl68lower, bins=gridres[1])
                axs[n, m].fill_between(log10.(centers), pctl_32, pctl_68, color="blue", alpha=0.2)
                axs[n, m].plot(log10.(centers), medians, "b-", title="Median", linewidth=5)
            elseif xscale == "linear" #tested
                medians, edges, bin_number     = stats.binned_statistic(xd, yd, statistic="median", bins=gridres[1])
                bin_width   = edges[2] - edges[1]
                centers     = edges[2:end] .- (bin_width/2)
                pctl_68, edges_up, bin_number_up     = stats.binned_statistic(xd, yd, statistic=pctl68upper, bins=gridres[1])
                pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(xd, yd, statistic=pctl68lower, bins=gridres[1])
                axs[n, m].fill_between(centers, pctl_32, pctl_68, color="blue", alpha=0.2)
                axs[n, m].plot(centers, medians, "b-", label="Median", linewidth=5)
            else
                error("Unknown xscale: $xscale")
            end
        end
    end
    

    #ax.set_yticks([1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6])
    #rc("axes", ticksize="40")
    fig.supxlabel(xlabel, fontsize="$(fontsize)")
    fig.supylabel(ylabel, fontsize="$(fontsize)")
    #axs[3,1].set_xscale(xscale)
    #axs[3,2].set_xscale(xscale)
    #axs[3,1].set_xlim([xmin, xmax])
    #axs[3,2].set_xlim([xmin, xmax])
    #axs[1,1].set_yscale(yscale)
    #axs[1,1].set_ylabel(ylabel, fontsize="$(2*scale)")
    #axs[1,1].set_ylim([ymin, ymax])

    #ticker   = pyimport("matplotlib.ticker")
    #nbins = length(axs[3,1].get_xticklabels()) # added 
    #axs[1,1].tick_params(axis="y", labelsize="$(scale)")
    #axs[2,1].tick_params(axis="y", labelsize="$(scale)")
    ##axs[2,1].yaxis.set_major_locator(ticker.MaxNLocator(prune="upper")) # added 
    #axs[3,1].tick_params(axis="both", labelsize="$(scale)")
    #axs[3,1].yaxis.set_major_locator(ticker.MaxNLocator(prune="lower")) # added 
    #axs[3,1].xaxis.set_major_locator(ticker.MaxNLocator(prune="upper")) # added 
    for i in 1:cols
        axs[rows,i].tick_params(axis="x", labelsize="$(scale)")
    end
    for i in 1:rows
        axs[i,1].tick_params(axis="y", labelsize="$(scale)")
    end
    if grid
        axs.grid()
    end
    colbar  = fig.colorbar(scatters[rows*cols], ax=axs[:,cols], anchor=(cbarshift, 0.5))
    colbar.ax.tick_params(axis="y", labelsize="$(1*scale)")
    
    #fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(scale, 0.7scale)
    fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
    close("all")

    PyPlot.matplotlib[:rc]("text", usetex=false)
    
    return scatters
end

print("'plot_scatter_nxm'   ")














@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_phasespace(; snapNR = "all", ymin=" ", ymax=" ",
    felixID=" ", rootID=" ", subID=" ", pointsize=2, opact=0.9, 
    parttype="STARS", r_by_rvir=0.1, weights=" ", 
    lognorm=lognorm, scale=15, gridres=(25, 15), colorlimits=(0.8,nothing), 
    indir=current_dir_stories, simbox=current_dir_simbox, 
    outdir="./OUT_plot_phasespace", outfile="PhSp", plottype="hexbin", 
    )
    #if colorlimits == ""

    
    if snapNR == "all"
        if typeof(felixID)==Int  # simple case since tree is provided
            mc = load(joinpath(current_dir_stories, find_rootID(felixID)[2]), "merger_collection_$parttype")
            rootID = find_rootID(felixID)[1]
        elseif typeof(rootID)==Int  # simple case since tree is provided
            mc = load(joinpath(current_dir_stories, "halo_$rootID.jld"), "merger_collection_$parttype")
        else
            error("Either felixID or rootID has to be provided.")
        end

        println("$(length(mc["snapNR"])) steps:")
        for i in 1:length(mc["snapNR"])
            print("$(mc["snapNR"][i]) ")
            flush(stdout)
            snapshot    = Snapshot(simbox, mc["snapNR"][i])
            galaxyID_list   = find_merging_progenitors(mc["snapNR"][i]; rootID=rootID, path_to_halostories=indir)
            g           = Galaxy(snapshot, galaxyID_list[1])
            rvir_group  = read_galaxy_prop(get_group(g), "RVIR", :physical)
            if parttype == "STARS"
                read_halo!(g, units=:physical, props=((:stars, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
                #velocities  = norm.(eachcol(g.stars.vel))
                radii       = norm.(eachcol(g.stars.pos))
                inside      = findcs(radii, leq=r_by_rvir*rvir_group)
                velocities  = get_radial_velocity.(eachcol(g.stars.pos),eachcol(g.stars.vel))
            elseif parttype == "GAS"
                read_halo!(g, units=:physical, props=((:gas, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
                #velocities  = norm.(eachcol(g.gas.vel))
                radii       = norm.(eachcol(g.gas.pos))
                inside      = findcs(radii, leq=r_by_rvir*rvir_group)
                velocities  = get_radial_velocity.(eachcol(g.gas.pos),eachcol(g.gas.vel))
            elseif parttype == "DM"
                read_halo!(g, units=:physical, props=((:dm, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
                #velocities  = norm.(eachcol(g.dm.vel))
                radii       = norm.(eachcol(g.dm.pos))
                inside      = findcs(radii, leq=r_by_rvir*rvir_group)
                velocities  = get_radial_velocity.(eachcol(g.dm.pos),eachcol(g.dm.vel))
            else
                error("You have failed me for the last time, parttype = $parttype")
            end

            if weights == "MASS"
                weights = g.stars.mass[inside]
            end
            if opact == "mass"
                opact = g.stars.mass[inside] ./ maximum(g.stars.mass[inside])
            end
            if plottype == "scatter"
                fig, ax = subplots()
                p = ax.scatter(pyplottable(radii[inside]), 
                        pyplottable(velocities[inside]) , 
                        zorder=2   , s=pointsize .* opact    , alpha=opact, 
                        cmap="rainbow_r", c=g.stars.mass[inside]) # brg_r  
                colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as["BVAL_0"    ]), b_ell, b_disk, maximum(as["BVAL_0"    ])])
                #colbar.ax.set_yticklabels(["$(round(minimum(as["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as["BVAL_0"    ]), digits=1))"])
                colbar.ax.tick_params(axis="y")
                #ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
                ax.set_title("Halo $(rootID), Snap $(mc["snapNR"][i])")
                ax.set_xlabel("r [kpc]")
                ax.set_ylabel("v [km/s]")
                ax.grid()
                #ax.autoscale()
                fig.set_size_inches(16scale, 9scale)
                fig.tight_layout()
                fig.savefig(joinpath(outdir,outfile*"_halo$(rootID)_snap$(mc["snapNR"][i])"*".png"), bbox_inches="tight", pad_inches=.1)
            elseif plottype == "hexbin"
                plot_hexbins(radii[inside], 
                            velocities[inside], ymin=ymin, ymax=ymax, 
                            colormod=weights, scale=scale, colorlimits=colorlimits, lognorm=lognorm, gridres=gridres, 
                            outfile=joinpath(outdir,outfile*"_halo$(rootID)_snap$(mc["snapNR"][i])"*".png"), 
                            xlabel="r [kpc]", 
                            ylabel="v [km/s]", 
                            title="Halo $(rootID), Snap $(mc["snapNR"][i])")#, R = $(round(r_by_rvir*rvir_group, digits=2))")
            else
                error("You have failed me for the last time, plottype = $(plottype)")
            end
        end
    else
        if typeof(felixID)==Int  # simple case since tree is provided
            galaxyID_list   = find_merging_progenitors(snapNR; felixID=felixID, path_to_halostories=indir)
            println(galaxyID_list)
            rootID = find_rootID(felixID)[1]
        elseif typeof(rootID)==Int  # simple case since tree is provided
            galaxyID_list   = find_merging_progenitors(snapNR; rootID=rootID, path_to_halostories=indir)
            println(galaxyID_list)
        elseif typeof(subID)==Int   # here we need to find the right tree first
            galaxyID_list   = find_merging_progenitors(snapNR; subID=subID, path_to_halostories=indir)
            println(galaxyID_list)
        else
            error("Either felixID, subID or rootID has to be provided.")
        end
    
        snapshot    = Snapshot(simbox, snapNR)
        g           = Galaxy(snapshot, galaxyID_list[1])
        rvir_group  = read_galaxy_prop(get_group(g), "RVIR", :physical)
        if parttype == "STARS"
            read_halo!(g, units=:physical, props=((:stars, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
            velocities  = get_radial_velocity.(eachcol(g.stars.pos),eachcol(g.stars.vel))
            radii       = norm.(eachcol(g.stars.pos))
            inside      = findcs(radii, leq=r_by_rvir*rvir_group)
        elseif parttype == "GAS"
            read_halo!(g, units=:physical, props=((:gas, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
            velocities  = get_radial_velocity.(eachcol(g.gas.pos),eachcol(g.gas.vel))
            radii       = norm.(eachcol(g.gas.pos))
            inside      = findcs(radii, leq=r_by_rvir*rvir_group)
        elseif parttype == "DM"
            read_halo!(g, units=:physical, props=((:dm, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
            velocities  = get_radial_velocity.(eachcol(g.dm.pos),eachcol(g.dm.vel))
            radii       = norm.(eachcol(g.dm.pos))
            inside      = findcs(radii, leq=r_by_rvir*rvir_group)
        else
            error("You have failed me for the last time, parttype = $parttype")
        end

        if weights == "MASS"
            weights = g.stars.mass[inside]
        end
        if opact == "mass"
            opact = g.stars.mass[inside] ./ maximum(g.stars.mass[inside])
        end
    
        if plottype == "scatter"
            fig, ax = subplots()
            p = ax.scatter(pyplottable(radii[inside]), 
                    pyplottable(velocities[inside]) , 
                    zorder=2   , s=pointsize .* opact   , alpha=opact, 
                    cmap="rainbow_r", c=g.stars.mass[inside]) # brg_r  
            colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(as["BVAL_0"    ]), b_ell, b_disk, maximum(as["BVAL_0"    ])])
            #colbar.ax.set_yticklabels(["$(round(minimum(as["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(as["BVAL_0"    ]), digits=1))"])
            colbar.ax.tick_params(axis="y")
            #ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
            ax.set_title("Halo $(rootID), Snap $(snapNR)")
            ax.set_xlabel("r [kpc]")
            ax.set_ylabel("v [km/s]")
            ax.grid()
            #ax.autoscale()
            fig.set_size_inches(16scale, 9scale)
            fig.tight_layout()
            fig.savefig(joinpath(outdir,outfile*"_halo$(rootID)_snap$(snapNR)"*".png"), bbox_inches="tight", pad_inches=.1)
        elseif plottype == "hexbin"
            plot_hexbins(radii[inside], 
                        velocities[inside], ymin=ymin, ymax=ymax, 
                        weights=weights, scale=scale, colorlimits=colorlimits, lognorm=lognorm, gridres=gridres, 
                        outfile=joinpath(outdir,outfile*"_halo$(rootID)_snap$(snapNR)"*".png"), 
                        xlabel="r [kpc]", 
                        ylabel="v [km/s]", 
                        title="Halo $(rootID), Snap $(snapNR)")#, R = $(round(r_by_rvir*rvir_group, digits=2))")
        else
            error("You have failed me for the last time, plottype = $(plottype)")
        end
    end
    
    return nothing
end

print("'plot_phasespace'   ")









include("/home/moon/sfortune/spinevo/pkg/plotting/plot_hist_nxm.jl")

println()