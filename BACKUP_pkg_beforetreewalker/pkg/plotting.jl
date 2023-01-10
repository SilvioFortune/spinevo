println("   Plotting...")
using PyPlot
using ColorSchemes
using JSServe
using WGLMakie
WGLMakie.activate!()

include("/home/moon/sfortune/spinevo/pkg/plot_group_collection.jl")


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
    return surface(-1..1, -1..1, z, colormap = :winter, shading = false)
end

print("'testconnection'   ")


@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_halostory_arrows(; outfile="./OUT_plot_halostory_arrows.png", felixID=" ", lastID=" "
    , parttype="STARS", snap=false
    , scale=1, plotwidth=16, plotheight=3, gridspace=0.6 #Gyr
    , indir=current_dir_jld
    )

    halo_filestring = find_halo_file(indir=indir,felixID=felixID,lastID=lastID)

    merger_collection = load(joinpath(indir, halo_filestring), string("merger_collection_", parttype))
    
    lbt     = Array{Float64}(undef, 1)
    lbt[1]  = merger_collection["LOOKBACKTIME"][1]
    rs      = Array{Float64}(undef, 1)
    rs[1]   = merger_collection["REDSHIFT"][1]
    snaps    = Array{Int64}(undef, 1)
    snaps[1] = merger_collection["SNAP"][1]
    for i in 2:length(merger_collection["SNAP"])-1
        if lbt[end] - merger_collection["LOOKBACKTIME"][i] > gridspace
            snaps    = vcat( snaps, merger_collection["SNAP"        ][i] )
            lbt     = vcat( lbt , merger_collection["LOOKBACKTIME"][i] )
            rs      = vcat( rs  , merger_collection["REDSHIFT"    ][i] )
        end
    end
    snaps    = vcat( snaps, merger_collection["SNAP"        ][end] )
    lbt     = vcat( lbt , merger_collection["LOOKBACKTIME"][end] )
    rs      = vcat( rs  , merger_collection["REDSHIFT"    ][end] )
    
    
    J_norm      = vcat( transpose(norm.(eachcol(merger_collection["J_main"]))), transpose(norm.(eachcol(merger_collection["J_main"]))), transpose(norm.(eachcol(merger_collection["J_main"]))) )
    
    proj_vec    = (merger_collection["J_main"]./J_norm)[1,:]
    x_vec       = (merger_collection["J_main"]./J_norm)[2,:]
    y_vec       = (merger_collection["J_main"]./J_norm)[3,:]
    
    pltcm       = pyimport("matplotlib.cm")
    pltcolors   = pyimport("matplotlib.colors")
    colormap    = pltcm.get_cmap(name="coolwarm")
    
    fig, ax = subplots()
    if typeof(lastID) == Int
        ax.set_title("Sub Halo $lastID")
    elseif typeof(felixID) == Int
        ax.set_title("Sub Halo $(find_lastID(felixID)[1])")
    end
    ax.set_ylim(-1, 1)
    ax.set_yticks([-1,0,1])
    ax.set_xlim( minimum(merger_collection["LOOKBACKTIME"])-0.2-2*maximum([0, x_vec[end]]), maximum(merger_collection["LOOKBACKTIME"])+0.2-2*minimum([0, x_vec[1]]) )
    ax.set_xticks( lbt )
    ax.set_xlabel("Lookback Time [Gyr]")
    
    ax.quiver(merger_collection["LOOKBACKTIME"], zeros(length(merger_collection["LOOKBACKTIME"])), 
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
function plot_halostory_data(; outfile="./OUT_plot_halostory_data.png", felixID=" ", lastID=" "
    , parttype="STARS", snap=false
    , scale=0.7, gridspace=0.6 #Gyr
    , indir=current_dir_jld
    )
    
    mcs = load(joinpath(indir, find_halo_file(indir=indir,felixID=felixID,lastID=lastID)), "merger_collection_$parttype")
    
    lbt     = Array{Float64}(undef, 1)
    lbt[1]  = mcs["LOOKBACKTIME"][1]
    rs      = Array{Float64}(undef, 1)
    rs[1]   = mcs["REDSHIFT"][1]
    snaps    = Array{Int64}(undef, 1)
    snaps[1] = mcs["SNAP"][1]
    for i in 2:length(mcs["SNAP"])-1
        if lbt[end] - mcs["LOOKBACKTIME"][i] > gridspace
            snaps    = vcat( snaps, mcs["SNAP"        ][i] )
            lbt     = vcat( lbt , mcs["LOOKBACKTIME"][i] )
            rs      = vcat( rs  , mcs["REDSHIFT"    ][i] )
        end
    end
    snaps    = vcat( snaps, mcs["SNAP"        ][end] )
    lbt     = vcat( lbt , mcs["LOOKBACKTIME"][end] )
    rs      = vcat( rs  , mcs["REDSHIFT"    ][end] )
    
    j_spherical = missings(Float64, 3, 0)
    for i in 1:length(mcs["j_main"][1,:])
            j_spherical = hcat( j_spherical, cartesian_to_spherical(mcs["j_main"][:,i]) )
    end
    
    ### Figure
    N_plots = 5
    fig, ax = subplots(5, sharex=true)
    fig.subplots_adjust(hspace=0)
    if typeof(lastID) == Int
        ax[1].set_title("Sub Halo $lastID")
    elseif typeof(felixID) == Int
        ax[1].set_title("Sub Halo $(find_lastID(felixID)[1])")
    end

    ax[1].set_xticks( lbt )
    
    ax[1].plot(mcs["LOOKBACKTIME"], replace(mcs["ϕ_flip"], missing => NaN), "b-", lw=5, label="ϕ_flip", alpha=1, zorder=3)
    ax[1].plot(mcs["LOOKBACKTIME"], replace(j_spherical[2,:], missing => NaN), "g-", lw=2, label="θ", alpha=1, zorder=1)
    ax[1].plot(mcs["LOOKBACKTIME"], replace(j_spherical[3,:], missing => NaN), "c-", lw=2, label="ϕ", alpha=1, zorder=2)
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
    
    ax[2].plot(mcs["LOOKBACKTIME"], mcs["BVAL"] .+ 0.5*log10.(1 .+ mcs["REDSHIFT"]), "-", lw=2, label="b-value converted to z=0", color="green", alpha=1, zorder=2)
    ax[2].grid()
    ax[2].axhline(-4.357,color="blue", label="Disks") # disks
    ax[2].axhline(-4.732,color="red", label="Ellipticals") # ellipticals
    ax[2].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
    
    ax[3].plot(mcs["LOOKBACKTIME"], norm.(eachcol(replace(mcs["j_main"], missing => 0.))), "-", lw=2, label="| j_$parttype |", color="navy", alpha=1, zorder=2)
    ax[3].set_ylabel("[kpc km/s]")
    ax[3].grid()
    ax[3].axhline(color="black")
    ax[3].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
    
    ax[4].plot(mcs["LOOKBACKTIME"], norm.(eachcol(replace(mcs["δj_main"], missing => 0.))), "-", lw=2, label="| Δj_$parttype |", color="navy", alpha=1, zorder=2)
    ax[4].set_ylabel("[kpc km/s]")
    ax[4].grid()
    ax[4].axhline(color="black")
    ax[4].legend(loc="upper left", frameon=true, borderpad=1, handlelength=1.8)
    
    ax[5].plot(mcs["LOOKBACKTIME"], mcs["M2_felix"], "-", lw=2, label="M2_$parttype", color="gold", alpha=1, zorder=4)
    ax[5].plot(mcs["LOOKBACKTIME"], mcs["M_felix"], "-", lw=2, label="M_$parttype", color="green", alpha=1, zorder=7)
    ax[5].plot(mcs["LOOKBACKTIME"], mcs["M2_felix"], ".", ms=7, color="gold", alpha=0.5, zorder=8)
    ax[5].plot(mcs["LOOKBACKTIME"], mcs["M_felix"],  ".", ms=7, color="green", alpha=0.5, zorder=9)
    ax[5].plot(mcs["Merger_Map"][7,:], mcs["Merger_Map"][2,:], ".", ms=7, color="gold", alpha=0.5, zorder=8)
    ax[5].plot(mcs["Merger_Map"][7,:], mcs["Merger_Map"][1,:], ".", ms=7, color="green", alpha=0.5, zorder=9)
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
function plot_hexbins(xdata, ydata; selection=:, outfile="./OUT_plot_hexbins.png"
    , xmin=" ", xmax=" ", ymin=" ", ymax=" ", title=" "
    , scale=15, gridres=" ", xscale="linear", yscale="linear", xlabel="", ylabel="", grid=false
    , lognorm=true, colorlimits=(0.8,nothing)
    , cmap="BuPu"
    , binned_median = " ", plot_bval=" "
    , weights=" "
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
    fig, ax = subplots()
    if weights == " "
        h = ax.hexbin( xd, yd
            , gridsize=gridres
            , extent=[xmin,xmax, ymin,ymax]
            , cmap=cmap, clim=colorlimits
            , zorder=1
            , norm=ifelse(lognorm, pltcolors.LogNorm(), nothing)
            )
    else
        h = ax.hexbin( xd, yd, C = weights
            , gridsize=gridres
            , extent=[xmin,xmax, ymin,ymax]
            , cmap=cmap#, clim=(0.8,nothing)
            , zorder=1
            , norm=ifelse(lognorm, pltcolors.LogNorm(), nothing)
            )
    end

    #println("Check")
    #a = h[:get_array]
    #println(typeof(a))
    #println(a)
    
    if plot_bval == "y"
        b_disk  = -4.357
        b_ell   = -4.732
        ax.hlines(b_disk, xmin=xmin, xmax=xmax, colors="blue")
        ax.hlines(b_ell, xmin=xmin, xmax=xmax, colors="red")
    elseif plot_bval == "x"
        b_disk  = -4.357
        b_ell   = -4.732
        ax.hlines(b_disk, ymin=ymin, ymax=ymax, colors="blue")
        ax.hlines(b_ell, ymin=ymin, ymax=ymax, colors="red")
    end


    if binned_median == "y"
        if yscale == "log"
            medians, edges, bin_number     = stats.binned_statistic(log10.(yd[allsane]), xd[allsane], statistic="median", bins=gridres[2])
            bin_width   = edges[2] - edges[1]
            centers     = edges[2:end] .- (bin_width/2)
            pctl_68, edges_up, bin_number_up     = stats.binned_statistic(log10.(yd[allsane]), xd[allsane], statistic=pctl68upper, bins=gridres[2])
            pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(log10.(yd[allsane]), xd[allsane], statistic=pctl68lower, bins=gridres[2])
            ax.fill_betweenx(log10.(centers), pctl_32, pctl_68, color="blue", alpha=0.2)#, transform=ax.get_xaxis_transform())
            ax.plot(medians, log10.(centers), "b-", title="Median")
        elseif yscale == "linear" # tested
            medians, edges, bin_number     = stats.binned_statistic(yd[allsane], xd[allsane], statistic="median", bins=gridres[2])
            bin_width   = edges[2] - edges[1]
            centers     = edges[2:end] .- (bin_width/2)
            pctl_68, edges_up, bin_number_up     = stats.binned_statistic(yd[allsane], xd[allsane], statistic=pctl68upper, bins=gridres[2])
            pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(yd[allsane], xd[allsane], statistic=pctl68lower, bins=gridres[2])
            ax.fill_betweenx(centers, pctl_32, pctl_68, color="blue", alpha=0.2)#, transform=ax.get_xaxis_transform())
            ax.plot(medians, centers, "b-", label="Median")
        else
            error("Unknown yscale: $yscale")
        end
    elseif binned_median == "x"
        if xscale == "log"
            medians, edges, bin_number     = stats.binned_statistic(log10.(xd[allsane]), yd[allsane], statistic="median", bins=gridres[1])
            bin_width   = edges[2] - edges[1]
            centers     = edges[2:end] .- (bin_width/2)
            pctl_68, edges_up, bin_number_up     = stats.binned_statistic(log10.(xd[allsane]), yd[allsane], statistic=pctl68upper, bins=gridres[1])
            pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(log10.(xd[allsane]), yd[allsane], statistic=pctl68lower, bins=gridres[1])
            ax.fill_between(log10.(centers), pctl_32, pctl_68, color="blue", alpha=0.2)
            ax.plot(log10.(centers), medians, "b-", title="Median", linewidth=5)
        elseif xscale == "linear" #tested
            medians, edges, bin_number     = stats.binned_statistic(xd[allsane], yd[allsane], statistic="median", bins=gridres[1])
            bin_width   = edges[2] - edges[1]
            centers     = edges[2:end] .- (bin_width/2)
            pctl_68, edges_up, bin_number_up     = stats.binned_statistic(xd[allsane], yd[allsane], statistic=pctl68upper, bins=gridres[1])
            pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(xd[allsane], yd[allsane], statistic=pctl68lower, bins=gridres[1])
            ax.fill_between(centers, pctl_32, pctl_68, color="blue", alpha=0.2)
            ax.plot(centers, medians, "b-", label="Median", linewidth=5)
        else
            error("Unknown yscale: $yscale")
        end
    end


    #ax.set_yticks([1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6])
    PyPlot.matplotlib[:rc]("text", usetex=true)
    PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
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
    if grid
        ax.grid()
    end
    colbar  = fig.colorbar(h, ax=ax)
    colbar.ax.tick_params(axis="y", labelsize="$(2*scale)")
    
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
function plot_hexbins_3x2(xdata, ydata; selection11=:, selection12=:, selection21=:, selection22=:, selection31=:, selection32=:
    , title_col1=" ", title_col2=" ", title_row1=" ", title_row2=" ", title_row3=" "
    , outfile="./OUT_plot_hexbins.png"
    , xmin=" ", xmax=" ", ymin=" ", ymax=" ", colormin=0.1, colormax=length(xdata)/6
    , scale=15, gridres=" ", xscale="linear", yscale="linear", xlabel="", ylabel="", grid=false
    , lognorm=true
    , cmap="BuPu"
    , binned_median = " "
    , calc_norm=true    # pyplottable
    )

    xd11 = pyplottable(xdata[selection11], calc_norm=calc_norm)
    yd11 = pyplottable(ydata[selection11], calc_norm=calc_norm)
    xd12 = pyplottable(xdata[selection12], calc_norm=calc_norm)
    yd12 = pyplottable(ydata[selection12], calc_norm=calc_norm)
    xd21 = pyplottable(xdata[selection21], calc_norm=calc_norm)
    yd21 = pyplottable(ydata[selection21], calc_norm=calc_norm)
    xd22 = pyplottable(xdata[selection22], calc_norm=calc_norm)
    yd22 = pyplottable(ydata[selection22], calc_norm=calc_norm)
    xd31 = pyplottable(xdata[selection31], calc_norm=calc_norm)
    yd31 = pyplottable(ydata[selection31], calc_norm=calc_norm)
    xd32 = pyplottable(xdata[selection32], calc_norm=calc_norm)
    yd32 = pyplottable(ydata[selection32], calc_norm=calc_norm)
    x11sane = findall(x-> !isnan.(x), xd11)
    y11sane = findall(x-> !isnan.(x), yd11)
    x12sane = findall(x-> !isnan.(x), xd12)
    y12sane = findall(x-> !isnan.(x), yd12)
    x21sane = findall(x-> !isnan.(x), xd21)
    y21sane = findall(x-> !isnan.(x), yd21)
    x22sane = findall(x-> !isnan.(x), xd22)
    y22sane = findall(x-> !isnan.(x), yd22)
    x31sane = findall(x-> !isnan.(x), xd31)
    y31sane = findall(x-> !isnan.(x), yd31)
    x32sane = findall(x-> !isnan.(x), xd32)
    y32sane = findall(x-> !isnan.(x), yd32)
    #x11sane = findall(x-> !isnan.(x), pyplottable(@view(xdata[selection11]), calc_norm=calc_norm))
    #x12sane = findall(x-> !isnan.(x), pyplottable(xdata[selection12], calc_norm=calc_norm))
    #x21sane = findall(x-> !isnan.(x), pyplottable(xdata[selection21], calc_norm=calc_norm))
    #x22sane = findall(x-> !isnan.(x), pyplottable(xdata[selection22], calc_norm=calc_norm))
    #x31sane = findall(x-> !isnan.(x), pyplottable(xdata[selection31], calc_norm=calc_norm))
    #x32sane = findall(x-> !isnan.(x), pyplottable(xdata[selection32], calc_norm=calc_norm))
    #y11sane = findall(x-> !isnan.(x), yd11)
    #y11sane = findall(x-> !isnan.(x), pyplottable(@view(ydata[selection11]), calc_norm=calc_norm))
    #y12sane = findall(x-> !isnan.(x), pyplottable(ydata[selection12], calc_norm=calc_norm))
    #y21sane = findall(x-> !isnan.(x), pyplottable(ydata[selection21], calc_norm=calc_norm))
    #y22sane = findall(x-> !isnan.(x), pyplottable(ydata[selection22], calc_norm=calc_norm))
    #y31sane = findall(x-> !isnan.(x), pyplottable(ydata[selection31], calc_norm=calc_norm))
    #y32sane = findall(x-> !isnan.(x), pyplottable(ydata[selection32], calc_norm=calc_norm))
    sane11  = x11sane[x11sane .∈ Ref(Set(y11sane))]
    sane12  = x12sane[x12sane .∈ Ref(Set(y12sane))]
    sane21  = x21sane[x21sane .∈ Ref(Set(y21sane))]
    sane22  = x22sane[x22sane .∈ Ref(Set(y22sane))]
    sane31  = x31sane[x31sane .∈ Ref(Set(y31sane))]
    sane32  = x32sane[x32sane .∈ Ref(Set(y32sane))]
    #xsane = findall(x-> !isnan.(x), pyplottable(@view(xdata[selection]), calc_norm=calc_norm))
    #ysane = findall(x-> !isnan.(x), pyplottable(@view(ydata[selection]), calc_norm=calc_norm))
    #sane = xsane[xsane .∈ Ref(Set(ysane))]
    if xmin == " "
        xmin = minimum(xdata[sane])
    end
    if xmax == " "
        xmax = maximum(xdata[sane])
    end
    if ymin == " "
        ymin = minimum(ydata[sane])
    end
    if ymax == " "
        ymax = maximum(ydata[sane])
    end
    if gridres == " "
        gridres = (50,25)
    elseif typeof(gridres) == Int
        gridres = ( gridres, Int(cld(gridres,2)) )
    end

    ### Figure
    #rc          = pyimport("matplotlib.rc")
    pltcolors   = pyimport("matplotlib.colors")
    fig, axs = subplots(3,2; sharex=true, sharey=true)

    #fig = plt.figure(figsize=(6.5,4.5), constrained_layout=true)

    #gs = fig.add_gridspec(3, 2, hspace=0, wspace=0)
    #ax11 = fig.add_subplot(gs[1, 1])
    #ax12 = fig.add_subplot(gs[1, 2])
    #ax21 = fig.add_subplot(gs[2, 1])
    #ax12 = fig.add_subplot(gs[2, 2])
    #ax31 = fig.add_subplot(gs[3, 1])
    #ax32 = fig.add_subplot(gs[3, 2])

    #(ax11, ax12), (ax21, ax22), (ax31, ax32) = gs.subplots(sharex="col", sharey="row")
    #fig, ax = subplots(3,2, sharex="col", sharey="row")
    #fig.subplots_adjust(hspace=0, wspace=0)

    #ax11 = fig.add_subplot(axs[1, 1])
    #ax12 = fig.add_subplot(axs[1, 2])
    #ax21 = fig.add_subplot(axs[2, 1])
    #ax12 = fig.add_subplot(axs[2, 2])
    #ax31 = fig.add_subplot(axs[3, 1])
    #ax32 = fig.add_subplot(axs[3, 2])



    h11 = axs[1, 1].hexbin( xd11[sane11], yd11[sane11]
        , gridsize=gridres
        , extent=[xmin,xmax, ymin,ymax]
        , cmap=cmap, zorder=1, norm=ifelse(lognorm, pltcolors.LogNorm(), nothing)
        )
    h12 = axs[1, 2].hexbin( xd12[sane12], yd12[sane12]
        , gridsize=gridres
        , extent=[xmin,xmax, ymin,ymax]
        , cmap=cmap, zorder=1, norm=ifelse(lognorm, pltcolors.LogNorm(), nothing)
        )
    h21 = axs[2, 1].hexbin( xd21[sane21], yd21[sane21]
        , gridsize=gridres
        , extent=[xmin,xmax, ymin,ymax]
        , cmap=cmap, zorder=1, norm=ifelse(lognorm, pltcolors.LogNorm(), nothing)
        )
    h22 = axs[2, 2].hexbin( xd22[sane22], yd22[sane22]
        , gridsize=gridres
        , extent=[xmin,xmax, ymin,ymax]
        , cmap=cmap, zorder=1, norm=ifelse(lognorm, pltcolors.LogNorm(), nothing)
        )
    h31 = axs[3, 1].hexbin( xd31[sane31], yd31[sane31]
        , gridsize=gridres
        , extent=[xmin,xmax, ymin,ymax]
        , cmap=cmap, zorder=1, norm=ifelse(lognorm, pltcolors.LogNorm(), nothing)
        )
    h32 = axs[3, 2].hexbin( xd32[sane32], yd32[sane32]
        , gridsize=gridres
        , extent=[xmin,xmax, ymin,ymax]
        , cmap=cmap, zorder=1, norm=ifelse(lognorm, pltcolors.LogNorm(), nothing)
        )
    
    if binned_median == "y"
        if yscale == "log"
            medians, edges, bin_number     = stats.binned_statistic(log10.(yd[allsane]), xd[allsane], statistic="median", bins=gridres[2])
            bin_width   = edges[2] - edges[1]
            centers     = edges[2:end] .- (bin_width/2)
            pctl_68, edges_up, bin_number_up     = stats.binned_statistic(log10.(yd[allsane]), xd[allsane], statistic=pctl68upper, bins=gridres[2])
            pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(log10.(yd[allsane]), xd[allsane], statistic=pctl68lower, bins=gridres[2])
            ax.fill_betweenx(log10.(centers), pctl_32, pctl_68, color="blue", alpha=0.2)#, transform=ax.get_xaxis_transform())
            ax.plot(medians, log10.(centers), "b-", title="Median")
        elseif yscale == "linear" # tested
            medians, edges, bin_number     = stats.binned_statistic(yd[allsane], xd[allsane], statistic="median", bins=gridres[2])
            bin_width   = edges[2] - edges[1]
            centers     = edges[2:end] .- (bin_width/2)
            pctl_68, edges_up, bin_number_up     = stats.binned_statistic(yd[allsane], xd[allsane], statistic=pctl68upper, bins=gridres[2])
            pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(yd[allsane], xd[allsane], statistic=pctl68lower, bins=gridres[2])
            ax.fill_betweenx(centers, pctl_32, pctl_68, color="blue", alpha=0.2)#, transform=ax.get_xaxis_transform())
            ax.plot(medians, centers, "b-", label="Median")
        else
            error("Unknown yscale: $yscale")
        end
    elseif binned_median == "x"
        if xscale == "log"
            medians, edges, bin_number     = stats.binned_statistic(log10.(xd[allsane]), yd[allsane], statistic="median", bins=gridres[1])
            bin_width   = edges[2] - edges[1]
            centers     = edges[2:end] .- (bin_width/2)
            pctl_68, edges_up, bin_number_up     = stats.binned_statistic(log10.(xd[allsane]), yd[allsane], statistic=pctl68upper, bins=gridres[1])
            pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(log10.(xd[allsane]), yd[allsane], statistic=pctl68lower, bins=gridres[1])
            ax.fill_between(log10.(centers), pctl_32, pctl_68, color="blue", alpha=0.2)
            ax.plot(log10.(centers), medians, "b-", title="Median", linewidth=5)
        elseif xscale == "linear" #tested
            medians, edges, bin_number     = stats.binned_statistic(xd[allsane], yd[allsane], statistic="median", bins=gridres[1])
            medians, edges, bin_number     = stats.binned_statistic(xd[allsane], yd[allsane], statistic="median", bins=gridres[1])
            medians, edges, bin_number     = stats.binned_statistic(xd[allsane], yd[allsane], statistic="median", bins=gridres[1])
            medians, edges, bin_number     = stats.binned_statistic(xd[allsane], yd[allsane], statistic="median", bins=gridres[1])
            medians, edges, bin_number     = stats.binned_statistic(xd[allsane], yd[allsane], statistic="median", bins=gridres[1])
            bin_width   = edges[2] - edges[1]
            centers     = edges[2:end] .- (bin_width/2)
            pctl_68, edges_up, bin_number_up     = stats.binned_statistic(xd[allsane], yd[allsane], statistic=pctl68upper, bins=gridres[1])
            pctl_32, edges_lo, bin_number_lo     = stats.binned_statistic(xd[allsane], yd[allsane], statistic=pctl68lower, bins=gridres[1])
            ax.fill_between(centers, pctl_32, pctl_68, color="blue", alpha=0.2)
            ax.plot(centers, medians, "b-", label="Median", linewidth=5)
        else
            error("Unknown xscale: $xscale")
        end
    end


    #ax.set_yticks([1e-2,1e-1,1e0,1e1,1e2,1e3,1e4,1e5,1e6])
    PyPlot.matplotlib[:rc]("text", usetex=true)
    PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
    #rc("axes", ticksize="40")
    axs[3,1].set_xlabel(xlabel, fontsize="$(2*scale)")
    axs[3,2].set_xlabel(xlabel, fontsize="$(2*scale)")
    axs[3,1].set_xscale(xscale)
    axs[3,2].set_xscale(xscale)
    axs[3,1].set_xlim([xmin, xmax])
    axs[3,2].set_xlim([xmin, xmax])
    axs[1,1].set_yscale(yscale)
    axs[1,1].set_ylabel(ylabel, fontsize="$(2*scale)")
    axs[1,1].set_ylim([ymin, ymax])

    ticker   = pyimport("matplotlib.ticker")
    #nbins = length(axs[3,1].get_xticklabels()) # added 
    axs[1,1].tick_params(axis="y", labelsize="$(scale)")
    axs[2,1].tick_params(axis="y", labelsize="$(scale)")
    #axs[2,1].yaxis.set_major_locator(ticker.MaxNLocator(prune="upper")) # added 
    axs[3,1].tick_params(axis="both", labelsize="$(scale)")
    axs[3,1].yaxis.set_major_locator(ticker.MaxNLocator(prune="lower")) # added 
    axs[3,1].xaxis.set_major_locator(ticker.MaxNLocator(prune="upper")) # added 
    axs[3,2].tick_params(axis="x", labelsize="$(scale)")
    #axs.tick_params(axis="both", labelsize="$(2.5*scale)")
    #axs.tick_params(axis="both", labelsize="$(2.5*scale)")
    #axs.tick_params(axis="both", labelsize="$(2.5*scale)")
    #axs.tick_params(axis="both", labelsize="$(2.5*scale)")
    if title_col1 != " "
        axs[1,1].set_title(title_col1, fontsize="$(2.5*scale)")
    end
    if title_col2 != " "
        axs[1,2].set_title(title_col2, fontsize="$(2.5*scale)")
    end
    if grid
        axs.grid()
    end
    #colbar  = fig.colorbar(h, ax=ax)
    #colbar.ax.tick_params(axis="y", labelsize="$(2*scale)")
    
    fig.tight_layout()
    fig.subplots_adjust(wspace=0, hspace=0)
    fig.set_size_inches(scale, 0.7scale)
    fig.savefig(outfile, bbox_inches="tight", pad_inches=.1)
    close("all")

    PyPlot.matplotlib[:rc]("text", usetex=false)
    
    return nothing
end

print("'plot_hexbins_3x2'   ")














@doc """
DESCRIPTION:\n
INPUT:\n
OUTPUT:\n
""" ->
function plot_phasespace(; snapNR = "all"
    , felixID=" ", lastID=" ", subID=" "
    , parttype="STARS", r_by_rvir=0.1, weights=" "
    , lognorm=lognorm, scale=15, gridres=(25, 15), colorlimits=(0.8,nothing)
    , indir=current_dir_jld, simbox=current_dir_simbox
    , outdir="./OUT_plot_phasespace", outfile="OUT"
    )

    
    if snapNR == "all"
        if typeof(felixID)==Int  # simple case since tree is provided
            mc = load(joinpath(current_dir_jld, find_lastID(felixID)[2]), "merger_collection_$parttype")
            lastID = find_lastID(felixID)[1]
        elseif typeof(lastID)==Int  # simple case since tree is provided
            mc = load(joinpath(current_dir_jld, find_felixID(lastID)[2]), "merger_collection_$parttype")
        else
            error("Either felixID or lastID has to be provided.")
        end

        println("$(length(mc["SNAP"])) steps:")
        for i in 1:length(mc["SNAP"])
            print("$i ")
            flush(stdout)
            snapshot    = Snapshot(simbox, mc["SNAP"][i])
            galaxyID_list   = find_merging_progenitors(mc["SNAP"][i]; lastID=lastID, path_to_halostories=indir)
            g           = Galaxy(snapshot, galaxyID_list[1])
            rvir_group  = read_galaxy_prop(get_group(g), "RVIR", :physical)
            if parttype == "STARS"
                read_halo!(g, units=:physical, props=((:stars, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
                velocities  = norm.(eachcol(g.stars.vel))
                radii       = norm.(eachcol(g.stars.pos))
            elseif parttype == "GAS"
                read_halo!(g, units=:physical, props=((:gas, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
                velocities  = norm.(eachcol(g.gas.vel))
                radii       = norm.(eachcol(g.gas.pos))
            elseif parttype == "DM"
                read_halo!(g, units=:physical, props=((:dm, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
                velocities  = norm.(eachcol(g.dm.vel))
                radii       = norm.(eachcol(g.dm.pos))
            else
                error("You have failed me for the last time, parttype = $parttype")
            end

            if weights == "MASS"
                weights = g.stars.mass
            end
        
            plot_hexbins(radii, 
                        velocities, 
                        weights=weights, scale=scale, colorlimits=colorlimits, lognorm=lognorm, gridres=gridres, 
                        outfile=joinpath(outdir,outfile*"_halo$(lastID)_snap$(mc["SNAP"][i])"*".png"), 
                        xlabel="r [kpc]", 
                        ylabel="v [km/s]", 
                        title="Halo $(lastID), Snap $(snapNR)")#, R = $(round(r_by_rvir*rvir_group, digits=2))")
        end
    else
        if typeof(felixID)==Int  # simple case since tree is provided
            galaxyID_list   = find_merging_progenitors(snapNR; felixID=felixID, path_to_halostories=indir)
            println(galaxyID_list)
            lastID = find_lastID(felixID)[1]
        elseif typeof(lastID)==Int  # simple case since tree is provided
            galaxyID_list   = find_merging_progenitors(snapNR; lastID=lastID, path_to_halostories=indir)
            println(galaxyID_list)
        elseif typeof(subID)==Int   # here we need to find the right tree first
            galaxyID_list   = find_merging_progenitors(snapNR; subID=subID, path_to_halostories=indir)
            println(galaxyID_list)
        else
            error("Either felixID, subID or lastID has to be provided.")
        end
    
        snapshot    = Snapshot(simbox, snapNR)
        g           = Galaxy(snapshot, galaxyID_list[1])
        rvir_group  = read_galaxy_prop(get_group(g), "RVIR", :physical)
        if parttype == "STARS"
            read_halo!(g, units=:physical, props=((:stars, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
            velocities  = norm.(eachcol(g.stars.vel))
            radii       = norm.(eachcol(g.stars.pos))
        elseif parttype == "GAS"
            read_halo!(g, units=:physical, props=((:gas, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
            velocities  = norm.(eachcol(g.gas.vel))
            radii       = norm.(eachcol(g.gas.pos))
        elseif parttype == "DM"
            read_halo!(g, units=:physical, props=((:dm, ["POS", "VEL", "MASS", "ID"]),), radius_units=:physical, radius=r_by_rvir*rvir_group)
            velocities  = norm.(eachcol(g.dm.vel))
            radii       = norm.(eachcol(g.dm.pos))
        else
            error("You have failed me for the last time, parttype = $parttype")
        end

        if weights == "MASS"
            weights = g.stars.mass
        end
    
        plot_hexbins(radii, 
                    velocities, 
                    weights=weights, scale=scale, colorlimits=colorlimits, lognorm=lognorm, gridres=gridres, 
                    outfile=joinpath(outdir,outfile*"_halo$(lastID)_snap$(snapNR)"*".png"), 
                    xlabel="r [kpc]", 
                    ylabel="v [km/s]", 
                    title="Halo $(lastID), Snap $(snapNR), R = $(round(r_by_rvir*rvir_group, digits=2))")
    end
    
    return nothing
end

print("'plot_phasespace'   ")

println()