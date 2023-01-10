
###############################################################################
# mergers

# 14, 15:17, 18:20, Immediate, imm
# 23, 24:26, 27:29, Earlier, earl
# 03, 0.3 Gyr, Instant
# 09, 1 Gyr, Long-term

as          = [as03, as09]
csw         = [csw03, csw09]
as_label1   = ["03", "09"]
as_label2   = ["0.3 Gyr", "1 Gyr"]
as_label3   = ["Instant", "Long-term"]

mt          = [Int.(LinRange(14,20,7)), Int.(LinRange(23,29,7))]
mt_label1   = ["imm", "earl"]
mt_label2   = ["Immediate", "Earlier"]

for (ias, ia1, ia2, ia3, ic) in zip(as, as_label1, as_label2, as_label3, csw)
    for (imt, im1, im2) in zip(mt, mt_label1, mt_label2)
        println("       Plotting $(ia3) + $(im2)")
        mmi = mergermap_indices(ias)
        delJ = convert( Array{Float64,2}, ias["ΔJ_main"][:,idxmatch(mmi["main"], ic[:transition] )] ) # all entries with mergers and matching ic[:transition]
        J_0  = convert( Array{Float64,2}, ias["J_main"][:,idxmatch(mmi["main"], ic[:transition] )] .- ias["ΔJ_main"][:,idxmatch(mmi["main"], ic[:transition] )] )
        J_1  = convert( Array{Float64,2}, ias["J_main"][:,idxmatch(mmi["main"], ic[:transition] )])
        idx = 1:length(mmi["main"])
        idx = idx[mmi["main"] .∈ Ref(Set(ic[:transition]))] # selection for merger map indices
        Jom  = missings(Float64, 3, 0)
        JoMM = missings(Float64, 3, 0)
        Jim  = missings(Float64, 3, 0)
        JiMM = missings(Float64, 3, 0)
        Jsm  = missings(Float64, 3, 0)
        JsMM = missings(Float64, 3, 0)
        idx_missed = Vector{Int64}(undef, 0) # in merger map
        idx_caught = Vector{Int64}(undef, 0) # in merger map
        for i in idx # iterate over y entries
                JoMM = hcat( JoMM,  ias["merger map"][imt[2:4],mmi["most_massive"][i]] .* ias["merger map"][imt[1],mmi["most_massive"][i]])
                JiMM = hcat( JiMM,  ias["merger map"][imt[5:7],mmi["most_massive"][i]] )
                if count(ismissing, JiMM[:,end]) == 0
                        JsMM = hcat( JsMM, JiMM[:,end] .+ JoMM[:,end] )
                else
                        JsMM = hcat( JsMM, JoMM[:,end] )
                end
                Jo_sum = zeros(3)
                Ji_sum = zeros(3)
                for ii in mmi["first"][i]:mmi["last"][i] # iterate over merger entries in merger map
                        if count(ismissing, ias["merger map"][imt[2:4],ii]) == 0
                                Jo_sum .+= ias["merger map"][imt[2:4],ii] .* ias["merger map"][imt[1],ii]
                        end
                        if count(ismissing, ias["merger map"][imt[5:7],ii]) == 0
                                Ji_sum .+= ias["merger map"][imt[5:7],ii]
                                idx_caught = vcat( idx_caught, ii )
                        else 
                                idx_missed = vcat( idx_missed, ii )
                        end
                end
                Jom = hcat( Jom,  ifelse(count(iszero, Jo_sum) == 3, missings(3), Jo_sum) )
                Jim = hcat( Jim,  ifelse(count(iszero, Ji_sum) == 3, missings(3), Ji_sum) ) # missings if zero vector
                Jsm = hcat( Jsm,  ifelse(count(iszero, Jo_sum) == 3, missings(3), Jo_sum) .+ Ji_sum )
        end
        
        angle_J0_Jom  = missings(Float64, length(J_0[1,:]))
        angle_J0_Jim  = missings(Float64, length(J_0[1,:]))
        angle_J0_Jsm  = missings(Float64, length(J_0[1,:]))
        angle_J0_JoMM = missings(Float64, length(J_0[1,:]))
        angle_J0_JiMM = missings(Float64, length(J_0[1,:]))
        angle_J0_JsMM = missings(Float64, length(J_0[1,:]))
        angle_J1_Jom  = missings(Float64, length(J_0[1,:]))
        angle_J1_Jim  = missings(Float64, length(J_0[1,:]))
        angle_J1_Jsm  = missings(Float64, length(J_0[1,:]))
        angle_J1_JoMM = missings(Float64, length(J_0[1,:]))
        angle_J1_JiMM = missings(Float64, length(J_0[1,:]))
        angle_J1_JsMM = missings(Float64, length(J_0[1,:]))
        angle_Jim_Jom = missings(Float64, length(J_0[1,:]))
        angle_JiMM_JoMM= missings(Float64, length(J_0[1,:]))
        for i in 1:length(J_0[1,:])
                if count(ismissing, Jom[:,i]) == 0
                        angle_J0_Jom[i] = 180/π * angle(J_0[:,i], J_0[:,i] .+ convert(Array{Float64,1},Jom[:,i]) ) 
                        angle_J1_Jom[i] = 180/π * angle(J_1[:,i], J_0[:,i] .+ convert(Array{Float64,1},Jom[:,i]) ) 
                        if count(ismissing, Jim[:,i]) == 0
                                angle_Jim_Jom[i] = 180/π * angle(convert(Array{Float64,1},Jim[:,i]), convert(Array{Float64,1},Jom[:,i]) ) 
                        end
                end
                if count(ismissing, JoMM[:,i]) == 0
                        angle_J0_JoMM[i] = 180/π * angle(J_0[:,i], J_0[:,i] .+ convert(Array{Float64,1},JoMM[:,i]) ) 
                        angle_J1_JoMM[i] = 180/π * angle(J_1[:,i], J_0[:,i] .+ convert(Array{Float64,1},JoMM[:,i]) ) 
                        if count(ismissing, JiMM[:,i]) == 0
                                angle_JiMM_JoMM[i] = 180/π * angle(convert(Array{Float64,1},JiMM[:,i]), convert(Array{Float64,1},JoMM[:,i]) ) 
                        end
                end
                if count(ismissing, Jsm[:,i]) == 0
                        angle_J0_Jsm[i] = 180/π * angle(J_0[:,i], J_0[:,i] .+ convert(Array{Float64,1},Jsm[:,i]) ) 
                        angle_J1_Jsm[i] = 180/π * angle(J_1[:,i], J_0[:,i] .+ convert(Array{Float64,1},Jsm[:,i]) ) 
                end
                if count(ismissing, JsMM[:,i]) == 0
                        angle_J0_JsMM[i] = 180/π * angle(J_0[:,i], J_0[:,i] .+ convert(Array{Float64,1},JsMM[:,i]) ) 
                        angle_J1_JsMM[i] = 180/π * angle(J_1[:,i], J_0[:,i] .+ convert(Array{Float64,1},JsMM[:,i]) ) 
                end
                if count(ismissing, Jim[:,i]) == 0
                        angle_J0_Jim[i] = 180/π * angle(J_0[:,i], J_0[:,i] .+ convert(Array{Float64,1},Jim[:,i]) ) 
                        angle_J1_Jim[i] = 180/π * angle(J_1[:,i], J_0[:,i] .+ convert(Array{Float64,1},Jim[:,i]) ) 
                end
                if count(ismissing, JiMM[:,i]) == 0
                        angle_J0_JiMM[i] = 180/π * angle(J_0[:,i], J_0[:,i] .+ convert(Array{Float64,1},JiMM[:,i]) ) 
                        angle_J1_JiMM[i] = 180/π * angle(J_1[:,i], J_0[:,i] .+ convert(Array{Float64,1},JiMM[:,i]) ) 
                end
        end
        
        
        ###################################################################
        # jim jom
        xmin = 8
        ymin = 8
        xmax = 16
        ymax = 16
        outnm     = joinpath(outputdir,"Jim_vs_Jom-$(im1)-$(ia1)")
        outfile = outnm*"_hxcount.pdf"
        x       = pyplottable(Jom)
        y       = pyplottable(Jim)
        plot_hexbins( log10.(x), log10.(y), 
            outfile=outfile,
            #title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) $(ia2), $(im2) Read",
            scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
            xlabel="log( J$(latexstring("\$_\\textrm{orb}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
            ylabel="log( J$(latexstring("\$_\\textrm{int}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
            grid=true,
            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
            lognorm=true,
            binned_median = "y", plot_bval=" ", plot_line=hcat([8,18],[8,18]),
            calc_norm=true,    # @pyplottable
            )
        #
        
        
        
        
        ###################################################################
        # predict delta J
        xmin = 11.5#" "#1e0
        ymin = 10#" "#xmin
        xmax = 15#" "#180
        ymax = 16#" "#xmax
        outnm     = joinpath(outputdir,"J0Jsm_v_Jdelta-$(im1)-$(ia1)")
        outfile = outnm*"_hxcount.pdf"
        x       = pyplottable(ias["ΔJ_main"][:,idxmatch(mmi["main"], ic[:transition])])
        y       = pyplottable(Jsm)
        #x       = pyplottable(ias["J_main"][:,idxmatch(mmi["main"], ic[:transition])])
        #y       = pyplottable(ias["J_main"][:,idxmatch(mmi["main"], ic[:transition])] .- ias["ΔJ_main"][:,idxmatch(mmi["main"], ic[:transition])] .+ Jsm)
        plot_hexbins( log10.(x), log10.(y), 
            outfile=outfile,
            #title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) $(ia2), $(im2) Read",
            scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
            xlabel="log( $(latexstring("\$\\Delta\$"))J [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
            ylabel="log( J$(latexstring("\$_\\textrm{Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
            grid=true,
            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
            lognorm=true,
            binned_median = "y", plot_bval=" ", plot_line=hcat([9,18],[9,18]),
            calc_norm=true,    # @pyplottable
            )
        #
        
        
        colr    = Vector{Float64}(undef, 0)
        for i in 1:length(idxmatch(mmi["main"], ic[:transition]))
                colr    = vcat( colr, (180/π) .* angle(convert(Array{Float64,1},replace(ias["ΔJ_main"][:,idxmatch(mmi["main"], ic[:transition])][:,i], missing => NaN)), convert(Array{Float64,1},replace(Jsm[:,i], missing => NaN)) ) )
        end
        
        outfile = outnm*"_hxalign.pdf"
        plot_hexbins( log10.(x), log10.(y), 
        outfile=outfile,
        #title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) $(ia2), $(im2) Read",
        scale=15, gridres=(15,7), xscale="linear", yscale="linear", 
        xlabel="log( $(latexstring("\$\\Delta\$"))J [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
        ylabel="log( J$(latexstring("\$_\\textrm{Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
        grid=true,
        xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
        lognorm=false,
        binned_median = " ", plot_bval=" ", plot_line=hcat([9,18],[9,18]),
        calc_norm=true,    # @pyplottable
        cmap = "seismic", colormod = colr, colorlimits=(0,180), cfunc=mean,#twilight_shifted
        )
        #
        
        
        order   = sortperm(colr)
        colr    = colr[order]
        maxcolr = maximum(colr)
        #slct    = idxmatch(slct2, n_switches)
        scale = 0.5
        PyPlot.matplotlib[:rc]("text", usetex=true)
        PyPlot.matplotlib[:rc]("font", family="serif", serif=["Computer Modern"])
        pointsize = 2
        fig, ax = subplots()
        ax.plot([11,15], [11,15], "-",color="black", linewidth=0.5)
        p = ax.scatter( log10.(x[order]), log10.(y[order]),
                zorder=2   , s= pointsize    , alpha=0.5,#colr ./ maxcolr, (pointsize/maxcolr) .* colr
                cmap="plasma_r", c=colr) # brg_r  , vmin=-6.1,vmax=-3.9
        #ax.plot([11,16], [11,16], "-",color="black", linewidth=0.5)
        colbar  = fig.colorbar(p, ax=ax)#, ticks=[minimum(ias["BVAL_0"    ]), b_ell, b_disk, maximum(ias["BVAL_0"    ])])
        #colbar.ax.set_yticklabels(["$(round(minimum(ias["BVAL_0"    ]), digits=1))", "ETG", "LTG", "$(round(maximum(ias["BVAL_0"    ]), digits=1))"])
        colbar.ax.tick_params(axis="y")
        colbar.ax.set_ylabel("Alignment Angle [°]")
        colbar.ax.yaxis.set_label_position("left")
        #ax.legend(loc="lower right", frameon=true, borderpad=1, handlelength=1.8)
        ax.set_ylim([ymin,ymax])
        ax.set_xlim([xmin,xmax])
        #ax.set_aspect("equal")
        ax.set_ylabel("log( J$(latexstring("\$_\\textrm{Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )")
        ax.set_xlabel("log( $(latexstring("\$\\Delta\$"))J [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )")
        ax.grid()
        #ax.autoscale()
        fig.set_size_inches(11*scale, 9*scale)
        fig.tight_layout()
        fig.savefig(outnm*"sctalign.pdf", bbox_inches="tight", pad_inches=.1)
        close("all")
        PyPlot.matplotlib[:rc]("text", usetex=false)
        
        
        #colr = ias["BVAL_0"][idxmatch(mmi["main"], ic[:transition])]
        outfile = outnm*"_hxbvalue.pdf"
        plot_hexbins( log10.(x), log10.(y), 
        outfile=outfile,
        #title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) $(ia2), $(im2) Read",
        scale=15, gridres=(30,13), xscale="linear", yscale="linear", 
        xlabel="log( $(latexstring("\$\\Delta\$"))J [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
        ylabel="log( J$(latexstring("\$_\\textrm{Mergers}\$")) [ M$(latexstring("\$_\\odot\$")) kpc km/s ] )", 
        grid=true,
        xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
        lognorm=false,
        binned_median = "y", plot_bval=" ", plot_line=hcat([9,18],[9,18]),
        calc_norm=true,    # @pyplottable
        cmap = ColorMap(get_bvalue_cmap().colors), colormod = ias["BVAL_0"][idxmatch(mmi["main"], ic[:transition])], colorlimits=(-6.1,-3.9), cfunc=mean,
        )
        
        
        
        
        
        ###################################################################
        # predict angle
        # log
        xmin = log10(1e-3)#" "#1e0
        ymin = log10(1e-3)#" "#xmin
        xmax = log10(200)#" "#180
        ymax = log10(200)#" "#xmax
        #y       = pyplottable(delJ[:,order] .- Jom[:,order])
        colr    = ias["BVAL_0"][idxmatch(mmi["main"], ic[:transition])] .- ias["ΔBVAL_0"][idxmatch(mmi["main"], ic[:transition])] # log10.(angle_J1_Jsm ./ x)
        x       = log10.(pyplottable(ias["ϕ_flip"][idxmatch(mmi["main"], ic[:transition])]))
        y       = log10.(pyplottable(angle_J0_Jsm))
        outnm     = joinpath(outputdir,"angmerger_v_ang-$(im1)-$(ia1)_log")
        outfile = outnm*"_hxcount.pdf"
        #x       = pyplottable(ias["J_main"][:,idxmatch(mmi["main"], ic[:transition])])
        #y       = pyplottable(ias["J_main"][:,idxmatch(mmi["main"], ic[:transition])] .- ias["ΔJ_main"][:,idxmatch(mmi["main"], ic[:transition])] .+ Jsm)
        plot_hexbins( x, y, 
            outfile=outfile,
            #title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) $(ia2), $(im2) Read",
            scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
            ylabel="log( Flip by Merger Contribution [°] )",
            xlabel="log( Flip Angle [°] )",
            grid=true,
            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
            lognorm=true,
            binned_median = "x", plot_bval=" ", plot_line=hcat([xmin,xmax],[ymin,ymax]),
            calc_norm=true,    # @pyplottable
            )
        #
        
        
        outfile = outnm*"_hxbvalue.pdf"
        plot_hexbins( x, y, 
            outfile=outfile,
            #title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) $(ia2), $(im2) Read",
            scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
            ylabel="log( Flip by Merger Contribution [°] )",
            xlabel="log( Flip Angle [°] )",
            grid=true,
            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
            lognorm=false,
            binned_median = "x", plot_bval=" ", plot_line=hcat([xmin,xmax],[ymin,ymax]),
            calc_norm=true,    # @pyplottable
            cmap=ColorMap(get_bvalue_cmap().colors), colormod = colr, cfunc=mean, colorlimits=(-6.1,-3.9),
        )
        #

        colr    = pyplottable( log10.(pyplottable(Jsm)) ./ log10.(pyplottable(ias["ΔJ_main"][:,idxmatch(mmi["main"], ic[:transition])])) )
        outfile = outnm*"_hxdJratio.pdf"
        plot_hexbins( x, y, 
            outfile=outfile,
            #title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) $(ia2), $(im2) Read",
            scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
            ylabel="log( Flip by Merger Contribution [°] )",
            xlabel="log( Flip Angle [°] )",
            grid=true,
            xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
            lognorm=false,
            binned_median = "x", plot_bval=" ", plot_line=hcat([xmin,xmax],[ymin,ymax]),
            calc_norm=true,    # @pyplottable
            cmap="seismic", colormod = colr, cfunc=mean, colorlimits=(0.8,1.2),
        )
        #
        
        for xmax in [5,30,180]
                xmin = 0#" "#1e0
                ymin = 0#" "#xmin
                ymax = xmax
                #y       = pyplottable(delJ[:,order] .- Jom[:,order])
                colr    = ias["BVAL_0"][idxmatch(mmi["main"], ic[:transition])] .- ias["ΔBVAL_0"][idxmatch(mmi["main"], ic[:transition])] # log10.(angle_J1_Jsm ./ x)
                x       = pyplottable(ias["ϕ_flip"][idxmatch(mmi["main"], ic[:transition])])
                y       = pyplottable(angle_J0_Jsm)
                outnm     = joinpath(outputdir,"angmerger_v_ang-$(im1)-$(ia1)_lin$(xmax)")
                outfile = outnm*"_hxcount.pdf"
                #x       = pyplottable(ias["J_main"][:,idxmatch(mmi["main"], ic[:transition])])
                #y       = pyplottable(ias["J_main"][:,idxmatch(mmi["main"], ic[:transition])] .- ias["ΔJ_main"][:,idxmatch(mmi["main"], ic[:transition])] .+ Jsm)
                plot_hexbins( x, y, 
                    outfile=outfile,
                    #title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) $(ia2), $(im2) Read",
                    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
                    ylabel="Flip by Merger Contribution [°]",
                    xlabel="Flip Angle [°]",
                    grid=true,
                    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                    lognorm=true,
                    binned_median = "y", plot_bval=" ", plot_line=hcat([xmin,xmax],[ymin,ymax]),
                    calc_norm=true,    # @pyplottable
                    )
                #
                
                
                outfile = outnm*"_hxbvalue.pdf"
                plot_hexbins( x, y, 
                    outfile=outfile,
                    #title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) $(ia2), $(im2) Read",
                    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
                    ylabel="Flip by Merger Contribution [°]",
                    xlabel="Flip Angle [°]",
                    grid=true,
                    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                    lognorm=false,
                    binned_median = "y", plot_bval=" ", plot_line=hcat([xmin,xmax],[ymin,ymax]),
                    calc_norm=true,    # @pyplottable
                    cmap=ColorMap(get_bvalue_cmap().colors), colormod = colr, colorlimits=(-6.1,-3.9), cfunc=mean,
                )
                
                
                colr    = pyplottable( log10.(pyplottable(Jsm)) ./ log10.(pyplottable(ias["ΔJ_main"][:,idxmatch(mmi["main"], ic[:transition])])) )
                outfile = outnm*"_hxdJratio.pdf"
                plot_hexbins( x, y, 
                    outfile=outfile,
                    #title="$(latexstring("\$\\Delta\$"))t $(latexstring("\$\\sim\$")) $(ia2), $(im2) Read",
                    scale=15, gridres=(20,8), xscale="linear", yscale="linear", 
                    ylabel="Flip by Merger Contribution [°]",
                    xlabel="Flip Angle [°]",
                    grid=true,
                    xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax,
                    lognorm=false,
                    binned_median = "y", plot_bval=" ", plot_line=hcat([xmin,xmax],[ymin,ymax]),
                    calc_norm=true,    # @pyplottable
                    cmap="seismic", colormod = colr, cfunc=mean, colorlimits=(0.8,1.2),
                )
        end
        
    end













    # flip limit multiplots
    outnm   = joinpath(outputdir,"flipangle_v_dJrelpre_$(ia1)Gyr")
    outfile = outnm*"_hxcount.pdf"
    y       = pyplottable(ias["ϕ_flip"])
    x       = pyplottable(log10.(pyplottable(ias["ΔJ_main"]) ./ pyplottable(ias["J_main"] .- ias["ΔJ_main"])))
    colr    = ias["BVAL_0"] .- ias["ΔBVAL_0"]
    slct    = ic[:transition]#idxclude(idxclude(ic[:start_thr], ic[:switches]), ic[:fakes])
    miny= 0#minimum(y[slct])
    maxy= 180#maximum(y[slct])
    minx= -2.5#minimum(x[slct])
    maxx= 1#maximum(x[slct])
    slcts       = Dict{Int64, Vector{Int64}}()
    slcts[1]    = slct
    slcts[4]    = idxmatch(slct, ic[:start_disk])
    slcts[2]    = idxmatch(slct, ic[:start_int])
    slcts[3]    = idxmatch(slct, ic[:start_ell])
    labels      = Dict{Int64, String}()
    labels[1]   = "Full Set"#Full Set, N = $(length(slcts[1]))"
    labels[4]   = "Disks"#$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
    labels[2]   = "Intermediates"#$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
    labels[3]   = "Ellipticals"#$(round(maximum(as["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
    labpos      = Array{String}(undef, 4)
    labpos     .= "upper left"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=true, gridres=(20,7), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="$(ia3) Flip Angle [°]",
                fontsize=30,
                binned_median="y",
                plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_vline=[0, 0.5, 1],
                )
    #
    
    outfile = outnm*"_hxbval.pdf"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=false, gridres=(20,7), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="$(ia3) Flip Angle [°]",
                fontsize=30,
                binned_median="y",
                cfunc=mean, colormod=colr, cmap=ColorMap(get_bvalue_cmap().colors), colorlimits=(-6.1,-3.9),
                plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_vline=[0, 0.5, 1],
                )
    #
    
    slcts       = Dict{Int64, Vector{Int64}}()
    slcts[1]    = slct
    slcts[3]    = idxmatch(slct, ic[:Mpre_hi])
    slcts[2]    = idxmatch(slct, ic[:Mpre_md])
    slcts[4]    = idxmatch(slct, ic[:Mpre_lo])
    labels      = Dict{Int64, String}()
    labels[1]   = "Full Set"#Full Set, N = $(length(slcts[1]))"
    labels[3]   = "High Mass"#$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
    labels[2]   = "Medium Mass"#$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
    labels[4]   = "Low Mass"#$(round(maximum(as["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
    colr    = log10.(ias["M"] .- ias["ΔM"])
    outfile = outnm*"_hxmass.pdf"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=false, gridres=(20,7), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="$(ia3) Flip Angle [°]",
                fontsize=30,
                binned_median="y",
                cfunc=mean, colormod=colr, cmap="plasma_r", colorlimits=(log10(2e10),log10(maximum(ias["M"][slct]))),
                plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_vline=[0, 0.5, 1],
                )
    #
    
    
    
    
    #outnm   = joinpath(outputdir,"flipangle_v_dJrelpost_$(ia1)Gyr")
    outnm   = joinpath(outputdir,"flipangle_v_dJrelpre_bpost_$(ia1)Gyr")
    outfile = outnm*"_hxcount.pdf"
    y       = pyplottable(ias["ϕ_flip"])
    #x       = pyplottable(log10.(pyplottable(ias["ΔJ_main"]) ./ pyplottable(ias["J_main"])))
    x       = pyplottable(log10.(pyplottable(ias["ΔJ_main"]) ./ pyplottable(ias["J_main"] .- ias["ΔJ_main"])))
    colr    = ias["BVAL_0"]
    slct    = ic[:transition]#idxclude(idxclude(ic[:start_thr], ic[:switches]), ic[:fakes])
    miny= 0#minimum(y[slct])
    maxy= 180#maximum(y[slct])
    minx= -2.5#minimum(x[slct])
    maxx= 1#maximum(x[slct])
    slcts       = Dict{Int64, Vector{Int64}}()
    slcts[1]    = slct
    slcts[4]    = idxmatch(slct, ic[:start_disk])
    slcts[2]    = idxmatch(slct, ic[:start_int])
    slcts[3]    = idxmatch(slct, ic[:start_ell])
    labels      = Dict{Int64, String}()
    labels[1]   = "Full Set"#Full Set, N = $(length(slcts[1]))"
    labels[4]   = "Disks"#$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
    labels[2]   = "Intermediates"#$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
    labels[3]   = "Ellipticals"#$(round(maximum(as["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
    labpos      = Array{String}(undef, 4)
    labpos     .= "upper left"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=true, gridres=(20,7), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                #xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,post}\$")) )",
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="$(ia3) Flip Angle [°]",
                fontsize=30,
                binned_median="y",
                plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_vline=[0, 0.5, 1],
                )
    #
    
    outfile = outnm*"_hxbval.pdf"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=false, gridres=(20,7), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                #xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,post}\$")) )",
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="$(ia3) Flip Angle [°]",
                fontsize=30,
                binned_median="y",
                cfunc=mean, colormod=colr, cmap=ColorMap(get_bvalue_cmap().colors), colorlimits=(-6.1,-3.9),
                plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_vline=[0, 0.5, 1],
                )
    #
    
    slcts       = Dict{Int64, Vector{Int64}}()
    slcts[1]    = slct
    slcts[3]    = idxmatch(slct, ic[:Mpre_hi])
    slcts[2]    = idxmatch(slct, ic[:Mpre_md])
    slcts[4]    = idxmatch(slct, ic[:Mpre_lo])
    labels      = Dict{Int64, String}()
    labels[1]   = "Full Set"#Full Set, N = $(length(slcts[1]))"
    labels[3]   = "High Mass"#$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
    labels[2]   = "Medium Mass"#$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
    labels[4]   = "Low Mass"#$(round(maximum(as["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
    colr    = log10.(ias["M"])
    outfile = outnm*"_hxmass.pdf"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=false, gridres=(20,7), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                #xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,post}\$")) )",
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="$(ia3) Flip Angle [°]",
                fontsize=30,
                binned_median="y",
                cfunc=mean, colormod=colr, cmap="plasma_r", colorlimits=(log10(2e10),log10(maximum(ias["M"][slct]))),
                plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_vline=[0, 0.5, 1],
                )
    #
    









    
    # Difference flip limit multiplots
    outnm   = joinpath(outputdir,"diff_flipangle_v_dJrelpre_$(ia1)Gyr")
    outfile = outnm*"_hxcount.pdf"
    x       = pyplottable(log10.(pyplottable(ias["ΔJ_main"]) ./ pyplottable(ias["J_main"] .- ias["ΔJ_main"])) )
    y       = pyplottable(ias["ϕ_flip"] .- avg_flipangled_v_dj.(10 .^ x) )
    colr    = ias["BVAL_0"] .- ias["ΔBVAL_0"]
    slct    = ic[:transition]#idxclude(idxclude(ic[:start_thr], ic[:switches]), ic[:fakes])
    miny= minimum(y[slct])
    maxy= maximum(y[slct])
    minx= -2.5#minimum(x[slct])
    maxx= 1#maximum(x[slct])
    slcts       = Dict{Int64, Vector{Int64}}()
    slcts[1]    = slct
    slcts[4]    = idxmatch(slct, ic[:start_disk])
    slcts[2]    = idxmatch(slct, ic[:start_int])
    slcts[3]    = idxmatch(slct, ic[:start_ell])
    labels      = Dict{Int64, String}()
    labels[1]   = "Full Set"#Full Set, N = $(length(slcts[1]))"
    labels[4]   = "Disks"#$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
    labels[2]   = "Intermediates"#$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
    labels[3]   = "Ellipticals"#$(round(maximum(as["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
    labpos      = Array{String}(undef, 4)
    labpos     .= "upper left"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=true, gridres=(20,7), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="$(ia3) Flip Angle Difference [°]",
                fontsize=30,
                binned_median="y",
                plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), zeros(1000)),
                #plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000)) .- avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                #plot_vline=[0, 0, 1],
                )
    #
    
    outfile = outnm*"_hxbval.pdf"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=false, gridres=(20,7), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="$(ia3) Flip Angle Difference [°]",
                fontsize=30,
                binned_median="y",
                cfunc=mean, colormod=colr, cmap=ColorMap(get_bvalue_cmap().colors), colorlimits=(-6.1,-3.9),
                plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), zeros(1000)),
                #plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000)) .- avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                #plot_vline=[0, 0, 1],
                )
    #
    
    slcts       = Dict{Int64, Vector{Int64}}()
    slcts[1]    = slct
    slcts[3]    = idxmatch(slct, ic[:Mpre_hi])
    slcts[2]    = idxmatch(slct, ic[:Mpre_md])
    slcts[4]    = idxmatch(slct, ic[:Mpre_lo])
    labels      = Dict{Int64, String}()
    labels[1]   = "Full Set"#Full Set, N = $(length(slcts[1]))"
    labels[3]   = "High Mass"#$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
    labels[2]   = "Medium Mass"#$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
    labels[4]   = "Low Mass"#$(round(maximum(as["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
    colr    = log10.(ias["M"] .- ias["ΔM"])
    outfile = outnm*"_hxmass.pdf"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=false, gridres=(20,7), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="$(ia3) Flip Angle Difference [°]",
                fontsize=30,
                binned_median="y",
                cfunc=mean, colormod=colr, cmap="plasma_r", colorlimits=(log10(2e10),log10(maximum(ias["M"][slct]))),
                plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), zeros(1000)),
                #plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000)) .- avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                #plot_vline=[0, 0, 1],
                )
    #
    
    
    
    
    outnm   = joinpath(outputdir,"diff_flipangle_v_dJrelpre_bpost_$(ia1)Gyr")
    outfile = outnm*"_hxcount.pdf"
    #x       = pyplottable(log10.(pyplottable(ias["ΔJ_main"]) ./ pyplottable(ias["J_main"])))
    x       = pyplottable(log10.(pyplottable(ias["ΔJ_main"]) ./ pyplottable(ias["J_main"] .- ias["ΔJ_main"])))
    y       = pyplottable(ias["ϕ_flip"] .- avg_flipangled_v_dj.(10 .^ x) )
    colr    = ias["BVAL_0"]
    slct    = ic[:transition]#idxclude(idxclude(ic[:start_thr], ic[:switches]), ic[:fakes])
    miny= minimum(y[slct])
    maxy= maximum(y[slct])
    minx= -2.5#minimum(x[slct])
    maxx= 1#maximum(x[slct])
    slcts       = Dict{Int64, Vector{Int64}}()
    slcts[1]    = slct
    slcts[4]    = idxmatch(slct, ic[:start_disk])
    slcts[2]    = idxmatch(slct, ic[:start_int])
    slcts[3]    = idxmatch(slct, ic[:start_ell])
    labels      = Dict{Int64, String}()
    labels[1]   = "Full Set"#Full Set, N = $(length(slcts[1]))"
    labels[4]   = "Disks"#$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
    labels[2]   = "Intermediates"#$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
    labels[3]   = "Ellipticals"#$(round(maximum(as["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
    labpos      = Array{String}(undef, 4)
    labpos     .= "upper left"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=true, gridres=(20,7), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                #xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,post}\$")) )",
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="$(ia3) Flip Angle Difference [°]",
                fontsize=30,
                binned_median="y",
                plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), zeros(1000)),
                #plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000)) .- avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                #plot_vline=[0, 0, 1],
                )
    #
    
    outfile = outnm*"_hxbval.pdf"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=false, gridres=(20,7), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                #xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,post}\$")) )",
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="$(ia3) Flip Angle Difference [°]",
                fontsize=30,
                binned_median="y",
                cfunc=mean, colormod=colr, cmap=ColorMap(get_bvalue_cmap().colors), colorlimits=(-6.1,-3.9),
                plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), zeros(1000)),
                #plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000)) .- avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                #plot_vline=[0, 0, 1],
                )
    #
    
    slcts       = Dict{Int64, Vector{Int64}}()
    slcts[1]    = slct
    slcts[3]    = idxmatch(slct, ic[:Mpre_hi])
    slcts[2]    = idxmatch(slct, ic[:Mpre_md])
    slcts[4]    = idxmatch(slct, ic[:Mpre_lo])
    labels      = Dict{Int64, String}()
    labels[1]   = "Full Set"#Full Set, N = $(length(slcts[1]))"
    labels[3]   = "High Mass"#$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
    labels[2]   = "Medium Mass"#$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
    labels[4]   = "Low Mass"#$(round(maximum(as["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
    colr    = log10.(ias["M"])
    outfile = outnm*"_hxmass.pdf"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=false, gridres=(20,7), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                #xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,post}\$")) )",
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="$(ia3) Flip Angle Difference [°]",
                fontsize=30,
                binned_median="y",
                cfunc=mean, colormod=colr, cmap="plasma_r", colorlimits=(log10(2e10),log10(maximum(ias["M"][slct]))),
                plot_line=hcat(log10.(LinRange(0.01, 10, 1000)), zeros(1000)),
                #plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), max_flipangled_v_dj.(LinRange(0.01, 10, 1000)) .- avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                #plot_vline=[0, 0, 1],
                )
    #
    
    



    





    # jrel vs Jrel multiplots
    outnm   = joinpath(outputdir,"djrelpre_v_dJrelpre_$(ia1)Gyr")
    outfile = outnm*"_hxcount.pdf"
    y       = pyplottable(log10.(pyplottable(ias["Δj_main"]) ./ pyplottable(ias["j_main"] .- ias["Δj_main"])))
    x       = pyplottable(log10.(pyplottable(ias["ΔJ_main"]) ./ pyplottable(ias["J_main"] .- ias["ΔJ_main"])))
    colr    = ias["ϕ_flip"]
    slct    = ic[:transition]#idxclude(idxclude(ic[:start_thr], ic[:switches]), ic[:fakes])
    miny= -2.5#minimum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))
    maxy= 1#maximum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))
    minx= -2.5#minimum(x[slct])
    maxx= 1#maximum(x[slct])
    slcts       = Dict{Int64, Vector{Int64}}()
    slcts[1]    = slct
    slcts[3]    = idxmatch(slct, ic[:start_disk])
    slcts[2]    = idxmatch(slct, ic[:start_int])
    slcts[4]    = idxmatch(slct, ic[:start_ell])
    labels      = Dict{Int64, String}()
    labels[1]   = "Full Set"#Full Set, N = $(length(slcts[1]))"
    labels[3]   = "Disks"#$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
    labels[2]   = "Intermediates"#$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
    labels[4]   = "Ellipticals"#$(round(maximum(as["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
    labpos      = Array{String}(undef, 4)
    labpos     .= "upper left"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=true, gridres=(40,14), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="log( $(latexstring("\$\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,pre}\$")) )",
                fontsize=30,
                binned_median="y",
                plot_line=hcat(LinRange(minx, maxx, 1000), LinRange(miny, maxy, 1000)),
                #plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                )
    #
    
    outfile = outnm*"_hxflip.pdf"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=false, gridres=(40,14), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,pre}\$")) )",
                ylabel="log( $(latexstring("\$\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,pre}\$")) )",
                fontsize=30,
                binned_median="y",
                cfunc=mean, colormod=colr, cmap="plasma_r", colorlimits=(0,180), #clabel="Flip Angle [°]",
                plot_line=hcat(LinRange(minx, maxx, 1000), LinRange(miny, maxy, 1000)),
                #plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                )
    #
    
    
    
    outnm   = joinpath(outputdir,"djrelpost_v_dJrelpost_$(ia1)Gyr")
    outfile = outnm*"_hxcount.pdf"
    y       = pyplottable(log10.(pyplottable(ias["Δj_main"]) ./ pyplottable(ias["j_main"])))
    x       = pyplottable(log10.(pyplottable(ias["ΔJ_main"]) ./ pyplottable(ias["J_main"])))
    colr    = ias["ϕ_flip"]
    slct    = ic[:transition]#idxclude(idxclude(ic[:start_thr], ic[:switches]), ic[:fakes])
    miny= -2.5#minimum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))
    maxy= 1#maximum(pyplottable(y[findall(x->x .!== NaN, pyplottable(y))]))
    minx= -2.5#minimum(x[slct])
    maxx= 1#maximum(x[slct])
    slcts       = Dict{Int64, Vector{Int64}}()
    slcts[1]    = slct
    slcts[4]    = idxmatch(slct, ic[:start_disk])
    slcts[2]    = idxmatch(slct, ic[:start_int])
    slcts[3]    = idxmatch(slct, ic[:start_ell])
    labels      = Dict{Int64, String}()
    labels[1]   = "Full Set"#Full Set, N = $(length(slcts[1]))"
    labels[4]   = "Disks"#$(round(minimum(as["redshift"][slcts[3]]), digits=1))  $(latexstring("\$>\$")) z  $(latexstring("\$>\$")) $(round(maximum(as["redshift"][slcts[4]]), digits=1)), N = $(length(slcts[2]))"
    labels[2]   = "Intermediates"#$(round(maximum(as["redshift"][slcts[3]]), digits=1)) $(latexstring("\$\\geq\$")) z $(latexstring("\$\\geq\$")) $(round(minimum(as["redshift"][slcts[3]]), digits=1)), N = $(length(slcts[3]))"
    labels[3]   = "Ellipticals"#$(round(maximum(as["redshift"][slcts[4]]), digits=1)) $(latexstring("\$\\geq\$")) z, N = $(length(slcts[4]))"
    labpos      = Array{String}(undef, 4)
    labpos     .= "upper left"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=true, gridres=(40,14), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,post}\$")) )",
                ylabel="log( $(latexstring("\$\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,post}\$")) )",
                fontsize=30,
                binned_median="y",
                plot_line=hcat(LinRange(minx, maxx, 1000), LinRange(miny, maxy, 1000)),
                #plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                )
    #
    
    outfile = outnm*"_hxflip.pdf"
    plot_hexbins_nxm( x, y, 
                slcts, labels, 2, 2,
                #scale=0.0,
                cbarshift=2.5,
                label_pos=labpos, 
                outfile=outfile, 
                lognorm=false, gridres=(40,14), 
                ymin=miny, ymax=maxy,
                xmin=minx, xmax=maxx, 
                grid=false,
                xlabel="log( $(latexstring("\$\\Delta\$"))J$(latexstring("\$_*\$")) / J$(latexstring("\$_\\textrm{*,post}\$")) )",
                ylabel="log( $(latexstring("\$\\Delta\$"))j$(latexstring("\$_*\$")) / j$(latexstring("\$_\\textrm{*,post}\$")) )",
                fontsize=30,
                binned_median="y",
                cfunc=mean, colormod=colr, cmap="plasma_r", colorlimits=(0,180), #clabel="Flip Angle [°]",
                plot_line=hcat(LinRange(minx, maxx, 1000), LinRange(miny, maxy, 1000)),
                #plot_line2=hcat(log10.(LinRange(0.01, 10, 1000)), avg_flipangled_v_dj.(LinRange(0.01, 10, 1000))),
                )
    #
    
    
end

