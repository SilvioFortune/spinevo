
include("/home/moon/sfortune/spinevo/pkg/meta.jl")


# Arrows 
halo = 3212
plot_halostory_arrows(; outfile="/home/moon/sfortune/spinevo/plots/halo$(halo)/halo_$(halo)_halostory_arrows.png", 
    felixID=" ", rootID=halo, 
    parttype="STARS", snap=true, 
    scale=1, plotwidth=16, plotheight=3, gridspace=0.6,
    indir=current_dir_stories, 
    )



# Data
halo = 3212
plot_halostory_data(; outfile="/home/moon/sfortune/spinevo/plots/halo$(halo)/halo_$(halo)_halostory_data.png", 
    felixID=" ", rootID=halo, 
    parttype="STARS", snap=true, 
    scale=0.7, gridspace=0.6,  #Gyr
    indir=current_dir_stories, 
    )


# Phase Space plot
res = 1
plot_phasespace(; #snapNR = 132, 
    felixID=" ", rootID=3212, subID=" " , ymin=-1000, ymax=1000, plottype="scatter", pointsize=5, opact="mass", 
    parttype="STARS", r_by_rvir=0.1, #weights="MASS" ,
    lognorm=true, scale=0.5, gridres=(25*res, 15*res), #colorlimits=(0.7,300), 
    outdir="/home/moon/sfortune/spinevo/plots/halo3212", outfile="phsp", 
    )





# Plot halo 3212
halo = 3212
snap = 36
prop = "vel"
plot_group(snap; rootID=halo, 
property=prop, stepsize=10, arsize=5, ptsize=5000, rad=200, res=(1920,1080), spin=false, 
bgcolor=:black, indir=current_dir_stories, boxfix=false, 
file="/home/moon/sfortune/spinevo/plots/halo3212/halo$(halo)_snap$(snap)_$(prop).html")
