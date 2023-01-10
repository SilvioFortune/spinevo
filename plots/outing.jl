
include("/home/moon/sfortune/spinevo/pkg/meta.jl")

outputdir = "/home/moon/sfortune/spinevo/plots/outing"

include("/home/moon/sfortune/spinevo/plots/load_data.jl")

root_id = 3968


rootfor136 = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=136),findcs(as03["subID"], eq=3968))]
rootfor132 = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=132),findcs(as03["subID"], eq=3735))]
rootfor128 = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=128),findcs(as03["subID"], eq=2931))]
rootfor124 = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=124),findcs(as03["subID"], eq=2807))]
rootfor120 = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=120),findcs(as03["subID"], eq=2663))]
rootfor116 = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=116),findcs(as03["subID"], eq=2639))]
rootfor112 = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=112),findcs(as03["subID"], eq=2479))]
rootfor108 = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=108),findcs(as03["subID"], eq=2165))]
rootfor106 = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=106),findcs(as03["subID"], eq=2219))]
rootfor104 = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=104),findcs(as03["subID"], eq=2192))]
rootfor102 = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=102),findcs(as03["subID"], eq=2203))]
rootfor100 = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=100),findcs(as03["subID"], eq=2249))]
rootfor96  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=96 ),findcs(as03["subID"], eq=2156))]
rootfor92  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=92 ),findcs(as03["subID"], eq=2095))]
rootfor88  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=88 ),findcs(as03["subID"], eq=2691))]
rootfor84  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=84 ),findcs(as03["subID"], eq=2677))]
rootfor80  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=80 ),findcs(as03["subID"], eq=3220))]
rootfor76  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=76 ),findcs(as03["subID"], eq=4841))]
rootfor72  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=72 ),findcs(as03["subID"], eq=4737))]
rootfor68  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=68 ),findcs(as03["subID"], eq=4117))]
rootfor64  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=64 ),findcs(as03["subID"], eq=3931))]
rootfor60  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=60 ),findcs(as03["subID"], eq=3525))]
rootfor58  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=58 ),findcs(as03["subID"], eq=3419))]
rootfor52  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=52 ),findcs(as03["subID"], eq=2875))]
rootfor48  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=48 ),findcs(as03["subID"], eq=2846))]
rootfor44  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=44 ),findcs(as03["subID"], eq=2250))]
rootfor40  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=40 ),findcs(as03["subID"], eq=2512))]
rootfor36  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=36 ),findcs(as03["subID"], eq=2133))]
rootfor32  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=32 ),findcs(as03["subID"], eq=1216))]
rootfor28  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=28 ),findcs(as03["subID"], eq=828 ))]
rootfor24  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=24 ),findcs(as03["subID"], eq=559 ))]
rootfor20  = as03["ID_ISUB"][idxmatch(findcs(as03["snapNR"], eq=20 ),findcs(as03["subID"], eq=187 ))]

# my halostory
plot_halostory_full(outfile=joinpath(outputdir, "halo_$(root_id)_halostory_full.pdf"), snap=true, rootID=root_id, title=true)
plot_halostory_full(outfile=joinpath(outputdir, "halo_$(24763)_halostory_full.pdf"), snap=true, rootID=24763, title=true)
plot_halostory_full(outfile=joinpath(outputdir, "halo_$(23709)_halostory_full.pdf"), snap=true, rootID=23709, title=true)

av = [Vec3f(v1,v2,v3) for (v1,v2,v3) in zip(LinRange(0, 10000, 10000) , LinRange(0, 20000, 10000) , LinRange(0, 30000, 10000))]
size(av)
#for (ias, ia1, ia2, ia3, ic) in zip(as, as_label1, as_label2, as_label3, csw)

resol = (1600,900)#(1920,1080)
h = 3968



s = 136
sts     = 10
plot_group(s; rootID=h, 
              property="vel", stepsize=sts, arsize=3, ptsize=5000, rad=300, res=resol, spin=true, 
              bgcolor=:white, indir=current_dir_stories, boxfix=false)
s



res = 5
plot_phasespace(; #snapNR = 128, 
    felixID=" ", rootID=3968, subID=" " ,
    parttype="STARS", r_by_rvir=0.1, #weights="MASS" ,
    lognorm=true, scale=15, gridres=(25*res, 15*res), colorlimits=(0.8,200), 
    outdir="/home/moon/sfortune/spinevo/plots/halo3968", outfile="PhSp", 
    )

#Vec3f.([0.2, 0.2], [0.2, 0.2], [0.3, 0.3])
