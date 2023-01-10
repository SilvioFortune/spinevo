println("   Starting import...")

current_dir_data    = "/home/moon/sfortune/spinevo/data"
current_dir_trees   = "/home/moon/sfortune/spinevo/data/newtrees_centrals_nobreak_physical"
current_dir_simbox  = "/HydroSims/Magneticum/Box4/uhr_test"
current_dir_stories = "/home/moon/sfortune/spinevo/data/silvio_stories_spinmap/"
current_root_list   = "/home/moon/sfortune/spinevo/data/root_list_20220306.jld"

little_h0   = 0.704
box4_fullsnaplist = 4:136
box4_snaplist = [12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 58, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 102, 104, 106, 108, 112, 116, 120, 124, 128, 132, 136]
b_disk  = -4.357
b_ell   = -4.732

include("/home/moon/sfortune/spinevo/pkg/packages.jl")
include("/home/moon/sfortune/spinevo/pkg/auxiliary.jl")
include("/home/moon/sfortune/spinevo/pkg/processing.jl")
include("/home/moon/sfortune/spinevo/pkg/physical.jl")
include("/home/moon/sfortune/spinevo/pkg/plotting.jl")

# processed variables
pv              = load("/home/moon/sfortune/spinevo/pkg/processed_variables.jld")
box4_fulllbt    = pv["box4_fulllbt"]
box4_lbt        = pv["box4_lbt"]

println("\nLoaded SpinEvo package at $(now())\n")

include("/home/moon/sfortune/spinevo/pkg/misc.jl")
println("\n")