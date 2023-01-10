println("   Starting import...")

current_dir_data   = "/home/moon/sfortune/spinevo/data"
current_dir_trees   = "/home/moon/sfortune/spinevo/data/newtrees_centrals_nobreak_physical"
current_dir_simbox  = "/HydroSims/Magneticum/Box4/uhr_test"
#current_dir_jld = "/home/moon/sfortune/spinevo/data/halostories_v20211127_min0.0Gyr"
current_dir_stories     = "/home/moon/sfortune/spinevo/data/centralstories_nobreak_physical_v20220127_min0.0Gyr"

include("/home/moon/sfortune/spinevo/pkg/packages.jl")
include("/home/moon/sfortune/spinevo/pkg/auxiliary.jl")
include("/home/moon/sfortune/spinevo/pkg/processing.jl")
include("/home/moon/sfortune/spinevo/pkg/physical.jl")
include("/home/moon/sfortune/spinevo/pkg/plotting.jl")

println("\n   DONE!")
