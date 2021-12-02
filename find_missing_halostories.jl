include("/home/moon/sfortune/spinevo/pkg/meta.jl")

trees       = readdir("/home/moon/sfortune/spinevo/newtrees")
storyfiles  = readdir("/home/moon/sfortune/spinevo/halostories_v20211127_min0.0Gyr")

for i in 1:length(trees)
    treefile    = replace( trees[i], ".dat" => "_" )
    included    = false
    for ii in 1:length(storyfiles)
        if occursin(treefile, storyfiles[ii])
            included    = true
        end
    end
    if !included
        println(trees[i])
    end
end