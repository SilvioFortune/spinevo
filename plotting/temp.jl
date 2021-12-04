
include("/home/moon/sfortune/spinevo/pkg/meta.jl")
assembly_STARS  = load("/home/moon/sfortune/spinevo/mergerimpact/assembly_Mstar_1.0e10.jld", "assembly_STARS")
centrals_136 = filter(assembly_STARS, condition="snap")


fig, ax = subplots()

ax.hist(assembly_STARS["M2_felix"][centrals_136["main"]], bins=100, label="M2_felix  at z ≈ $(@sprintf("%.2f", assembly_STARS["REDSHIFT"][centrals_136["main"]][1]))", rwidth=0.9, color="navy")


ax.set_xlabel("M_⊙")
ax.set_xscale("log")
ax.set_yscale("log")
ax.legend(loc="upper right", frameon=true, borderpad=1, handlelength=1.8)
ax.grid()

scale=0.5

fig.set_size_inches(16scale, 9scale)
fig.tight_layout()

println(minimum(assembly_STARS["M2_felix"][centrals_136["main"]]))
println(maximum(assembly_STARS["M2_felix"][centrals_136["main"]]))
println(Statistics.median(assembly_STARS["M2_felix"][centrals_136["main"]]))
#fig.savefig(joinpath(@__DIR__, "masses_z$(@sprintf("%.2f", head.z)).png"), bbox_inches="tight", pad_inches=.1)

m_smalls    = 10^10.5
m_intermediates = 10^11.5
smalls          = missings(Int64,0)
intermediates   = missings(Int64,0)
bigs            = missings(Int64,0)

for i in centrals_136["main"]
    if count(ismissing, assembly_STARS["J_SUMorbital"][:,i]) > 0 || count(ismissing, assembly_STARS["J_main"][:,i]) > 0
    elseif assembly_STARS["M2_felix"][i] < m_smalls
        smalls  = vcat(smalls, i)
    elseif assembly_STARS["M2_felix"][i] < m_intermediates
        intermediates  = vcat(intermediates, i)
    else
        bigs    = vcat(bigs,i)
    end
end
