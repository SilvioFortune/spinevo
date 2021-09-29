### Settings

haloID      = 1
flipangle   = 180 #[°]

box         = "/HydroSims/Magneticum/Box4/uhr_test/"


### Packages

using Printf
using LinearAlgebra
using JLD
using QuadGK
using Statistics
using PyCall
using PyPlot
using LaTeXStrings
using GadgetIO
using GadgetUnits
using GadgetGalaxies
using Unitful
using UnitfulAstro
using Missings
using HypothesisTests
using Distributions
using CSV
using DataFrames
using Cosmology

np      = pyimport("numpy")
pltcm   = pyimport("matplotlib.cm")
plt     = pyimport("matplotlib.pyplot")
plt.rcParams["figure.figsize"] = [7.00, 3.50]
plt.rcParams["figure.autolayout"] = true
fig = figure()
ax = fig.add_subplot(projection="3d")
r = 0.05
#u, v = np.mgrid[0:2*π:30, 0:π:20]
#u, v = np.mgrid[LinRange(0,2*π,30), LinRange(0,π,20)]

x = cos.(u) .* sin.(v)
y = sin.(u) .* sin.(v)
z = cos.(v)
ax.plot_surface(x, y, z, cmap=pltcm.YlGnBu_r)
show()

