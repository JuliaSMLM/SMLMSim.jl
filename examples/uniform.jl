# Generate uniform random data

using SMLMSim
using SMLMData
using GLMakie


## create a blinking fluorophore 
γ = 1e5
q = [0 50
   1e-2 0]
f = SMLMSim.GenericFluor(γ, q)

# pattern is a single molecule 
pattern = SMLMSim.Point2D()
ρ = 1.0
xsize = 25.6
ysize = 25.6
smd_true = SMLMSim.uniform2D(ρ, pattern, xsize, ysize)

fig = Figure()
ax = Axis(fig[1, 1],
   xlabel="x",
   ylabel="y",
   aspect=AxisAspect(1)
)
scatter!(smd_true.x, smd_true.y; color=:red)

# generate the kinetic model
nframes = 2000
framerate = 50.0
smd_model = SMLMSim.kineticmodel(smd_true, f, nframes, framerate; ndatasets=10)

ax2 = Axis(fig[1, 2],
   xlabel="Photons",
   ylabel="Counts"
)
hist!(smd_model.photons)

# make noisy coordinates
σ_psf = 1.3
smd_noisy = SMLMSim.noise(smd_model, σ_psf)
scatter!(ax, smd_noisy.x, smd_noisy.y; color=:green)

#show figure
fig







