# Generate uniform random data

using SMLMSim
using SMLMData
using GLMakie

## create a blinking fluorophore 
γ=1e5
q=[0 50
   1e-2 0]
f=SMLMSim.GenericFluor(γ,q)

# pattern is a single molecule 
pattern=SMLMSim.Point3D()
ρ=1.0 
xsize=25.6
ysize=25.6
smd_true=SMLMSim.uniform3D(ρ,pattern,xsize,ysize)

fig=Figure()
ax=Axis3(fig[1, 1], xlabel = "x", ylabel = "y",zlabel = "z", title = "3D Localizations")
scatter!(smd_true.x, smd_true.y,smd_true.z, color=:black)

# generate the kinetic model
nframes=2000
framerate=50.0
smd_model=SMLMSim.kineticmodel(smd_true,f,nframes,framerate;ndatasets=10)

# make noisy coordinates
σ_psf=.13 
smd_noisy=SMLMSim.noise(smd_model,[σ_psf,σ_psf,2*σ_psf])
scatter!(smd_noisy.x, smd_noisy.y,smd_noisy.z, color=:green)
fig




