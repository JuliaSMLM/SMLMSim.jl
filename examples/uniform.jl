# Generate uniform random data

using Revise
using SMLMSim
using SMLMData
using PlotlyJS
using Plots 

## create a blinking fluorophore 
γ=1e5
q=[0 50
   1e-2 0]
f=SMLMSim.GenericFluor(γ,q)

# pattern is a single molecule 
pattern=SMLMSim.Point2D()
ρ=1.0 
xsize=25.6
ysize=25.6
smd_true=SMLMSim.uniform2D(ρ,pattern,xsize,ysize)
plt=PlotlyJS.plot(scattergl(x=smd_true.x, y=smd_true.y, mode="markers"))
display(plt)

# generate the kinetic model
nframes=2000
framerate=50.0
smd_model=SMLMSim.kineticmodel(smd_true,f,nframes,framerate;ndatasets=10)
plt=Plots.histogram(smd_model.photons;xlabel="photons",ylabel="counts")
display(plt)

# make noisy coordinates
σ_psf=1.3
smd_noisy=SMLMSim.noise(smd_model,σ_psf)
plt=PlotlyJS.plot(scattergl(x=smd_noisy.x, y=smd_noisy.y, mode="markers"))
display(plt)





