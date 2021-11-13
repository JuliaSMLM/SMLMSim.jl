# Generate noisy coorindate data of Nmers

using Revise
using SMLMSim
using SMLMData
using PlotlyJS
using Plots 

## Simulation setup 

γ=1e5 # Fluorophore emission rate
q=[0 50
   1e-2 0] # Fluorophore blinking rates
n=6 # Nmer rank
d=.1 # Nmer diameter
ρ=0.1 # density of Nmers 
xsize=25.6 # image size
ysize=25.6
nframes=2000 # number of frames
framerate=50.0 # framerate
σ_psf=1.3 # psf sigma used for uncertainty calcs
minphotons=500 # minimum number of photons per frame accepted

# Simulation sequence
f=SMLMSim.GenericFluor(γ,q)
pattern=SMLMSim.Nmer2D(n,d)
smd_true=SMLMSim.uniform2D(ρ,pattern,xsize,ysize)
smd_model=SMLMSim.kineticmodel(smd_true,f,nframes,framerate;ndatasets=10,minphotons=minphotons)
smd_noisy=SMLMSim.noise(smd_model,σ_psf)

# Plotting
plt=PlotlyJS.plot(scattergl(x=smd_noisy.x, y=smd_noisy.y, mode="markers"))
display(plt)





