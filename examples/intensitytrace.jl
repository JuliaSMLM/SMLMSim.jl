# generate an intensity trace from a blinking fluorophore

using SMLMSim
using GLMakie

## create a blinking fluorophore 
γ = 10000.0
q = [0 10
   1e-1 0]
f = SMLMSim.GenericFluor(γ, q)

## Simulate intensity trace 

# setup sim
nframes = 1000
framerate = 10.0

## generate CTMC
starttime = 0.0
endtime = (nframes) / framerate
state1 = 1
ctmc = SMLMSim.CTMC(q, endtime, state1)

##  generate integrated photons 
photons = SMLMSim.intensitytrace(f, nframes, framerate)
fig = Figure()
ax = Axis(fig[1, 1],
   xlabel="time (frames)",
   ylabel="Intensity (photons)"
)
lines!((1:nframes), photons)
fig

