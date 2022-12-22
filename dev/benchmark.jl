using Revise
using SMLMSim
using Plots  
using BenchmarkTools

## create a blinking fluorophore 
γ=10000.0
q=[0 10
   1e-1 0]
f=SMLMSim.GenericFluor(γ,q)

## Simulate intensity trace 

# setup sim
nframes = 100000
framerate = 10.0

## generate CTMC
starttime=0.0
endtime = (nframes) / framerate
state1=1
ctmc=SMLMSim.CTMC(q,endtime,state1)



@btime SMLMSim.intensitytrace(f,nframes,framerate)

@profview SMLMSim.intensitytrace(f,nframes,framerate)


