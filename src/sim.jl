
"""
CTMC

Continous Time Markov Chain    
"""
mutable struct CTMC    
    τ::AbstractFloat
    transitiontimes::Vector{AbstractFloat}
    states::Vector{Int}
    CTMC() = new()
end


function CTMC(q::Array, τ::AbstractFloat,state1::Int)

# create empty CTMC structure
    ctmc = CTMC()   
    ctmc.τ=τ

    lastchange = 0.0 
    currentstate = state1

    ctmc.states = [state1]
    ctmc.transitiontimes = [lastchange]

    sidx = Vector((1:size(q, 1)))

    while lastchange < τ
    
# get time for state change
        k_tot = sum(q[currentstate,:])
        Δt = rand(Exponential(1 / k_tot))

# get the new state
        ps = q[currentstate,:]
        ps = ps ./ sum(ps)
        deleteat!(ps,currentstate)
        xs = sidx[:]
        deleteat!(xs,currentstate)
        newstate = rand(DiscreteNonParametric(xs, ps))

# update CTMC
        lastchange +=Δt
        push!(ctmc.states,newstate)
        push!(ctmc.transitiontimes,lastchange)
        currentstate = newstate
    end
    return ctmc
end


"""
    getstate(ctmc::CTMC,t::AbstractFloat)

 return the state at time t   
"""
function getstate(ctmc::CTMC,t::AbstractFloat)   
    nstates=length(ctmc.states)
    for nn=2:nstates
        if t<ctmc.transitiontimes[nn]
            return ctmc.states[nn-1]
        end
    end
end

"""
    getnext(ctmc::CTMC,t::AbstractFloat)

 return the time and state of next transision
"""
function getnext(ctmc::CTMC,t::AbstractFloat)   
    nstates=length(ctmc.states)
    for nn=1:nstates
        if t<ctmc.transitiontimes[nn]
            return ctmc.states[nn],ctmc.transitiontimes[nn]
        end
    end
end



"""
    intensitytrace(f::GenericFluor, nframes::Int, framerate::AbstractFloat;state1=1)    

Calculate an intensity trace.     

"""
function intensitytrace(f::GenericFluor, nframes::Int, framerate::AbstractFloat;state1=1)

    endtime = (nframes) / framerate

    # generate CTMC
    ctmc=SMLMSim.CTMC(f.q,endtime,state1)

    # generate integrated photons 
    exptime=1/framerate
    photons=zeros(nframes)
    for nn=1:nframes
        t=(nn-1)*exptime
        frameend=(nn*exptime)
        while t<frameend
            currentstate=SMLMSim.getstate(ctmc,t)
            (~,nexttime)=SMLMSim.getnext(ctmc,t)
            tend=min(frameend,nexttime)
            if currentstate==1 # the fluorescent state 
                photons[nn]+=f.γ*(tend-t)
            end
            t=tend    
        end
    end
    return photons
end

function genroistack(psf::MicroscopePSFs.PSF, xsize, ysize, nframes, framerate)


end

"""
    updatestate(state::Int,q::Matrix)

Update the state and return the elapsed time.    


q has format of i,j where i is current state, j is rate to another state. 
"""
function updatestate(state::Int, q::Matrix)
       
    
    

end
