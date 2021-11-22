
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
    function kineticmodel(smd_true::SMLMData.SMLD2D,f::Molecule,nframes::Int,framerate::AbstractFloat;ndatasets::Int=1,minphotons=50.0)

generate noise-free blinking model from smd_true
"""
function kineticmodel(smd_true::SMLMData.SMLD2D,f::Molecule,nframes::Int,framerate::AbstractFloat;ndatasets::Int=1,minphotons=50.0)

    state1=2;

    smd=SMLMData.SMLD2D(0)
    smd.ndatasets=ndatasets
    smd.nframes=nframes
    for dd=1:ndatasets, ll=1:length(smd_true.x)
        photons=SMLMSim.intensitytrace(f,nframes,framerate;state1=state1)    
        framenum=findall(photons.>minphotons)
        n=length(framenum)
        push!(smd.photons,photons[framenum]...)
        push!(smd.framenum,framenum...)
        for nn=1:n
            push!(smd.x,smd_true.x[ll])
            push!(smd.y,smd_true.y[ll])
            push!(smd.datasetnum,dd)
        end
    end
    return smd
end

function noise(smd_model::SMLMData.SMLD2D,σ_psf::AbstractFloat)
    n=length(smd_model.x)
    smd=deepcopy(smd_model)
    smd.σ_x=zeros(n)
    smd.σ_y=zeros(n)

    for nn=1:n
        σ=σ_psf/sqrt(smd_model.photons[nn])
        smd.x[nn]=smd_model.x[nn]+randn()*σ
        smd.y[nn]=smd_model.y[nn]+randn()*σ
        smd.σ_x[nn]=σ
        smd.σ_y[nn]=σ
    end
    return smd
end

