
"""
CTMC

Continous Time Markov Chain    
"""
mutable struct CTMC{T<:AbstractFloat,U<:Int}
    τ::T
    transitiontimes::Vector{T}
    states::Vector{U}
    # CTMC() = new()
end


function CTMC(q::Array{T}, τ::T, state1::Int) where {T<:AbstractFloat}

    # create empty CTMC structure
    # ctmc = CTMC()
    # ctmc.τ = τ

    lastchange = 0.0
    currentstate = state1

    # ctmc.states = [state1]
    # ctmc.transitiontimes = [lastchange]

    states = [state1]
    transitiontimes = [lastchange]


    sidx = Vector((1:size(q, 1)))

    while lastchange < τ

        # get time for state change
        k_tot = sum(q[currentstate, :])
        Δt = rand(Exponential(1 / k_tot))

        # get the new state
        ps = q[currentstate, :]
        ps = ps ./ sum(ps)
        deleteat!(ps, currentstate)
        xs = sidx[:]
        deleteat!(xs, currentstate)
        newstate = rand(DiscreteNonParametric(xs, ps))

        # update CTMC
        lastchange += Δt
        push!(states, newstate)
        push!(transitiontimes, lastchange)
        currentstate = newstate
    end
    return CTMC(τ, transitiontimes, states)
end


"""
    getstate(ctmc::CTMC,t::AbstractFloat)

 return the state at time t   
"""
function getstate(ctmc::CTMC, t::AbstractFloat)
    nstates = length(ctmc.states)
    for nn = 2:nstates
        if t < ctmc.transitiontimes[nn]
            return ctmc.states[nn-1]
        end
    end
end

"""
    getnext(ctmc::CTMC,t::AbstractFloat)

 return the time and state of next transision
"""
function getnext(ctmc::CTMC, t::AbstractFloat)
    nstates = length(ctmc.states)
    for nn = 1:nstates
        if t < ctmc.transitiontimes[nn]
            return ctmc.states[nn], ctmc.transitiontimes[nn]
        end
    end
end



"""
    intensitytrace(f::GenericFluor, nframes::Int, framerate::AbstractFloat;state1=1)    

Calculate an intensity trace.     

"""
function intensitytrace(f::GenericFluor, nframes::Int, framerate::Real; state1=1)

    endtime = (nframes) / framerate

    # generate CTMC
    ctmc = SMLMSim.CTMC(f.q, endtime, state1)

    # generate integrated photons 
    exptime = 1 / framerate
    photons = zeros(nframes)
    for nn = 1:nframes
        t = (nn - 1) * exptime
        frameend = (nn * exptime)
        while t < frameend
            currentstate = SMLMSim.getstate(ctmc, t)
            (_, nexttime) = SMLMSim.getnext(ctmc, t)
            tend = min(frameend, nexttime)
            if currentstate == 1 # the fluorescent state 
                photons[nn] += f.γ * (tend - t)
            end
            t = tend
        end
    end
    return photons
end

"""
    function kineticmodel(smd_true::SMLMData.SMLD,f::Molecule,nframes::Int,framerate::AbstractFloat;ndatasets::Int=1,minphotons=50.0)

generate noise-free blinking model from smd_true
"""
function kineticmodel end

function kineticmodel(y_true::Vector{<:Real}, x_true::Vector{<:Real}, f::Molecule, nframes::Int, framerate::Real; ndatasets::Int=1, minphotons=50.0)

    state1 = 2

    smd = SMLMData.SMLD2D(0)
    smd.ndatasets = ndatasets
    smd.nframes = nframes
    
    for dd = 1:ndatasets, ll = 1:length(y_true)
        photons = SMLMSim.intensitytrace(f, nframes, framerate; state1=state1)
        framenum = findall(photons .> minphotons)
        n = length(framenum)
        push!(smd.photons, photons[framenum]...)
        push!(smd.framenum, framenum...)
        for nn = 1:n
            push!(smd.x, x_true[ll])
            push!(smd.y, y_true[ll])
            push!(smd.connectID, ll)
            push!(smd.datasetnum, dd)
        end
    end
    return smd.y, smd.x, smd.photons, smd.framenum, smd.datasetnum, smd.connectID
end

function kineticmodel(y_true::Vector{<:Real}, x_true::Vector{<:Real}, z_true::Vector{<:Real},  
    f::Molecule, nframes::Int, framerate::Real; ndatasets::Int=1, minphotons=50.0)

    state1 = 2

    smd = SMLMData.SMLD3D(0)
    smd.ndatasets = ndatasets
    smd.nframes = nframes
    
    for dd = 1:ndatasets, ll = 1:length(x_true)
        photons = SMLMSim.intensitytrace(f, nframes, framerate; state1=state1)
        framenum = findall(photons .> minphotons)
        n = length(framenum)
        push!(smd.photons, photons[framenum]...)
        push!(smd.framenum, framenum...)
        for nn = 1:n
            push!(smd.x, x_true[ll])
            push!(smd.y, y_true[ll])
            push!(smd.z, z_true[ll])
            push!(smd.connectID, ll)
            push!(smd.datasetnum, dd)
        end
    end
    return smd.y, smd.x, smd.z, smd.photons, smd.framenum, smd.datasetnum, smd.connectID
end


"""
    noise(smd_model::SMLMData.SMLD2D,σ_psf::AbstractFloat) 

Add zero mean Gaussian noise to coordinates with σ = σ_pdf/sqrt(photons) 
"""
function noise end

function noise(smd_model::SMLMData.SMLD2D, σ_psf::AbstractFloat)
    n = length(smd_model.x)
    smd = deepcopy(smd_model)
    smd.σ_x = zeros(n)
    smd.σ_y = zeros(n)

    for nn = 1:n
        σ = σ_psf / sqrt(smd_model.photons[nn])
        smd.x[nn] = smd_model.x[nn] + randn() * σ
        smd.y[nn] = smd_model.y[nn] + randn() * σ
        smd.σ_x[nn] = σ
        smd.σ_y[nn] = σ
    end
    return smd
end

""" 
    noise(smd_model::SMLMData.SMLD3D,σ_psf::Vector{<:AbstractFloat})
3D data requries `σ_psf = [σ_x,σ_y,σ_z]`
"""
function noise(smd_model::SMLMData.SMLD3D, σ_psf::Vector{<:AbstractFloat})
    n = length(smd_model.x)
    smd = deepcopy(smd_model)
    smd.σ_x = zeros(n)
    smd.σ_y = zeros(n)
    smd.σ_z = zeros(n)

    for nn = 1:n
        σ = σ_psf ./ sqrt(smd_model.photons[nn])
        smd.x[nn] = smd_model.x[nn] + randn() * σ[1]
        smd.y[nn] = smd_model.y[nn] + randn() * σ[2]
        smd.z[nn] = smd_model.z[nn] + randn() * σ[3]

        smd.σ_x[nn] = σ[1]
        smd.σ_y[nn] = σ[2]
        smd.σ_y[nn] = σ[3]
    end
    return smd
end
