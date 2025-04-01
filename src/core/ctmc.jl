"""
    CTMC - Continuous Time Markov Chain Module

This module provides a general implementation of Continuous Time Markov Chains
that can be used for simulating stochastic processes in various contexts including
fluorophore blinking, conformational dynamics, and other state transitions.
"""

# Using statements moved to module level (Core.jl)

"""
    CTMC{T<:AbstractFloat, U<:Int}

A Continuous Time Markov Chain representation storing the full trajectory of state transitions.

# Fields
- `simulation_time::T`: Total simulation time span
- `transition_times::Vector{T}`: Time points at which state changes occurred, starting at 0.0 
- `states::Vector{U}`: Sequence of states entered at each transition time, starting with initial state

# Type Parameters
- `T`: Floating point type for time values
- `U`: Integer type for state indices

# Note
The `states` and `transition_times` vectors have the same length, with each entry in `states[i]`
representing the state entered at time `transition_times[i]`. The system remains in `states[i]` 
until time `transition_times[i+1]`.
"""
mutable struct CTMC{T<:AbstractFloat,U<:Int}
    simulation_time::T
    transition_times::Vector{T}
    states::Vector{U}
end

"""
    CTMC(q::Array{T}, simulation_time::T, state1::Int) where {T<:AbstractFloat}

Construct a Continuous Time Markov Chain simulation from a rate matrix.

# Arguments
- `q::Array{T}`: Rate matrix where q[i,j] is the transition rate from state i to j
- `simulation_time::T`: Total time to simulate
- `state1::Int`: Initial state

# Returns
- `CTMC{T,Int}`: Simulated CTMC with transition times and states

# Details
Simulates a CTMC using the Gillespie algorithm:
1. Start in state1 at time 0
2. For current state i:
   - Calculate total exit rate k_tot = Σ_j q[i,j]
   - Sample time until next transition from Exp(k_tot)
   - Sample next state j with probability q[i,j]/k_tot
3. Repeat until exceeding simulation_time
"""
function CTMC(q::Array{T}, simulation_time::T, state1::Int) where {T<:AbstractFloat}
    # Validate inputs
    if size(q, 1) != size(q, 2)
        throw(ArgumentError("Rate matrix q must be square"))
    end
    if any(diag(q) .!= 0)
        @warn "Diagonal elements of rate matrix should be zero"
    end
    if state1 < 1 || state1 > size(q, 1)
        throw(ArgumentError("Initial state must be between 1 and size(q,1)"))
    end
    
    lastchange = zero(T)
    currentstate = state1

    states = [state1]
    transition_times = [lastchange]

    sidx = Vector((1:size(q, 1)))

    while lastchange < simulation_time
        # get time for state change
        k_tot = sum(q[currentstate, :])
        if k_tot ≈ 0
            # If total rate is zero, we're in an absorbing state
            break
        end
        
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
        push!(transition_times, lastchange)
        currentstate = newstate
    end
    return CTMC(simulation_time, transition_times, states)
end

"""
    get_state(ctmc::CTMC, t::AbstractFloat)

Get the state of the CTMC at a specific time point.

# Arguments
- `ctmc::CTMC`: The CTMC to query
- `t::AbstractFloat`: Time point of interest

# Returns
- `Int`: State of the chain at time t

# Note
Searches through transition times to find the state active at time t.
Returns the state that was entered at the last transition before t.
"""
function get_state(ctmc::CTMC, t::AbstractFloat)
    if t < 0 || t > ctmc.simulation_time
        throw(ArgumentError("Time point must be between 0 and simulation_time"))
    end
    
    nstates = length(ctmc.states)
    for nn = 2:nstates
        if t < ctmc.transition_times[nn]
            return ctmc.states[nn-1]
        end
    end
    
    # If we get here, t is after the last transition
    return ctmc.states[end]
end

"""
    get_next(ctmc::CTMC, t::AbstractFloat)

Get the next state transition after a specific time point.

# Arguments
- `ctmc::CTMC`: The CTMC to query
- `t::AbstractFloat`: Current time point

# Returns
- `Tuple{Int,AbstractFloat}`: (next_state, transition_time)

# Note
Returns the next state that will be entered and when it will be entered,
searching from the current time point forward.
"""
function get_next(ctmc::CTMC, t::AbstractFloat)
    if t < 0 || t > ctmc.simulation_time
        throw(ArgumentError("Time point must be between 0 and simulation_time"))
    end
    
    nstates = length(ctmc.states)
    for nn = 1:nstates
        if t < ctmc.transition_times[nn]
            return ctmc.states[nn], ctmc.transition_times[nn]
        end
    end
    
    # If we get here, t is after the last transition
    return ctmc.states[end], ctmc.simulation_time
end
