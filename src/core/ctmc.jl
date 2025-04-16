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
    Base.show(io::IO, ctmc::CTMC)

Custom display method for CTMC showing basic properties.
"""
function Base.show(io::IO, ctmc::CTMC)
    num_transitions = length(ctmc.states) - 1
    print(io, "CTMC(time=$(round(ctmc.simulation_time, digits=2)), $(num_transitions) transitions)")
end

"""
    Base.show(io::IO, ::MIME"text/plain", ctmc::CTMC)

Extended display method for CTMC in REPL and other text contexts.
"""
function Base.show(io::IO, ::MIME"text/plain", ctmc::CTMC)
    num_transitions = length(ctmc.states) - 1
    unique_states = length(Set(ctmc.states))
    
    println(io, "CTMC (Continuous Time Markov Chain):")
    println(io, "  Simulation time: $(ctmc.simulation_time)")
    println(io, "  Number of transitions: $(num_transitions)")
    println(io, "  Number of unique states: $(unique_states)")
    println(io, "  Initial state: $(ctmc.states[1])")
    println(io, "  Final state: $(ctmc.states[end])")
    
    # Show detailed transition information if there are not too many
    if num_transitions <= 10
        println(io, "  Transitions:")
        for i in 1:num_transitions
            from_state = ctmc.states[i]
            to_state = ctmc.states[i+1]
            time = ctmc.transition_times[i+1]
            println(io, "    [$(round(time, digits=3))s] $(from_state) → $(to_state)")
        end
    end
end

"""
    CTMC(q::Array{T}, simulation_time::T, state1::Int) where {T<:AbstractFloat}

Construct a Continuous Time Markov Chain simulation from a rate matrix.

# Arguments
- `q::Array{T}`: Rate matrix where q[i,j] for i≠j is the transition rate from state i to j,
  and q[i,i] is the negative exit rate from state i
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
    if any(diag(q) .>= 0)
        @warn "Diagonal elements of rate matrix should be negative (representing exit rates)"
    end
    
    # Check for row sums approximately zero (rate matrix property)
    row_sums = sum(q, dims=2)
    if !all(abs.(row_sums) .< 1e-10)
        @warn "Rate matrix rows should sum to zero; matrix may not represent a valid CTMC"
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
        # get time for state change (negative of diagonal element)
        k_tot = -q[currentstate, currentstate]
        if k_tot ≈ 0
            # If total rate is zero, we're in an absorbing state
            break
        end
        
        Δt = rand(Exponential(1 / k_tot))

        # get the new state
        ps = copy(q[currentstate, :])
        ps[currentstate] = 0.0  # Zero out the diagonal element
        sum_ps = sum(ps)
        
        if sum_ps ≈ 0
            # No valid transitions
            break
        end
        
        ps = ps ./ sum_ps  # Normalize probabilities
        
        # Create a vector without the current state
        valid_indices = [i for i in 1:length(ps) if i != currentstate && ps[i] > 0]
        
        if isempty(valid_indices)
            break  # No valid transitions
        end
        
        valid_probs = ps[valid_indices]
        valid_probs = valid_probs ./ sum(valid_probs)  # Re-normalize
        
        # Sample the new state
        newstate = sample_discrete_with_probs(valid_indices, valid_probs)

        # update CTMC
        lastchange += Δt
        push!(states, newstate)
        push!(transition_times, lastchange)
        currentstate = newstate
    end
    return CTMC(simulation_time, transition_times, states)
end

"""
    sample_discrete_with_probs(indices, probs)

Sample a discrete value from indices with corresponding probabilities.

# Arguments
- `indices::Vector{Int}`: Vector of possible indices to sample
- `probs::Vector{<:AbstractFloat}`: Corresponding probabilities

# Returns
- `Int`: Sampled index
"""
function sample_discrete_with_probs(indices, probs)
    r = rand()
    cdf = 0.0
    for i in 1:length(indices)
        cdf += probs[i]
        if r <= cdf
            return indices[i]
        end
    end
    return indices[end]  # Fallback for numerical precision issues
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
function get_state(ctmc::CTMC{T}, t::T) where T <: AbstractFloat
    if t < zero(T) || t > ctmc.simulation_time
        throw(ArgumentError("Time point must be between 0 and simulation_time"))
    end
    
    # Use binary search to find the correct time interval
    idx = searchsortedlast(ctmc.transition_times, t)
    return ctmc.states[max(1, idx)]
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
function get_next(ctmc::CTMC{T}, t::T) where T <: AbstractFloat
    if t < zero(T) || t > ctmc.simulation_time
        throw(ArgumentError("Time point must be between 0 and simulation_time"))
    end
    
    # Use binary search to find the next transition
    idx = searchsortedfirst(ctmc.transition_times, t)
    
    # Adjust index if t is exactly at a transition time
    if idx <= length(ctmc.transition_times) && ctmc.transition_times[idx] == t
        idx += 1
    end

    if idx > length(ctmc.transition_times)
        # No more transitions after t
        return ctmc.states[end], ctmc.simulation_time
    else
        return ctmc.states[idx], ctmc.transition_times[idx]
    end
end
