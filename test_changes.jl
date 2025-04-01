# Test script to verify our changes
println("Loading packages...")
using SMLMSim
using SMLMData

# Check if StaticSMLMParams is exported correctly
println("Checking if StaticSMLMParams is available from main module...")
params = StaticSMLMParams(ρ=1.0)
println("StaticSMLMParams successfully created with ρ = ", params.ρ)

# Check if get_state and get_next are available
println("Checking if get_state and get_next are available...")
q = [0.0 1.0; 0.1 0.0]
ctmc = CTMC(q, 10.0, 1)
state = get_state(ctmc, 0.5)
next_state, next_time = get_next(ctmc, 0.5)
println("CTMC state at t=0.5: ", state)
println("Next state after t=0.5: ", next_state, " at time ", next_time)

println("All checks passed!")