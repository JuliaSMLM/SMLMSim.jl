# test_photophysics.jl - Demonstrate fluorophore photophysics and plot intensity trace

using Pkg
Pkg.activate("dev")  # Activate the dev project environment

using SMLMSim
using GLMakie
using SMLMSim.Core: CTMC, GenericFluor, intensity_trace
# using LinearAlgebra   # For matrix operations
# using Distributions   # For sampling from exponential distribution

println("Simulating fluorophore photophysics...")

# Simulation parameters
nframes = 1000           # Total frames to simulate
framerate = 10.0         # 10 frames per second
total_time = nframes / framerate  # Total time in seconds

# Kinetics parameters
k_off = framerate        # Off state rate (Hz) - corresponds to 1 frame duration
n_blinks_per_fluor = 10  # Desired number of events per 1000 frames
k_on = n_blinks_per_fluor / total_time  # On state rate (Hz)

println("Fluorophore parameters:")
println("- Off rate: $k_off Hz (average on time: $(1/k_off) seconds)")
println("- On rate: $k_on Hz (average off time: $(1/k_on) seconds)")
println("- Expected events per $nframes frames: $n_blinks_per_fluor")

# Create fluorophore model with on/off blinking
# Rate matrix: [on→off on←off; off→on off←on]
fluor = GenericFluor(
    1e3,                        # Photon emission rate (Hz)
    [-k_off k_off; k_on -k_on]  # Rate matrix
)

# Directly simulate intensity trace using the photophysics functions
# Start in off state (state=2) so we see clear on transitions
photons = intensity_trace(fluor, nframes, framerate; state1=2)

# Convert to binary state values for plotting
# State is 'on' when photons > 0
states = zeros(nframes)
for i in 1:nframes
    states[i] = photons[i] > 0 ? 1.0 : 0.0
end

# Calculate time points in seconds
times = (1:nframes) ./ framerate

# Count number of on events (transitions from off to on)
on_events = count(diff(states) .== 1)
println("Number of observed on events: $on_events")

# Plot the intensity trace
fig = Figure(size=(900, 500))
ax = Axis(fig[1, 1], 
    title="Fluorophore State vs Time", 
    xlabel="Time (seconds)", 
    ylabel="State (0=off, 1=on)",
    yticklabelsize=14)

# Plot the binary states
lines!(ax, times, states, linewidth=2, color=:blue)

# Add annotation for rates
text!(ax, 0.05, 0.95, 
    text="k_off = $(round(k_off, digits=2)) Hz\nk_on = $(round(k_on, digits=2)) Hz\nEvents: $on_events",
    align=(:left, :top), 
    space=:relative, 
    fontsize=14)

display(fig)

println("\nIntensity trace displayed. Close the window to exit.")

# Optionally, save the figure
# save("fluorophore_trace.png", fig)