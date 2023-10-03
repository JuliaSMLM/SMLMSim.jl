using Revise
using SMLMSim
using CairoMakie

# run the smoluchowski simulation
box_size = 10
dt = .005
state_history, args = SMLMSim.InteractionDiffusion.smoluchowski(; 
    dt=dt, box_size=box_size, t_max = 5.0, density = 10, d_dimer = 0.05);

SMLMSim.show_frame(state_history, 1, args)
    
SMLMSim.gen_movie(state_history, args)


box_size = 1.0
x = rand(10)
y = rand(10)
    f = Figure()
    ax = Axis(f[1,1]; 
    title="Positions at frame 1")
    scatter!(ax, x, y, markersize=10)
    xlims!(ax, (0, box_size))
    ylims!(ax, ( box_size, 0))
    display(f)


