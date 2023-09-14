using Revise
using SMLMSim
using CairoMakie

# run the smoluchowski simulation
box_size = 10
dt = .005
state_history = SMLMSim.InteractionDiffusion.smoluchowski(; 
    dt, box_size, t_max = 5.0, density = 10, d_dimer = 0.05);


function display_positions(states, time_step, box_size)
    x = [mol.x for mol in states.States[time_step]]
    y = [mol.y for mol in states.States[time_step]]
    colors = [mol.state == 2 ? :red : :blue for mol in states.States[time_step]]
    f = Figure()
    ax = Axis(f[1,1]; title="Positions at time step $time_step")
    # f[1,1] = Axis(f[1,1]; title="Positions at time step $time_step")
    scatter!(ax, x, y, markersize=10, color=colors)
    xlims!(ax, (0, box_size))
    ylims!(ax, (0, box_size))
    display(f)
end

display_positions(state_history, 1, box_size)



function generate_animation(states, box_size, filename, dt)
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="x", ylabel="y", 
    limits=(0, box_size, 0, box_size), 
    title="Smoluchowski Simulation")
    sc = scatter!(ax, 
        [mol.x for mol in states.States[1]], 
        [mol.y for mol in states.States[1]], 
        color=[mol.state == 2 ? :red : :blue for mol in states.States[1]], 
        markersize=10)
        display(sc)
    framerate = 1/dt
    timestamps = 1:length(states.States)
    record(fig, filename, timestamps; framerate=framerate) do i
        ax.title = "Positions at time step $i"
        sc[1] = [mol.x for mol in states.States[i]]
        sc[2] = [mol.y for mol in states.States[i]]
        sc.color = [mol.state == 2 ? :red : :blue for mol in states.States[i]]
        # scatter!(ax, sc.x, sc.y, color=sc.color, markersize=10)
    end
end

generate_animation(state_history, box_size, "smoluchowski.mp4", 0.01)


