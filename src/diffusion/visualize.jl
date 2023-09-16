

"""
    show_frame(states::MoleculeStates, framenum::Int64, args::ArgsSmol)

Display a scatter plot of the positions of all molecules in the system at a given time step.

# Arguments
- `states::MoleculeHistory`: a `MoleculeHistory` object containing the states of the system at each time step
- `framenum::Int64`: the time step to display
- `args::ArgsSmol`: a struct containing the simulation parameters

# Returns
- `nothing`

"""
function show_frame(states::MoleculeHistory, framenum, args)
    box_size = args.box_size
    x = [mol.x for mol in states.frames[framenum].molecules]
    y = [mol.y for mol in states.frames[framenum].molecules]
    colors = [mol.state == 2 ? :red : :blue for mol in states.frames[framenum].molecules]
    f = Figure()
    ax = Axis(f[1,1]; title="Positions at frame $framenum")
    scatter!(ax, x, y, markersize=10, color=colors)
    xlims!(ax, (0, box_size))
    ylims!(ax, (0, box_size))
    display(f)
end


"""
gen_movie(states::MoleculeStates, args::ArgsSmol; filename::String="smoluchowski.mp4")

Generate an animation of the positions of all molecules in the system over time.

# Arguments
- `states::MoleculeHistory`: a `MoleculeHistory` object containing the states of the system at each time step
- `args::ArgsSmol`: a struct containing the simulation parameters
- `filename::String`: the name of the output file (default: "smoluchowski.mp4")

# Returns
- `nothing`

"""
function gen_movie(states::MoleculeHistory, args; filename::String="smoluchowski.mp4")
    dt = args.dt
    box_size = args.box_size
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="x", ylabel="y", 
    limits=(0, box_size, 0, box_size), 
    title="Smoluchowski Simulation")
    sc = scatter!(ax, 
        [mol.x for mol in states.frames[1].molecules],
        [mol.y for mol in states.frames[1].molecules], 
        color=[mol.state == 2 ? :red : :blue for mol in states.frames[1].molecules], 
        markersize=10)
        display(sc)
    framerate = Int(round(1/dt))
    timestamps = 1:length(states.frames)
    record(fig, filename, timestamps; framerate=framerate) do i
        ax.title = "Positions at time step $i"
        sc[1] = [mol.x for mol in states.frames[i].molecules]
        sc[2] = [mol.y for mol in states.frames[i].molecules]
        sc.color = [mol.state == 2 ? :red : :blue for mol in states.frames[i].molecules]
    end
    return nothing
end