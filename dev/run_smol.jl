using Revise
using SMLMSim
using CairoMakie

# run the smoluchowski simulation

state_history = SMLMSim.InteractionDiffusion.smoluchowski()


function display_positions(states, time_step)
    x = [mol.x for mol in states.States[time_step]]
    y = [mol.y for mol in states.States[time_step]]
    colors = [mol.state == 2 ? :red : :blue for mol in states.States[time_step]]
    f = Figure()
    f[1,1]= Axis(f[1,1]; title="Positions at time step $time_step")
    scatter!(f[1,1], x, y, markersize=10, color=colors)
    display(f)
end

display_positions(state_history, 200)



# function animate_positions(states::SMLMSim.InteractionDiffusion.MoleculeStates, filename::AbstractString)
#     n_time_steps = length(states.States)
#     positions = [[mol.x mol.y mol.z] for mol in states.States[1]]
#     colors = [mol.state == 2 ? :red : :blue for mol in states.States[1]]
#     fig = Figure(resolution=(800, 800))
#     scatter!(fig[1, 1], positions[:, 1], positions[:, 2], positions[:, 3], markersize=10, color=colors)
#     xlabel!(fig[1, 1], "X")
#     ylabel!(fig[1, 1], "Y")
#     zlabel!(fig[1, 1], "Z")
#     title!(fig[1, 1], "Positions at time step 1")
#     for i in 2:n_time_steps
#         positions = [mol.x mol.y mol.z for mol in states.States[i]]
#         colors = [mol.state == 2 ? :red : :blue for mol in states.States[i]]
#         scatter!(fig[1, 1], positions[:, 1], positions[:, 2], positions[:, 3], markersize=10, color=colors)
#         title!(fig[1, 1], "Positions at time step $i")
#         Makie.save(filename * "_$(i-1).png", fig)
#     end
# end