module SMLMSim

using SMLMData
using Distributions
using LinearAlgebra

include("typedefs.jl")
include("molecules.jl")
include("cameras.jl")
include("patterns.jl")
include("sim.jl")
include("smld.jl")
include("interface.jl")

# Submodules
include("diffusion/InteractionDiffusion.jl")
include("cell_detection.jl")
include("neural_network.jl")  # Include the NeuralNetwork module

using SMLMSim.InteractionDiffusion
using SMLMSim.CellDetection  # Use the CellDetection module
using SMLMSim.NeuralNetwork  # Use the NeuralNetwork module

end