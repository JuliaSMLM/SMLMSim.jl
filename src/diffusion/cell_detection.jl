include("neural_network.jl")

module CellDetection

export detect_cells, detect_monomers, detect_dimers

using .NeuralNetwork: create_model, detect_with_model

model = create_model()

function detect_cells(image)
    println("Detecting cells using neural network...")
    detect_with_model(model, image)
end

function detect_monomers(image)
    println("Detecting monomers using neural network...")
    detect_with_model(model, image)
end

function detect_dimers(image)
    println("Detecting dimers using neural network...")
    detect_with_model(model, image)
end

end