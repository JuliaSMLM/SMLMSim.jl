module CellDetection

include("neural_network.jl")

export CellDetector, detect_cells, detect_monomers, detect_dimers

using .NeuralNetwork: create_model, detect_with_model

struct CellDetector
    model
    state
end

function CellDetector()
    model, state = create_model()
    new(model, state)
end

function detect_cells(detector::CellDetector, image)
    println("Detecting cells using neural network...")
    detect_with_model(detector.model, detector.state, image)
end

function detect_monomers(detector::CellDetector, image)
    println("Detecting monomers using neural network...")
    detect_with_model(detector.model, detector.state, image)
end

function detect_dimers(detector::CellDetector, image)
    println("Detecting dimers using neural network...")
    detect_with_model(detector.model, detector.state, image)
end

end