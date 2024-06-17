# Tools for analyzing dimers

#get molecules that are in linked states
"""

    get_dimers(states::MoleculeHistory)

create a array of vector of type MoleculeFrame that takes the frames in a the input states and pass it to the MoleculeHistory function

# Arguments
    - states::MoleculeHistory  : the estates of the molecule

# Returns
    - MoleculeHistory(states.dt, dimer_frames)  : return the function MoleculeHistory

"""
function get_dimers(states::MoleculeHistory)
    dimer_frames = Vector{MoleculeFrame}()
    for frame in states.frames
        push!(dimer_frames, get_dimers(frame))
    end
    return MoleculeHistory(states.dt, dimer_frames)
end


function get_dimers(frame::MoleculeFrame)
    linked_molecules = Vector{Monomer}()
    for molecule in frame.molecules
        if molecule.state == 2
            push!(linked_molecules, molecule)
        end
    end
    return MoleculeFrame(frame.framenum, linked_molecules)
end
