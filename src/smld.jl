# tools for working with SMLMData SMLD structures

import SMLMData.SMLD2D 
import SMLMData.SMLD3D

# new constructors that use a camera structure
""" 
    function SMLD2D(cam::Camera, y_microns::Vector{T}, x_microns::Vector{T}, photons::Vector{T},
    framenum::Vector{Int}, datasetnum::Vector{Int}, connectID::Vector{Int} ) where {T<:Real}

Generate SMLD2D using camera structure and physical units.
"""
function SMLD2D(cam::Camera, y_microns::Vector{T}, x_microns::Vector{T}, photons::Vector{T},
    framenum::Vector{Int}, datasetnum::Vector{Int}, connectID::Vector{Int} ) where {T<:Real}
    
    smld = SMLD2D(length(y_microns))
    smld.x = x_microns ./ cam.pixelsize
    smld.y = y_microns ./ cam.pixelsize
    smld.photons = photons
    smld.framenum = framenum
    smld.datasetnum = datasetnum
    smld.connectID = connectID
    smld.datasize = [cam.ypixels, cam.xpixels]

    return smld
end


""" 
    function SMLD3D(cam::Camera, y_microns::Vector{T}, x_microns::Vector{T}, z_microns::Vector{T}, photons::Vector{T},
    framenum::Vector{Int}, datasetnum::Vector{Int}, connectID::Vector{Int} ) where {T<:Real}

Generate SMLD3D using camera structure and physical units.
"""
function SMLMData.SMLD3D(cam::Camera, y_microns::Vector{T}, x_microns::Vector{T}, z_microns::Vector{T}, photons::Vector{T},
    framenum::Vector{Int}, datasetnum::Vector{Int}, connectID::Vector{Int} ) where {T<:Real}
    
    smld = SMLD3D(length(y_microns))
    smld.x = x_microns ./ cam.pixelsize
    smld.y = y_microns ./ cam.pixelsize
    smld.z = z_microns 
    smld.photons = photons
    smld.framenum = framenum
    smld.datasetnum = datasetnum
    smld.connectID = connectID
    smld.datasize = [cam.ypixels, cam.xpixels]

    return smld
end

