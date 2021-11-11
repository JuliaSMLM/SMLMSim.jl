
"""
    Camera 

"""
abstract type Camera end 



"""
    IdealCamera



"""
mutable struct IdealCamera <: Camera
    pixelsize
    xpixels
    ypixels
    gain
    offset
end




