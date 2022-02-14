
"""
    Camera 

"""
abstract type Camera end 



"""
    IdealCamera <: Camera

A camera with no added noise. 

#Fields
- pixelsize 
- xpixels
- ypixels
- gain
- offset

```
IdealCamera(;
pixelsize=0.1,
xpixels::Int=256,
ypixels::Int=256,
gain=1.0,
offset=0.0,
)

```
"""
mutable struct IdealCamera <: Camera
    pixelsize
    xpixels
    ypixels
    gain
    offset
end
function IdealCamera(;
    pixelsize=0.1,
    xpixels::Int=256,
    ypixels::Int=256,
    gain=1.0,
    offset=0.0,
    )
    return IdealCamera(pixelsize,xpixels,ypixels,gain,offset)
end





