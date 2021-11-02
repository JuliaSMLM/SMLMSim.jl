
```
    Camera 

```
abstract type Camera end 



```
    IdealCamera



```

mutable struct IdealCamera <: Camera
    pixelsize
    xpixels
    ypixels
    gain=1
    offset=0
end




