#define patterned stuctures

"""
    Pattern

Abstract type for structured patterns of molecules    
"""
abstract type Pattern end


"""
    Nmer <: Pattern

"""
mutable struct Nmer <: Pattern
    d::AbstractFloat
    
end


"""
    Uniform <: Pattern

Generate data with uniform random distribution.    

# Fields
- `ρ`: molecule density

"""
mutable struct Point2D <: Pattern
    n::Int
    x::Vector{AbstractFloat}
    y::Vector{AbstractFloat}
    Point2D()=new(1,[0.0],[0.0])
end

"""
    function uniformPattern2D(ρ,xsize::AbstractFloat,ysize::AbstractFloat, p::Pattern)

create true positions of molecules from uniformly randomly placed patterns
"""
function uniform2D(ρ,p::Pattern,xsize::AbstractFloat,ysize::AbstractFloat)

npatterns=rand(Poisson(xsize*ysize*ρ))    
ntotal=npatterns*p.n

θ=2*pi*rand()

#make smd 
smd=SMLMData.SMLD2D(ntotal)

for nn=1:npatterns, mm=1:p.n
    idx=(p.n)*(nn-1)+mm
    smd.x[idx]=p.x[mm]*cos(θ)-p.y[mm]*sin(θ)+rand()*xsize
    smd.y[idx]=p.x[mm]*sin(θ)-p.y[mm]*cos(θ)+rand()*ysize
end

return smd
end

