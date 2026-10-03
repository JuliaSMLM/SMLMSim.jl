# Gaussian PSF renderer (exact erf pixel integrals) and out-of-focus stamp tables

"""
    RenderBuffer(rmax::Integer)

Scratch space for [`render_gaussian!`](@ref): the 1D x and y pixel weights of one
emitter, sized for windows up to `2rmax+1` pixels wide. One buffer per task; a buffer
is not safe to share between concurrently running tasks.
"""
struct RenderBuffer
    wx::Vector{Float64}
    wy::Vector{Float64}
end

RenderBuffer(rmax::Integer) = RenderBuffer(zeros(2 * rmax + 1), zeros(2 * rmax + 1))

# Pixel weights along one axis. Pixel k spans [k-1, k) in pixel units and the centre is
# at c. Each edge costs one erfc of |t|/s (tail-accurate on both sides of the centre).
function _axis_weights!(w::Vector{Float64}, c::Float64, s::Float64, rng::UnitRange{Int})
    n = length(rng)
    length(w) < n && resize!(w, n)
    isempty(rng) && return w
    a = (first(rng) - 1 - c) / s
    ea = 0.5 * erfc(abs(a))
    @inbounds for (m, k) in enumerate(rng)
        b = (k - c) / s
        eb = 0.5 * erfc(abs(b))
        w[m] = a >= 0 ? ea - eb : (b <= 0 ? eb - ea : 1.0 - eb - ea)
        a, ea = b, eb
    end
    return w
end

"""
    render_gaussian!(img, buf::RenderBuffer, u, v, σ, photons, r::Integer) -> img
    render_gaussian!(img, buf::RenderBuffer, u, v, σ, photons, cols::UnitRange{Int}, rows::UnitRange{Int}) -> img

Add the exact pixel integrals of a 2D Gaussian with `photons` photons to `img`.

`u` and `v` are the emitter position in pixel units: pixel (row `i`, column `j`) spans
`u ∈ [j-1, j)` and `v ∈ [i-1, i)`. `σ` is in pixels. With `u = (x - x0)/pixel_size`
this is the erf-based integral over each pixel, the model erf-based fitters assume.

The `r` method uses the window `(floor(u)+1-r):(floor(u)+1+r)` clipped to `img`; the
range method uses the given (clipped) `cols` and `rows`. Weights are not renormalised,
so mass outside the window is lost. Allocation-free while the window fits `buf`.
"""
function render_gaussian!(img::AbstractMatrix, buf::RenderBuffer, u::Real, v::Real,
                          σ::Real, photons::Real, cols::UnitRange{Int}, rows::UnitRange{Int})
    cols = intersect(cols, axes(img, 2))
    rows = intersect(rows, axes(img, 1))
    (isempty(cols) || isempty(rows)) && return img
    s = sqrt(2.0) * Float64(σ)
    wx = _axis_weights!(buf.wx, Float64(u), s, cols)
    wy = _axis_weights!(buf.wy, Float64(v), s, rows)
    T = eltype(img)
    p = Float64(photons)
    @inbounds for (mj, j) in enumerate(cols)
        pw = p * wx[mj]
        for (mi, i) in enumerate(rows)
            img[i, j] += T(pw * wy[mi])
        end
    end
    return img
end

function render_gaussian!(img::AbstractMatrix, buf::RenderBuffer, u::Real, v::Real,
                          σ::Real, photons::Real, r::Integer)
    cu = floor(Int, u) + 1
    cv = floor(Int, v) + 1
    return render_gaussian!(img, buf, u, v, σ, photons, (cu - r):(cu + r), (cv - r):(cv + r))
end

# z extent a PSF declares, or nothing
_psf_z_extent(::AbstractPSF) = nothing
_psf_z_extent(psf::SplinePSF) = psf.z_range === nothing ? nothing : extrema(psf.z_range)

"""
    StampTable(psf::AbstractPSF, pixel_size::Real, zs::AbstractRange{<:Real}; radius::Integer,
               oversample::Integer=4, zinterp::Symbol=:linear)

Precomputed pixel stamps of `psf` (any 3D PSF, rings included) on the z planes `zs`, for
rendering out-of-focus emitters with [`render_stamp!`](@ref). Each plane holds
`oversample^2` sub-pixel phases of `(2radius+1)^2` pixels, each normalised to unit sum.
`pixel_size` is in μm. The table is immutable and can be shared by every task.

`zinterp = :linear` blends the two nearest planes; `:nearest` uses the nearest plane.
Recommended plane spacing is 0.2 μm or less with `:linear` and 0.1 μm or less with
`:nearest` for astigmatic or other rapidly varying PSFs. Rendering picks the nearest
phase, a position error of at most `1/(2oversample)` pixel.

Throws `ArgumentError` when `zs` is not strictly increasing, when it extends beyond the z
range the PSF declares (a `SplinePSF`'s `z_range`, no tolerance), or when a stamp's mass is
not finite and positive; rendering never clamps or draws empty stamps.
"""
struct StampTable{R<:AbstractRange{Float64}}
    stamps::Array{Float64,4}    # (2radius+1, 2radius+1, oversample^2, length(zs)); phase = a + os*b + 1
    zs::R
    radius::Int
    oversample::Int
    zinterp::Symbol
    pixel_size::Float64         # μm; the pixel size the stamps were built for
end

function StampTable(psf::AbstractPSF, pixel_size::Real, zs::AbstractRange{<:Real};
                    radius::Integer, oversample::Integer=4, zinterp::Symbol=:linear)
    zinterp in (:linear, :nearest) || throw(ArgumentError("zinterp must be :linear or :nearest, got :$zinterp"))
    radius >= 0 || throw(ArgumentError("radius must be >= 0"))
    oversample >= 1 || throw(ArgumentError("oversample must be >= 1"))
    isempty(zs) && throw(ArgumentError("zs must not be empty"))
    (length(zs) == 1 || step(zs) > 0) ||
        throw(ArgumentError("zs must be strictly increasing, got step $(step(zs))"))
    ext = _psf_z_extent(psf)
    if ext !== nothing
        zlo, zhi = extrema(zs)
        (zlo < ext[1] || zhi > ext[2]) &&
            throw(ArgumentError("z plane $(zlo < ext[1] ? zlo : zhi) of zs = $(zlo)..$(zhi) is outside the PSF's z range $(ext[1])..$(ext[2])"))
    end
    zr = convert(AbstractRange{Float64}, zs)
    r, os = Int(radius), Int(oversample)
    n = 2r + 1
    nfine = (2r + 2) * os - 1            # emitter sits at the centre of fine cell e0
    e0 = r * os + os - 1
    h = Float64(pixel_size) / os
    edges = ((0:nfine) .- (e0 + 0.5)) .* h
    stamps = zeros(n, n, os^2, length(zr))
    fine = zeros(nfine, nfine)
    for (iz, z) in enumerate(zr)
        fine .= 0.0
        integrate_pixels!(fine, psf, edges, edges, Emitter3D(0.0, 0.0, Float64(z), 1.0); threaded=false)
        for b in 0:os-1, a in 0:os-1
            sx, sy = os - 1 - a, os - 1 - b
            st = view(stamps, :, :, a + os * b + 1, iz)
            for j in 1:n, i in 1:n
                acc = 0.0
                for q in 1:os, p in 1:os
                    acc += fine[sy + (i - 1) * os + p, sx + (j - 1) * os + q]
                end
                st[i, j] = acc
            end
            mass = sum(st)
            (isfinite(mass) && mass > 0) ||
                throw(ArgumentError("stamp at z = $z has mass $mass; it must be finite and positive"))
            st ./= mass
        end
    end
    return StampTable{typeof(zr)}(stamps, zr, r, os, zinterp, Float64(pixel_size))
end

"""
    render_stamp!(img, t::StampTable, u, v, z, photons) -> img

Add `photons` times the stamp of `t` nearest to sub-pixel position `(u, v)` (pixel units,
as in [`render_gaussian!`](@ref)) and height `z` to `img`, clipped to the image.
Throws `ArgumentError` if `z` is outside `extrema(t.zs)`. Allocation-free.
"""
function render_stamp!(img::AbstractMatrix, t::StampTable, u::Real, v::Real, z::Real, photons::Real)
    zs = t.zs
    zlo, zhi = first(zs), last(zs)
    (z < zlo || z > zhi) && _stamp_z_error(z, zlo, zhi)
    nz = length(zs)
    os, r = t.oversample, t.radius
    fu, fv = floor(u), floor(v)
    cu, cv = Int(fu) + 1, Int(fv) + 1
    a = min(floor(Int, (u - fu) * os), os - 1)
    b = min(floor(Int, (v - fv) * os), os - 1)
    ph = a + os * b + 1
    if nz == 1
        _add_stamp!(img, t, ph, 1, cu, cv, Float64(photons))
    elseif t.zinterp === :nearest
        _add_stamp!(img, t, ph, clamp(round(Int, (z - zlo) / step(zs)) + 1, 1, nz), cu, cv, Float64(photons))
    else
        x = (z - zlo) / step(zs)
        i0 = clamp(floor(Int, x), 0, nz - 2)
        w = x - i0
        _add_stamp!(img, t, ph, i0 + 1, cu, cv, (1 - w) * photons)
        _add_stamp!(img, t, ph, i0 + 2, cu, cv, w * photons)
    end
    return img
end

@noinline _stamp_z_error(z, zlo, zhi) =
    throw(ArgumentError("emitter z = $z is outside the StampTable's z range $zlo..$zhi"))

function _add_stamp!(img::AbstractMatrix, t::StampTable, ph::Int, iz::Int, cu::Int, cv::Int, p::Float64)
    r = t.radius
    T = eltype(img)
    st = t.stamps
    @inbounds for j in max(1, cu - r):min(size(img, 2), cu + r)
        sj = j - cu + r + 1
        for i in max(1, cv - r):min(size(img, 1), cv + r)
            img[i, j] += T(p * st[i - cv + r + 1, sj, ph, iz])
        end
    end
    return nothing
end
