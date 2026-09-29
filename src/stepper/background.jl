# Structured background: level per stretch, per-exposure jitter, lattice-OU pattern and
# illumination; gen_background builds a background movie from a SimWorld

# Cubic B-spline weights of the four nodes i0-1 .. i0+2 at fractional position f in [0, 1)
@inline function _bspline_weights(f::Float64)
    f2 = f * f
    f3 = f2 * f
    return ((1 - f)^3 / 6, (3 * f3 - 6 * f2 + 4) / 6, (-3 * f3 + 3 * f2 + 3 * f + 1) / 6, f3 / 6)
end

# Per-axis lattice: first node index and normalised weights of every pixel, and the node count
function _lattice_axis(n::Int, px::Float64, h::Float64, off::Float64)
    first = zeros(Int, n)
    wts = zeros(4, n)
    for k in 1:n
        s = 2 + off + (k - 0.5) * px / h
        i0 = floor(Int, s)
        w = _bspline_weights(s - i0)
        rn = sqrt(sum(abs2, w))
        first[k] = i0 - 1
        for a in 1:4
            wts[a, k] = w[a] / rn
        end
    end
    return first, wts, maximum(first) + 3
end

# Gaussian illumination profile along one axis, normalised to mean 1 over the pixels
function _illumination_axis(n::Int, px::Float64, centre::Float64, width::Float64)
    v = [exp(-((k - 0.5) * px - centre)^2 / (2 * width^2)) for k in 1:n]
    return v ./ (sum(v) / n)
end

_draw_level!(out::Vector{Float64}, rng::AbstractRNG, level::Real) = (out[1] = Float64(level); nothing)
_draw_level!(out::Vector{Float64}, rng::AbstractRNG, level) = (out[1] = Float64(rand(rng, level)); nothing)

function BackgroundState(rng::AbstractRNG, m::BackgroundModel, ny::Int, nx::Int, px::Float64, t0::Float64)
    (m.level isa Real ? m.level >= 0 : true) || throw(ArgumentError("level must be >= 0"))
    m.stretch > 0 || throw(ArgumentError("stretch must be > 0"))
    m.jitter >= 0 || throw(ArgumentError("jitter must be >= 0"))
    m.feature_size > 0 || throw(ArgumentError("feature_size must be > 0"))
    m.contrast >= 0 || throw(ArgumentError("contrast must be >= 0"))
    m.correlation_time > 0 || throw(ArgumentError("correlation_time must be > 0"))
    m.illumination_width > 0 || throw(ArgumentError("illumination_width must be > 0"))
    P = ones(ny, nx)
    if isfinite(m.illumination_width)
        cx = nx * px * rand(rng)
        cy = ny * px * rand(rng)
        vx = _illumination_axis(nx, px, cx, m.illumination_width)
        vy = _illumination_axis(ny, px, cy, m.illumination_width)
        P .= vy .* vx'
    end
    if m.contrast > 0
        h = m.feature_size / sqrt(2 / 3)
        ox, oy = rand(rng), rand(rng)
        ix, wx, nnx = _lattice_axis(nx, px, h, ox)
        iy, wy, nny = _lattice_axis(ny, px, h, oy)
        g = randn(rng, nny, nnx)
    else
        ix, wx = zeros(Int, nx), zeros(4, nx)
        iy, wy = zeros(Int, ny), zeros(4, ny)
        g = zeros(0, 0)
    end
    bs = BackgroundState(m.level, m.stretch, m.jitter, m.contrast, m.correlation_time, t0, t0, 0, 0.0, 0.0, zeros(1), P, g, zeros(ny, size(g, 2)), iy, ix, wy, wx)
    _draw_level!(bs.draw, rng, bs.level_src)
    bs.level = bs.draw[1]
    return bs
end

# Advance the background to the exposure [t_a, t_b) and fill w.structured
function _update_background!(w, bs::BackgroundState, t_a::Float64, t_b::Float64)
    rng = w.rng
    if isfinite(bs.stretch)
        idx = floor(Int, (t_a - bs.t0 + 1e-9 * max(1.0, abs(w.t))) / bs.stretch)
        if idx > bs.stretch_idx
            bs.stretch_idx = idx
            _draw_level!(bs.draw, rng, bs.level_src)
            bs.level = bs.draw[1]
        end
    end
    J = bs.jitter > 0 ? max(0.0, 1 + bs.jitter * randn(rng)) : 1.0
    lvl = bs.level * J * (t_b - t_a)
    bs.frame_level = lvl
    S = w.structured
    P = bs.P
    if bs.contrast == 0
        @inbounds @. S = lvl * P
    else
        g = bs.g
        ρ = exp(-(t_a - bs.t_prev) / bs.tau)
        if ρ < 1
            s = sqrt(1 - ρ^2)
            @inbounds for k in eachindex(g)
                g[k] = ρ * g[k] + s * randn(rng)
            end
        end
        tmp, iy, ix, wy, wx = bs.tmp, bs.iy, bs.ix, bs.wy, bs.wx
        c = bs.contrast
        shift = c^2 / 2
        @inbounds for jl in 1:size(g, 2), i in eachindex(iy)
            k = iy[i]
            tmp[i, jl] = wy[1, i] * g[k, jl] + wy[2, i] * g[k+1, jl] + wy[3, i] * g[k+2, jl] + wy[4, i] * g[k+3, jl]
        end
        @inbounds for j in eachindex(ix), i in eachindex(iy)
            k = ix[j]
            pat = wx[1, j] * tmp[i, k] + wx[2, j] * tmp[i, k+1] + wx[3, j] * tmp[i, k+2] + wx[4, j] * tmp[i, k+3]
            S[i, j] = lvl * P[i, j] * exp(c * pat - shift)
        end
    end
    bs.t_prev = t_a
    return nothing
end

"""
    gen_background(rng, camera, bg, n_frames; oof = Population[], frame_time, exposure = frame_time,
                   n_sub = 1, t_burn = 0.0) -> (structured, oof, level)

Noise-free background movie: the expected photons per pixel of every frame, as two arrays of
size `(ny, nx, n_frames)` and the per-frame level.

- `bg`: a [`BackgroundModel`](@ref) (model A: level, jitter, pattern, illumination) or
  `nothing`.
- `oof`: [`Population`](@ref)s with `layer = :oof` (model B: blinking, bleaching
  out-of-focus emitters). The models combine: A alone, B alone, or A and B.
- `frame_time`: s between exposure starts (required); `exposure` (at most `frame_time`) is
  the exposure length; `n_sub` is sub-steps per exposure.
- `t_burn`: s of unrecorded kinetics before frame 1. The stretch grid starts at the world's `t0`, so when `t_burn` is not a multiple of `stretch` the first recorded stretch is short. Use at least 5 mean residence times of
  the `oof` populations, since their initial ensemble is stationary only for an infinite budget.

`structured` is the pattern map (zeros for `bg = nothing`); `oof` is the out-of-focus light
(zeros without populations); `level[k]` is `L J exposure`, the true mean of `structured` in
frame `k`. The two maps are kept apart so the caller decides the background target: `structured`,
`oof` or their sum. Camera noise is a separate call, for example
`gen_images(smld, psf; bg = structured .+ oof, camera_noise = true, rng)`.

`rng` drives the background and the out-of-focus kinetics.
"""
function gen_background(rng::AbstractRNG, camera::Union{IdealCamera,SCMOSCamera},
                        bg::Union{Nothing,BackgroundModel}, n_frames::Integer;
                        oof::AbstractVector{Population}=Population[], frame_time::Real,
                        exposure::Real=frame_time, n_sub::Integer=1, t_burn::Real=0.0)
    n_frames >= 1 || throw(ArgumentError("n_frames must be >= 1"))
    frame_time > 0 || throw(ArgumentError("frame_time must be > 0"))
    0 < exposure <= frame_time || throw(ArgumentError("exposure must satisfy 0 < exposure <= frame_time"))
    t_burn >= 0 || throw(ArgumentError("t_burn must be >= 0"))
    all(p -> p.layer === :oof, oof) || throw(ArgumentError("every oof population must have layer = :oof"))
    w = SimWorld(rng, camera, Population[oof...]; n_sub, background=bg)
    ny, nx = size(w.expected)
    structured = zeros(ny, nx, n_frames)
    oofmap = zeros(ny, nx, n_frames)
    level = zeros(n_frames)
    for k in 1:n_frames
        t_a = t_burn + (k - 1) * frame_time
        step!(w, t_a, t_a + exposure)
        structured[:, :, k] .= w.structured
        oofmap[:, :, k] .= w.oof
        level[k] = w.bg === nothing ? 0.0 : w.bg.frame_level
    end
    return structured, oofmap, level
end
