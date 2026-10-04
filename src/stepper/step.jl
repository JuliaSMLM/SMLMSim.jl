# SimWorld construction, step! and layers

# Kernel half-width in μm of one population (5σ for a Gaussian, the stamp radius for a stamp table)
_half_width(p::Population) =
    p.psf isa StampTable ? p.psf.radius * p.psf.pixel_size : 5 * p.psf.σ

default_margin(pops::AbstractVector{Population}) =
    isempty(pops) ? 0.0 : maximum(_half_width, pops)

# Pixel geometry of a camera: (ny, nx, px, x0, y0); pixels must be uniform and square
function _pixel_geometry(camera)
    ex, ey = camera.pixel_edges_x, camera.pixel_edges_y
    (length(ex) >= 2 && length(ey) >= 2) || throw(ArgumentError("camera has no pixels"))
    px = Float64(ex[2] - ex[1])
    tol = 1e-6 * px
    (isfinite(_uniform_pitch(ex)) && isfinite(_uniform_pitch(ey)) && abs(Float64(ey[2] - ey[1]) - px) <= tol) ||
        throw(ArgumentError("SimWorld needs uniform square pixels"))
    return length(ey) - 1, length(ex) - 1, px, Float64(ex[1]), Float64(ey[1])
end

# Stamp table with a concrete range type, so rendering does not dispatch dynamically
function _concrete_stamp(t::StampTable)
    zs = t.zs
    zr = range(Float64(first(zs)), Float64(last(zs)); length=length(zs))
    return StampTable{typeof(zr)}(t.stamps, zr, t.radius, t.oversample, t.zinterp, t.pixel_size)
end

function _pop_state(rng::AbstractRNG, p::Population, px::Float64, box::NTuple{4,Float64}, t0::Float64)
    q = Matrix{Float64}(p.fluor.q)
    ns = size(q, 1)
    exitrate = [sum(q[s, j] for j in 1:ns if j != s; init=0.0) for s in 1:ns]
    isg = p.psf isa GaussianPSF
    σpx = isg ? p.psf.σ / px : 0.0
    A = (box[2] - box[1]) * (box[4] - box[3])
    (isfinite(p.density * A) && isfinite(p.birth_rate * A)) ||
        throw(ArgumentError("population :$(p.name): density or birth_rate times the world-box area $A μm² is not finite"))
    n0 = _poisson_count(rng, p.density * A)
    cap = _capacity(n0)
    ps = PopState(p, Float64(p.fluor.γ), q, exitrate, _stationary(q), σpx,
                  isg ? ceil(Int, 5 * σpx) + 1 : 0,
                  isg ? nothing : _concrete_stamp(p.psf),
                  0, zeros(Int, cap), zeros(cap), zeros(cap), zeros(cap), zeros(cap), zeros(cap),
                  zeros(Int32, cap), zeros(UInt8, cap), zeros(cap), zeros(cap), zeros(cap), zeros(cap),
                  zeros(cap), zeros(cap), zeros(cap), zeros(cap), zeros(cap), zeros(cap), zeros(cap),
                  zeros(cap), zeros(cap), zeros(cap), zeros(cap), zeros(cap),
                  zeros(Int32, cap), zeros(Int32, cap), zeros(Int, cap), zeros(cap), zeros(cap),
                  zeros(cap), zeros(cap), zeros(cap), zeros(Bool, cap), zeros(Int, cap), Inf)
    return ps, n0, A
end

function SimWorld(rng::AbstractRNG, camera::Union{IdealCamera,SCMOSCamera}, pops::AbstractVector{Population};
                  background::Union{Nothing,BackgroundModel}=nothing, n_sub::Integer,
                  boundary::Symbol=:reflecting,
                  margin::Real=default_margin(pops), t0::Real=0.0, merge_radius::Real=0.0,
                  dimers::Union{Nothing,DimerKinetics}=nothing)
    boundary in (:reflecting, :periodic) || throw(ArgumentError("boundary must be :reflecting or :periodic"))
    n_sub >= 1 || throw(ArgumentError("n_sub must be >= 1"))
    margin, t0, merge_radius = Float64(margin), Float64(t0), Float64(merge_radius)
    (margin >= 0 && isfinite(margin)) || throw(ArgumentError("margin must be finite and >= 0"))
    isfinite(t0) || throw(ArgumentError("t0 must be finite"))
    merge_radius >= 0 || throw(ArgumentError("merge_radius must be >= 0"))
    ny, nx, px, x0, y0 = _pixel_geometry(camera)
    for p in pops
        p.psf isa StampTable && !isapprox(p.psf.pixel_size, px; rtol=1e-6) &&
            throw(ArgumentError("population :$(p.name) has a StampTable built at pixel size $(p.psf.pixel_size) μm, the camera's is $px μm"))
    end
    if dimers !== nothing
        for p in pops
            p.binds && p.multiplicity > 1 &&
                throw(ArgumentError("population :$(p.name) binds with multiplicity $(p.multiplicity); dimers need multiplicity <= 1 (set binds = false)"))
            p.binds && p.birth_rate > 0 && !isfinite(p.lifetime) &&
                throw(ArgumentError("population :$(p.name) binds and has births with an infinite lifetime; with dimers its bleached molecules stay until they depart, so births need a finite lifetime"))
        end
    end
    box = (x0 - margin, x0 + nx * px + margin, y0 - margin, y0 + ny * px + margin)
    if dimers !== nothing
        need = 2 * max(_split_sep(box, dimers.r_react), dimers.d_dimer * (1 + 1e-9))
        (box[2] - box[1] > need && box[4] - box[3] > need) ||
            throw(ArgumentError("every side of the world box must exceed 2 max(s, d_dimer) = $need μm, s the split separation (r_react plus a rounding margin)"))
    end
    ncx = dimers === nothing ? 0 : _cells_per_axis(box[2] - box[1], dimers.r_react)
    ncy = dimers === nothing ? 0 : _cells_per_axis(box[4] - box[3], dimers.r_react)
    w = SimWorld(rng, camera, px, x0, y0, box, boundary, Int(n_sub), t0, 0, PopState[], nothing,
                 zeros(ny, nx), zeros(ny, nx), zeros(ny, nx), zeros(ny, nx), 0,
                 RenderBuffer(maximum((p.psf isa StampTable ? p.psf.radius : ceil(Int, 5 * p.psf.σ / px) + 1
                                       for p in pops); init=0)),
                 0, FrameTruth[], 0, merge_radius, Bool[], Int[], t0, t0,
                 dimers, zeros(Int32, ncx * ncy), Int32[], Int32[], Int32[], ncx, ncy)
    for p in pops
        ps, n0, A = _pop_state(rng, p, px, box, t0)
        push!(w.pops, ps)
        for _ in 1:n0
            _add_emitter!(w, ps, t0)
        end
        ps.t_next_birth = p.birth_rate > 0 ? t0 + randexp(rng) / (p.birth_rate * A) : Inf
    end
    background === nothing || (w.bg = BackgroundState(rng, background, ny, nx, px, t0))
    resize!(w.truth, max(1, 2 * sum(ps -> length(ps.x), w.pops; init=0)))
    return w
end

"""
    SMLMSim.layers(world) -> (signal, oof, structured, expected)

Public but not exported (the name is generic): call it as `SMLMSim.layers`. The noise-free photon maps of the last [`step!`](@ref), without allocating: the `:signal`
populations, the `:oof` populations, the structured background and their sum.
"""
layers(w::SimWorld) = (signal=w.signal, oof=w.oof, structured=w.structured, expected=w.expected)

"""
    SMLMSim.step!(world, t_a, t_b, excitation = UniformExcitation()) -> Matrix{Float64}

Public but not exported (SciML and Agents.jl also export a `step!`): call it as `SMLMSim.step!`.
Advance `world` over the exposure `[t_a, t_b)` in `world.n_sub` sub-steps and return
`world.expected`, the noise-free photons per pixel (the same buffer every call). The world's
RNG drives kinetics only; add camera noise with a separate call on your own RNG.

`excitation(x, y, z, t)::Float64` gives the relative intensity (at least 0; 1 is the
baseline) at an emitter's position in μm, its own height `z` in μm from focus, and time in
s. It scales emission, photon-budget draw-down and exits from the emitting state. It must be
type-stable and non-allocating for `step!` to be allocation-free. [`next_switch`](@ref)
declares discontinuities so they act at their exact time.

Time never runs backward. A call with `t_a` within `1e-9 max(1, |world.t|)` of `world.t` is
contiguous (`world.t` is set to `t_a`). A later `t_a` advances the gap unrecorded in
`ceil(gap/h)` equal sub-steps, `h = (t_b - t_a)/n_sub`. A nonfinite time, `t_b <= t_a` or an earlier `t_a`
throws `ArgumentError`.

Every sub-step is a half-open window `[t0, t1)` and the last one ends exactly at `t_b`. An event exactly at `t_b`
(a departure, bleach, blink, switch, birth, formation or break) belongs to the next exposure.
"""
function step!(w::SimWorld, t_a::Real, t_b::Real, excitation::E=UniformExcitation()) where {E}
    t_a, t_b = Float64(t_a), Float64(t_b)
    (isfinite(t_a) && isfinite(t_b)) || throw(ArgumentError("t_a and t_b must be finite"))
    t_b > t_a || throw(ArgumentError("t_b must be greater than t_a"))
    ε = 1e-9 * max(1.0, abs(w.t))
    t_a < w.t - ε && throw(ArgumentError("t_a = $t_a is before the world time $(w.t); time does not run backward"))
    h = (t_b - t_a) / w.n_sub
    if abs(t_a - w.t) <= ε
        w.t = t_a
    else
        g = t_a - w.t
        n_g = ceil(Int, g / h)
        hg = g / n_g
        tg = w.t
        for k in 1:n_g
            _substep!(w, tg + (k - 1) * hg, k == n_g ? t_a : tg + k * hg, excitation, false)
            w.n_gap += 1
        end
        w.t = t_a
    end
    fill!(w.signal, 0.0)
    fill!(w.oof, 0.0)
    w.bg === nothing || _update_background!(w, w.bg, t_a, t_b)
    w.t_a, w.t_b = t_a, t_b
    _begin_truth!(w)
    for k in 1:w.n_sub
        _substep!(w, t_a + (k - 1) * h, k == w.n_sub ? t_b : t_a + k * h, excitation, true)
    end
    _finish_truth!(w)
    w.t = t_b
    w.frame += 1
    @. w.expected = w.signal + w.oof + w.structured
    return w.expected
end
