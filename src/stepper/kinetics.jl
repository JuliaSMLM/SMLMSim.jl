# Per-emitter kinetics: creation, clocks, bleaching, departure, removal, motion and walls

# Poisson(mean) count as the number of partial sums of Exp(1) draws that stay <= mean
function _poisson_count(rng::AbstractRNG, mean::Float64)
    n = 0
    s = randexp(rng)
    while s <= mean
        n += 1
        s += randexp(rng)
    end
    return n
end

_capacity(n0::Integer) = n0 + ceil(Int, 10 * sqrt(n0)) + 64

function _grow!(ps::PopState, cap::Int)
    for v in (ps.id, ps.x, ps.y, ps.z, ps.D, ps.γ, ps.m, ps.state, ps.clock, ps.budget,
              ps.t_depart, ps.t_birth)
        resize!(v, cap)
    end
    return ps
end

function _draw_state(rng::AbstractRNG, π0::Vector{Float64})
    length(π0) == 1 && return UInt8(1)
    u = rand(rng)
    acc = 0.0
    for s in 1:length(π0)-1
        acc += π0[s]
        u < acc && return UInt8(s)
    end
    return UInt8(length(π0))
end

function _draw_D(rng::AbstractRNG, mobility::Vector{Tuple{Float64,Float64}})
    length(mobility) == 1 && return mobility[1][2]
    u = rand(rng)
    acc = 0.0
    for k in 1:length(mobility)-1
        acc += mobility[k][1]
        u < acc && return mobility[k][2]
    end
    return mobility[end][2]
end

# Append one emitter born at t_birth with a fresh state; returns its index.
function _add_emitter!(w::SimWorld, ps::PopState, t_birth::Float64)
    p = ps.p
    rng = w.rng
    ps.n == length(ps.x) && _grow!(ps, 2 * ps.n)
    i = ps.n += 1
    w.id_counter += 1
    ps.id[i] = w.id_counter
    xmin, xmax, ymin, ymax = w.box
    ps.x[i] = xmin + (xmax - xmin) * rand(rng)
    ps.y[i] = ymin + (ymax - ymin) * rand(rng)
    ps.D[i] = _draw_D(rng, p.mobility)
    ps.γ[i] = p.brightness_sigma > 0 ? ps.γ0 * exp(p.brightness_sigma * randn(rng)) : ps.γ0
    ps.z[i] = p.z[1] == p.z[2] ? p.z[1] : p.z[1] + (p.z[2] - p.z[1]) * rand(rng)
    ps.m[i] = p.multiplicity
    ps.state[i] = _draw_state(rng, ps.π0)
    ps.clock[i] = randexp(rng)
    ps.budget[i] = isfinite(p.budget) ? p.budget * randexp(rng) : Inf
    ps.t_depart[i] = isfinite(p.lifetime) ? t_birth + p.lifetime * randexp(rng) : Inf
    ps.t_birth[i] = t_birth
    return i
end

# Swap-remove emitter i (the last emitter takes its slot)
function _remove!(ps::PopState, i::Int)
    j = ps.n
    if i != j
        ps.id[i] = ps.id[j]; ps.x[i] = ps.x[j]; ps.y[i] = ps.y[j]; ps.z[i] = ps.z[j]
        ps.D[i] = ps.D[j]; ps.γ[i] = ps.γ[j]; ps.m[i] = ps.m[j]; ps.state[i] = ps.state[j]
        ps.clock[i] = ps.clock[j]; ps.budget[i] = ps.budget[j]; ps.t_depart[i] = ps.t_depart[j]
        ps.t_birth[i] = ps.t_birth[j]
    end
    ps.n = j - 1
    return nothing
end

@noinline _negative_excitation(I) =
    throw(DomainError(I, "excitation must return a relative intensity >= 0"))

@noinline _bad_switch(ts, t) =
    throw(ArgumentError("next_switch returned $ts, which is not after t = $t"))

@inline function _next_switch(excitation, t::Float64)
    ts = Float64(next_switch(excitation, t))
    ts > t || _bad_switch(ts, t)
    return ts
end

@inline function _excite(excitation, x::Float64, y::Float64, z::Float64, t::Float64)
    I = Float64(excitation(x, y, z, t))
    I < 0 && _negative_excitation(I)
    return I
end

# Destination of a CTMC exit from state s (a rand is drawn only if more than one rate is positive)
function _exit_to(rng::AbstractRNG, ps::PopState, s::Int)
    q = ps.q
    n = size(q, 1)
    npos = 0
    last = 0
    for j in 1:n
        if j != s && q[s, j] > 0
            npos += 1
            last = j
        end
    end
    npos == 1 && return UInt8(last)
    u = rand(rng) * ps.exitrate[s]
    acc = 0.0
    for j in 1:n
        (j != s && q[s, j] > 0) || continue
        acc += q[s, j]
        u < acc && return UInt8(j)
        last = j
    end
    return UInt8(last)
end

# Advance emitter i over [t0 + τ0, t0 + h): the event loop of one sub-step. Returns the photons
# emitted and whether the emitter is still present.
function _advance!(w::SimWorld, ps::PopState, i::Int, t0::Float64, h::Float64, τ0::Float64, excitation)
    rng = w.rng
    x, y, z = ps.x[i], ps.y[i], ps.z[i]
    γi = ps.γ[i]
    bmean = ps.p.budget
    tdep = ps.t_depart[i]
    s = Int(ps.state[i])
    m = ps.m[i]
    clock = ps.clock[i]
    budget = ps.budget[i]
    τ = τ0
    e = 0.0
    I = _excite(excitation, x, y, z, t0 + τ0)
    ts = _next_switch(excitation, t0 + τ0)
    alive = true
    while true
        tnow = t0 + τ
        λx = s == 1 ? ps.exitrate[1] * I : ps.exitrate[s]
        tx = λx > 0 ? clock / λx : Inf
        ρe = s == 1 ? m * γi * I : 0.0
        tbl = ρe > 0 ? budget / ρe : Inf
        tdp = tdep - tnow
        trem = h - τ
        tsw = ts - tnow
        Δ0 = min(tx, tbl, tdp, trem, tsw)
        Δ = max(Δ0, 0.0)
        if Δ0 == tdp
            e += ρe * Δ
            alive = false
            break
        elseif Δ0 == tbl
            e += ρe * Δ
            τ += Δ
            clock -= λx * Δ
            m -= Int32(1)
            if m == 0
                alive = false
                break
            end
            budget = bmean * randexp(rng)
        elseif Δ0 == tx
            e += ρe * Δ
            τ += Δ
            budget -= ρe * Δ
            s = Int(_exit_to(rng, ps, s))
            clock = randexp(rng)
        elseif Δ0 == tsw && Δ0 < trem
            e += ρe * Δ
            budget -= ρe * Δ
            clock -= λx * Δ
            τ = ts - t0
            I = _excite(excitation, x, y, z, ts)
            ts = _next_switch(excitation, ts)
        else
            e += ρe * Δ
            budget -= ρe * Δ
            clock -= λx * Δ
            break
        end
    end
    if alive
        ps.state[i] = UInt8(s)
        ps.m[i] = m
        ps.clock[i] = clock
        ps.budget[i] = budget
    end
    return e, alive
end

# Fold a coordinate into [lo, hi] by reflection (any step length) or wrap it periodically
@inline function _reflect(x::Float64, lo::Float64, hi::Float64)
    (x >= lo && x <= hi) && return x
    L = hi - lo
    d = mod(x - lo, 2L)
    d > L && (d = 2L - d)
    return lo + d
end

@inline function _wrap(x::Float64, lo::Float64, hi::Float64)
    (x >= lo && x <= hi) && return x
    return lo + mod(x - lo, hi - lo)
end

@inline _fold(x::Float64, lo::Float64, hi::Float64, periodic::Bool) =
    periodic ? _wrap(x, lo, hi) : _reflect(x, lo, hi)

function _move!(w::SimWorld, ps::PopState, h::Float64)
    rng = w.rng
    xmin, xmax, ymin, ymax = w.box
    periodic = w.boundary === :periodic
    @inbounds for i in 1:ps.n
        D = ps.D[i]
        D > 0 || continue
        sd = sqrt(2 * D * h)
        ps.x[i] = _fold(ps.x[i] + sd * randn(rng), xmin, xmax, periodic)
        ps.y[i] = _fold(ps.y[i] + sd * randn(rng), ymin, ymax, periodic)
    end
    return nothing
end

# Add e photons of emitter i to its population's layer
@inline function _render!(w::SimWorld, ps::PopState, i::Int, e::Float64)
    img = ps.p.layer === :signal ? w.signal : w.oof
    u = (ps.x[i] - w.x0) / w.px
    v = (ps.y[i] - w.y0) / w.px
    st = ps.stamp
    if st === nothing
        render_gaussian!(img, w.buf, u, v, ps.sigma_px, e, ps.rwin)
    else
        render_stamp!(img, st, u, v, ps.z[i], e)
    end
    return nothing
end

# One sub-step [t0, t0 + h) of every population: kinetics and rendering, births, removal, motion
function _substep!(w::SimWorld, t0::Float64, h::Float64, excitation, record::Bool)
    for ps in w.pops
        i = 1
        while i <= ps.n
            e, alive = _advance!(w, ps, i, t0, h, 0.0, excitation)
            record && e > 0 && _render!(w, ps, i, e)
            if alive
                i += 1
            else
                _remove!(ps, i)
            end
        end
        if ps.t_next_birth < t0 + h
            xmin, xmax, ymin, ymax = w.box
            gap_scale = 1 / (ps.p.birth_rate * (xmax - xmin) * (ymax - ymin))
            while ps.t_next_birth < t0 + h
                tb = ps.t_next_birth
                j = _add_emitter!(w, ps, tb)
                e, alive = _advance!(w, ps, j, t0, h, tb - t0, excitation)
                record && e > 0 && _render!(w, ps, j, e)
                alive || _remove!(ps, j)
                ps.t_next_birth += randexp(w.rng) * gap_scale
            end
        end
        _move!(w, ps, h)
    end
    return nothing
end
