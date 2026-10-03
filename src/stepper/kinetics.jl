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

# Every per-emitter vector of a PopState: the one list that growth and swap-removal walk. Kept as two
# tuples because inference gives up on tuples longer than 32 elements and the loop would allocate.
_vectors(ps::PopState) = ((ps.id, ps.x, ps.y, ps.z, ps.D, ps.γ, ps.m, ps.state, ps.clock, ps.budget,
                          ps.t_depart, ps.t_birth, ps.lj, ps.xr, ps.yr, ps.sx, ps.sy, ps.sxp, ps.syp,
                          ps.sph, ps.t_present, ps.t_lit, ps.sI, ps.t_bleach_f, ps.partner, ps.partner_pop),
                          (ps.partner_id, ps.θ, ps.t_form, ps.t_break_due, ps.t_bound, ps.t_break_f, ps.t_litb, ps.vis_form,
                           ps.partner_f))

@inline function _each_vector(f, ps::PopState)
    a, b = _vectors(ps)
    foreach(f, a)
    foreach(f, b)
    return nothing
end

function _grow!(ps::PopState, cap::Int)
    _each_vector(v -> resize!(v, cap), ps)
    return ps
end

# Start a frame's accumulators for emitter i at its current position
@inline function _reset_acc!(ps::PopState, i::Int)
    ps.xr[i] = ps.x[i]; ps.yr[i] = ps.y[i]
    ps.sx[i] = 0.0; ps.sy[i] = 0.0; ps.sxp[i] = 0.0; ps.syp[i] = 0.0
    ps.sph[i] = 0.0; ps.t_present[i] = 0.0; ps.t_lit[i] = 0.0; ps.sI[i] = 0.0
    ps.t_bleach_f[i] = NaN
    ps.t_bound[i] = 0.0; ps.t_break_f[i] = NaN; ps.t_litb[i] = 0.0
    return nothing
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
    ps.lj[i] = p.brightness_jitter > 0 ? p.brightness_jitter * randn(rng) : 0.0
    ps.partner[i] = 0; ps.partner_pop[i] = 0; ps.partner_id[i] = 0
    ps.θ[i] = 0.0; ps.t_form[i] = -Inf; ps.t_break_due[i] = -Inf; ps.vis_form[i] = false; ps.partner_f[i] = 0
    _reset_acc!(ps, i)
    return i
end

# Swap-remove emitter i (the last emitter takes its slot). A bound partner is unlinked (it keeps
# partner_pop and partner_id for its truth row); a bound emitter moved into slot i re-points its partner.
function _remove!(w::SimWorld, ps::PopState, i::Int)
    pi = ps.partner[i]
    pi != 0 && (w.pops[ps.partner_pop[i]].partner[pi] = 0)
    j = ps.n
    i != j && _each_vector(v -> (v[i] = v[j]), ps)
    ps.n = j - 1
    if i != j && ps.partner[i] != 0
        w.pops[ps.partner_pop[i]].partner[ps.partner[i]] = i
    end
    return nothing
end

@noinline _negative_excitation(I) =
    throw(DomainError(I, "excitation must return a relative intensity >= 0"))

@noinline _bad_switch(ts, t) =
    throw(ArgumentError("next_switch returned $ts, which is not after t = $t"))

@inline function _next_switch(excitation::E, t::Float64) where {E}
    ts = Float64(next_switch(excitation, t))
    ts > t || _bad_switch(ts, t)
    return ts
end

@inline function _excite(excitation::E, x::Float64, y::Float64, z::Float64, t::Float64) where {E}
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
# emitted, whether the emitter is still present and whether it left by departure.
function _advance!(w::SimWorld, ps::PopState, i::Int, t0::Float64, h::Float64, τ0::Float64, excitation::E) where {E}
    rng = w.rng
    x, y, z = ps.x[i], ps.y[i], ps.z[i]
    γi = ps.p.brightness_jitter > 0 ? ps.γ[i] * exp(ps.lj[i]) : ps.γ[i]
    bmean = ps.p.budget
    tdep = ps.t_depart[i]
    s = Int(ps.state[i])
    m = ps.m[i]
    clock = ps.clock[i]
    budget = ps.budget[i]
    τ = τ0
    e = 0.0
    tp = ps.t_present[i]
    tl = ps.t_lit[i]
    sI = ps.sI[i]
    tbf = ps.t_bleach_f[i]
    tb = ps.t_bound[i]
    tlb = ps.t_litb[i]
    tform = ps.t_form[i]
    tbrk = ps.t_break_due[i]
    dimers = w.dimers !== nothing
    phase = t0 + τ0 < tform ? 0 : (t0 + τ0 < tbrk ? 1 : 2)   # before, inside, after the bound interval
    I = _excite(excitation, x, y, z, t0 + τ0)
    ts = _next_switch(excitation, t0 + τ0)
    alive = true
    departed = false
    while true
        lit = s == 1 && m > 0
        tnow = t0 + τ
        λx = m == 0 ? 0.0 : (s == 1 ? ps.exitrate[1] * I : ps.exitrate[s])
        tx = λx > 0 ? clock / λx : Inf
        ρe = s == 1 ? m * γi * I : 0.0
        tbl = ρe > 0 ? budget / ρe : Inf
        tdp = tdep - tnow
        trem = h - τ
        tsw = ts - tnow
        tbk = phase == 0 ? tform - tnow : (phase == 1 ? tbrk - tnow : Inf)
        bnd = phase == 1
        Δ0 = min(tx, tbl, tdp, trem, tsw, tbk)
        Δ = max(Δ0, 0.0)
        if Δ0 == tdp
            e += ρe * Δ
            tp += Δ; sI += I * Δ; lit && (tl += Δ); bnd && (tb += Δ; lit && (tlb += Δ))
            alive = false
            departed = true
            break
        elseif Δ0 == tbl
            e += ρe * Δ
            tp += Δ; sI += I * Δ; lit && (tl += Δ); bnd && (tb += Δ; lit && (tlb += Δ))
            τ += Δ
            clock -= λx * Δ
            m -= Int32(1)
            if m == 0
                tbf = t0 + τ
                dimers || (alive = false; break)   # with dimers a bleached emitter stays, dark
            else
                budget = bmean * randexp(rng)
            end
        elseif Δ0 == tx
            e += ρe * Δ
            tp += Δ; sI += I * Δ; lit && (tl += Δ); bnd && (tb += Δ; lit && (tlb += Δ))
            τ += Δ
            budget -= ρe * Δ
            s = Int(_exit_to(rng, ps, s))
            clock = randexp(rng)
        elseif Δ0 == tsw && Δ0 < trem
            e += ρe * Δ
            tp += Δ; sI += I * Δ; lit && (tl += Δ); bnd && (tb += Δ; lit && (tlb += Δ))
            budget -= ρe * Δ
            clock -= λx * Δ
            τ = ts - t0
            I = _excite(excitation, x, y, z, ts)
            ts = _next_switch(excitation, ts)
        elseif Δ0 == tbk && Δ0 < trem
            e += ρe * Δ
            tp += Δ; sI += I * Δ; lit && (tl += Δ); bnd && (tb += Δ; lit && (tlb += Δ))
            budget -= ρe * Δ
            clock -= λx * Δ
            τ += Δ
            phase += 1
        else
            e += ρe * Δ
            tp += Δ; sI += I * Δ; lit && (tl += Δ); bnd && (tb += Δ; lit && (tlb += Δ))
            budget -= ρe * Δ
            clock -= λx * Δ
            break
        end
    end
    ps.m[i] = m
    ps.t_present[i] = tp
    ps.t_lit[i] = tl
    ps.sI[i] = sI
    ps.t_bleach_f[i] = tbf
    ps.t_bound[i] = tb
    ps.t_litb[i] = tlb
    if alive
        ps.state[i] = UInt8(s)
        ps.clock[i] = clock
        ps.budget[i] = budget
    end
    return e, alive, departed
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

function _move!(w::SimWorld, k::Int, ps::PopState, h::Float64)
    rng = w.rng
    xmin, xmax, ymin, ymax = w.box
    periodic = w.boundary === :periodic
    @inbounds for i in 1:ps.n
        if ps.partner[i] != 0
            _move_pair!(w, k, ps, i, h)
            continue
        end
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

# Advance emitter i and fold its photons and presence into the frame accumulators; returns
# (e, alive, departed)
@inline function _advance_acc!(w::SimWorld, ps::PopState, i::Int, t0::Float64, h::Float64, τ0::Float64,
                               excitation::E) where {E}
    tp0 = ps.t_present[i]
    e, alive, departed = _advance!(w, ps, i, t0, h, τ0, excitation)
    Δp = ps.t_present[i] - tp0
    dx = ps.x[i] - ps.xr[i]
    dy = ps.y[i] - ps.yr[i]
    if w.boundary === :periodic
        Lx = w.box[2] - w.box[1]
        Ly = w.box[4] - w.box[3]
        dx -= Lx * round(dx / Lx)
        dy -= Ly * round(dy / Ly)
    end
    ps.sx[i] += dx * e; ps.sy[i] += dy * e; ps.sph[i] += e
    ps.sxp[i] += dx * Δp; ps.syp[i] += dy * Δp
    return e, alive, departed
end

# One sub-step [t0, t0 + h) of every population: kinetics and rendering, births, removal, motion
function _substep!(w::SimWorld, t0::Float64, h::Float64, excitation::E, record::Bool) where {E}
    w.dimers === nothing || _dimer_events!(w, t0, h)
    for (k, ps) in enumerate(w.pops)
        i = 1
        while i <= ps.n
            e, alive, departed = _advance_acc!(w, ps, i, t0, h, 0.0, excitation)
            record && e > 0 && _render!(w, ps, i, e)
            if alive
                i += 1
            else
                record && _write_row!(w, k, ps, i, departed)
                _remove!(w, ps, i)
            end
        end
        if ps.t_next_birth < t0 + h
            xmin, xmax, ymin, ymax = w.box
            gap_scale = 1 / (ps.p.birth_rate * (xmax - xmin) * (ymax - ymin))
            while ps.t_next_birth < t0 + h
                tb = ps.t_next_birth
                j = _add_emitter!(w, ps, tb)
                e, alive, departed = _advance_acc!(w, ps, j, t0, h, tb - t0, excitation)
                record && e > 0 && _render!(w, ps, j, e)
                if !alive
                    record && _write_row!(w, k, ps, j, departed)
                    _remove!(w, ps, j)
                end
                ps.t_next_birth += randexp(w.rng) * gap_scale
            end
        end
        _move!(w, k, ps, h)
        ps.p.brightness_jitter > 0 && _jitter!(w, ps, h)
    end
    return nothing
end

# Advance every emitter's log-brightness multiplier by h: the exact OU (AR(1)) step
function _jitter!(w::SimWorld, ps::PopState, h::Float64)
    s = ps.p.brightness_jitter
    τ = ps.p.jitter_time
    a, b = _ou_coeffs(s, τ, h)
    rng = w.rng
    @inbounds for i in 1:ps.n
        ps.lj[i] = a * ps.lj[i] + b * randn(rng)
    end
    return nothing
end
