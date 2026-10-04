# Cross-population dimers: cell list, formation, dissociation with the unbinding rule, pair motion

# Pair geometry follows the diffusion path (docs/src/diffusion/rules.md): the numeric primitives _disp, _dir, _angles and
# _norm are the diffusion path's; the placement is the stepper's own, in its rectangular box with a nonzero origin

# Cells per axis: the side is at least r_react; an axis with fewer than 3 cells gets 1 (no neighbour visits). The
# quotient saturates at 128 before conversion, so a tiny r_react cannot overflow Int
function _cells_per_axis(L::Float64, r::Float64)
    n = floor(Int, min(L / r, 128.0))
    return n < 3 ? 1 : n
end

# Separation of split partners closer than r_react: r_react plus a margin that survives rounding coordinates as
# large as the box's largest |coordinate| (rules.md, Definitions), and at least 1e-9 r_react
_split_sep(box::NTuple{4,Float64}, r::Float64) = max(r * (1 + 1e-9), r + 8 * eps(Float64) * (maximum(abs, box) + r))

# In-plane displacement from (A, ia) to (B, ib): each axis in Float64, rounded once, the minimum image under
# :periodic (rules.md, Distance)
@inline function _pair_disp(w::SimWorld, A::PopState, ia::Int, B::PopState, ib::Int)
    if w.boundary === :periodic
        return (_disp(A.x[ia], B.x[ib], w.box[2] - w.box[1]), _disp(A.y[ia], B.y[ib], w.box[4] - w.box[3]))
    end
    return (_disp(A.x[ia], B.x[ib], nothing), _disp(A.y[ia], B.y[ib], nothing))
end

# The one contact distance of formation and split: the scaled norm of (dx, dy, dz), each emitter at its own z
@inline function _contact(w::SimWorld, A::PopState, ia::Int, B::PopState, ib::Int)
    d = _pair_disp(w, A, ia, B, ib)
    return _norm((d[1], d[2], B.z[ib] - A.z[ia]))
end

# (A, ia) to c - (s/2)u and (B, ib) to c + (s/2)u. :periodic wraps the centre, then each end. :reflecting keeps the
# centre (s/2)|u_k| from each wall: clamped when placing (formation, split), folded as often as it crosses when
# moving (reflect = true: the rigid-body reflection of rules.md); each end is then clamped into the box
@inline function _put_pair!(w::SimWorld, A::PopState, ia::Int, B::PopState, ib::Int,
                            cx::Float64, cy::Float64, ux::Float64, uy::Float64, s::Float64, reflect::Bool)
    xmin, xmax, ymin, ymax = w.box
    hs = 0.5 * s
    if w.boundary === :periodic
        cx = _wrap(cx, xmin, xmax); cy = _wrap(cy, ymin, ymax)
        A.x[ia] = _wrap(cx - hs * ux, xmin, xmax); A.y[ia] = _wrap(cy - hs * uy, ymin, ymax)
        B.x[ib] = _wrap(cx + hs * ux, xmin, xmax); B.y[ib] = _wrap(cy + hs * uy, ymin, ymax)
    else
        hx = hs * abs(ux); hy = hs * abs(uy)
        cx = reflect ? _reflect(cx, xmin + hx, xmax - hx) : clamp(cx, xmin + hx, xmax - hx)
        cy = reflect ? _reflect(cy, ymin + hy, ymax - hy) : clamp(cy, ymin + hy, ymax - hy)
        A.x[ia] = clamp(cx - hs * ux, xmin, xmax); A.y[ia] = clamp(cy - hs * uy, ymin, ymax)
        B.x[ib] = clamp(cx + hs * ux, xmin, xmax); B.y[ib] = clamp(cy + hs * uy, ymin, ymax)
    end
    return nothing
end

# Place the pair (A, ia) < (B, ib) at in-plane separation s along its current axis (the displacement from A to B).
# An anchored pair (D_dimer = :min with an immobile member; of two immobile members the smaller id) keeps the anchor
function _place_pair!(w::SimWorld, A::PopState, ia::Int, B::PopState, ib::Int, s::Float64)
    if w.dimers.D_dimer === :min && min(A.D[ia], B.D[ib]) == 0
        xmin, xmax, ymin, ymax = w.box
        a_anchor = A.D[ia] == 0 && (B.D[ib] != 0 || A.id[ia] < B.id[ib])
        P, ip, Q, iq = a_anchor ? (A, ia, B, ib) : (B, ib, A, ia)
        ux, uy = _dir(_pair_disp(w, P, ip, Q, iq))
        nx = P.x[ip] + s * ux
        ny = P.y[ip] + s * uy
        if w.boundary === :periodic
            nx = _wrap(nx, xmin, xmax)
            ny = _wrap(ny, ymin, ymax)
        else
            (nx < xmin || nx > xmax) && (nx = P.x[ip] - s * ux)
            (ny < ymin || ny > ymax) && (ny = P.y[ip] - s * uy)
            nx = clamp(nx, xmin, xmax)
            ny = clamp(ny, ymin, ymax)
        end
        Q.x[iq] = nx
        Q.y[iq] = ny
        return nothing
    end
    d = _pair_disp(w, A, ia, B, ib)
    ux, uy = _dir(d)
    _put_pair!(w, A, ia, B, ib, A.x[ia] + 0.5 * d[1], A.y[ia] + 0.5 * d[2], ux, uy, s, false)
    return nothing
end

# A bound pair moves once, at its member with the larger (population, index); an anchored pair does not move. The
# orientation is that of the current displacement from the partner, plus the rotational step
function _move_pair!(w::SimWorld, k::Int, ps::PopState, i::Int, h::Float64)
    kp = Int(ps.partner_pop[i])
    ip = Int(ps.partner[i])
    (k > kp || (k == kp && i > ip)) || return nothing
    dk = w.dimers
    Q = w.pops[kp]
    dk.D_dimer === :min && min(ps.D[i], Q.D[ip]) == 0 && return nothing
    Dc = dk.D_dimer === :min ? min(ps.D[i], Q.D[ip]) : dk.D_dimer
    rng = w.rng
    ξ1 = randn(rng); ξ2 = randn(rng); ξ3 = randn(rng)
    d = _pair_disp(w, Q, ip, ps, i)
    sd = sqrt(2 * Dc * h)
    φ = _angles(d) + sqrt(2 * dk.D_rot * h) * ξ3
    _put_pair!(w, Q, ip, ps, i, Q.x[ip] + 0.5 * d[1] + sd * ξ1, Q.y[ip] + 0.5 * d[2] + sd * ξ2,
               cos(φ), sin(φ), dk.d_dimer, true)
    return nothing
end

# 1a: cell list over present, unbound, binding emitters
function _build_cells!(w::SimWorld)
    total = 0
    for ps in w.pops
        total += length(ps.x)
    end
    if length(w.next) < total
        resize!(w.next, total); resize!(w.ent_pop, total); resize!(w.ent_idx, total)
    end
    fill!(w.head, Int32(0))
    xmin, xmax, ymin, ymax = w.box
    ncx, ncy = w.ncx, w.ncy
    cwx = (xmax - xmin) / ncx
    cwy = (ymax - ymin) / ncy
    ne = 0
    for (k, ps) in enumerate(w.pops)
        ps.p.binds || continue
        for i in 1:ps.n
            ps.partner[i] == 0 || continue
            cx = clamp(floor(Int, (ps.x[i] - xmin) / cwx), 0, ncx - 1)
            cy = clamp(floor(Int, (ps.y[i] - ymin) / cwy), 0, ncy - 1)
            c = cy * ncx + cx + 1
            ne += 1
            w.ent_pop[ne] = k
            w.ent_idx[ne] = i
            w.next[ne] = w.head[c]
            w.head[c] = ne
        end
    end
    return nothing
end

# 1b, one candidate: form the pair of entries e1, e2 if they are in contact and neither has bound
function _try_pair!(w::SimWorld, e1::Int, e2::Int, t0::Float64, t1::Float64)
    ka, ia = Int(w.ent_pop[e1]), Int(w.ent_idx[e1])
    kb, ib = Int(w.ent_pop[e2]), Int(w.ent_idx[e2])
    A = w.pops[ka]
    B = w.pops[kb]
    (A.partner[ia] != 0 || B.partner[ib] != 0) && return nothing
    dk = w.dimers
    _contact(w, A, ia, B, ib) < dk.r_react || return nothing
    if kb < ka || (kb == ka && ib < ia)
        ka, ia, kb, ib = kb, ib, ka, ia
        A, B = B, A
    end
    rng = w.rng
    E = dk.k_on == Inf ? 0.0 : randexp(rng) / dk.k_on
    tform = t0 + E
    # Formation is a Poisson hazard while in contact, drawn fresh at each sub-step's contact test, not stored state:
    # a draw at or after t1 is discarded, and the next sub-step tests contact at t1 and draws again, which is exact
    # in law by memorylessness (and a draw exactly at t1 has probability 0).
    tform < t1 || return nothing
    (A.t_depart[ia] <= tform || B.t_depart[ib] <= tform) && return nothing
    _place_pair!(w, A, ia, B, ib, dk.d_dimer)
    tbd = tform + (dk.k_off > 0 ? randexp(rng) / dk.k_off : Inf)
    τc = max(A.p.lifetime, B.p.lifetime)
    tdep = isfinite(τc) ? tform + τc * randexp(rng) : Inf
    A.partner[ia] = ib; A.partner_pop[ia] = kb; A.partner_id[ia] = B.id[ib]
    B.partner[ib] = ia; B.partner_pop[ib] = ka; B.partner_id[ib] = A.id[ia]
    A.t_form[ia] = tform; B.t_form[ib] = tform
    vis = A.state[ia] == 1 && A.m[ia] > 0 && B.state[ib] == 1 && B.m[ib] > 0
    A.vis_form[ia] = vis; B.vis_form[ib] = vis
    A.partner_f[ia] = B.id[ib]; B.partner_f[ib] = A.id[ia]
    A.t_break_due[ia] = tbd; B.t_break_due[ib] = tbd
    A.t_depart[ia] = tdep; B.t_depart[ib] = tdep
    return nothing
end

# 1b: visit candidates cell by cell (x fastest): own list, then the forward cells (+x), (-x,+y), (+y), (+x,+y)
function _form_pairs!(w::SimWorld, t0::Float64, t1::Float64)
    ncx, ncy = w.ncx, w.ncy
    periodic = w.boundary === :periodic
    for cy in 1:ncy, cx in 1:ncx
        e = Int(w.head[(cy - 1) * ncx + cx])
        while e != 0
            e2 = Int(w.next[e])
            while e2 != 0
                _try_pair!(w, e, e2, t0, t1)
                e2 = Int(w.next[e2])
            end
            for (ox, oy) in ((1, 0), (-1, 1), (0, 1), (1, 1))
                (ox != 0 && ncx == 1) && continue
                (oy != 0 && ncy == 1) && continue
                nx = cx + ox
                ny = cy + oy
                if periodic
                    nx = nx < 1 ? ncx : (nx > ncx ? 1 : nx)
                    ny = ny > ncy ? 1 : ny
                elseif nx < 1 || nx > ncx || ny > ncy
                    continue
                end
                e2 = Int(w.head[(ny - 1) * ncx + nx])
                while e2 != 0
                    _try_pair!(w, e, e2, t0, t1)
                    e2 = Int(w.next[e2])
                end
            end
            e = Int(w.next[e])
        end
    end
    return nothing
end

# 1c: split every pair whose break falls in [t0, t1) before its departure; the unbinding rule keeps a
# split pair from re-forming in place; each member draws a fresh departure time
function _split_pairs!(w::SimWorld, t0::Float64, t1::Float64)
    dk = w.dimers
    rng = w.rng
    s = _split_sep(w.box, dk.r_react)
    for (k, ps) in enumerate(w.pops)
        for i in 1:ps.n
            ip = Int(ps.partner[i])
            ip == 0 && continue
            kp = Int(ps.partner_pop[i])
            (k < kp || (k == kp && i < ip)) || continue
            tbd = ps.t_break_due[i]
            (tbd < t1 && tbd < ps.t_depart[i]) || continue
            Q = w.pops[kp]
            ps.t_break_f[i] = tbd; Q.t_break_f[ip] = tbd
            ps.partner[i] = 0; ps.partner_pop[i] = 0; ps.partner_id[i] = 0
            Q.partner[ip] = 0; Q.partner_pop[ip] = 0; Q.partner_id[ip] = 0
            _contact(w, ps, i, Q, ip) < dk.r_react && _place_pair!(w, ps, i, Q, ip, s)
            ps.t_depart[i] = isfinite(ps.p.lifetime) ? tbd + ps.p.lifetime * randexp(rng) : Inf
            Q.t_depart[ip] = isfinite(Q.p.lifetime) ? tbd + Q.p.lifetime * randexp(rng) : Inf
        end
    end
    return nothing
end

# Dimer events of one sub-step [t0, t1), at t0 positions: 1a cell list, 1b formation, 1c dissociation
function _dimer_events!(w::SimWorld, t0::Float64, t1::Float64)
    _build_cells!(w)
    _form_pairs!(w, t0, t1)
    _split_pairs!(w, t0, t1)
    return nothing
end
