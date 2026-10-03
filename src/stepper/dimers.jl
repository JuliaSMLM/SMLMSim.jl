# Cross-population dimers: cell list, formation, dissociation with the unbinding rule, pair motion

# Cells per axis: the side is at least r_react; an axis with fewer than 3 cells gets 1 (no neighbour visits)
function _cells_per_axis(L::Float64, r::Float64)
    n = clamp(floor(Int, L / r), 1, 128)
    return n < 3 ? 1 : n
end

# Put the members of a pair (A, ia) < (B, ib) at centre ± (s/2) u; reflecting walls clamp the centre, periodic wrap
@inline function _put_centered!(w::SimWorld, A::PopState, ia::Int, B::PopState, ib::Int,
                                cx::Float64, cy::Float64, ux::Float64, uy::Float64, s::Float64)
    xmin, xmax, ymin, ymax = w.box
    periodic = w.boundary === :periodic
    if !periodic
        hx = 0.5 * s * abs(ux)
        hy = 0.5 * s * abs(uy)
        cx = clamp(cx, xmin + hx, xmax - hx)
        cy = clamp(cy, ymin + hy, ymax - hy)
    end
    A.x[ia] = _fold(cx - 0.5 * s * ux, xmin, xmax, periodic)
    A.y[ia] = _fold(cy - 0.5 * s * uy, ymin, ymax, periodic)
    B.x[ib] = _fold(cx + 0.5 * s * ux, xmin, xmax, periodic)
    B.y[ib] = _fold(cy + 0.5 * s * uy, ymin, ymax, periodic)
    return nothing
end

# Place the pair (A, ia) < (B, ib) at in-plane separation s along its current axis (the 0.7.3 rule);
# returns the unit axis from A to B. An anchored pair (D_dimer = :min with an immobile member) keeps the anchor.
function _place_pair!(w::SimWorld, A::PopState, ia::Int, B::PopState, ib::Int, s::Float64)
    xmin, xmax, ymin, ymax = w.box
    periodic = w.boundary === :periodic
    Lx, Ly = xmax - xmin, ymax - ymin
    if w.dimers.D_dimer === :min && min(A.D[ia], B.D[ib]) == 0
        a_anchor = A.D[ia] == 0 && (B.D[ib] != 0 || A.id[ia] < B.id[ib])
        P, ip, Q, iq = a_anchor ? (A, ia, B, ib) : (B, ib, A, ia)
        dx = _sep(Q.x[iq] - P.x[ip], Lx, periodic)
        dy = _sep(Q.y[iq] - P.y[ip], Ly, periodic)
        n = hypot(dx, dy)
        ux, uy = n > 0 ? (dx / n, dy / n) : (1.0, 0.0)
        nx = P.x[ip] + s * ux
        ny = P.y[ip] + s * uy
        if periodic
            nx = _wrap(nx, xmin, xmax)
            ny = _wrap(ny, ymin, ymax)
        else
            (nx < xmin || nx > xmax) && (nx = P.x[ip] - s * ux)
            (ny < ymin || ny > ymax) && (ny = P.y[ip] - s * uy)
        end
        Q.x[iq] = nx
        Q.y[iq] = ny
        return a_anchor ? (ux, uy) : (-ux, -uy)
    end
    dx = _sep(B.x[ib] - A.x[ia], Lx, periodic)
    dy = _sep(B.y[ib] - A.y[ia], Ly, periodic)
    n = hypot(dx, dy)
    ux, uy = n > 0 ? (dx / n, dy / n) : (1.0, 0.0)
    _put_centered!(w, A, ia, B, ib, A.x[ia] + 0.5 * dx, A.y[ia] + 0.5 * dy, ux, uy, s)
    return ux, uy
end

# A bound pair moves once, at its member with the larger (population, index); an anchored pair does not move
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
    xmin, xmax, ymin, ymax = w.box
    periodic = w.boundary === :periodic
    dx = _sep(ps.x[i] - Q.x[ip], xmax - xmin, periodic)
    dy = _sep(ps.y[i] - Q.y[ip], ymax - ymin, periodic)
    sd = sqrt(2 * Dc * h)
    cx = _fold(Q.x[ip] + 0.5 * dx + sd * ξ1, xmin, xmax, periodic)
    cy = _fold(Q.y[ip] + 0.5 * dy + sd * ξ2, ymin, ymax, periodic)
    θ = Q.θ[ip] + sqrt(2 * dk.D_rot * h) * ξ3
    Q.θ[ip] = θ
    ps.θ[i] = θ
    _put_centered!(w, Q, ip, ps, i, cx, cy, cos(θ), sin(θ), dk.d_dimer)
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
    periodic = w.boundary === :periodic
    dx = _sep(A.x[ia] - B.x[ib], w.box[2] - w.box[1], periodic)
    dy = _sep(A.y[ia] - B.y[ib], w.box[4] - w.box[3], periodic)
    dz = A.z[ia] - B.z[ib]
    dx * dx + dy * dy + dz * dz < dk.r_react^2 || return nothing
    if kb < ka || (kb == ka && ib < ia)
        ka, ia, kb, ib = kb, ib, ka, ia
        A, B = B, A
    end
    rng = w.rng
    E = dk.k_on == Inf ? 0.0 : randexp(rng) / dk.k_on
    tform = t0 + E
    tform < t1 || return nothing
    (A.t_depart[ia] <= tform || B.t_depart[ib] <= tform) && return nothing
    ux, uy = _place_pair!(w, A, ia, B, ib, dk.d_dimer)
    θ = atan(uy, ux)
    tbd = tform + (dk.k_off > 0 ? randexp(rng) / dk.k_off : Inf)
    τc = max(A.p.lifetime, B.p.lifetime)
    tdep = isfinite(τc) ? tform + τc * randexp(rng) : Inf
    A.partner[ia] = ib; A.partner_pop[ia] = kb; A.partner_id[ia] = B.id[ib]
    B.partner[ib] = ia; B.partner_pop[ib] = ka; B.partner_id[ib] = A.id[ia]
    A.θ[ia] = θ; B.θ[ib] = θ
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
    periodic = w.boundary === :periodic
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
            dx = _sep(ps.x[i] - Q.x[ip], w.box[2] - w.box[1], periodic)
            dy = _sep(ps.y[i] - Q.y[ip], w.box[4] - w.box[3], periodic)
            dz = ps.z[i] - Q.z[ip]
            dx * dx + dy * dy + dz * dz < dk.r_react^2 &&
                _place_pair!(w, ps, i, Q, ip, dk.r_react * (1 + 1e-9))
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
