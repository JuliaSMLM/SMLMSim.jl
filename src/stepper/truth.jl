# Per-exposure truth rows: one FrameTruth per emitter present in the recorded exposure

"""
    SMLMSim.frame_truth(world) -> AbstractVector{FrameTruth}

The [`FrameTruth`](@ref) rows of the last [`step!`](@ref): one per emitter present at any time in
the exposure. A view of the world's buffer, valid until the next `step!`.
"""
frame_truth(w::SimWorld) = view(w.truth, 1:w.n_truth)

# Start the recorded part of a step!: empty buffer, every present emitter's accumulators restarted
function _begin_truth!(w::SimWorld)
    w.n_truth = 0
    for ps in w.pops
        for i in 1:ps.n
            _reset_acc!(ps, i)
        end
    end
    return nothing
end

# Append the row of emitter i of population k; `departed` marks a removal by departure at t_depart
function _write_row!(w::SimWorld, k::Int, ps::PopState, i::Int, departed::Bool)
    w.n_truth == length(w.truth) && resize!(w.truth, 2 * max(1, length(w.truth)))
    T = w.t_b - w.t_a
    sph = ps.sph[i]
    tp = ps.t_present[i]
    x = ps.xr[i] + (sph > 0 ? ps.sx[i] / sph : tp > 0 ? ps.sxp[i] / tp : 0.0)
    y = ps.yr[i] + (sph > 0 ? ps.sy[i] / sph : tp > 0 ? ps.syp[i] / tp : 0.0)
    if w.boundary === :periodic
        x = _wrap(x, w.box[1], w.box[2])
        y = _wrap(y, w.box[3], w.box[4])
    end
    tb = ps.t_birth[i]
    w.n_truth += 1
    pid = ps.partner_id[i]
    tf = ps.t_form[i]
    w.truth[w.n_truth] = FrameTruth(w.frame + 1, ps.id[i], Int32(k), ps.m[i], x, y, ps.z[i], sph,
                                    ps.t_lit[i] / T, tp > 0 ? ps.sI[i] / tp : NaN,
                                    tb >= w.t_a ? tb : NaN, ps.t_bleach_f[i],
                                    departed ? ps.t_depart[i] : NaN, pid, pid == 0 ? Int32(0) : ps.partner_pop[i],
                                    ps.t_bound[i] / T, tf >= w.t_a ? tf : NaN, ps.t_break_f[i], 0.0,
                                    false, false, false)
    return nothing
end

# End of a step!: rows for every emitter still present, then the overlap pass
function _finish_truth!(w::SimWorld)
    for (k, ps) in enumerate(w.pops)
        for i in 1:ps.n
            _write_row!(w, k, ps, i, false)
        end
    end
    w.merge_radius > 0 && _mark_overlaps!(w)
    return nothing
end

# Minimum-image separation of one coordinate under :periodic
@inline _sep(d::Float64, L::Float64, periodic::Bool) = periodic ? d - L * round(d / L) : d

# Flag both rows of every pair of :signal rows within merge_radius that are not partners
function _mark_overlaps!(w::SimWorld)
    n = w.n_truth
    length(w.overlap_scratch) < n && resize!(w.overlap_scratch, length(w.truth))
    flags = w.overlap_scratch
    fill!(flags, false)
    rows = w.truth
    periodic = w.boundary === :periodic
    Lx = w.box[2] - w.box[1]
    Ly = w.box[4] - w.box[3]
    r2 = w.merge_radius^2
    for a in 1:n-1
        ra = rows[a]
        w.pops[ra.pop].p.layer === :signal || continue
        for b in a+1:n
            rb = rows[b]
            w.pops[rb.pop].p.layer === :signal || continue
            (ra.partner == rb.id || rb.partner == ra.id) && continue
            dx = _sep(ra.x - rb.x, Lx, periodic)
            dy = _sep(ra.y - rb.y, Ly, periodic)
            if dx * dx + dy * dy <= r2
                flags[a] = true
                flags[b] = true
            end
        end
    end
    for a in 1:n
        flags[a] || continue
        rows[a] = _with_overlap(rows[a])
    end
    return nothing
end

@inline function _with_overlap(r::FrameTruth)
    N = fieldcount(FrameTruth)
    return FrameTruth(ntuple(k -> k == N ? true : getfield(r, k), Val(N))...)
end
