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
    length(w.vis_partner) < length(w.truth) && resize!(w.vis_partner, length(w.truth))
    w.vis_partner[w.n_truth] = ps.partner_f[i]
    pid = ps.partner_id[i]
    tf = ps.t_form[i]
    w.truth[w.n_truth] = FrameTruth(w.frame + 1, ps.id[i], Int32(k), ps.m[i], x, y, ps.z[i], sph,
                                    ps.t_lit[i] / T, tp > 0 ? ps.sI[i] / tp : NaN,
                                    tb >= w.t_a ? tb : NaN, ps.t_bleach_f[i],
                                    departed ? ps.t_depart[i] : NaN, pid, pid == 0 ? Int32(0) : ps.partner_pop[i],
                                    ps.t_bound[i] / T, tf >= w.t_a ? tf : NaN, ps.t_break_f[i], ps.t_litb[i] / T,
                                    tf >= w.t_a && ps.vis_form[i], false, false)
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
    w.dimers === nothing || _mark_visible_pairs!(w)
    return nothing
end

# Minimum-image separation of one coordinate under :periodic
@inline _sep(d::Float64, L::Float64, periodic::Bool) = periodic ? d - L * round(d / L) : d

# Rows `r` and `o` were a pair in this exposure: partners at its end, or bound during it with `o` as r's last
# partner (`vpr`, its vis_partner entry), so a pair that broke in the exposure is not an overlap
@inline _paired(r::FrameTruth, vpr::Int, o::FrameTruth) = r.partner == o.id || (vpr == o.id && r.bound > 0)

# Flag both rows of every pair of :signal rows within merge_radius that were not a pair in this exposure
function _mark_overlaps!(w::SimWorld)
    n = w.n_truth
    length(w.overlap_scratch) < n && resize!(w.overlap_scratch, length(w.truth))
    flags = w.overlap_scratch
    fill!(flags, false)
    rows = w.truth
    vp = w.vis_partner
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
            (_paired(ra, vp[a], rb) || _paired(rb, vp[b], ra)) && continue
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
        rows[a] = _with_flag(rows[a], fieldcount(FrameTruth))
    end
    return nothing
end

@inline function _with_flag(r::FrameTruth, j::Int)
    N = fieldcount(FrameTruth)
    return FrameTruth(ntuple(k -> k == j ? true : getfield(r, k), Val(N))...)
end

# vis_bound: a row bound at some point in this frame is visible as a pair when it and its last partner's row of
# this frame both have lit_bound > 0 (vis_partner, so a pair that broke counts). Runs after all rows are written:
# one pass over the N rows, and each of the B bound, lit rows scans the N rows for its partner's row, O(N + B·N).
function _mark_visible_pairs!(w::SimWorld)
    n = w.n_truth
    rows = w.truth
    vp = w.vis_partner
    for a in 1:n
        ra = rows[a]
        (vp[a] != 0 && ra.lit_bound > 0) || continue
        for b in 1:n
            rb = rows[b]
            if rb.id == vp[a] && rb.lit_bound > 0
                rows[a] = _with_flag(ra, fieldcount(FrameTruth) - 1)
                break
            end
        end
    end
    return nothing
end

"""
    params_dict(world) -> Dict{String,Any}

The configuration of `world` as a flat dictionary whose values are only `String`, `Bool`, `Int`,
`Float64` or `Vector{Float64}`, so each can be written as an HDF5 attribute. Keys:

- `"smlmsim.version"`, `"rng.type"`;
- `"world.n_sub"`, `"world.boundary"`, `"world.box_um"` (xmin, xmax, ymin, ymax), `"world.merge_radius"`;
- `"camera.type"`, `"camera.nx"`, `"camera.ny"`, `"camera.pixel_size_um"` and `"camera.offset"`, `"camera.gain"`,
  `"camera.readnoise"`, `"camera.qe"` for an sCMOS camera: the scalar, or `".mean"` appended to the key for a
  per-pixel map;
- `"pop<k>.<field>"` for every field of the `k`th [`Population`](@ref): `mobility` as `.mobility.fraction` and
  `.mobility.D`, `fluor` as `.fluor.gamma` and the rate matrix `.fluor.q` (row-major) with its size `.fluor.q.n`,
  `z` as a two-element vector, and `psf` as `.psf.sigma_um` or `.psf.stamp.z_min`, `.z_max`, `.z_step`,
  `.radius` and `.oversample`;
- `"dimers.k_on"`, `.r_react`, `.k_off`, `.D_rot`, `.d_dimer` and `"dimers.D_dimer"` (a number, or the String
  `"min"`), only when the world has `dimers`;
- `"bg.<field>"` for every field of the [`BackgroundModel`](@ref) (the level as a number, or
  `string(distribution)`), only when the world has a background.

The commit of the code is not recoverable from an installed package: a caller that needs it adds its own
entries with `merge`, for example the seed.
"""
function params_dict(w::SimWorld)
    d = Dict{String,Any}()
    d["smlmsim.version"] = string(pkgversion(parentmodule(@__MODULE__)))
    d["rng.type"] = string(typeof(w.rng))
    d["world.n_sub"] = w.n_sub
    d["world.boundary"] = string(w.boundary)
    d["world.box_um"] = collect(Float64, w.box)
    d["world.merge_radius"] = w.merge_radius
    c = w.camera
    d["camera.type"] = string(nameof(typeof(c)))
    d["camera.nx"] = size(w.signal, 2)
    d["camera.ny"] = size(w.signal, 1)
    d["camera.pixel_size_um"] = w.px
    if c isa SCMOSCamera
        for f in (:offset, :gain, :readnoise, :qe)
            v = getfield(c, f)
            if v isa AbstractMatrix
                d["camera.$f.mean"] = Float64(sum(v) / length(v))
            else
                d["camera.$f"] = Float64(v)
            end
        end
    end
    for (k, ps) in enumerate(w.pops)
        p = ps.p
        pre = "pop$k."
        d[pre * "name"] = string(p.name)
        d[pre * "layer"] = string(p.layer)
        d[pre * "density"] = p.density
        d[pre * "lifetime"] = p.lifetime
        d[pre * "birth_rate"] = p.birth_rate
        d[pre * "mobility.fraction"] = Float64[f for (f, _) in p.mobility]
        d[pre * "mobility.D"] = Float64[D for (_, D) in p.mobility]
        d[pre * "fluor.gamma"] = Float64(p.fluor.γ)
        q = Matrix{Float64}(p.fluor.q)
        d[pre * "fluor.q"] = vec(permutedims(q))
        d[pre * "fluor.q.n"] = size(q, 1)
        d[pre * "brightness_sigma"] = p.brightness_sigma
        d[pre * "budget"] = p.budget
        d[pre * "multiplicity"] = p.multiplicity
        d[pre * "z"] = Float64[p.z[1], p.z[2]]
        if p.psf isa StampTable
            zs = p.psf.zs
            d[pre * "psf.stamp.z_min"] = Float64(first(zs))
            d[pre * "psf.stamp.z_max"] = Float64(last(zs))
            d[pre * "psf.stamp.z_step"] = length(zs) > 1 ? Float64(zs[2] - zs[1]) : 0.0
            d[pre * "psf.stamp.radius"] = p.psf.radius
            d[pre * "psf.stamp.oversample"] = p.psf.oversample
        else
            d[pre * "psf.sigma_um"] = Float64(p.psf.σ)
        end
        d[pre * "brightness_jitter"] = p.brightness_jitter
        d[pre * "jitter_time"] = p.jitter_time
        d[pre * "binds"] = p.binds
    end
    dk = w.dimers
    if dk !== nothing
        d["dimers.k_on"] = dk.k_on
        d["dimers.r_react"] = dk.r_react
        d["dimers.k_off"] = dk.k_off
        d["dimers.D_rot"] = dk.D_rot
        d["dimers.d_dimer"] = dk.d_dimer
        d["dimers.D_dimer"] = dk.D_dimer === :min ? "min" : Float64(dk.D_dimer)
    end
    if w.bg !== nothing
        m = w.bg.model
        for f in fieldnames(typeof(m))
            v = getfield(m, f)
            d["bg.$f"] = v isa Real ? Float64(v) : string(v)
        end
    end
    return d
end
