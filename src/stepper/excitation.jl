# Focused-spot excitation: Gaussian-beam spots with optional evanescent part, tilt and per-spot switching

"""
    Spot(; x, y, σ, gain, z_R, t_on = -Inf, t_off = Inf)

One focused spot of a [`SpotExcitation`](@ref). Built by keyword; `x`, `y`, `σ` (μm), `gain`
(relative intensity at the centre in focus) and `z_R` (μm) are required. `σ` is the Gaussian σ
of the focal intensity, so the waist is `w0 = 2σ`; `z_R` is the Rayleigh range (`Inf`: a
z-independent column). The spot contributes only for `t_on <= t < t_off`, in seconds on the world
clock, the same `t` the excitation receives in [`step!`](@ref); the defaults keep it always on.

Throws `ArgumentError` unless `σ > 0`, `gain >= 0`, `z_R > 0`, `t_on < t_off` and neither time is NaN.
"""
struct Spot
    x::Float64
    y::Float64
    σ::Float64
    gain::Float64
    z_R::Float64
    t_on::Float64
    t_off::Float64

    function Spot(; x::Real, y::Real, σ::Real, gain::Real, z_R::Real, t_on::Real=-Inf, t_off::Real=Inf)
        σ > 0 || throw(ArgumentError("σ must be > 0, got $σ"))
        gain >= 0 || throw(ArgumentError("gain must be >= 0, got $gain"))
        z_R > 0 || throw(ArgumentError("z_R must be > 0 (Inf allowed), got $z_R"))
        (isnan(t_on) || isnan(t_off)) && throw(ArgumentError("t_on and t_off must not be NaN"))
        t_on < t_off || throw(ArgumentError("t_on must be < t_off, got t_on = $t_on, t_off = $t_off"))
        return new(Float64(x), Float64(y), Float64(σ), Float64(gain), Float64(z_R), Float64(t_on), Float64(t_off))
    end
end

"""
    SpotExcitation(; spots, base = 1.0, λ = 0.642, n = 1.33, f_evan = 0.0, d_evan = 0.1, tilt = (0.0, 0.0))

Relative excitation of focused spots, for [`step!`](@ref): a callable `e(x, y, z, t)` in μm and s.
Built by keyword. Each [`Spot`](@ref) follows the Gaussian-beam law, with `s = σ √(1 + (z/z_R)²)`,

    I = base + Σ_k [(1 - f_evan) g_k + evanescent_k],   g_k = gain (σ/s)² exp(-r²/(2s²)),

over the spots with `t_on <= t < t_off`; `r` is the distance to the spot centre
`(x_k, y_k) + z tan θ (cos φ, sin φ)` for `tilt = (tan θ, φ)`, the oblique mean ray. The
evanescent term, `gain f_evan exp(-z/d_evan) exp(-r0²/(2σ²))` for `z >= 0` (`r0` the distance to the
untilted centre), is the fraction `f_evan` of the gain that lives within `d_evan` μm of the
interface. With `f_evan = 0` it is skipped and the result is the plain law.

- `base >= 0`: the baseline intensity (1 is the calibrated widefield/TIRF level).
- `λ` (μm), `n`: the wavelength and refractive index, used only to check each finite `z_R` against
  `π (2σ)² n/λ`: a warning names a spot whose `z_R` is more than 2x away.
- `f_evan` in [0, 1], `d_evan > 0` (μm; the rig's range is 0.08-0.3).
- `tilt`: `tan θ` (finite, `>= 0`; the rig measured about 2) and the azimuth `φ` (rad).

[`next_switch`](@ref) is the smallest `t_on` or `t_off` of the spots after `t`, so a switch inside
an exposure acts at its exact time. The call allocates nothing.

Throws `ArgumentError` unless `base >= 0`, `0 <= f_evan <= 1`, `d_evan > 0`, `λ > 0`, `n > 0`,
`tan θ` finite and `>= 0`, and `φ` finite.
"""
struct SpotExcitation
    base::Float64
    spots::Vector{Spot}
    λ::Float64
    n::Float64
    f_evan::Float64
    d_evan::Float64
    tilt::Tuple{Float64,Float64}
    tan_cos::Float64                      # tan θ cos φ and tan θ sin φ, so the call does no trigonometry
    tan_sin::Float64

    function SpotExcitation(; spots::Vector{Spot}, base::Real=1.0, λ::Real=0.642, n::Real=1.33,
                            f_evan::Real=0.0, d_evan::Real=0.1, tilt=(0.0, 0.0))
        base >= 0 || throw(ArgumentError("base must be >= 0, got $base"))
        0 <= f_evan <= 1 || throw(ArgumentError("f_evan must be in [0, 1], got $f_evan"))
        d_evan > 0 || throw(ArgumentError("d_evan must be > 0, got $d_evan"))
        λ > 0 || throw(ArgumentError("λ must be > 0, got $λ"))
        n > 0 || throw(ArgumentError("n must be > 0, got $n"))
        (isfinite(tilt[1]) && tilt[1] >= 0) || throw(ArgumentError("tilt[1] (tan θ) must be finite and >= 0, got $(tilt[1])"))
        isfinite(tilt[2]) || throw(ArgumentError("tilt[2] (the azimuth) must be finite, got $(tilt[2])"))
        for (k, sp) in enumerate(spots)
            isfinite(sp.z_R) || continue
            z_exp = π * (2 * sp.σ)^2 * n / λ
            (sp.z_R > 2 * z_exp || sp.z_R < z_exp / 2) &&
                @warn "spot $k: z_R = $(sp.z_R) μm differs by more than 2x from π (2σ)² n/λ = $z_exp μm"
        end
        tt, φ = Float64(tilt[1]), Float64(tilt[2])
        return new(Float64(base), copy(spots), Float64(λ), Float64(n), Float64(f_evan), Float64(d_evan),
                   (tt, φ), tt * cos(φ), tt * sin(φ))
    end
end

function (e::SpotExcitation)(x::Float64, y::Float64, z::Float64, t::Float64)
    I = e.base
    f = e.f_evan
    for sp in e.spots
        (sp.t_on <= t && t < sp.t_off) || continue
        s = sp.σ * sqrt(1 + (z / sp.z_R)^2)
        dx = x - (sp.x + z * e.tan_cos)
        dy = y - (sp.y + z * e.tan_sin)
        r2 = dx * dx + dy * dy
        g = sp.gain * (sp.σ / s)^2 * exp(-r2 / (2 * s^2))
        I += (1 - f) * g
        if f > 0 && z >= 0
            ex = x - sp.x
            ey = y - sp.y
            I += sp.gain * f * exp(-z / e.d_evan) * exp(-(ex * ex + ey * ey) / (2 * sp.σ^2))
        end
    end
    return I
end

function next_switch(e::SpotExcitation, t::Float64)
    ts = Inf
    for sp in e.spots
        sp.t_on > t && sp.t_on < ts && (ts = sp.t_on)
        sp.t_off > t && sp.t_off < ts && (ts = sp.t_off)
    end
    return ts
end
