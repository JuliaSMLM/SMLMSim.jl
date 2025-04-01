using Revise
using SMLMSim
using SMLMData

# 2D Example
camera = IdealCamera(1:512, 1:512, 0.1)
pattern = SMLMSim.Nmer2D(n=6, d=0.2)
params = StaticSMLMParams(ρ=1.0)
smld_true, smld_model, smld_noisy = simulate(
    params,
    pattern=pattern,
    camera=camera
)

# 3D Example
pattern3d = SMLMSim.Nmer3D(n=8, d=0.3)
params3d = StaticSMLMParams(
    ρ=1.0,
    ndims=3,
    zrange=[-2.0, 2.0]
)
smld_true, smld_model, smld_noisy = simulate(
    params3d,
    pattern=pattern3d,
    camera=camera
)

# Complex Pattern Example
line3d = SMLMSim.Line3D(λ=5.0, endpoints=[(-1.0, 0.0, -0.5), (1.0, 0.0, 0.5)])
params_line = StaticSMLMParams(
    ρ=0.5,
    σ_psf=0.15,
    ndims=3
)
smld_true, smld_model, smld_noisy = simulate(
    params_line,
    pattern=line3d,
    camera=camera
)