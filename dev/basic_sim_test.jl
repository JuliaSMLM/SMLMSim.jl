using Pkg
Pkg.activate("dev")

using Revise
using SMLMSim
using SMLMData
using MicroscopePSFs
using CairoMakie

# 2D Example
camera = IdealCamera(1:512, 1:512, 0.1)
pattern = SMLMSim.Nmer2D(n=6, d=0.2)
params = StaticSMLMParams(ρ=1.0, nframes = 100)
smld_true, smld_model, smld_noisy = simulate(
    params,
    pattern=pattern,
    camera=camera
)

# generate images
psf = AiryPSF(1.2, .555)
psf_spline = SplinePSF(psf, -1:.1:1, -1:.1:1)
img = gen_images(smld_model, psf_spline; dataset=1,  support = 1.0)

image(img[:, :, 1])




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

