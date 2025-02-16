using Revise
using SMLMSim
using SMLMData 



# 2D Example
camera = IdealCamera(1:512, 1:512, 0.1)
pattern = SMLMSim.Nmer2D(n=6, d=0.2)
smld_true, smld_model, smld_noisy = sim(
    ρ=1.0,
    pattern=pattern,
    camera=camera
)

# 3D Example
pattern3d = SMLMSim.Nmer3D(n=8, d=0.3)
smld_true, smld_model, smld_noisy = sim(
    ρ=1.0,
    pattern=pattern3d,
    camera=camera,
    zrange=[-2.0, 2.0]
)

# Complex Pattern Example
line3d = SMLMSim.Line3D(λ=5.0, endpoints=[(-1.0, 0.0, -0.5), (1.0, 0.0, 0.5)])
smld_true, smld_model, smld_noisy = sim(
    ρ=0.5,
    pattern=line3d,
    σ_psf=0.15,
    camera=camera
)