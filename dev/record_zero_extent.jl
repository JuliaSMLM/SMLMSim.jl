# Record simulate(::StaticSMLMConfig) output for zero-extent patterns, as TOML, for
# test/test_field_extent.jl ("zero-extent patterns as at 64e88f7"). SMLMDeepFit's retrain is pinned
# to commit 64e88f7; zero-extent patterns must give it bitwise the same output for the same seed.
#
# The reference test/field_extent_64e88f7.toml was recorded at 64e88f7 on Julia 1.10.11, 1.11.9,
# 1.12.7 and 1.13.0 (identical on all four). If that test fails only because the toolchain changed
# its random numbers (a new Julia or Distributions release), re-record it from 64e88f7 with the new
# toolchain:
#
#   git archive 64e88f7 | tar -x -C /path/to/smlmsim-64e88f7
#   julia dev/record_zero_extent.jl /path/to/smlmsim-64e88f7 test/field_extent_64e88f7.toml 64e88f7
#
# A failure that re-recording at 64e88f7 does not clear is a real change in placement.
# The cases below must stay identical to ZERO_EXTENT_CASES in test/test_field_extent.jl.
using Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.develop(path=ARGS[1]; io=devnull)
Pkg.add(["Random", "TOML"]; io=devnull)
using SMLMSim, Random, TOML

const ZERO_EXTENT_CASES = Dict(
    "point2d" => (seed=64881, pattern=Nmer2D(n=1, d=0.0), camera=IdealCamera(8, 8, 0.1), ndims=2),
    "point3d" => (seed=64882, pattern=Nmer3D(n=1, d=0.0), camera=IdealCamera(8, 8, 0.1), ndims=3),
    "point2d_offset" => (seed=64883, pattern=Nmer2D(n=1, d=0.0),
                         camera=IdealCamera(collect(5.0:0.1:5.8), collect(-2.0:0.1:-1.2)), ndims=2),
)

function run_zero_extent_case(case)
    Random.seed!(case.seed)
    cfg = StaticSMLMConfig(density=30.0, nframes=5, ndatasets=1, minphotons=0, ndims=case.ndims,
                           zrange=[-0.5, 0.5])
    noisy, info = SMLMSim.simulate(cfg; pattern=case.pattern, camera=case.camera,
                                   molecule=GenericFluor(photons=1e4, k_off=50.0, k_on=20.0))
    t, e = info.smld_true.emitters, noisy.emitters
    out = Dict{String,Any}(
        "n_patterns" => info.n_patterns,
        "true_x" => [v.x for v in t], "true_y" => [v.y for v in t],
        "noisy_x" => [v.x for v in e], "noisy_y" => [v.y for v in e],
        "noisy_photons" => [v.photons for v in e], "noisy_frame" => [v.frame for v in e],
        "noisy_track" => [v.track_id for v in e], "noisy_id" => [v.id for v in e])
    if case.ndims == 3
        out["true_z"] = [v.z for v in t]
        out["noisy_z"] = [v.z for v in e]
    end
    return out
end

out = Dict{String,Any}("julia" => "recorded on Julia $(VERSION)", "source" => ARGS[3])
for (name, case) in ZERO_EXTENT_CASES
    out[name] = run_zero_extent_case(case)
end
open(ARGS[2], "w") do io
    TOML.print(io, out; sorted=true)
end
println("wrote ", ARGS[2], ": ", join(["$k $(length(out[k]["noisy_x"])) records" for k in keys(ZERO_EXTENT_CASES)], ", "))
