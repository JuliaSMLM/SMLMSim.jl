# Records the golden values of test/test_brightness.jl ("P1: absent keyword is main's output and
# random stream") from SMLMSim at commit b016876, main before the brightness keyword:
#
#   git archive b016876 | tar -x -C <dir>
#   julia dev/record_golden_brightness.jl <dir>
#
# The output was identical on Julia 1.10.11, 1.11.9, 1.12.7 and 1.13.0, with Distributions
# 0.25.134, SMLMData 0.7.1 and MicroscopePSFs 0.5.7. The two case functions below must stay
# identical to br_golden_simulate and br_golden_kinetic in test/test_brightness.jl (there they
# also pass the keywords under test).
using Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.develop(path=ARGS[1]; io=devnull)
Pkg.add(["Random"]; io=devnull)
using SMLMSim, Random

const GOLDEN_FLUOR = GenericFluor(photons=1e4, k_off=50.0, k_on=5.0)

function br_golden_simulate(; kw...)
    Random.seed!(20261108)
    params = StaticSMLMConfig(density=2.0, nframes=40, ndatasets=3, minphotons=50)
    noisy, info = simulate(params; pattern=Nmer2D(n=6, d=0.2), camera=IdealCamera(16, 16, 0.1),
                           molecule=GOLDEN_FLUOR, kw...)
    return noisy.emitters, info, rand()
end

function br_golden_kinetic(; kw...)
    emitters = [Emitter2DFit{Float64}(0.1 * k, 0.2 * k, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0;
                                      σ_xy=0.0, frame=1, dataset=1,
                                      track_id=isodd(k) ? k : 0, id=k) for k in 1:8]
    smld = BasicSMLD(emitters, IdealCamera(32, 32, 0.1), 1, 1, Dict{String,Any}())
    Random.seed!(20261109)
    out = kinetic_model(smld, GenericFluor(1e4, [-30.0 30.0; 5.0 -5.0]), 40, 20.0;
                        ndatasets=2, minphotons=0.0, kw...)
    return out.emitters, out, rand()
end

# What the test compares: the record count, every record's photons and x summed in order
# (foldl, so no reassociation), the first, middle and last records, and the next rand().
function fingerprint(em, next)
    n = length(em)
    rec(e) = (x=e.x, y=e.y, photons=e.photons, frame=e.frame, dataset=e.dataset, track_id=e.track_id)
    return (n=n, sum_photons=foldl(+, (e.photons for e in em)), sum_x=foldl(+, (e.x for e in em)),
            first=rec(em[1]), middle=rec(em[cld(n, 2)]), last=rec(em[n]), next_rand=next)
end

println("# julia ", VERSION, ", SMLMSim at ", ARGS[1])
for name in ("Distributions", "SMLMData", "MicroscopePSFs")
    v = Pkg.dependencies()[findfirst(p -> p.name == name, Pkg.dependencies())].version
    println("# ", name, " ", v)
end
em, _, next = br_golden_simulate()
println("simulate = ", repr(fingerprint(em, next)))
em, _, next = br_golden_kinetic()
println("kinetic = ", repr(fingerprint(em, next)))
