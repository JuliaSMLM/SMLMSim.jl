# Records the golden values of test/test_directed_motion.jl ("no velocity: main's output and random stream") from
# SMLMSim at commit b016876, main before directed motion:
#
#   git archive b016876 | tar -x -C <dir>
#   julia dev/record_golden_diffusion.jl <dir>
#
# The output was identical on Julia 1.10.11, 1.11.9, 1.12.7 and 1.13.0. Julia does not promise its random streams
# across versions, and a few rare randn draws differ between them (seed 4712 of the 3D case differed on 1.10), so the
# test keeps cases that agree on all four. The cases below must stay identical to the test's.
using Pkg
Pkg.activate(; temp=true, io=devnull)
Pkg.develop(path=ARGS[1]; io=devnull)
Pkg.add(["Random"]; io=devnull)
using SMLMSim, Random
const DM = SMLMSim.InteractionDiffusion

coords(e) = hasproperty(e, :z) ? (e.x, e.y, e.z) : (e.x, e.y)
fingerprint(s) = (n = length(s.emitters);
    (n, count(e -> e.state == :dimer, s.emitters), [coords(s.emitters[i]) for i in (1, n ÷ 2, n)],
     sum(Float64(e.x) for e in s.emitters), rand()))

Random.seed!(4711)
s, _ = simulate(DiffusionSMLMConfig(density=20.0, box_size=1.0, t_max=0.3, r_react=0.03, d_dimer=0.04, k_off=5.0,
                                    monomer_mobility=[(0.5, 0.0), (0.5, 0.1)]); γ=1000.0)
println("periodic 2D: ", repr(fingerprint(s)))

Random.seed!(4713)
s, _ = simulate(DiffusionSMLMConfig(ndims=3, boundary="reflecting", density=5.0, box_size=1.0, t_max=0.2,
                                    r_react=0.05, d_dimer=0.06, k_off=3.0, pair_mobility=:min,
                                    monomer_mobility=[(0.3, 0.0), (0.7, 0.2)]); γ=1000.0)
println("reflecting 3D: ", repr(fingerprint(s)))

Random.seed!(4713)
start = [DM.DiffusingEmitter2D{Float32}(Float32(0.1 * i), Float32(0.05 * i), 10.0, 0.0, 1, 1, i, :monomer, nothing)
         for i in 1:9]
s, _ = simulate(DiffusionSMLMConfig(boundary="reflecting", box_size=1.0, t_max=0.3, r_react=0.04, d_dimer=0.05,
                                    diff_monomer=0.2, k_off=4.0); starting_conditions=start, γ=1000.0)
println("Float32: ", repr(fingerprint(s)))
