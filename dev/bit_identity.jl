# Bit identity of diffusion `simulate` output between two versions of SMLMSim on one Julia.
# RNG streams are not promised across Julia versions, so this is a dev check, not a test.
# Run it once per version (each worktree's own project), then compare the two files:
#
#   julia --project=<worktree A> dev/bit_identity.jl run A.jls
#   julia --project=<worktree B> dev/bit_identity.jl run B.jls
#   julia dev/bit_identity.jl compare A.jls B.jls
#
# Every case uses only keywords that exist in both versions. Compared: every field of every
# emitter, in order, plus the next rand() after the run (the RNG stream consumed).
using Serialization, Random

if ARGS[1] == "run"
    using SMLMSim
    fields(e) = ntuple(i -> getfield(e, i), fieldcount(typeof(e)))
    dimers(; kw...) = DiffusionSMLMConfig(; density=40.0, box_size=2.0, diff_monomer=0.3, diff_dimer=0.1,
        diff_dimer_rot=0.5, k_off=5.0, r_react=0.05, d_dimer=0.03, dt=0.001, t_max=0.1,
        camera_framerate=100.0, camera_exposure=0.01, kw...)
    function case(name, cfg; kw...)
        Random.seed!(20260929)
        smld, _ = simulate(cfg; kw...)
        return (name, [fields(e) for e in smld.emitters], rand()), smld
    end
    out = Any[]
    r, _ = case("defaults", DiffusionSMLMConfig(); γ=1e5)
    push!(out, r)
    r, s = case("dimers periodic", dimers(); γ=500.0)
    push!(out, r)
    r, _ = case("dimers reflecting 3D", dimers(ndims=3, density=8.0, boundary="reflecting"); γ=500.0)
    push!(out, r)
    r, _ = case("mobility mixture", dimers(monomer_mobility=[(0.5, 0.3), (0.5, 0.0)]); γ=500.0)
    push!(out, r)
    r, _ = case("continuation", dimers(); starting_conditions=s)
    push!(out, r)
    serialize(ARGS[2], out)
    println("wrote ", length(out), " cases to ", ARGS[2], " (SMLMSim ", pkgversion(SMLMSim), ", Julia ", VERSION, ")")
else
    function compare(a, b)
        same = length(a) == length(b)
        for ((n, ea, ra), (nb, eb, rb)) in zip(a, b)
            ok = n == nb && isequal(ea, eb) && ra === rb
            same &= ok
            println(rpad(n, 22), length(ea), " records: ", ok ? "identical" : "DIFFERENT")
        end
        println(same ? "BIT-IDENTICAL: all $(length(a)) cases" : "NOT IDENTICAL")
    end
    compare(deserialize(ARGS[2]), deserialize(ARGS[3]))
end
