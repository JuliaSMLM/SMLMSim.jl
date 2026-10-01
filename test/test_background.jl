using SMLMSim, Test, Statistics, Random, StableRNGs, MicroscopePSFs, Distributions, LinearAlgebra

const T = 0.01

cam(n, px) = IdealCamera(1:n, 1:n, px)
one_state(γ) = GenericFluor(; γ, q=zeros(1, 1))
two_state(γ, k_off, k_on) = GenericFluor(; γ, q=[-k_off k_off; k_on -k_on])

# log pattern of frame k of a movie, with the mean level and lognormal shift removed
function log_pattern(S, level, c)
    return [(log(S[i, j, k] / level[k]) + c^2 / 2) / c for i in axes(S, 1), j in axes(S, 2), k in axes(S, 3)]
end

@testset "background/level_mean" begin
    # contrast 0, flat illumination: the frame mean is the level exactly
    bg = BackgroundModel(level=87.7, jitter=0.03)
    S, O, L = gen_background(StableRNG(1), cam(32, 0.1), bg, 200; frame_time=T)
    @test all(isapprox(mean(S[:, :, k]), L[k]; rtol=1e-12) for k in 1:200)
    @test all(L .> 0)
    S0, _, L0 = gen_background(StableRNG(1), cam(32, 0.1), BackgroundModel(level=87.7), 5; frame_time=T)
    @test all(isapprox.(L0, 87.7 * T; rtol=1e-12))
    @test all(iszero, O)
    # contrast 0.5 with tau << T: the mean of the frame-mean ratio is 1
    bg = BackgroundModel(level=50.0, contrast=0.5, feature_size=0.3, correlation_time=1e-9)
    S, _, L = gen_background(StableRNG(2), cam(64, 0.1), bg, 400; frame_time=T)
    r = [mean(S[:, :, k]) / L[k] for k in 1:400]
    @test abs(mean(r) - 1) < 3 * std(r) / sqrt(400)
end

@testset "background/contrast" begin
    for c in (0.3, 0.6)
        bg = BackgroundModel(level=20.0, contrast=c, feature_size=0.3, correlation_time=1e-9)
        S, _, L = gen_background(StableRNG(3), cam(64, 0.1), bg, 100; frame_time=T)
        z = [log(S[i, j, k] / L[k]) for i in 1:64, j in 1:64, k in 1:100]
        @test isapprox(std(z), c; rtol=0.05)
    end
end

# Gaussian fit of the autocorrelation along both axes: slope of log rho against d^2 for rho > 0.5
function fit_sigma(z, px)
    zc = z .- mean(z)
    v = mean(abs2, zc)
    ds, lr = Float64[], Float64[]
    for d in 1:40
        ρ = 0.5 * (mean(zc[:, 1:end-d, :] .* zc[:, 1+d:end, :]) + mean(zc[1:end-d, :, :] .* zc[1+d:end, :, :])) / v
        ρ > 0.5 || break
        push!(ds, (d * px)^2)
        push!(lr, log(ρ))
    end
    slope = sum(ds .* lr) / sum(abs2, ds)
    return sqrt(-1 / (2 * slope))
end

@testset "background/spatial_corr" begin
    for (fs, px, n) in ((0.4, 0.078, 128), (1.6, 0.1, 192))
        bg = BackgroundModel(level=10.0, contrast=0.3, feature_size=fs, correlation_time=1e-9)
        S, _, L = gen_background(StableRNG(4), cam(n, px), bg, 40; frame_time=T)
        z = [log(S[i, j, k] / L[k]) for i in 1:n, j in 1:n, k in 1:40]
        @test isapprox(fit_sigma(z, px), fs; rtol=0.15)
    end
end

@testset "background/temporal_corr" begin
    tau = 0.05
    c = 0.4
    bg = BackgroundModel(level=10.0, contrast=c, feature_size=0.3, correlation_time=tau)
    n = 800
    S, _, L = gen_background(StableRNG(5), cam(48, 0.1), bg, n; frame_time=T)
    z = log_pattern(S, L, c)
    v = mean(abs2, z)
    for k in (1, 4, 16)
        ρ = mean(z[:, :, 1:n-k] .* z[:, :, 1+k:n]) / v
        @test abs(ρ - exp(-k * T / tau)) < 0.03
    end
end

@testset "background/stretch_jitter" begin
    # constant level within a stretch; levels follow the distribution
    lu = LogUniform(2.0, 100.0)
    bg = BackgroundModel(level=lu, stretch=10T)
    n = 5000
    _, _, L = gen_background(StableRNG(6), cam(4, 0.1), bg, n; frame_time=T)
    @test all(all(isapprox.(L[10j+1:10j+10], L[10j+1]; rtol=1e-12)) for j in 0:n÷10-1)
    @test all(L[10j+11] != L[10j+1] for j in 0:n÷10-2)
    lev = sort(L[1:10:n] ./ T)
    m = length(lev)
    D = maximum(max(abs(i / m - cdf(lu, lev[i])), abs((i - 1) / m - cdf(lu, lev[i]))) for i in 1:m)
    @test D < 1.628 / sqrt(m)                # KS, p > 0.01
    # the sd of the jitter multiplier
    bg = BackgroundModel(level=10.0, jitter=0.1)
    _, _, L = gen_background(StableRNG(7), cam(4, 0.1), bg, 4000; frame_time=T)
    @test isapprox(std(L ./ (10.0 * T)), 0.1; rtol=0.1)
end

@testset "background/illumination" begin
    px, w = 0.1, 2.0
    bg = BackgroundModel(level=1.0 / T, illumination_width=w)
    S, _, L = gen_background(StableRNG(8), cam(64, px), bg, 3; frame_time=T)
    @test isapprox(L[1], 1.0; rtol=1e-12)
    @test isapprox(mean(S[:, :, 1]), 1.0; rtol=1e-12)
    logS = log.(S[:, :, 1])
    d2 = logS[:, 1:end-2] .- 2 .* logS[:, 2:end-1] .+ logS[:, 3:end]
    @test all(isapprox.(d2, -(px / w)^2; rtol=1e-6))
    # separable: rank one
    @test isapprox(svdvals(S[:, :, 1])[2], 0.0; atol=1e-12)
end

# uniform-density static OOF population, one frame
function oof_pop(; density, lifetime=Inf, budget=Inf, birth_rate=nothing, γ=5000.0, σ=0.39, fluor=one_state(γ))
    kw = birth_rate === nothing ? (;) : (; birth_rate)
    return Population(; name=:oof, layer=:oof, density, lifetime, budget, fluor, z=(0.5, 1.0),
                      psf=GaussianPSF(σ), kw...)
end

@testset "background/oof_campbell" begin
    ρ, γ, σ, px = 0.3, 5000.0, 0.39, 0.1
    a = px^2
    nw = 40
    means, vars = Float64[], Float64[]
    for s in 1:nw
        _, O, _ = gen_background(StableRNG(100 + s), cam(128, px), nothing, 1;
                                 oof=[oof_pop(density=ρ)], frame_time=T)
        push!(means, mean(O))
        push!(vars, var(O))
    end
    A = (12.8 + 2 * 5 * σ)^2
    @test isapprox(mean(means), ρ * γ * T * a; rtol=3 / sqrt(ρ * A * nw))
    @test isapprox(mean(vars), ρ * (γ * T)^2 * a^2 / (4π * σ^2); rtol=0.1)
    # lifetime: the lag-k correlation decays as exp(-kT/tau)
    tau = 0.2
    n = 400
    _, O, _ = gen_background(StableRNG(9), cam(128, 0.2), nothing, n;
                             oof=[oof_pop(density=0.3, lifetime=tau)], frame_time=T)
    z = O .- mean(O; dims=(1, 2))
    v = mean(abs2, z)
    for k in (1, 4, 16)
        ρk = mean(z[:, :, 1:n-k] .* z[:, :, 1+k:n]) / v
        @test abs(ρk - exp(-k * T / tau)) < 0.05
    end
end

@testset "background/level_sweep" begin
    scam = SCMOSCamera(64, 64, 0.1, 1.6; offset=100.0, gain=2.0, qe=1.0)
    K = 50
    npx = 64 * 64
    pop = oof_pop(density=0.05, fluor=two_state(3000.0, 200.0, 50.0), lifetime=0.5, budget=500.0)
    for (i, lvl) in enumerate((3.2, 8.1, 16.0, 32.0, 81.0))
        bg = BackgroundModel(level=lvl / T)
        S, O, L = gen_background(StableRNG(200 + i), scam, bg, K; oof=[pop], frame_time=T)
        @test all(isapprox(mean(S[:, :, k]), L[k]; rtol=1e-12) for k in 1:K)
        @test all(isapprox.(L, lvl; rtol=1e-12))
        nrng = StableRNG(300 + i)
        d = 0.0
        λ = 0.0
        for k in 1:K
            lam = S[:, :, k] .+ O[:, :, k]
            img = copy(lam)
            SMLMSim.scmos_noise!(nrng, img, scam)
            d += (mean(img) - 100.0) * 2.0 - mean(lam)
            λ += mean(lam)
        end
        d /= K
        λ /= K
        @test abs(d) < 3 * sqrt((λ + 1.6^2) / (npx * K))
    end
end

@testset "background/threads_reproducible" begin
    pops = [oof_pop(density=0.2, lifetime=0.1, budget=300.0, fluor=two_state(4000.0, 100.0, 20.0))]
    bg = BackgroundModel(level=30.0 / T, jitter=0.05, contrast=0.3, correlation_time=0.05, illumination_width=3.0)
    master = StableRNG(11)
    seeds = [rand(master, UInt32) for _ in 1:4]
    run(seed) = gen_background(StableRNG(seed), cam(32, 0.1), bg, 30; oof=pops, frame_time=T, n_sub=2, t_burn=0.1)
    serial = [run(s) for s in seeds]
    tasks = [Threads.@spawn run(s) for s in seeds]
    @test fetch.(tasks) == serial
end

@testset "background/t_burn" begin
    # start-up transient: 10 emitters/um^2 at t = 0, but births hold only 0.4/um^2 (birth_rate x residence)
    pop = oof_pop(density=10.0, birth_rate=20.0, budget=100.0)     # residence 100/5000 = 0.02 s
    trans(t_burn, s) = begin
        _, O, _ = gen_background(StableRNG(s), cam(32, 0.1), nothing, 221; oof=[pop], frame_time=T, t_burn)
        return mean(O[:, :, 1:20]), mean(O[:, :, 200:220])
    end
    nw = 30
    burned = [trans(0.2, 500 + s) for s in 1:nw]            # 10 residence times
    d = first.(burned) .- last.(burned)
    @test abs(mean(d)) < 3 * std(d) / sqrt(nw)
    raw = [trans(0.0, 600 + s) for s in 1:nw]
    @test mean(first.(raw)) > 1.5 * mean(last.(raw))
    @test_throws ArgumentError gen_background(StableRNG(1), cam(8, 0.1), nothing, 2; frame_time=T, t_burn=-1.0)
end

@testset "background/models" begin
    pop = oof_pop(density=0.2, lifetime=0.1, budget=300.0, fluor=two_state(4000.0, 100.0, 20.0))
    bg = BackgroundModel(level=20.0 / T, contrast=0.2)
    S, O, L = gen_background(StableRNG(12), cam(32, 0.1), bg, 5; frame_time=T)
    @test all(iszero, O) && sum(S) > 0
    S, O, L = gen_background(StableRNG(12), cam(32, 0.1), nothing, 5; oof=[pop], frame_time=T)
    @test all(iszero, S) && all(iszero, L) && sum(O) > 0
    S, O, L = gen_background(StableRNG(12), cam(32, 0.1), bg, 5; oof=[pop], frame_time=T)
    @test sum(S) > 0 && sum(O) > 0
    sig = Population(density=1.0, fluor=one_state(1000.0), psf=GaussianPSF(0.13))
    @test_throws ArgumentError gen_background(StableRNG(1), cam(8, 0.1), nothing, 2; oof=[sig], frame_time=T)
    @test_throws ArgumentError BackgroundModel(level=-1.0) |> b -> SimWorld(StableRNG(1), cam(8, 0.1), Population[]; n_sub=1, background=b)
    # exposure shorter than the frame: level scales with the exposure
    _, _, L = gen_background(StableRNG(13), cam(8, 0.1), BackgroundModel(level=10.0), 3; frame_time=T, exposure=T / 2)
    @test all(isapprox.(L, 10.0 * T / 2; rtol=1e-12))
end
