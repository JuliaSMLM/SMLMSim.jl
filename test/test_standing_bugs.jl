using SMLMSim, Test, Random, Statistics, Distributions, MicroscopePSFs

fit2(x, y; kwargs...) = Emitter2DFit{Float64}(x, y, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; kwargs...)
fit3(x, y, z; kwargs...) = Emitter3DFit{Float64}(x, y, z, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; kwargs...)
const BLINKING = GenericFluor(photons=1e4, k_off=20.0, k_on=20.0)

# Dimer analysis needs only these fields, without timestamps or partner information.
struct StateRecord <: AbstractEmitter
    frame::Int
    track_id::Int
    state::Symbol
end

@testset "Must 1: absorbing states finish the final frame" begin
    # p02: 3 * (1/10) rounds past 3/10; 4 frames is the terminating control.
    fluor = GenericFluor(1e4, [-60.0 10.0 50.0; 1.0 -1.0 0.0; 0.0 0.0 0.0])
    for nframes in (4, 3)
        Random.seed!(2)
        photons = intensity_trace(fluor, nframes, 10.0; state1=1)
        @test length(photons) == nframes
        @test photons[1] > 0
        @test all(iszero, photons[2:end])  # bleached after its first frame
        @test all(iszero, intensity_trace(fluor, nframes, 10.0; state1=3))
    end

    # p04: the same rounding hang through the static simulation at 50 Hz.
    Random.seed!(2)
    fluor = GenericFluor(1e4, [-15.0 10.0 5.0; 1.0 -1.0 0.0; 0.0 0.0 0.0])
    smld, info = simulate(StaticSMLMConfig(nframes=35, framerate=50.0, ndatasets=1, density=1.0);
                          molecule=fluor, state1=1, camera=IdealCamera(1:16, 1:16, 0.1))
    @test smld.n_frames == 35
    @test !isempty(info.smld_model.emitters)
    @test all(e -> isfinite(e.photons) && e.photons > 0, info.smld_model.emitters)
end

@testset "Must 2: images and totals select one dataset" begin
    Random.seed!(1)
    camera = IdealCamera(1:32, 1:32, 0.1)
    e1, e2 = fit2(0.8, 0.8; frame=1, dataset=1), fit2(2.4, 2.4; frame=1, dataset=2)
    e1.photons = e2.photons = 1000.0
    smld = BasicSMLD([e1, e2], camera, 1, 2)
    psf = GaussianPSF(0.13)
    for dataset in (1, 2)
        images, info = gen_images(smld, psf; dataset=dataset, threaded=false)
        own, other = dataset == 1 ? (8, 24) : (24, 8)
        @test sum(images) ≈ 1000.0 rtol=1e-6
        @test images[own, own, 1] > 50.0
        @test images[other, other, 1] ≈ 0.0 atol=1e-40
        @test info.n_photons_total == 1000.0
    end
    images, info = gen_images(smld, psf; dataset=3, threaded=false)
    @test all(iszero, images)
    @test info.n_photons_total == 0.0
end

@testset "Must 3: dimer fraction averages molecules over exposure" begin
    Random.seed!(10)
    params = DiffusionSMLMConfig(density=2.0, box_size=5.0, r_react=0.05, k_off=1.0,
                               t_max=3.0, dt=0.01, camera_framerate=10.0, camera_exposure=0.1)
    smld, _ = simulate(params)
    frames, fractions = analyze_dimer_fraction(smld)
    @test frames == collect(1:30)
    @test any(>(0), fractions)
    @test all(f -> 0 <= f <= 1, fractions)
    for (frame, fraction) in zip(frames, fractions)
        records = filter(e -> e.frame == frame, smld.emitters)
        @test fraction ≈ count(e -> e.state == :dimer, records) / length(records)
    end
    Random.seed!(7403)
    shuffled = BasicSMLD(shuffle(smld.emitters), smld.camera, smld.n_frames, smld.n_datasets, copy(smld.metadata))
    @test isequal(analyze_dimer_fraction(shuffled), (frames, fractions))

    camera = IdealCamera(1:32, 1:32, 0.1)
    # Eight molecules, four bound: one record per molecule must give 0.5.
    records = [DiffusingEmitter2D{Float64}(0.1i, 1.0, 1000.0, 0.01, 1, 1, i,
               i <= 4 ? :dimer : :monomer, i <= 4 ? (isodd(i) ? i + 1 : i - 1) : nothing) for i in 1:8]
    @test analyze_dimer_fraction(BasicSMLD(records, camera, 1, 1)) == ([1], [0.5])

    # Partial exposure occupancy, with unequal record counts: weight molecules equally.
    partial = [DiffusingEmitter2D{Float64}(0.1i, 1.0, 1000.0, 0.01k, 1, 1, i,
               k == 1 ? :dimer : :monomer, k == 1 ? 3 - i : nothing)
               for i in 1:2 for k in 1:(2i)]
    @test analyze_dimer_fraction(BasicSMLD(partial, camera, 1, 1)) == ([1], [0.375])
    empty_smld = BasicSMLD(DiffusingEmitter2D{Float64}[], camera, 1, 1)
    @test analyze_dimer_fraction(empty_smld) == (Int[], Float64[])

    minimal = [StateRecord(e.frame, e.track_id, e.state) for e in partial]
    append!(minimal, [StateRecord(3, 1, :dimer), StateRecord(3, 1, :dimer),
                      StateRecord(3, 2, :monomer), StateRecord(3, 2, :dimer), StateRecord(3, 2, :monomer)])
    expected = ([1, 3], [0.375, 2 / 3])
    @test isequal(analyze_dimer_fraction(BasicSMLD(minimal, camera, 3, 1)), expected)
    Random.seed!(7403)
    @test isequal(analyze_dimer_fraction(BasicSMLD(shuffle(minimal), camera, 3, 1)), expected)
    @test analyze_dimer_fraction(BasicSMLD(StateRecord[], camera, 3, 1)) == (Int[], Float64[])
end

@testset "Must 4: conflicting track positions warn without changing output" begin
    Random.seed!(3)
    camera = IdealCamera(1:64, 1:64, 0.1)
    params = StaticSMLMConfig(nframes=200, framerate=50.0, ndatasets=1)
    emitters = [fit2(x, 3.0) for x in range(0.5, 5.5; length=10)]
    actual, info = @test_logs (:warn, r"repeated track_id 0.*only the first record's position is simulated.*each molecule needs its own track_id") simulate(
        params; starting_conditions=emitters, camera=camera, molecule=BLINKING)
    next_draw = rand()
    @test info.n_emitters == 10
    @test unique(e.x for e in info.smld_model.emitters) == [first(emitters).x]
    @test unique(e.track_id for e in info.smld_model.emitters) == [0]
    @test all(e -> e.track_id == 0, emitters)  # no mutation of starting conditions

    # The original grouping simulates just the first of these zero-id records, bitwise.
    Random.seed!(3)
    expected, expected_info = @test_logs simulate(params; starting_conditions=emitters[1:1], camera=camera, molecule=BLINKING)
    fields(smld) = [Tuple(getfield(e, f) for f in fieldnames(typeof(e))) for e in smld.emitters]
    @test isequal(fields(actual), fields(expected))
    @test isequal(fields(info.smld_model), fields(expected_info.smld_model))
    @test rand() === next_draw

    for dimension in (2, 3)
        makefit(x, id) = dimension == 2 ? fit2(x, 1.0; track_id=id) : fit3(x, 1.0, x; track_id=id)
        inputs = [makefit(1.0, 7), makefit(0.5, 0), makefit(1.5, 2)]
        source = BasicSMLD(inputs, camera, 1, 1)
        Random.seed!(74)
        model = @test_logs kinetic_model(source, BLINKING, 40, 50.0; ndatasets=2)
        @test unique(e.track_id for e in model.emitters) == [0, 2, 7]
        @test unique((e.track_id, e.x) for e in model.emitters) == [(0, 0.5), (2, 1.5), (7, 1.0)]
        @test all(d -> Set(e.track_id for e in model.emitters if e.dataset == d) == Set([0, 2, 7]), 1:2)

        for id in (0, 7)
            fixed = BasicSMLD([makefit(0.5, id), makefit(0.5, id)], camera, 1, 1)
            @test_logs kinetic_model(fixed, BLINKING, 40, 50.0)
        end

        repeated = BasicSMLD(vcat(inputs, [makefit(0.5, 0), makefit(1.0, 7)]), camera, 1, 1)
        Random.seed!(74)
        same = @test_logs kinetic_model(repeated, BLINKING, 40, 50.0; ndatasets=2)
        @test isequal(fields(same), fields(model))

        conflict = dimension == 2 ? makefit(2.5, 7) : fit3(1.0, 1.0, 2.5; track_id=7)
        conflicted = BasicSMLD(vcat(inputs, [conflict, conflict]), camera, 1, 1)
        # Both repeated inputs and repeated calls at this logger emit only one warning.
        Random.seed!(74)
        models = @test_logs (:warn, r"repeated track_id 7.*only the first record's position is simulated.*each molecule needs its own track_id") respect_maxlog=true begin
            [kinetic_model(conflicted, BLINKING, 40, 50.0; ndatasets=2) for _ in 1:2]
        end
        @test isequal(fields(models[1]), fields(model))
    end
end

@testset "Must 5: line copies draw independent sites" begin
    # p14, without edge clipping: keep every copy, including zero-site copies.
    for dimension in (2, 3)
        Random.seed!(9)
        line = dimension == 2 ? Line2D(λ=10.0, endpoints=[(-0.5, 0.0), (0.5, 0.0)]) :
                               Line3D(λ=10.0, endpoints=[(-0.5, 0.0, 0.0), (0.5, 0.0, 0.0)])
        Random.seed!(9)
        n_copies = rand(Poisson(0.8 * 25.6 * 25.6))
        Random.seed!(9)
        output = dimension == 2 ? uniform2D(0.8, line, 25.6, 25.6) : uniform3D(0.8, line, 25.6, 25.6)
        coords, ids = output[1:end-1], output[end]
        counts = zeros(Int, n_copies)
        for id in ids
            counts[id] += 1
        end
        μ, N = 10.0, length(counts)
        @test abs(mean(counts) - μ) <= 4 * sqrt(μ / N)
        # Var(sample variance) for N independent Poisson(μ) draws.
        @test abs(var(counts) - μ) <= 4 * sqrt(μ / N + 2μ^2 / (N - 1))
        @test length(unique(counts)) > 1
        candidates = findall(>=(2), counts)
        function offsets(id)
            idx = findall(==(id), ids)
            sort(round.([sqrt(sum((c[i] - c[first(idx)])^2 for c in coords)) for i in idx]; digits=10))
        end
        first_offsets, second_offsets = offsets(candidates[1]), offsets(candidates[2])
        common_length = min(length(first_offsets), length(second_offsets))
        @test common_length >= 2
        @test first_offsets[1:common_length] != second_offsets[1:common_length]
        @test SMLMSim.StaticSMLM._pattern_radius(line) == 0.5

        # The radius also covers lines whose constructor happened to draw no sites.
        empty_line = dimension == 2 ? Line2D(λ=0.0, endpoints=[(-0.5, 0.0), (0.5, 0.0)]) :
                                     Line3D(λ=0.0, endpoints=[(-0.5, 0.0, 0.0), (0.5, 0.0, 0.0)])
        @test SMLMSim.StaticSMLM._pattern_radius(empty_line) == 0.5
        Random.seed!(9)
        _, info = simulate(StaticSMLMConfig(density=0.05, nframes=1, ndatasets=1, ndims=dimension);
                           pattern=line, camera=IdealCamera(1:256, 1:256, 0.1))
        @test all(e -> 0 <= e.x <= 25.6 && 0 <= e.y <= 25.6, info.smld_true.emitters)
        @test length(unique(count(e -> e.id == id, info.smld_true.emitters)
                            for id in unique(e.id for e in info.smld_true.emitters))) > 1
    end

    # Direct pre-0.7.4 computation: no stored coordinates or RNG-specific fixtures.
    pattern, ρ, fx, fy = Nmer2D(n=6, d=0.2), 0.5, 10.0, 8.0
    Random.seed!(7405)
    x, y, ids = Float64[], Float64[], Int[]
    for id in 1:rand(Poisson(fx * fy * ρ))
        θ, x0, y0 = 2 * pi * rand(), rand() * fx, rand() * fy
        for i in 1:pattern.n
            push!(x, pattern.x[i] * cos(θ) - pattern.y[i] * sin(θ) + x0)
            push!(y, pattern.x[i] * sin(θ) + pattern.y[i] * cos(θ) + y0)
            push!(ids, id)
        end
    end
    next_draw = rand()
    Random.seed!(7405)
    @test isequal(uniform2D(ρ, pattern, fx, fy), (x, y, ids))
    @test rand() === next_draw
end

@testset "Should 1: sCMOS draws shot noise once with both flags" begin
    camera = SCMOSCamera(64, 64, 0.1, 0.0; offset=0.0, gain=1.0, qe=1.0)
    smld = BasicSMLD([fit2(10.0, 10.0; frame=2)], camera, 2, 1)
    psf = GaussianPSF(0.13)
    for (poisson_noise, camera_noise) in ((true, false), (false, true), (true, true))
        Random.seed!(4)
        images, _ = gen_images(smld, psf; frames=[1], bg=100.0, poisson_noise=poisson_noise, camera_noise=camera_noise)
        μ, N = 100.0, length(images)
        @test abs(mean(images) - μ) <= 4 * sqrt(μ / N)
        @test abs(var(images) - μ) <= 4 * sqrt(μ / N + 2μ^2 / (N - 1))
    end
    Random.seed!(4)
    both, _ = gen_images(smld, psf; frames=[1], bg=100.0, poisson_noise=true, camera_noise=true)
    next_draw = rand()
    Random.seed!(4)
    camera_only, _ = gen_images(smld, psf; frames=[1], bg=100.0, camera_noise=true)
    @test isequal(both, camera_only)
    @test rand() === next_draw
end

@testset "Should 4: diffusion density documents area and volume units" begin
    # This item changes documentation only: the probe's counts already obey these units.
    @test occursin(r"density.*molecules/um² in 2D.*molecules/um³ in 3D", string(@doc DiffusionSMLMConfig))
    for (dimension, expected) in ((2, 25), (3, 125))
        Random.seed!(12)
        _, info = simulate(DiffusionSMLMConfig(density=1.0, box_size=5.0, ndims=dimension, t_max=0.1))
        @test info.n_emitters == expected
    end
end
