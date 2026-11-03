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

struct CameraFreeSMLD{E<:AbstractEmitter} <: SMLMSim.SMLD
    emitters::Vector{E}
end

struct OptionalCameraSMLD{E<:AbstractEmitter} <: SMLMSim.SMLD
    emitters::Vector{E}
    camera::Union{Nothing,AbstractCamera}
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

@testset "Must 4: unlinked static emitters are independent molecules" begin
    Random.seed!(3)
    camera = IdealCamera(1:64, 1:64, 0.1)
    params = StaticSMLMConfig(nframes=200, framerate=50.0, ndatasets=1)
    emitters = [fit2(x, 3.0) for x in range(0.5, 5.5; length=10)]
    _, info = simulate(params; starting_conditions=emitters, camera=camera, molecule=BLINKING)
    @test info.n_emitters == 10
    @test length(unique(e.x for e in info.smld_model.emitters)) == 10
    @test sort(unique(e.track_id for e in info.smld_model.emitters)) == collect(1:10)
    @test all(e -> e.track_id == 0, emitters)  # no mutation of starting conditions

    # Linked ids remain sorted and keep their first record, followed by unlinked inputs.
    for dimension in (2, 3)
        makefit(x, id) = dimension == 2 ? fit2(x, 1.0; track_id=id) : fit3(x, 1.0, x; track_id=id)
        unique_source = BasicSMLD([makefit(1.0, 7), makefit(1.5, 2)], camera, 1, 1)
        @test_logs kinetic_model(unique_source, BLINKING, 40, 50.0)
        inputs = [makefit(0.5, 0), makefit(1.0, 7), makefit(1.5, 2), makefit(2.0, 0), makefit(1.0, 7)]
        source = BasicSMLD(inputs, camera, 1, 1)
        Random.seed!(74)
        model = @test_logs kinetic_model(source, BLINKING, 40, 50.0; ndatasets=2)
        @test unique(e.track_id for e in model.emitters) == [2, 7, 8, 9]
        @test unique((e.track_id, e.x) for e in model.emitters) == [(2, 1.5), (7, 1.0), (8, 0.5), (9, 2.0)]
        @test all(d -> Set(e.track_id for e in model.emitters if e.dataset == d) == Set([2, 7, 8, 9]), 1:2)

        conflict = dimension == 2 ? makefit(2.5, 7) : fit3(1.0, 1.0, 2.5; track_id=7)
        conflicted = BasicSMLD(vcat(inputs, [conflict, conflict]), camera, 1, 1)
        # Both repeated inputs and repeated calls at this logger emit only one warning.
        Random.seed!(74)
        models = @test_logs (:warn, r"repeated non-zero track_id 7.*different positions.*first record") respect_maxlog=true begin
            [kinetic_model(conflicted, BLINKING, 40, 50.0; ndatasets=2) for _ in 1:2]
        end
        @test [e.x for e in models[1].emitters] == [e.x for e in model.emitters]
        @test [e.photons for e in models[1].emitters] == [e.photons for e in model.emitters]
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

@testset "Should 2: diffusion warns once about ignored keywords" begin
    params = DiffusionSMLMConfig(density=0.5, box_size=2.0, t_max=0.2)
    Random.seed!(7)
    baseline, _ = @test_logs simulate(params)
    next_draw = rand()
    Random.seed!(7)
    results = @test_logs (:warn, r"ignored keywords: gamma, Gamma, photon, camera_exposure.*deprecated.*0.8.0") respect_maxlog=true begin
        first = simulate(params; gamma=2e4, Gamma=2e4, photon=10.0, camera_exposure=0.05)[1]
        next = rand()
        simulate(params; gamma=2e4)
        (first, next)
    end
    @test all(e -> e.photons == 1000.0, results[1].emitters)
    @test [(e.x, e.y, e.photons) for e in results[1].emitters] == [(e.x, e.y, e.photons) for e in baseline.emitters]
    @test results[2] === next_draw
    bright, _ = @test_logs simulate(params; γ=2e4)
    @test all(e -> e.photons == 200.0, bright.emitters)
    @test_logs (:warn, r"ignored keywords: gamma;.*0.8.0") simulate(params; gamma=2e4)
end

@testset "Should 3: static starting SMLD supplies its camera" begin
    Random.seed!(6)
    camera = SCMOSCamera(64, 64, 0.1, 1.6; offset=100.0, gain=0.5, qe=0.8)
    emitters = [fit2(1.0 + 0.1i, 3.0; track_id=i) for i in 1:5]
    source = BasicSMLD(emitters, camera, 1, 1)
    params = StaticSMLMConfig(nframes=50, ndatasets=1)
    for kw in ((;), (; camera=nothing))
        smld, info = simulate(params; starting_conditions=source, molecule=BLINKING, kw...)
        @test smld.camera === camera
        @test info.smld_true.camera === camera
        @test info.smld_model.camera === camera
    end
    override = IdealCamera(1:32, 1:32, 0.1)
    smld, info = simulate(params; starting_conditions=source, molecule=BLINKING, camera=override)
    @test smld.camera === override
    @test info.smld_true.camera === override
    smld, _ = simulate(params; starting_conditions=emitters, molecule=BLINKING)
    @test smld.camera isa IdealCamera
    @test length(smld.camera.pixel_edges_x) - 1 == 128

    # No starting conditions: implicit and explicit default cameras consume the same draws.
    Random.seed!(7403)
    implicit, _ = simulate(params; molecule=BLINKING)
    @test length(implicit.camera.pixel_edges_x) - 1 == 128
    @test length(implicit.camera.pixel_edges_y) - 1 == 128
    next_draw = rand()
    Random.seed!(7403)
    explicit, _ = simulate(params; molecule=BLINKING, camera=IdealCamera(1:128, 1:128, 0.1))
    @test [(e.x, e.y, e.photons, e.track_id) for e in implicit.emitters] ==
          [(e.x, e.y, e.photons, e.track_id) for e in explicit.emitters]
    @test rand() === next_draw
end

@testset "Lesser 1: Float32 rate matrices simulate with Float64 time" begin
    Random.seed!(7401)
    fluor = GenericFluor(1f4, Float32[-50 50; 1f-2 -1f-2])
    smld, _ = simulate(StaticSMLMConfig(density=2.0, nframes=20, ndatasets=1, minphotons=0);
                       molecule=fluor, state1=1, camera=IdealCamera(16, 16, 0.1))
    @test !isempty(smld.emitters)
    @test all(e -> isfinite(e.photons) && e.photons > 0, smld.emitters)

    # Promotion gives the same trace and random stream as an explicit Float64 matrix.
    Random.seed!(7401)
    trace = intensity_trace(fluor, 20, 50.0; state1=1)
    next_draw = rand()
    Random.seed!(7401)
    reference = GenericFluor(Float64(fluor.γ), Float64.(fluor.q))
    @test isequal(trace, intensity_trace(reference, 20, 50.0; state1=1))
    @test rand() === next_draw

    # Other floating point matrix types also work with the usual Float64 time arguments.
    for T in (Float16, BigFloat)
        Random.seed!(7401)
        trace = intensity_trace(GenericFluor(T(1e4), T[-50 50; 0.01 -0.01]), 20, 50.0; state1=1)
        @test length(trace) == 20
        @test all(isfinite, trace)
        @test sum(trace) > 0
    end
end

@testset "Lesser 2: GenericFluor accepts integer keyword values" begin
    integer = GenericFluor(; photons=10000, k_off=50, k_on=1)
    reference = GenericFluor(; photons=10000.0, k_off=50.0, k_on=1.0)
    @test integer.γ === reference.γ
    @test isequal(integer.q, reference.q)
    @test eltype(integer.q) === Float64
    for fluor in (GenericFluor(; photons=10000), GenericFluor(; k_off=50), GenericFluor(; k_on=1))
        @test fluor.γ isa Float64
        @test eltype(fluor.q) === Float64
    end
    rational = GenericFluor(; photons=10000//1, k_off=100//2, k_on=1//1)
    @test rational.γ === reference.γ
    @test isequal(rational.q, reference.q)
    Random.seed!(7402)
    trace = intensity_trace(integer, 20, 50.0; state1=1)
    Random.seed!(7402)
    @test isequal(trace, intensity_trace(reference, 20, 50.0; state1=1))
end

@testset "Lesser 3: saturated dissociation probability warns once" begin
    params = DiffusionSMLMConfig(density=0.5, box_size=2.0, t_max=0.2, k_off=200.0, dt=0.01)
    Random.seed!(7403)
    @test_logs (:warn, r"k_off=200.0.*dt=0.01.*> 1.*mean lifetime is dt=0.01.*1/k_off=0.005") respect_maxlog=true begin
        simulate(params)
        simulate(params)
    end
    @test_logs simulate(DiffusionSMLMConfig(density=0.5, box_size=2.0, t_max=0.2))
    @test_logs simulate(DiffusionSMLMConfig(density=0.5, box_size=2.0, t_max=0.2, k_off=100.0, dt=0.01))
end

@testset "Lesser 4: lifetime averages only intervals with an observed break" begin
    Random.seed!(7408)
    params = DiffusionSMLMConfig(density=5.0, box_size=5.0, r_react=0.3, k_off=0.5,
                               dt=0.01, t_max=0.6, camera_framerate=100.0, camera_exposure=0.01)
    smld, _ = simulate(params)
    # One recorded step per frame: independently count observed bound-to-unbound intervals.
    completed = Float64[]
    unfinished = 0
    for id in unique(e.track_id for e in smld.emitters)
        records = sort(filter(e -> e.track_id == id, smld.emitters), by=e -> e.frame)
        first_bound = nothing
        for e in records
            if e.state == :dimer
                first_bound === nothing && (first_bound = e.frame)
            elseif first_bound !== nothing
                push!(completed, (e.frame - 1) / params.camera_framerate - (first_bound - 1) / params.camera_framerate)
                first_bound = nothing
            end
        end
        unfinished += first_bound !== nothing
    end
    @test !isempty(completed)
    @test unfinished > length(completed)
    lifetime = analyze_dimer_lifetime(smld)
    @test lifetime ≈ mean(completed)
    @test lifetime > sum(completed) / (length(completed) + unfinished)

    # Real simulator output with bound dimers but no breaks, and with no dimers at all.
    Random.seed!(7408)
    bound, _ = simulate(DiffusionSMLMConfig(density=5.0, box_size=5.0, r_react=0.3, k_off=0.0,
                                           dt=0.01, t_max=0.2, camera_framerate=100.0, camera_exposure=0.01))
    @test any(e -> e.state == :dimer, bound.emitters)
    @test analyze_dimer_lifetime(bound) == 0.0
    monomer, _ = simulate(DiffusionSMLMConfig(box_size=2.0, t_max=0.1); override_count=1)
    @test analyze_dimer_lifetime(monomer) == 0.0
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

@testset "Should 5: plain static emitters fail on entry with a useful error" begin
    Random.seed!(3)
    camera = IdealCamera(1:64, 1:64, 0.1)
    params = StaticSMLMConfig(nframes=200, ndatasets=1)
    # Abstractly typed containers of valid fits are accepted, as through 0.7.3:
    # deepcopy.(starting_conditions) narrows them to the concrete emitter type.
    concrete = [fit2(1.0, 3.0), fit2(2.0, 3.0)]
    abstractly_typed = AbstractEmitter[concrete...]
    Random.seed!(7405)
    expected, expected_info = simulate(params; starting_conditions=concrete, camera=camera, molecule=BLINKING)
    next_draw = rand()
    Random.seed!(7405)
    actual, actual_info = simulate(params; starting_conditions=abstractly_typed, camera=camera, molecule=BLINKING)
    fields(smld) = [Tuple(getfield(e, f) for f in fieldnames(typeof(e))) for e in smld.emitters]
    @test actual_info.n_emitters == expected_info.n_emitters == 2
    @test isequal(fields(actual), fields(expected))
    @test isequal(fields(actual_info.smld_true), fields(expected_info.smld_true))
    @test isequal(fields(actual_info.smld_model), fields(expected_info.smld_model))
    @test rand() === next_draw

    for emitters in ([Emitter2D{Float64}(1.0, 3.0, 1.0)], [Emitter3D{Float64}(1.0, 3.0, 0.5, 1.0)])
        mixed = AbstractEmitter[fit2(1.0, 3.0), only(emitters)]
        for source in (emitters, BasicSMLD(emitters, camera, 1, 1), mixed, BasicSMLD(mixed, camera, 1, 1))
            err = try
                simulate(params; starting_conditions=source, camera=camera, molecule=BLINKING)
                nothing
            catch e
                e
            end
            @test err isa ArgumentError
            @test occursin(string(eltype(emitters)), sprint(showerror, err))
            @test occursin("Emitter2DFit or Emitter3DFit is needed", sprint(showerror, err))
        end
    end
end

@testset "Round 5: camera-free SMLD uses the default camera" begin
    emitters = [fit2(1.0, 2.0; track_id=1), fit2(2.0, 3.0; track_id=2)]
    params = StaticSMLMConfig(nframes=20, ndatasets=1, minphotons=0)
    for source in (CameraFreeSMLD(emitters), OptionalCameraSMLD(emitters, nothing))
        Random.seed!(7501)
        smld, info = simulate(params; starting_conditions=source, molecule=BLINKING, state1=1)
        @test smld.camera isa IdealCamera
        @test length(smld.camera.pixel_edges_x) - 1 == 128
        @test length(smld.camera.pixel_edges_y) - 1 == 128
        @test info.n_emitters == 2
        @test smld.camera === info.smld_true.camera === info.smld_model.camera
        @test Set(e.track_id for e in info.smld_model.emitters) == Set([1, 2])
        override = IdealCamera(16, 32, 0.1)
        smld, info = simulate(params; starting_conditions=source, camera=override, molecule=BLINKING, state1=1)
        @test smld.camera === override
        @test info.smld_true.camera === override
    end
end

@testset "Round 5: floating point keywords retain their precision" begin
    precise = big"10000.0000000000007"
    fluor = GenericFluor(; photons=precise, k_off=0.0, k_on=1.0)
    @test fluor.γ === precise
    @test isequal(intensity_trace(fluor, 1, 10), [1000.0000000000001])
    for T in (Float16, Float32, Float64, BigFloat)
        photons, k_off, k_on = T(10000), T(50), T(1)
        fluor = GenericFluor(; photons=photons, k_off=k_off, k_on=k_on)
        @test fluor.γ === photons
        @test eltype(fluor.q) === T
        @test isequal(fluor.q, [-k_off k_off; k_on -k_on])
    end
    integer = GenericFluor(; photons=10000, k_off=50, k_on=1)
    reference = GenericFluor(; photons=10000.0, k_off=50.0, k_on=1.0)
    @test integer.γ === reference.γ
    @test isequal(integer.q, reference.q)
end

@testset "Round 5: unlinked molecules accept abstract and mixed precision containers" begin
    camera = IdealCamera(32, 32, 0.1)
    params = StaticSMLMConfig(nframes=20, ndatasets=2, minphotons=0)
    for dimension in (2, 3)
        concrete = dimension == 2 ? [fit2(0.5, 1.0), fit2(1.5, 2.0)] :
                                    [fit3(0.5, 1.0, 0.0), fit3(1.5, 2.0, 0.25)]
        abstractly_typed = AbstractEmitter[concrete...]
        Random.seed!(7503)
        model = kinetic_model(BasicSMLD(abstractly_typed, camera, 1, 1), BLINKING, 20, 50.0;
                              ndatasets=2, minphotons=0.0, state1=1)
        @test unique(e.track_id for e in model.emitters) == [1, 2]
        @test unique((e.x, e.y) for e in model.emitters) == [(0.5, 1.0), (1.5, 2.0)]
        @test eltype(model.emitters) === (dimension == 2 ? Emitter2DFit{Float64} : Emitter3DFit{Float64})

        single = dimension == 2 ?
            Emitter2DFit{Float32}(1.5f0, 2f0, 1f0, 0f0, 0f0, 0f0, 0f0, 0f0) :
            Emitter3DFit{Float32}(1.5f0, 2f0, 0.25f0, 1f0, 0f0, 0f0, 0f0, 0f0, 0f0, 0f0)
        family = dimension == 2 ? Emitter2DFit : Emitter3DFit
        mixed = family[concrete[1], single]
        for inputs in (mixed, reverse(mixed), AbstractEmitter[mixed...])
            Random.seed!(7503)
            smld, info = simulate(params; starting_conditions=inputs, camera=camera, molecule=BLINKING, state1=1)
            @test info.n_emitters == 2
            @test unique(e.track_id for e in info.smld_model.emitters) == [1, 2]
            @test unique((e.x, e.y) for e in info.smld_model.emitters) == [(Float64(e.x), Float64(e.y)) for e in inputs]
            @test eltype(smld.emitters) === (dimension == 2 ? Emitter2DFit{Float64} : Emitter3DFit{Float64})
            @test all(d -> Set(e.track_id for e in info.smld_model.emitters if e.dataset == d) == Set([1, 2]), 1:2)
        end
    end
    mixed_dimensions = BasicSMLD(AbstractEmitter[fit2(1.0, 2.0), fit3(1.0, 2.0, 0.0)], camera, 1, 1)
    @test_throws ArgumentError kinetic_model(mixed_dimensions, BLINKING, 20, 50.0; state1=1)
end

@testset "Round 6: unlinked id allocation rejects overflow" begin
    camera = IdealCamera(32, 32, 0.1)
    largest_id = typemax(Int)
    source = BasicSMLD([fit2(0.5, 1.0; track_id=largest_id), fit2(1.5, 2.0)], camera, 1, 1)
    Random.seed!(7602)
    expected_next = rand()
    Random.seed!(7602)
    err = try
        kinetic_model(source, BLINKING, 1, 50.0)
        nothing
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("largest non-zero id $largest_id", sprint(showerror, err))
    @test rand() == expected_next

    # The final available id is valid; only the addition of the unlinked offset counts.
    source = BasicSMLD([fit2(0.5, 1.0; track_id=largest_id - 1), fit2(1.5, 2.0)], camera, 1, 1)
    model = kinetic_model(source, BLINKING, 1, 50.0; state1=1, minphotons=0.0)
    @test [e.track_id for e in model.emitters] == [largest_id - 1, largest_id]
    linked_only = BasicSMLD([fit2(0.5, 1.0; track_id=largest_id)], camera, 1, 1)
    model = kinetic_model(linked_only, BLINKING, 1, 50.0; state1=1, minphotons=0.0)
    @test only(model.emitters).track_id == largest_id
end

@testset "Round 7: empty starting conditions are valid" begin
    camera = IdealCamera(32, 24, 0.1)
    inputs = ((AbstractEmitter[], nothing), (Emitter2DFit[], 2), (Emitter3DFit[], 3),
              (Emitter2DFit{Float64}[], 2), (Emitter3DFit{Float64}[], 3))
    for (emitters, declared_dimension) in inputs
        source = BasicSMLD(emitters, camera, 1, 1)
        model = kinetic_model(source, BLINKING, 5, 50.0; ndatasets=2)
        @test isempty(model.emitters)
        @test eltype(model.emitters) === eltype(emitters)
        @test model.emitters !== emitters
        @test model.camera === camera
        @test model.n_frames == 5
        @test model.n_datasets == 2

        for dimension in (2, 3)
            params = StaticSMLMConfig(nframes=5, ndatasets=2, ndims=dimension)
            expected_dimension = something(declared_dimension, dimension)
            expected_family = expected_dimension == 2 ? Emitter2DFit : Emitter3DFit
            for (starting_conditions, kwargs) in ((emitters, (; camera=camera)), (source, (;)))
                smld, info = simulate(params; starting_conditions=starting_conditions, kwargs...)
                @test isempty(smld.emitters)
                @test isempty(info.smld_true.emitters)
                @test isempty(info.smld_model.emitters)
                @test info.n_emitters == info.n_localizations == info.n_patterns == 0
                @test smld.camera === info.smld_true.camera === info.smld_model.camera === camera
                @test smld.n_frames == info.smld_model.n_frames == params.nframes
                @test smld.n_datasets == info.smld_model.n_datasets == params.ndatasets
                @test eltype(info.smld_true.emitters) === eltype(emitters)
                @test info.smld_true.emitters !== emitters
                @test eltype(smld.emitters) <: expected_family
                @test eltype(info.smld_model.emitters) <: expected_family
                expected_psf = expected_dimension == 2 ? params.σ_psf :
                               [params.σ_psf, params.σ_psf, 3 * params.σ_psf]
                @test smld.metadata["psf_width"] == expected_psf
            end
        end
    end

    smld, info = simulate(StaticSMLMConfig(nframes=5, ndatasets=2); starting_conditions=AbstractEmitter[])
    @test isempty(smld.emitters)
    @test smld.camera === info.smld_true.camera === info.smld_model.camera
    @test length(smld.camera.pixel_edges_x) - 1 == 128
    @test length(smld.camera.pixel_edges_y) - 1 == 128
end

@testset "Round 8: diffusion analysis ignores record order" begin
    Random.seed!(7408)
    params = DiffusionSMLMConfig(density=5.0, box_size=5.0, r_react=0.3, k_off=0.5,
                               dt=0.01, t_max=0.6, camera_framerate=10.0, camera_exposure=0.1)
    smld, _ = simulate(params)
    track_changes = SMLMSim.InteractionDiffusion.track_state_changes
    history = track_changes(smld)
    lifetime = analyze_dimer_lifetime(smld)
    fractions = analyze_dimer_fraction(smld)
    truth = frame_dimer_truth(smld)
    @test lifetime ≈ 0.23636363636363636

    Random.seed!(7808)
    for records in (reverse(smld.emitters), shuffle(smld.emitters), shuffle(smld.emitters))
        reordered = BasicSMLD(records, smld.camera, smld.n_frames, smld.n_datasets, copy(smld.metadata))
        original_order = copy(records)
        @test isequal(track_changes(reordered), history)
        @test isequal(analyze_dimer_lifetime(reordered), lifetime)
        @test isequal(analyze_dimer_fraction(reordered), fractions)
        @test isequal(frame_dimer_truth(reordered), truth)
        @test records == original_order
    end

    # Records without timestamps retain frame-only support.
    minimal = [StateRecord(1, 1, :monomer), StateRecord(2, 1, :dimer), StateRecord(4, 1, :monomer)]
    source = BasicSMLD(minimal, smld.camera, 4, 1, copy(smld.metadata))
    expected = Dict(1 => [(1, :monomer), (2, :dimer), (4, :monomer)])
    expected_lifetime = analyze_dimer_lifetime(source)
    @test expected_lifetime ≈ 0.2
    for records in (minimal, reverse(minimal), shuffle(minimal))
        reordered = BasicSMLD(records, source.camera, 4, 1, copy(source.metadata))
        @test track_changes(reordered) == expected
        @test isequal(analyze_dimer_lifetime(reordered), expected_lifetime)
    end

    # Stable sorting retains the existing sequence when temporal keys tie.
    tied = [DiffusingEmitter2D{Float64}(1.0, 1.0, 100.0, 0.01, 1, 1, 1, state, partner)
            for (state, partner) in ((:monomer, nothing), (:dimer, 2), (:monomer, nothing))]
    tied_source = BasicSMLD(tied, smld.camera, 1, 1, copy(smld.metadata))
    @test track_changes(tied_source)[1] == [(1, :monomer), (1, :dimer), (1, :monomer)]
    tied_minimal = BasicSMLD([StateRecord(e.frame, e.track_id, e.state) for e in tied], smld.camera, 1, 1)
    @test track_changes(tied_minimal)[1] == [(1, :monomer), (1, :dimer), (1, :monomer)]
end
