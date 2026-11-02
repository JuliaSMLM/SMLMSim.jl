using SMLMSim, Test, Random, TOML

# Molecules cover the whole camera, its last pixel row and column included (T27).
#
# Static: the field is the camera's full extent, from its first to its last pixel edge
# (SMLMData's pixel-edge convention). Through 0.7.2 it was the span of the pixel centers,
# one pixel short, starting at 0: the last pixel row and column never held an emitter and
# the density was ((n-1)/n)^2 of the request (0.879 at 16 px, 0.951 at 40 px). Pattern
# centers are drawn past the camera by the pattern's radius and only sites on the camera
# kept, so patterns larger than a point are flat to the edge too.
#
# Diffusion: the default camera is the fewest 100 nm pixels that cover the box. Through 0.7.2
# a box of a whole number of pixels up to rounding (24 * 0.1 is 24.000000000000004 pixels)
# got one more pixel, which no molecule could reach.
#
# Point patterns (one emitter per pattern, d = 0) so that every emitter is a pattern center
# and nothing overhangs the field. The count of a static run is Poisson with mean
# density·area, and each emitter falls in a pixel column with probability width/extent;
# counts are checked at 4σ, at a fixed seed.

const POINT2D = Nmer2D(n=1, d=0.0)
const POINT3D = Nmer3D(n=1, d=0.0)
const FLUOR = GenericFluor(photons=1e4, k_off=50.0, k_on=1.0)

true_emitters(density, pattern, camera; ndims=2) = SMLMSim.simulate(
    StaticSMLMConfig(density=density, nframes=1, ndatasets=1, minphotons=0, ndims=ndims,
                     zrange=[-0.5, 0.5]);
    pattern=pattern, molecule=FLUOR, camera=camera)[2].smld_true.emitters

# Count within 4σ of a Poisson mean λ
within_poisson(n, λ) = abs(n - λ) <= 4 * sqrt(λ)
# Count within 4σ of a binomial(N, p) mean
within_binomial(k, N, p) = abs(k - N * p) <= 4 * sqrt(N * p * (1 - p))
# A count of emitters from Poisson-many patterns of at most m sites each has variance at most
# m times its mean; within 4σ of its expected value μ
within_clustered(k, μ, m) = abs(k - μ) <= 4 * sqrt(m * μ)

# Zero-extent cases recorded at 64e88f7, the commit SMLMDeepFit's retrain is pinned to. The
# recorder, dev/record_zero_extent.jl, uses these same cases (keep them identical) and says how
# to re-record when only the toolchain's random numbers changed.
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

# Coordinates `c` in pixel `i` of the edges `edges`, and that pixel's share of the extent
in_pixel(c, edges, i) = count(v -> edges[i] <= v < edges[i+1], c)
share(edges, i) = (edges[i+1] - edges[i]) / (last(edges) - first(edges))

# Every coordinate on the camera; the first and last pixels hold their share, and the
# coordinates reach the outer edge: the largest of N uniform draws is within δ = 16L/N of
# it except with probability about e^-16 (L the extent), where 0.7.2 stopped a pixel short
function covers(c, edges)
    N, n = length(c), length(edges) - 1
    @test all(v -> first(edges) <= v <= last(edges), c)
    @test within_binomial(in_pixel(c, edges, 1), N, share(edges, 1))
    @test within_binomial(in_pixel(c, edges, n), N, share(edges, n))
    δ = 16 * (last(edges) - first(edges)) / N
    @test δ < edges[n+1] - edges[n]  # so the check below separates the fix from 0.7.2
    @test maximum(c) > last(edges) - δ
end

@testset "molecules cover the whole camera" begin
    static_cases = (
        ("16 px, 2D", IdealCamera(16, 16, 0.1), 2, 4000.0),
        ("40 px, 2D", IdealCamera(40, 40, 0.1), 2, 1500.0),
        ("16 px, 3D", IdealCamera(16, 16, 0.1), 3, 4000.0),
        ("40 px, 3D", IdealCamera(40, 40, 0.1), 3, 1500.0),
        ("16 × 40 px", IdealCamera(16, 40, 0.1), 2, 1500.0),
        ("edges that do not start at 0",
         IdealCamera(collect(5.0:0.1:6.6), collect(-2.0:0.1:-0.4)), 2, 4000.0),
        ("sCMOS, Float32 rectangular pixels", SCMOSCamera(16, 24, (0.1f0, 0.15f0), 1.6f0), 2, 4000.0),
    )
    @testset "static: $name" for (k, (name, camera, dim, density)) in enumerate(static_cases)
        Random.seed!(2700 + k)
        ex, ey = camera.pixel_edges_x, camera.pixel_edges_y
        es = true_emitters(density, dim == 2 ? POINT2D : POINT3D, camera; ndims=dim)
        area = (Float64(last(ex)) - first(ex)) * (Float64(last(ey)) - first(ey))
        # The measured density is the request
        @test within_poisson(length(es), density * area)
        covers([e.x for e in es], ex)
        covers([e.y for e in es], ey)
        dim == 3 && @test all(e -> -0.5 <= e.z <= 0.5, es)
    end

    # Patterns larger than a point are flat to the edge: every emitter on the camera, the total
    # at density × sites × area, and the edge pixel columns and rows as full as the middle ones.
    # At 64e88f7 the default ring at 16 px gave 0.965 of the request and edge columns about
    # 0.87 of a middle one.
    @testset "static: flat to the edge, $name" for (k, (name, pattern, dim, density)) in enumerate((
            ("ring of 8, 0.1 µm, 2D", Nmer2D(), 2, 10000.0),
            ("ring of 8, 0.2 µm, 3D", Nmer3D(n=8, d=0.2), 3, 3000.0)))
        Random.seed!(2770 + k)
        camera = IdealCamera(16, 16, 0.1)
        ex, ey = camera.pixel_edges_x, camera.pixel_edges_y
        es = true_emitters(density, pattern, camera; ndims=dim)
        m = pattern.n
        μ = density * m * (last(ex) - first(ex)) * (last(ey) - first(ey))  # expected on the camera
        x, y = [e.x for e in es], [e.y for e in es]
        @test all(v -> first(ex) <= v <= last(ex), x)
        @test all(v -> first(ey) <= v <= last(ey), y)
        @test within_clustered(length(es), μ, m)
        for (c, edges) in ((x, ex), (y, ey)), i in (1, 2, 8, 15, 16)
            @test within_clustered(in_pixel(c, edges, i), μ / 16, m)
        end
    end

    # Zero-extent patterns are drawn over the camera extent itself: for the same seed, the true
    # positions are exactly a direct uniform2D/uniform3D draw over the extent, shifted to its first
    # edges (no stored numbers, so this holds on any toolchain)
    @testset "static: zero-extent patterns are a direct draw over the camera, $name" for (name, case) in ZERO_EXTENT_CASES
        ex, ey = case.camera.pixel_edges_x, case.camera.pixel_edges_y
        x0, y0 = Float64(first(ex)), Float64(first(ey))
        fx, fy = Float64(last(ex)) - x0, Float64(last(ey)) - y0
        Random.seed!(case.seed)
        direct = case.ndims == 2 ? uniform2D(30.0, case.pattern, fx, fy) :
                 uniform3D(30.0, case.pattern, fx, fy; zrange=[-0.5, 0.5])
        got = run_zero_extent_case(case)
        @test got["true_x"] == direct[1] .+ x0
        @test got["true_y"] == direct[2] .+ y0
        case.ndims == 3 && @test got["true_z"] == direct[3]
    end

    # Zero-extent patterns are placed exactly as at 64e88f7: the same seed gives bitwise the same
    # emitters, hence the same rendered frames (recorded on Julia 1.10.11, 1.11.9, 1.12.7 and
    # 1.13.0, identical on all four; dev/record_zero_extent.jl re-records it)
    reference = TOML.parsefile(joinpath(@__DIR__, "field_extent_64e88f7.toml"))
    @testset "static: zero-extent patterns as at 64e88f7, $name" for (name, case) in ZERO_EXTENT_CASES
        ref, got = reference[name], run_zero_extent_case(case)
        @test sort(collect(keys(got))) == sort(collect(keys(ref)))
        for key in sort(collect(keys(ref)))
            @test isequal(got[key], ref[key])
        end
    end

    @testset "diffusion default camera: box of $n pixels" for n in (3, 24, 48)
        Random.seed!(2750 + n)
        params = DiffusionSMLMConfig(box_size=n * 0.1, density=100.0, t_max=0.1)
        smld, _ = SMLMSim.simulate(params)
        ex, ey = smld.camera.pixel_edges_x, smld.camera.pixel_edges_y
        @test length(ex) - 1 == n
        @test length(ey) - 1 == n
        # One record per molecule: its first in frame 1
        first_record = Dict{Int,Any}()
        for e in smld.emitters
            e.frame == 1 && !haskey(first_record, e.track_id) && (first_record[e.track_id] = e)
        end
        es = collect(values(first_record))
        if n > 3  # too few molecules at 3 px for the share test
            covers([e.x for e in es], ex)
            covers([e.y for e in es], ey)
        end
    end

    @testset "diffusion default camera: box from a Float32 pixel size" begin
        smld, _ = SMLMSim.simulate(DiffusionSMLMConfig(box_size=16 * 0.1f0, density=10.0, t_max=0.1))
        @test length(smld.camera.pixel_edges_x) - 1 == 16
    end

    # The box is the sample and a passed camera a view into it: a camera reaching outside the box
    # warns once, naming the pixels outside; a camera inside it, the default one included, does not
    @testset "diffusion with a passed camera" begin
        p = DiffusionSMLMConfig(box_size=2.4, density=1.0, t_max=0.2)
        @test_logs (:warn, r"columns 25-30 and rows 25-30 lie wholly outside it") SMLMSim.simulate(
            p; camera=IdealCamera(30, 30, 0.1))
        @test_logs (:warn, r"columns 1-2 lie wholly outside it.*columns 3 lie partly outside it") SMLMSim.simulate(
            p; camera=IdealCamera(collect(-0.25:0.1:1.05), collect(0.0:0.1:1.0)))
        quiet = Base.CoreLogging.Warn
        @test_logs min_level=quiet SMLMSim.simulate(p; camera=IdealCamera(16, 16, 0.1))
        @test_logs min_level=quiet SMLMSim.simulate(p; camera=IdealCamera(collect(1.0:0.1:2.4), collect(0.5:0.1:1.5)))
        own = SMLMSim.simulate(p)[1].camera
        @test_logs min_level=quiet SMLMSim.simulate(p; camera=own)
    end

    @testset "diffusion default camera: box between whole pixels" begin
        smld, _ = SMLMSim.simulate(DiffusionSMLMConfig(box_size=2.45, density=10.0, t_max=0.1))
        @test length(smld.camera.pixel_edges_x) - 1 == 25
        @test last(smld.camera.pixel_edges_x) >= 2.45
    end
end
