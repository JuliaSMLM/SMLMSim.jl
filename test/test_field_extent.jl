using SMLMSim, Test, Random

# Molecules cover the whole camera, its last pixel row and column included (T27).
#
# Static: the field is the camera's full extent, from its first to its last pixel edge
# (SMLMData's pixel-edge convention). Through 0.7.2 it was the span of the pixel centers,
# one pixel short, starting at 0: the last pixel row and column never held an emitter and
# the density was ((n-1)/n)^2 of the request (0.879 at 16 px, 0.951 at 40 px).
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

    @testset "diffusion default camera: box between whole pixels" begin
        smld, _ = SMLMSim.simulate(DiffusionSMLMConfig(box_size=2.45, density=10.0, t_max=0.1))
        @test length(smld.camera.pixel_edges_x) - 1 == 25
        @test last(smld.camera.pixel_edges_x) >= 2.45
    end
end
