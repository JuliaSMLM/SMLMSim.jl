using SMLMSim, Test, Statistics, Random, StableRNGs, MicroscopePSFs, Distributions
using SpecialFunctions: erf, erfc
using SMLMSim.CameraImages: RenderBuffer, render_gaussian!, StampTable, render_stamp!,
                            scmos_noise!, poisson_noise!

em2(x, y, ph, frame) = Emitter2DFit{Float64}(x, y, ph, 0.0, 0.0, 0.0, 0.0, 0.0; frame=frame)
em3(x, y, z, ph, frame) = Emitter3DFit{Float64}(x, y, z, ph, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; frame=frame)

# analytic 1D box integral of a unit Gaussian, pixel k spanning [k-1, k)
box(k, c, σ) = 0.5 * (erf((k - c) / (sqrt(2) * σ)) - erf((k - 1 - c) / (sqrt(2) * σ)))

@testset "render/erf_exact" begin
    buf = RenderBuffer(8)
    for (u, v, σ) in ((7.3, 9.9, 1.2), (5.0, 5.5, 0.8), (10.01, 3.2, 1.6))
        img = zeros(20, 20)
        render_gaussian!(img, buf, u, v, σ, 1000.0, 8)
        cu, cv = floor(Int, u) + 1, floor(Int, v) + 1
        for i in 1:20, j in 1:20
            inwin = abs(j - cu) <= 8 && abs(i - cv) <= 8
            ref = inwin ? 1000.0 * box(j, u, σ) * box(i, v, σ) : 0.0
            @test isapprox(img[i, j], ref; atol=1e-14 * 1000.0, rtol=0)
        end
        wx = sum(box(j, u, σ) for j in (cu - 8):(cu + 8) if 1 <= j <= 20)
        wy = sum(box(i, v, σ) for i in (cv - 8):(cv + 8) if 1 <= i <= 20)
        @test isapprox(sum(img), 1000.0 * wx * wy; rtol=1e-13)
    end
    # tail accuracy far from the centre
    img = zeros(1, 40)
    render_gaussian!(img, buf, 0.5, 0.5, 1.0, 1.0, 1:40, 1:1)
    @test img[1, 12] > 0 && isapprox(img[1, 12], 0.5 * (erfc(10.5 / sqrt(2)) - erfc(11.5 / sqrt(2))) * box(1, 0.5, 1.0); rtol=1e-8)
end

@testset "render/zero_alloc" begin
    img = zeros(64, 64)
    buf = RenderBuffer(8)
    render_gaussian!(img, buf, 20.3, 30.7, 1.2, 1000.0, 7)
    render_gaussian!(img, buf, 20.3, 30.7, 1.2, 1000.0, 1:64, 1:64)
    @test (@allocated render_gaussian!(img, buf, 20.3, 30.7, 1.2, 1000.0, 7)) == 0
    @test (@allocated render_gaussian!(img, buf, 20.3, 30.7, 1.2, 1000.0, 10:30, 20:40)) == 0
    tbl = StampTable(GaussianPSF(0.13), 0.1, range(-0.3, 0.3, length=5); radius=6)
    render_stamp!(img, tbl, 20.3, 30.7, 0.1, 500.0)
    @test (@allocated render_stamp!(img, tbl, 20.3, 30.7, 0.1, 500.0)) == 0
    tbl_n = StampTable(GaussianPSF(0.13), 0.1, range(-0.3, 0.3, length=5); radius=6, zinterp=:nearest)
    render_stamp!(img, tbl_n, 20.3, 30.7, 0.1, 500.0)
    @test (@allocated render_stamp!(img, tbl_n, 20.3, 30.7, 0.1, 500.0)) == 0
end

@testset "gen_images/gaussian_parity" begin
    cam = IdealCamera(32, 32, 0.1)
    es = [em2(1.23, 1.71, 1000.0, 1), em2(2.2, 0.9, 700.0, 1), em2(1.6, 1.6, 900.0, 2)]
    smld = BasicSMLD(es, cam, 2, 1)
    for σpx in (1.0, 1.2, 1.3)
        psf = GaussianPSF(σpx * 0.1)
        new, _ = gen_images(smld, psf)
        s16, _ = gen_images(smld, psf; sampling=16)
        s2, _ = gen_images(smld, psf; sampling=2)
        pk = maximum(new)
        @test maximum(abs.(new .- s16)) <= 5e-4 * pk
        @test maximum(abs.(new .- s2)) <= 2.5e-2 * pk
        @test isapprox(sum(new), sum(s16); rtol=1e-9)
        sup, _ = gen_images(smld, psf; support=0.4)
        supold, _ = gen_images(smld, psf; support=0.4, sampling=2)
        ref = zeros(32, 32, 2)
        for (k, f) in ((1, 1), (2, 2))
            ref[:, :, k] = integrate_pixels(psf, cam, filter(e -> e.frame == f, es); support=0.4)
        end
        @test (sup .!= 0) == (ref .!= 0)
        tup, _ = gen_images(smld, psf; support=(0.5, 2.5, 0.5, 2.5))
        @test all(isfinite, tup) && sum(tup) > 0
    end
    # frame totals for the default support agree to round-off with the full image
    psf = GaussianPSF(0.12)
    new, _ = gen_images(smld, psf)
    @test isapprox(sum(new[:, :, 1]), 1700.0; rtol=1e-9)
end

@testset "gen_images/fallback_bitwise" begin
    cam = IdealCamera(24, 24, 0.1)
    zr = range(-0.6, 0.6, length=7)
    xr = range(-0.8, 0.8, length=17)
    sp = SplinePSF(ScalarPSF(1.4, 0.532, 1.518), xr, xr, zr)
    es3 = [em3(1.1, 1.3, 0.1, 1000.0, 1), em3(1.15, 1.25, -0.2, 800.0, 1), em3(0.7, 1.9, 0.3, 500.0, 2)]
    es2 = [em2(1.1, 1.3, 1000.0, 1), em2(1.15, 1.25, 800.0, 1), em2(0.7, 1.9, 500.0, 2)]
    for (psf, es) in ((AiryPSF(1.4, 0.532), es2), (sp, es3))
        smld = BasicSMLD(es, cam, 3, 1)
        for support in (Inf, 0.5, (0.3, 1.6, 0.4, 1.7)), threaded in (true, false)
            imgs, _ = gen_images(smld, psf; support=support, threaded=threaded, bg=2.5)
            @test eltype(imgs) == Float64
            for f in 1:3
                fe = filter(e -> e.frame == f, es)
                old = isempty(fe) ? zeros(24, 24) : integrate_pixels(psf, cam, fe; support=support, sampling=2)
                @test imgs[:, :, f] == 2.5 .+ old
            end
        end
    end
end

@testset "gen_images/bg_stack" begin
    cam = IdealCamera(16, 16, 0.1)
    smld = BasicSMLD([em2(0.8, 0.8, 500.0, 1), em2(0.9, 0.7, 300.0, 2)], cam, 3, 1)
    psf = GaussianPSF(0.12)
    stack = rand(StableRNG(1), 16, 16, 3) .* 10
    a, _ = gen_images(smld, psf; bg=stack)
    b, _ = gen_images(smld, psf; bg=0.0)
    @test a == b .+ stack
    c, _ = gen_images(smld, psf; bg=4.0)
    @test c == b .+ 4.0
    @test_throws DimensionMismatch gen_images(smld, psf; bg=zeros(16, 16, 2))
    @test_throws DimensionMismatch gen_images(smld, psf; bg=zeros(16, 15, 3))
end

@testset "gen_images/rng" begin
    cam = SCMOSCamera(16, 16, 0.1, 1.6)
    smld = BasicSMLD([em2(0.8, 0.8, 500.0, 1), em2(0.9, 0.7, 300.0, 2)], cam, 2, 1)
    psf = GaussianPSF(0.12)
    r1, _ = gen_images(smld, psf; bg=3.0, camera_noise=true, rng=Xoshiro(7))
    r2, _ = gen_images(smld, psf; bg=3.0, camera_noise=true, rng=Xoshiro(7))
    @test r1 == r2
    r3, _ = gen_images(smld, psf; bg=3.0, camera_noise=true, rng=Xoshiro(8))
    @test r1 != r3
    p1, _ = gen_images(smld, psf; bg=3.0, poisson_noise=true, rng=StableRNG(3))
    p2, _ = gen_images(smld, psf; bg=3.0, poisson_noise=true, rng=StableRNG(3))
    @test p1 == p2
    # default RNG path equals the old inline algorithm
    clean, _ = gen_images(smld, psf; bg=3.0)
    Random.seed!(1)
    d, _ = gen_images(smld, psf; bg=3.0, poisson_noise=true)
    Random.seed!(1)
    old = copy(clean)
    for i in eachindex(old)
        λ = max(old[i], 0.0)
        old[i] = λ > 0.0 ? Float64(rand(Poisson(λ))) : 0.0
    end
    @test d == old
    Random.seed!(2)
    d, _ = gen_images(smld, psf; bg=3.0, camera_noise=true)
    Random.seed!(2)
    old = copy(clean)
    for f in 1:2
        fr = @view old[:, :, f]
        for j in 1:16, i in 1:16
            pe = max(fr[i, j], 0.0) * cam.qe
            pe = pe > 0.0 ? Float64(rand(Poisson(pe))) : 0.0
            cam.readnoise > 0.0 && (pe += randn() * cam.readnoise)
            fr[i, j] = pe / cam.gain + cam.offset
        end
    end
    @test d == old
end

@testset "noise/scmos_moments" begin
    n = 200
    μ = 40.0
    for cam in (SCMOSCamera(n, n, 0.1, 2.0; offset=100.0, gain=2.0, qe=0.8),
                SCMOSCamera(n, n, 0.1, fill(2.0, n, n); offset=fill(100.0, n, n), gain=fill(2.0, n, n), qe=fill(0.8, n, n)))
        img = fill(μ, n, n)
        scmos_noise!(StableRNG(11), img, cam)
        @test abs(mean(img) - (100.0 + 0.8 * μ / 2.0)) < 3 * sqrt((0.8 * μ + 4.0) / 4.0 / n^2)
        @test isapprox(var(img), (0.8 * μ + 4.0) / 4.0; rtol=0.03)
        Random.seed!(5)
        a = fill(μ, n, n); scmos_noise!(a, cam)
        Random.seed!(5)
        b = fill(μ, n, n); scmos_noise!(Random.default_rng(), b, cam)
        @test a == b
    end
    # the old algorithm's stream
    cam = SCMOSCamera(8, 8, 0.1, 1.6)
    Random.seed!(9)
    a = fill(30.0, 8, 8); scmos_noise!(a, cam)
    Random.seed!(9)
    b = fill(30.0, 8, 8)
    for j in 1:8, i in 1:8
        pe = Float64(rand(Poisson(b[i, j] * cam.qe))) + randn() * cam.readnoise
        b[i, j] = pe / cam.gain + cam.offset
    end
    @test a == b
end

@testset "background/stamp_matches_psf" begin
    px = 0.1
    psf = ScalarPSF(1.4, 0.532, 1.518)
    zs = range(-0.4, 0.4, length=5)
    os, r = 4, 6
    tbl = StampTable(psf, px, zs; radius=r, oversample=os)
    n = 2r + 1
    edges = collect(0:n) .* px
    for (iz, z) in enumerate(zs)
        # phase 0: the emitter at 0.5/os of the centre pixel
        e = Emitter3D(r * px + 0.5 * px / os, r * px + 0.5 * px / os, z, 1.0)
        ref = integrate_pixels(psf, edges, edges, e; sampling=2 * os)
        ref ./= sum(ref)
        @test isapprox(tbl.stamps[:, :, 1, iz], ref; atol=1e-12, rtol=0)
        @test isapprox(sum(tbl.stamps[:, :, 1, iz]), 1.0; atol=1e-12)
    end
    # rendering at a phase-0 position reproduces the stamp
    img = zeros(20, 20)
    render_stamp!(img, tbl, 10.0 + 0.5 / os, 10.0 + 0.5 / os, zs[3], 100.0)
    @test isapprox(img[10-r+1:10+r+1, 10-r+1:10+r+1] ./ 100, tbl.stamps[:, :, 1, 3]; atol=1e-12) ||
          isapprox(img[(11-r):(11+r), (11-r):(11+r)] ./ 100, tbl.stamps[:, :, 1, 3]; atol=1e-12)
end

@testset "StampTable/z_range_and_interp" begin
    px = 0.1
    xr = range(-0.8, 0.8, length=17)
    sp = SplinePSF(ScalarPSF(1.4, 0.532, 1.518), xr, xr, range(-0.6, 0.6, length=7))
    # (R-a) construction beyond the PSF's z range throws
    @test_throws ArgumentError StampTable(sp, px, range(-0.8, 0.8, length=9); radius=5)
    tbl = StampTable(sp, px, range(-0.6, 0.6, length=7); radius=5)
    img = zeros(20, 20)
    # render beyond the table's z range throws
    @test_throws ArgumentError render_stamp!(img, tbl, 10.0, 10.0, 0.61, 1.0)
    @test_throws ArgumentError render_stamp!(img, tbl, 10.0, 10.0, -0.7, 1.0)
    @test_throws ArgumentError StampTable(sp, px, range(0.0, 0.0, length=1); radius=5, zinterp=:cubic)
    # (R-b) linear midway between planes equals the mean of the two stamps
    zs = tbl.zs
    zm = (zs[3] + zs[4]) / 2
    a = zeros(20, 20); render_stamp!(a, tbl, 10.0, 10.0, zm, 1.0)
    b = zeros(20, 20); render_stamp!(b, tbl, 10.0, 10.0, zs[3], 1.0)
    c = zeros(20, 20); render_stamp!(c, tbl, 10.0, 10.0, zs[4], 1.0)
    @test isapprox(a, (b .+ c) ./ 2; atol=1e-14)
    tn = StampTable(sp, px, range(-0.6, 0.6, length=7); radius=5, zinterp=:nearest)
    d = zeros(20, 20); render_stamp!(d, tn, 10.0, 10.0, zs[3] + 0.01, 1.0)
    @test isapprox(d, b; atol=1e-14)
end
