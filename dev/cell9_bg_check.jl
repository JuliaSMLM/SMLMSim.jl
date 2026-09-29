# Cell9-like background draw (design section 4, "Cell9 mapping").
# Prints the statistics of the noise-free maps next to the Cell9 measurements; asserts nothing.
# The Cell9 match itself is PPIDetect's Gate A comparison.
#
# Usage: julia --project=dev dev/cell9_bg_check.jl

using SMLMSim, MicroscopePSFs, Random, Statistics

T = 0.01                      # frame time, s
px = 0.078                    # pixel size, μm
n = 128
nframes = 400
camera = IdealCamera(1:n, 1:n, px)

# far out-of-focus blobs: rho 0.04 um^-2, gamma 5000/s, brightness sigma 0.3, sigma 0.39 um, tau 0.1 s, D ~ U(0, 0.5)
blobs = Population(name = :oof, layer = :oof, density = 0.04, lifetime = 0.1,
    mobility = [(0.25, 0.0625), (0.25, 0.1875), (0.25, 0.3125), (0.25, 0.4375)],
    brightness_sigma = 0.3, fluor = GenericFluor(; γ = 5000.0, q = zeros(1, 1)),
    z = (0.5, 1.0), psf = GaussianPSF(0.39))

# haze, the one-frame part: rho_h 0.07 -> birth_rate rho_h/(2T) = 3.5 um^-2 s^-1, lifetime 2 ms
haze = Population(name = :haze, layer = :oof, density = 3.5 * 0.002, birth_rate = 3.5, lifetime = 0.002,
    fluor = GenericFluor(; γ = 25000.0, q = zeros(1, 1)), z = (0.5, 1.0), psf = GaussianPSF(0.195))

# level 87.7 ph/px/s, jitter 0.025, illumination width 30 um
bg = BackgroundModel(level = 87.7, jitter = 0.025, contrast = 0.1, feature_size = 0.8,
                     correlation_time = 0.3, illumination_width = 30.0)

structured, oof, level = gen_background(Random.Xoshiro(1), camera, bg, nframes;
                                        oof = [blobs, haze], frame_time = T, t_burn = 1.0)
total = structured .+ oof

gmean = vec(mean(total; dims = (1, 2)))
relstep(k) = sqrt(mean(((gmean[1+k:end] .- gmean[1:end-k]) ./ mean(gmean)) .^ 2) / 2)
blockmeans(A, b) = [mean(A[i:i+b-1, j:j+b-1, k]) for i in 1:b:size(A, 1)-b+1, j in 1:b:size(A, 2)-b+1, k in axes(A, 3)]

println("statistic                                   simulated    Cell9")
println("masked mean, ph/px/frame                    ", round(mean(total); digits = 3), "        0.91")
println("global step lag 1, relative                 ", round(relstep(1); digits = 4), "       0.020")
println("global step lag 4, relative                 ", round(relstep(4); digits = 4), "       0.025")
for b in (8, 16, 32)
    bm = blockmeans(oof, b)
    println("OOF block rms, ", lpad(b, 2), " px blocks, ph/px          ", round(sqrt(mean((bm .- mean(bm)) .^ 2)); digits = 3), "        0.06 (8-32 px)")
end
bm8 = blockmeans(total, 8)
for lag in (1, 4, 10, 30)
    d = bm8[:, :, 1+lag:end] .- bm8[:, :, 1:end-lag]
    println("D8(", lpad(lag, 2), "), rms of 8x8 block-mean differences  ", round(sqrt(mean(d .^ 2)); digits = 3))
end
println("mean level, ph/px/frame                     ", round(mean(level); digits = 3), "        ", round(87.7 * T; digits = 3))
