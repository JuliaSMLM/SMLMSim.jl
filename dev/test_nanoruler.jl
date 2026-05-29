# test_nanoruler.jl - Simulate GATTAquant-style Nanoruler structures and verify
# the geometry round-trips through the simulation pipeline.
#
# SCOPE: SMLMSim is a *simulation* package. This script exercises the simulation
# pipeline (pattern -> random placement/rotation -> blinking -> localization
# noise) and quantitatively verifies the ground-truth geometry of the new
# Nanoruler patterns. The full *analysis* round-trip (camera images ->
# localization fitting -> blink grouping -> clustering -> recovered spacing)
# belongs in SMLMAnalysis, not here -- SMLMSim has no fitting code.
#
# API: v0.7.0 -- StaticSMLMConfig + 2-tuple `simulate` return (smld, info).
# Runs headless by default; set SHOW_PLOTS = true for interactive GLMakie figures.

using Pkg
Pkg.activate(@__DIR__)
using Revise

using SMLMSim
using Statistics
using LinearAlgebra

const SHOW_PLOTS = false   # set true for interactive GLMakie figures

# --- helpers ---------------------------------------------------------------

"Group emitter indices by pattern id (one group per ruler instance)."
function group_by_id(emitters)
    groups = Dict{Int,Vector{Int}}()
    for (i, e) in enumerate(emitters)
        push!(get!(groups, e.id, Int[]), i)
    end
    return groups
end

"""
    recovered_gaps(pts) -> Vector

Recover the adjacent-mark spacings of (near-)collinear points by projecting
them onto their principal axis (PCA) and taking sorted gaps. `pts` is an
`(n, d)` matrix. Rotation-invariant; exact for noise-free collinear marks.
"""
function recovered_gaps(pts::AbstractMatrix)
    μ = mean(pts, dims=1)
    centered = pts .- μ
    C = Symmetric(centered' * centered)
    axis = eigvecs(C)[:, end]              # principal direction (largest eigenvalue)
    proj = sort(vec(centered * axis))
    return diff(proj)
end

# --- simulation ------------------------------------------------------------

spacing = 0.02                              # 20 nm between adjacent marks
camera  = IdealCamera(1:128, 1:128, 0.1)    # 12.8 μm field, 100 nm pixels
ruler   = Nanoruler2D(spacing=spacing)      # default: 3-mark GATTAquant ruler

# High photon budget -> localization precision (σ_psf / sqrt(photons)) well below
# the 20 nm spacing, so the sub-diffraction ruler is resolvable in the noisy data.
fluor = GenericFluor(photons=1e5, k_off=50.0, k_on=1e-2)

params = StaticSMLMConfig(
    density    = 1.0,    # rulers per μm^2
    σ_psf      = 0.13,   # 130 nm PSF; localization precision = σ_psf / sqrt(photons)
    minphotons = 50,
    ndatasets  = 1,
    nframes    = 1000,
    framerate  = 50.0,
    ndims      = 2,
)

println("Simulating Nanoruler2D (n=$(ruler.n), spacing=$(spacing*1e3) nm)...")
smld, info = simulate(params; pattern=ruler, molecule=fluor, camera=camera)

println("\nSimulation summary")
println("  rulers placed      : ", info.n_patterns)
println("  ground-truth marks : ", info.n_emitters,
        "  (expected ", info.n_patterns * ruler.n, ")")
println("  noisy localizations: ", info.n_localizations)

# --- 1. ground-truth geometry check (rigorous) -----------------------------

println("\n[1] Ground-truth geometry (info.smld_true)")
gt = info.smld_true.emitters
groups = group_by_id(gt)

all_gaps = Float64[]
nmark_ok = true
for (id, idxs) in groups
    global nmark_ok &= length(idxs) == ruler.n
    pts = reduce(vcat, ([gt[i].x gt[i].y] for i in idxs))
    append!(all_gaps, recovered_gaps(pts))
end

println("  rulers checked          : ", length(groups))
println("  all have n=$(ruler.n) marks    : ", nmark_ok)
println("  recovered spacing       : ", round(mean(all_gaps)*1e3, digits=4), " ± ",
        round(std(all_gaps)*1e3, digits=4), " nm  (true = ", spacing*1e3, " nm)")
pass = nmark_ok && isapprox(mean(all_gaps), spacing; atol=1e-6)
println("  RESULT: ", pass ? "PASS" : "FAIL")

# --- 2. localization precision (is the structure resolvable after imaging?) -

println("\n[2] Localization precision in noisy data")
photons = [e.photons for e in smld.emitters]
σ_loc = params.σ_psf ./ sqrt.(photons)      # SMLM precision = σ_psf / sqrt(N)
println("  median photons/loc : ", round(median(photons), digits=0))
println("  median σ_loc       : ", round(median(σ_loc)*1e3, digits=2), " nm")
println("  spacing / σ_loc    : ", round(spacing / median(σ_loc), digits=1),
        "   (> ~2 means the 20 nm marks are resolvable)")
println("\nNote: recovering spacing from the *noisy* localizations (grouping blinks,")
println("      fitting, clustering) is an analysis task -> SMLMAnalysis.")

# --- 3. optional visualization ---------------------------------------------

if SHOW_PLOTS
    using GLMakie
    noisy_groups = group_by_id(smld.emitters)
    # pick the ruler with the most noisy localizations for a clear close-up
    id_best = argmax(id -> length(noisy_groups[id]), keys(noisy_groups))
    gi, ni = groups[id_best], noisy_groups[id_best]

    fig = Figure(size=(950, 460))
    ax1 = Axis(fig[1, 1], title="All rulers (ground truth)",
               xlabel="X (μm)", ylabel="Y (μm)", aspect=DataAspect())
    scatter!(ax1, [e.x for e in gt], [e.y for e in gt],
             color=[e.id for e in gt], colormap=:viridis, markersize=4)

    ax2 = Axis(fig[1, 2], title="One ruler: truth (×) vs noisy (·)",
               xlabel="X (μm)", ylabel="Y (μm)", aspect=DataAspect())
    scatter!(ax2, [smld.emitters[i].x for i in ni], [smld.emitters[i].y for i in ni],
             color=(:red, 0.3), markersize=5, label="noisy")
    scatter!(ax2, [gt[i].x for i in gi], [gt[i].y for i in gi],
             color=:black, marker=:xcross, markersize=16, label="true marks")
    axislegend(ax2)
    display(fig)
    println("\nClose the window to exit.")
end
