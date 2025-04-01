using Revise 
using SMLMSim
using SMLMData
using GLMakie

# Create camera with physical pixel size
camera = IdealCamera(1:128, 1:256, 0.1)  # 128x256 pixels, 100nm pixels

# Create StaticSMLMParams with simulation parameters
params = StaticSMLMParams(
    ρ=100.0,              # emitters per μm²
    σ_psf=0.13,           # PSF width in μm
    minphotons=50,        # minimum photons for detection
    ndatasets=10,         # number of independent datasets
    nframes=1000,         # frames per dataset
    framerate=50.0,       # frames per second
    ndims=3,              # 3D simulation
    zrange=[-1.0, 1.0]    # Set z range for 3D patterns
)

# Run simulation with type-based interface
pattern = SMLMSim.Nmer3D()  # Using 3D pattern
molecule = SMLMSim.GenericFluor(; q=[0 50; 1e-2 0])

smld_true, smld_model, smld_noisy = simulate(
    params,
    pattern=pattern,
    molecule=molecule,
    camera=camera
)

# Extract coordinates from emitters
x_noisy = [e.x for e in smld_noisy.emitters]
y_noisy = [e.y for e in smld_noisy.emitters]
z_noisy = [e.z for e in smld_noisy.emitters]
photons = [e.photons for e in smld_noisy.emitters]

# Create multi-panel figure
fig = Figure(size=(1200, 800))

# Create the 3D axes
ax1 = Axis3(fig[1, 1],
    title="True Positions",
    xlabel="x (μm)",
    ylabel="y (μm)",
    zlabel="z (μm)",
    aspect=:data)

ax2 = Axis3(fig[1, 2],
    title="Model with Blinking",
    xlabel="x (μm)",
    ylabel="y (μm)",
    zlabel="z (μm)",
    aspect=:data)

ax3 = Axis3(fig[2, 1:2],
    title="Noisy Localizations",
    xlabel="x (μm)",
    ylabel="y (μm)",
    zlabel="z (μm)",
    aspect=:data)

# Plot true positions
x_true = [e.x for e in smld_true.emitters]
y_true = [e.y for e in smld_true.emitters]
z_true = [e.z for e in smld_true.emitters]
scatter!(ax1, x_true, y_true, z_true, 
    color=:blue, 
    markersize=4)

# Plot model with blinking
x_model = [e.x for e in smld_model.emitters]
y_model = [e.y for e in smld_model.emitters]
z_model = [e.z for e in smld_model.emitters]
scatter!(ax2, x_model, y_model, z_model, 
    color=:green, 
    markersize=4, 
    alpha=0.5)

# Plot noisy localizations
scatter!(ax3, x_noisy, y_noisy, z_noisy,
    color=photons,
    colormap=:viridis,
    markersize=4,
    alpha=0.5)
Colorbar(fig[2, 3], label="Photons")

# Add 2D projections
ax_xy = Axis(fig[3, 1],
    title="XY Projection",
    xlabel="x (μm)",
    ylabel="y (μm)",
    aspect=DataAspect())

ax_xz = Axis(fig[3, 2],
    title="XZ Projection",
    xlabel="x (μm)",
    ylabel="z (μm)",
    aspect=DataAspect())

scatter!(ax_xy, x_noisy, y_noisy,
    color=photons,
    colormap=:viridis,
    markersize=2,
    alpha=0.3)

scatter!(ax_xz, x_noisy, z_noisy,
    color=photons,
    colormap=:viridis,
    markersize=2,
    alpha=0.3)

# Add visualization controls
ctrl_panel = fig[1:3, 4] = GridLayout()
Label(ctrl_panel[1, 1], "Visualization Controls")

# Point size slider
point_size = Slider(ctrl_panel[2, 1], range=0.5:0.5:10, startvalue=4)
Label(ctrl_panel[2, 2], "Point Size")
on(point_size.value) do val
    for axis in [ax1, ax2, ax3]
        for plt in axis.scene.plots
            if plt isa Scatter
                plt.markersize = val
            end
        end
    end
end

# Alpha slider
alpha_slider = Slider(ctrl_panel[3, 1], range=0.1:0.1:1.0, startvalue=0.5)
Label(ctrl_panel[3, 2], "Opacity")
on(alpha_slider.value) do val
    for axis in [ax2, ax3]
        for plt in axis.scene.plots
            if plt isa Scatter
                plt.alpha = val
            end
        end
    end
end

# Reset view button
reset_btn = Button(ctrl_panel[4, 1:2], label="Reset View")
on(reset_btn.clicks) do n
    # Reset to default view
    zoom!(ax1, 1.0)
    zoom!(ax2, 1.0)
    zoom!(ax3, 1.0)
end

# Add statistics
supertitle = "3D Simulation Results\n" *
    "$(length(x_true)) true positions, " *
    "$(length(x_noisy)) localizations\n" *
    "Z range: [$(round(minimum(z_true), digits=2)), $(round(maximum(z_true), digits=2))] μm"
Label(fig[0, :], supertitle)

# Link axes and camera
# linkaxes!(ax1, ax2, ax3)

fig