# stack_viewer.jl - Simple image stack viewer using GLMakie
using GLMakie

"""
    view_stack(img_stack; title="Image Stack Viewer")

Create a minimalist viewer for a 3D image stack with interactive browsing
and contrast adjustment.

# Arguments
- `img_stack`: 3D array of images (height × width × frames)
- `title`: Window title (optional)

# Example
```julia
# Example with a random image stack
img_stack = rand(Float32, 512, 512, 100)
view_stack(img_stack)

# Example with SMLM image data
using SMLMSim
psf = MicroscopePSFs.Gaussian2D(0.15)  # 150nm PSF width
images = gen_image_sequence(psf, systems)
view_stack(images)
```
"""
function view_stack(img_stack; title="Image Stack Viewer")
    # Get stack dimensions
    height, width, n_frames = size(img_stack)
    
    # Create the figure
    fig = Figure(size=(800, 650), title=title)
    
    # Create layout for image and controls
    image_panel = fig[1, 1]
    controls_panel = fig[2, 1]
    
    # Main image display
    ax = Axis(image_panel[1, 1], aspect=DataAspect(), yreversed=true)
    hidedecorations!(ax)  # Hide axis decorations for cleaner look
    
    # Initial frame to display
    current_frame = Observable(1)
    
    # State for contrast mode
    global_contrast = Observable(false)
    
    # Calculate global min/max for global contrast
    global_min = minimum(img_stack)
    global_max = maximum(img_stack)
    
    # Function to get current frame data
    current_image = @lift(img_stack[:, :, $current_frame]')
    
    # Function to determine color range based on contrast mode
    color_range = @lift begin
        if $global_contrast
            (global_min, global_max)
        else
            frame_data = img_stack[:, :, $current_frame]
            (minimum(frame_data), maximum(frame_data))
        end
    end
    
    # Display the image with dynamic color range
    heatmap!(ax, current_image, colormap=:grays, colorrange=color_range)
    
    # Create controls
    
    # Slider for browsing frames with label
    slider_label = Label(controls_panel[1, 1], "Frame:")
    slider = Slider(controls_panel[1, 2:10], range=1:n_frames, startvalue=1)
    frame_counter = Label(controls_panel[1, 11], @lift("$($current_frame) / $n_frames"))
    
    # Connect slider to current frame
    connect!(current_frame, slider.value)
    
    # Add toggle button for contrast mode
    contrast_btn = Button(controls_panel[2, 1:11], label="Local Contrast")
    
    # Button click handler
    on(contrast_btn.clicks) do _
        global_contrast[] = !global_contrast[]
        contrast_btn.label = global_contrast[] ? "Global Contrast" : "Local Contrast"
    end
    
    # Add a small gap between controls and image
    rowgap!(fig.layout, 10)
    
    # Return the figure
    display(fig)
    return fig
end