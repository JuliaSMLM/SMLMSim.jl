# stack_viewer.jl - Simple image/stack viewer using GLMakie
using GLMakie

"""
    view_stack(img; title="Image Viewer")

Create a minimalist viewer for a 2D image or 3D image stack with interactive browsing
and contrast adjustment.

# Arguments
- `img`: 2D image or 3D array of images (height × width × frames)
- `title`: Window title (optional)

# Example
```julia
# Example with a single 2D image
img = rand(Float32, 512, 512)
view_stack(img)

# Example with a 3D image stack
img_stack = rand(Float32, 512, 512, 100)
view_stack(img_stack)

# Example with SMLM image data
using SMLMSim
psf = MicroscopePSFs.Gaussian2D(0.15)  # 150nm PSF width
images = gen_image_sequence(psf, systems)
view_stack(images)
```
"""
function view_stack(img; title="Image Viewer")
    # Check dimensions and handle either 2D or 3D input
    dims = ndims(img)
    
    if dims == 2
        # Convert 2D image to 3D stack with a single frame
        img_stack = reshape(img, size(img,1), size(img,2), 1)
        is_stack = false
    elseif dims == 3
        # Use 3D stack directly
        img_stack = img
        is_stack = true
    else
        error("Input must be a 2D image or 3D image stack")
    end
    
    # Get stack dimensions
    height, width, n_frames = size(img_stack)
    
    # Create the figure - adjust title based on input type
    auto_title = is_stack ? "Image Stack Viewer" : "Image Viewer"
    fig = Figure(size=(800, 650), title=title == "Image Viewer" ? auto_title : title)
    
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
    # Only show frame navigation controls for 3D stacks
    if is_stack
        # Slider for browsing frames with label
        slider_label = Label(controls_panel[1, 1], "Frame:")
        slider = Slider(controls_panel[1, 2:10], range=1:n_frames, startvalue=1)
        frame_counter = Label(controls_panel[1, 11], @lift("$($current_frame) / $n_frames"))
        
        # Connect slider to current frame
        connect!(current_frame, slider.value)
        
        # Button row is 2 for stacks
        button_row = 2
        shortcuts_row = 3
    else
        # Skip frame navigation for single images
        button_row = 1
        shortcuts_row = 2
    end
    
    # Add toggle button for contrast mode
    contrast_btn = Button(controls_panel[button_row, 1:11], label="Local Contrast")
    
    # Button click handler
    on(contrast_btn.clicks) do _
        global_contrast[] = !global_contrast[]
        contrast_btn.label = global_contrast[] ? "Global Contrast" : "Local Contrast"
    end
    
    # Add a small gap between controls and image
    rowgap!(fig.layout, 10)
    
    # Add an info label about keyboard shortcuts
    if is_stack
        shortcuts_text = "Keyboard: [←/→] Change frame, [r] Reset view"
    else
        shortcuts_text = "Keyboard: [r] Reset view"
    end
    shortcuts_label = Label(controls_panel[shortcuts_row, 1:11], shortcuts_text)
    
    # Setup key bindings
    on(events(fig).keyboardbutton) do event
        if event.action == Keyboard.press || event.action == Keyboard.repeat
            # Left/right arrow keys for frame navigation (only for stacks)
            if is_stack && event.key == Keyboard.left && current_frame[] > 1
                # Just update the slider value - this handles both the current_frame and the slider position
                # as they are connected with the connect! function
                current_frame[] = current_frame[] - 1
                return true
            elseif is_stack && event.key == Keyboard.right && current_frame[] < n_frames
                # Just update the slider value - this handles both the current_frame and the slider position
                # as they are connected with the connect! function
                current_frame[] = current_frame[] + 1
                return true
            # 'r' key to reset the axis view (works for both 2D and 3D)
            elseif event.key == Keyboard.r
                # Reset to full view showing entire image
                autolimits!(ax)
                return true
            end
        end
        return false
    end
    
    # Return the figure
    display(fig)
    return fig
end