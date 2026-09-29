# Generate camera images from SMLD and PSF

"""
    gen_images(smld::SMLD, psf::AbstractPSF; kwargs...) -> Tuple{Array{T, 3}, ImageInfo} where T<:Real

Generate camera images from SMLD data using the specified PSF model.

# Arguments
- `smld::SMLD`: Single molecule localization data container
- `psf::AbstractPSF`: Point spread function model

# Keyword arguments
- `dataset::Int=1`: Dataset number to use from SMLD
- `frames=nothing`: Specific frames to generate (default: all frames in smld.n_frames)
- `support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}}=Inf`: PSF support region size:
  - `Inf` (default): Calculate PSF over the entire image (most accurate but slowest). On the
    `GaussianPSF` fast path the kernel stops at 8σ + 1 px; the dropped mass is the Gaussian at 8σ,
    exp(-32) ≈ 1.3e-14 of the peak value
  - `Real`: Circular region with specified radius (in microns) around each emitter
  - `Tuple{<:Real,<:Real,<:Real,<:Real}`: Explicit region as (xmin, xmax, ymin, ymax) in microns
- `sampling::Int=2`: Supersampling factor for PSF integration (ignored for `GaussianPSF` on
  uniform square pixels, where the pixel integral is exact)
- `threaded::Bool=true`: Render frames in parallel
- `bg::Union{Real,AbstractArray{<:Real,3}}=0.0`: Background signal (photons per pixel): a scalar,
  or a stack of size `(height, width, length(frames))` added before noise. Any other size throws.
- `poisson_noise::Bool=false`: Apply Poisson noise only (for simple shot noise)
- `camera_noise::Bool=false`: Apply full camera noise model (requires SCMOSCamera)
  - For SCMOSCamera: applies QE, Poisson, read noise, gain, and offset
  - For IdealCamera: ignored (use poisson_noise instead)
- `rng::AbstractRNG=Random.default_rng()`: Random number generator for the noise draws

# Returns
- `Tuple{Array{T,3}, ImageInfo}`: (images, info)
    - images: 3D array of camera images with dimensions [height, width, num_frames]
    - info: ImageInfo containing timing and image statistics

# Rendering
A `GaussianPSF` on uniform square pixels is rendered with exact erf pixel integrals
(`render_gaussian!`), so the pixel values are the true integrals over each pixel. Any other
PSF is integrated with `MicroscopePSFs.integrate_pixels!` on the support window of each
emitter, giving the same values as `integrate_pixels` without allocating a full image per
emitter. Frames are independent of each other and are rendered in parallel when `threaded`;
rendering draws no random numbers, so the output does not depend on the thread count.

# Performance Note
For the `support` parameter, using a finite radius (typically 3-5× the PSF width)
provides a good balance between accuracy and performance. For example, with a PSF
width of 0.15μm, a support radius of 0.5-1.0μm is usually sufficient.
"""
function gen_images(smld::SMLD, psf::AbstractPSF;
          dataset::Int=1,
          frames=nothing,
          support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}}=Inf,
          sampling::Int=2,
          threaded::Bool=true,
          bg::Union{Real,AbstractArray{<:Real,3}}=0.0,
          poisson_noise::Bool=false,
          camera_noise::Bool=false,
          rng::AbstractRNG=Random.default_rng())

    start_time = time_ns()

    # Filter for the specified dataset
    dataset_smld = @filter(smld, dataset == dataset)

    # Determine frames to process
    if frames === nothing
        # Use all frames from 1 to n_frames
        frames = 1:smld.n_frames
    end

    # Get camera dimensions from SMLD
    camera = smld.camera
    width = length(camera.pixel_edges_x) - 1
    height = length(camera.pixel_edges_y) - 1
    nframes = length(frames)

    if !(bg isa Real) && size(bg) != (height, width, nframes)
        throw(DimensionMismatch("bg stack has size $(size(bg)), expected $((height, width, nframes))"))
    end

    # Determine the element type from emitters (if any)
    if isempty(dataset_smld.emitters)
        T = Float64 # default if no emitters
    else
        T = typeof(dataset_smld.emitters[1].photons)
    end

    # Pre-allocate output array
    images = zeros(T, height, width, nframes)

    # Bucket emitters by frame in one stable pass (within-frame order preserved)
    emitters = dataset_smld.emitters
    buckets = Dict{Int,Vector{eltype(emitters)}}()
    for e in emitters
        push!(get!(() -> eltype(emitters)[], buckets, e.frame), e)
    end
    empty_bucket = eltype(emitters)[]

    renderer = _frame_renderer(psf, camera, support, sampling, T)
    parallel = threaded && nframes > 1 && Threads.nthreads() > 1
    inner_threaded = threaded && !parallel

    function render_range(range)
        ws = _Workspace(T, height, width)
        for i in range
            frame_emitters = get(buckets, frames[i], empty_bucket)
            frame = view(images, :, :, i)
            isempty(frame_emitters) || renderer(frame, frame_emitters, ws, inner_threaded)
            # Add background to the accumulated frame
            if bg isa Real
                frame .+= bg
            else
                frame .+= view(bg, :, :, i)
            end
        end
    end

    if parallel
        chunk = cld(nframes, Threads.nthreads())
        tasks = [Threads.@spawn(render_range(r)) for r in Iterators.partition(1:nframes, chunk)]
        foreach(wait, tasks)
    else
        render_range(1:nframes)
    end

    # Track total photons
    n_photons_total = 0.0
    for frame_num in frames
        frame_emitters = get(buckets, frame_num, empty_bucket)
        isempty(frame_emitters) || (n_photons_total += sum(e.photons for e in frame_emitters))
    end

    # Apply Poisson noise if requested
    if poisson_noise
        poisson_noise!(rng, images)
    end

    # Apply camera noise if requested
    if camera_noise
        if camera isa SCMOSCamera
            # Apply full sCMOS noise model (QE, Poisson, read noise, gain, offset) to each frame
            for frame_idx in 1:size(images, 3)
                frame = @view images[:, :, frame_idx]
                scmos_noise!(rng, frame, camera)
            end
        else
            # For IdealCamera, camera_noise flag is ignored (IdealCamera is Poisson-only)
            @warn "Camera noise only supported for SCMOSCamera. Use poisson_noise=true for IdealCamera." maxlog=1
        end
    end

    elapsed_s = (time_ns() - start_time) / 1e9

    # Build ImageInfo
    info = ImageInfo(
        elapsed_s=elapsed_s,
        backend=:cpu,
        device_id=-1,
        frames_generated=nframes,
        n_photons_total=n_photons_total,
        output_size=(height, width, nframes)
    )

    return images, info
end

# Per-task scratch: Gaussian weights and the fallback's full-frame emitter image
mutable struct _Workspace{T}
    buf::RenderBuffer
    scratch::Matrix{T}
    height::Int
    width::Int
end

_Workspace(::Type{T}, height, width) where T = _Workspace{T}(RenderBuffer(8), Matrix{T}(undef, 0, 0), height, width)

function _scratch!(ws::_Workspace{T}) where T
    size(ws.scratch) == (ws.height, ws.width) || (ws.scratch = zeros(T, ws.height, ws.width))
    return ws.scratch
end

# Private copy of MicroscopePSFs' get_pixel_indices rule (integration_core.jl), for an
# emitter at (x, y): the (i_range, j_range) of pixels overlapping the support region.
function _support_window(edges_x::AbstractVector, edges_y::AbstractVector, x::Real, y::Real,
                         support::Union{Real,Tuple{<:Real,<:Real,<:Real,<:Real}})
    if support isa Real
        if isinf(support)
            return 1:length(edges_x)-1, 1:length(edges_y)-1
        end
        x_min = x - support
        x_max = x + support
        y_min = y - support
        y_max = y + support
    else
        x_min, x_max, y_min, y_max = support
    end

    i_min = searchsortedlast(edges_x, x_min)
    i_max = searchsortedfirst(edges_x, x_max)
    j_min = searchsortedlast(edges_y, y_min)
    j_max = searchsortedfirst(edges_y, y_max)

    if i_min >= length(edges_x) || i_max <= 1 ||
       j_min >= length(edges_y) || j_max <= 1 ||
       i_min >= i_max || j_min >= j_max
        return UnitRange{Int}(1:0), UnitRange{Int}(1:0)
    end

    i_min = max(1, i_min)
    i_max = min(length(edges_x)-1, max(i_min, i_max))
    j_min = max(1, j_min)
    j_max = min(length(edges_y)-1, max(j_min, j_max))

    return i_min:i_max, j_min:j_max
end

# Pixel pitch if the edges are uniform, else NaN
function _uniform_pitch(edges::AbstractVector)
    length(edges) < 2 && return NaN
    d = (edges[end] - edges[1]) / (length(edges) - 1)
    ok = all(k -> abs((edges[k+1] - edges[k]) - d) <= 1e-6 * abs(d), 1:length(edges)-1)
    return ok ? Float64(d) : NaN
end

# Returns render!(frame, frame_emitters, ws, inner_threaded) that accumulates one frame
function _frame_renderer(psf::AbstractPSF, camera, support, sampling, ::Type{T}) where T
    ex, ey = camera.pixel_edges_x, camera.pixel_edges_y
    if psf isa GaussianPSF
        px, py = _uniform_pitch(ex), _uniform_pitch(ey)
        if isfinite(px) && isfinite(py) && abs(px - py) <= 1e-9 * px
            σpx = psf.σ / px
            x0, y0 = Float64(ex[1]), Float64(ey[1])
            rinf = ceil(Int, 8σpx) + 1
            return function (frame, es, ws, inner_threaded)
                for e in es
                    u = (e.x - x0) / px
                    v = (e.y - y0) / px
                    if support isa Real && isinf(support)
                        render_gaussian!(frame, ws.buf, u, v, σpx, e.photons, rinf)
                    else
                        cols, rows = _support_window(ex, ey, e.x, e.y, support)
                        render_gaussian!(frame, ws.buf, u, v, σpx, e.photons, cols, rows)
                    end
                end
            end
        end
    end
    # Any other PSF: integrate_pixels! on the support window, same operations as integrate_pixels
    return function (frame, es, ws, inner_threaded)
        scratch = _scratch!(ws)
        for e in es
            i_range, j_range = _support_window(ex, ey, e.x, e.y, support)
            (isempty(i_range) || isempty(j_range)) && continue
            edges_x = view(ex, first(i_range):last(i_range)+1)
            edges_y = view(ey, first(j_range):last(j_range)+1)
            win = view(scratch, j_range, i_range)
            integrate_pixels!(win, psf, edges_x, edges_y, e; sampling=sampling, threaded=inner_threaded)
            view(frame, j_range, i_range) .+= win
        end
    end
end

"""
    gen_image(smld::SMLD, psf::AbstractPSF, frame::Int; kwargs...) -> Tuple{Matrix{T}, ImageInfo} where T<:Real

Generate a single camera image for a specific frame from SMLD data.
See `gen_images` for full documentation of parameters.

# Returns
- `Tuple{Matrix{T}, ImageInfo}`: (image, info)
    - image: 2D camera image as Matrix{T}
    - info: ImageInfo containing timing and image statistics
"""
function gen_image(smld::SMLD, psf::AbstractPSF, frame::Int; kwargs...)
    # Call gen_images with a single frame
    images, info = gen_images(smld, psf; frames=[frame], kwargs...)
    return images[:,:,1], info
end
