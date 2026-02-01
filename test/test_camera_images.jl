@testset "Camera Images" begin
    # Create some emitters for testing (using Emitter2DFit which has frame and other parameters)
    emitters = Vector{Emitter2DFit{Float64}}()

    # Create emitters with explicit frame numbers, IDs, and datasets
    push!(emitters, Emitter2DFit{Float64}(5.0, 5.0, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0; frame=1, dataset=1, track_id=1, id=1))
    push!(emitters, Emitter2DFit{Float64}(15.0, 15.0, 2000.0, 0.0, 0.0, 0.0, 0.0, 0.0; frame=1, dataset=1, track_id=2, id=2))
    push!(emitters, Emitter2DFit{Float64}(5.0, 5.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0; frame=2, dataset=1, track_id=3, id=3))

    # Create a camera for testing
    camera = IdealCamera(32, 32, 0.1)  # 32x32 pixels, 100 nm pixel size

    # Create a PSF model
    psf = GaussianPSF(0.13)  # 130 nm PSF width

    # Create basic SMLD container
    smld = BasicSMLD(emitters, camera, 2, 1)  # 2 frames, 1 dataset

    # Test gen_images for all frames - now returns (images, info) tuple
    images, info = gen_images(smld, psf)

    # Check that we got ImageInfo
    @test isa(info, ImageInfo)
    @test info.elapsed_ns > 0
    @test info.backend == :cpu
    @test info.device_id == -1
    @test info.frames_generated == 2
    @test info.output_size == (32, 32, 2)

    # Check dimensions
    @test size(images) == (32, 32, 2)  # 2 frames

    # First locate which pixels actually have intensity
    # For frame 1
    max_val_frame1, max_idx_frame1 = findmax(images[:,:,1])
    # For frame 2
    max_val_frame2, max_idx_frame2 = findmax(images[:,:,2])

    # Check that we have intensity at the maximum positions
    @test max_val_frame1 > 0  # First frame has signal
    @test max_val_frame2 > 0  # Second frame has signal

    # Check that the frames have some signal
    @test max_val_frame1 > 0
    @test max_val_frame2 > 0

    # Test single frame image generation - now returns (image, info) tuple
    image1, img_info = gen_image(smld, psf, 1)

    # Check ImageInfo for single frame
    @test isa(img_info, ImageInfo)
    @test img_info.frames_generated == 1
    @test img_info.output_size == (32, 32, 1)

    # Check single frame dimensions
    @test size(image1) == (32, 32)

    # Check that image has signal
    @test maximum(image1) > 0

    # Check that the image from gen_images matches the image from gen_image
    @test image1 ≈ images[:, :, 1]

    # Test with background
    images_bg, info_bg = gen_images(smld, psf, bg=10.0)

    # Check that background is added
    @test all(images_bg .>= images)
    @test images_bg[1, 1, 1] ≈ 10.0  # Only background at this position

    # Test with Poisson noise
    images_noisy, info_noisy = gen_images(smld, psf, bg=10.0, poisson_noise=true)

    # Check that noise is applied (images should differ)
    @test images_noisy != images_bg

    # The mean should be approximately preserved despite noise
    @test abs(mean(images_noisy) - mean(images_bg)) < 1.0

    @testset "ImageInfo struct" begin
        # Test ImageInfo fields
        images, info = gen_images(smld, psf)

        @test info.elapsed_ns > 0
        @test info.backend == :cpu
        @test info.device_id == -1
        @test info.frames_generated == 2
        @test info.n_photons_total > 0  # Should have some photons
        @test info.output_size == (32, 32, 2)
    end
end

@testset "SCMOSCamera Integration" begin
    # Create emitters for testing
    emitters = Vector{Emitter2DFit{Float64}}()
    push!(emitters, Emitter2DFit{Float64}(5.0, 5.0, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0; frame=1, dataset=1, track_id=1, id=1))
    push!(emitters, Emitter2DFit{Float64}(15.0, 15.0, 2000.0, 0.0, 0.0, 0.0, 0.0, 0.0; frame=1, dataset=1, track_id=2, id=2))

    # Create an sCMOS camera (32x32 pixels, 100 nm pixel size, 1.6 e⁻ read noise)
    camera_scmos = SCMOSCamera(32, 32, 0.1, 1.6)

    # Create a PSF model
    psf = GaussianPSF(0.13)

    # Create SMLD container with sCMOS camera
    smld = BasicSMLD(emitters, camera_scmos, 1, 1)

    # Test gen_images with camera_noise - now returns (images, info) tuple
    images_camera_noise, info = gen_images(smld, psf, camera_noise=true)

    # Check dimensions
    @test size(images_camera_noise) == (32, 32, 1)

    # Check that we have signal
    @test maximum(images_camera_noise) > 0

    # Check ImageInfo
    @test isa(info, ImageInfo)
    @test info.frames_generated == 1

    # Test that camera noise produces different results than clean image
    images_clean, _ = gen_images(smld, psf)
    @test images_camera_noise != images_clean

    # Test scmos_noise function directly
    clean_image = ones(32, 32) * 100.0
    noisy_scmos = scmos_noise(clean_image, camera_scmos)

    # Check dimensions preserved
    @test size(noisy_scmos) == (32, 32)

    # Check that noise was applied (should be different from input)
    @test noisy_scmos != clean_image

    # Test scmos_noise! in-place version
    test_image = ones(32, 32) * 100.0
    original_copy = copy(test_image)
    scmos_noise!(test_image, camera_scmos)

    # Check that image was modified in-place
    @test test_image != original_copy

    # Test with background and camera noise
    images_bg_camera, _ = gen_images(smld, psf, bg=10.0, camera_noise=true)
    @test size(images_bg_camera) == (32, 32, 1)
    # Note: With default offset=0 and readnoise, some pixels can be slightly negative
    # This is realistic - real cameras can have negative ADU values before dark frame subtraction
    @test minimum(images_bg_camera) > -10.0  # Check for reasonable range, not strictly >= 0

    # Test warning for IdealCamera with camera_noise
    camera_ideal = IdealCamera(32, 32, 0.1)
    smld_ideal = BasicSMLD(emitters, camera_ideal, 1, 1)

    # This should warn but not error
    images_ideal_warn, _ = gen_images(smld_ideal, psf, camera_noise=true)
    @test size(images_ideal_warn) == (32, 32, 1)
end

@testset "Diffusion with SCMOSCamera" begin
    # Test that diffusion simulation accepts camera parameter
    params = DiffusionSMLMParams(
        density = 0.5,
        box_size = 5.0,
        dt = 0.01,
        t_max = 0.1,
        camera_framerate = 10.0,
        camera_exposure = 0.05
    )

    # Create an sCMOS camera
    pixel_size = 0.1
    n_pixels = ceil(Int, params.box_size / pixel_size)
    camera_scmos = SCMOSCamera(n_pixels, n_pixels, pixel_size, 1.6)

    # Run simulation with sCMOS camera - now returns (smld, info) tuple
    smld_scmos, info = simulate(params; camera=camera_scmos, override_count=5)

    # Check that the camera is correctly set
    @test smld_scmos.camera isa SCMOSCamera
    @test !isempty(smld_scmos.emitters)

    # Check SimInfo
    @test isa(info, SimInfo)
    @test info.elapsed_ns > 0

    # Test default behavior (should create IdealCamera) - now returns (smld, info) tuple
    smld_default, _ = simulate(params; override_count=5)
    @test smld_default.camera isa IdealCamera
end
