@testset "Camera Images" begin
    # Create some emitters for testing (using Emitter2DFit which has frame and other parameters)
    emitters = Vector{Emitter2DFit{Float64}}()
    
    # Create emitters with explicit frame numbers, IDs, and datasets
    push!(emitters, Emitter2DFit(5.0, 5.0, 1000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, 1, 1, 1))
    push!(emitters, Emitter2DFit(15.0, 15.0, 2000.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1, 1, 2, 2))
    push!(emitters, Emitter2DFit(5.0, 5.0, 800.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2, 1, 3, 3))
    
    # Create a camera for testing
    camera = IdealCamera(32, 32, 0.1)  # 32x32 pixels, 100 nm pixel size
    
    # Create a PSF model
    psf = GaussianPSF(0.13)  # 130 nm PSF width
    
    # Create basic SMLD container
    smld = BasicSMLD(emitters, camera, 2, 1)  # 2 frames, 1 dataset
    
    # Test gen_images for all frames
    images = gen_images(smld, psf)
    
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
    
    # Test single frame image generation
    image1 = gen_image(smld, psf, 1)
    
    # Check single frame dimensions
    @test size(image1) == (32, 32)
    
    # Check that image has signal
    @test maximum(image1) > 0
    
    # Test that the image from gen_images matches the image from gen_image
    @test image1 ≈ images[:, :, 1]
    
    # Test with background
    images_bg = gen_images(smld, psf, bg=10.0)
    
    # Check that background is added
    @test all(images_bg .>= images)
    @test images_bg[1, 1, 1] ≈ 10.0  # Only background at this position
    
    # Test with Poisson noise
    images_noisy = gen_images(smld, psf, bg=10.0, poisson_noise=true)
    
    # Check that noise is applied (images should differ)
    @test images_noisy != images_bg
    
    # The mean should be approximately preserved despite noise
    @test abs(mean(images_noisy) - mean(images_bg)) < 1.0
end