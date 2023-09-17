using SMLMSim
using Test

@testset "SMLMSim.jl" begin
    # Write your tests here.

    @testset "Interaction Diffusion" begin
        using MicroscopePSFs

        box_size = 1
        dt = 0.005
        density = 2
        state_history, args = SMLMSim.InteractionDiffusion.smoluchowski(;
            dt=dt, box_size=box_size, t_max=5.0, density=density, d_dimer=0.05)
        @test length(state_history.frames) == 1000
        
        pixelsize = 0.1
        pixels = Int64(round(box_size / pixelsize))
        psf = MicroscopePSFs.Airy2D(1.3, 0.6, pixelsize)
        camera = SMLMSim.IdealCamera(xpixels=pixels, ypixels=pixels, pixelsize=pixelsize)
        dimer_history = SMLMSim.InteractionDiffusion.get_dimers(state_history)
        dimer_stack = SMLMSim.gen_image_stack(psf, dimer_history, camera;
            photons=1000.0, bg=5.0, poissonnoise=true, frame_integration=10)
        @test size(dimer_stack) == (pixels, pixels, 100)
    end


end
