using Revise
using SMLMSim
using CairoMakie
using MicroscopePSFs
using Images
using ImageView

# run the smoluchowski simulation
box_size = 10
dt = .005
density = 2
state_history, args = SMLMSim.InteractionDiffusion.smoluchowski(; 
    dt=dt, box_size=box_size, t_max = 5.0, density = density, d_dimer = 0.05);

SMLMSim.show_frame(state_history, 1, args)
pixelsize = 0.1
pixels = Int64(round(box_size/pixelsize))
psf = MicroscopePSFs.Airy2D(1.3,0.6,pixelsize)
camera = SMLMSim.IdealCamera(xpixels=pixels, ypixels=pixels,pixelsize=pixelsize)

image_frame = SMLMSim.gen_image(psf, state_history, camera, 200; photons=1000.0, bg=5.0, poissonnoise=false)
Gray.(image_frame./maximum(image_frame))

image_stack = SMLMSim.gen_image_stack(psf, state_history, camera; 
    photons=1000.0, bg=5.0, poissonnoise=true, frame_integration=10);

imshow(image_stack)

SMLMSim.gen_movie(state_history, args; filename="smoluchowski.mp4")

# dimers 
dimer_history = SMLMSim.InteractionDiffusion.get_dimers(state_history)
dimer_stack = SMLMSim.gen_image_stack(psf, dimer_history, camera; 
    photons=1000.0, bg=5.0, poissonnoise=true, frame_integration=10);

imshow(dimer_stack)


