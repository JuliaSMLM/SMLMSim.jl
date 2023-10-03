var documenterSearchIndex = {"docs":
[{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"CurrentModule = SMLMSim.InteractionDiffusion\nDocTestSetup = quote\n    using SMLMSim\nend","category":"page"},{"location":"diffusion/#SMLMSim.InteractionDiffusion","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion","text":"","category":"section"},{"location":"diffusion/#Overview","page":"Interaction-Diffusion","title":"Overview","text":"","category":"section"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"This module simulates interacting particles within a box. Monomers are consdiered point particles, and dimers and monomers kept at fixed separation. At each time step of the simulation, the following actions are taken:","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"If two free monomers are within the reaction radius, they are linked to form dimers.  \nPre-existing dimers are broken with a probabilty of p = k_mathrmoffdt\nMonomer position and dimer center-of-mass positions are updated with diffusion simulated as isotropic Brownian motion using their respective diffusion constant.  \nDimers undergo rotational diffusion. \nBoundary conditions are applied to the monomer position and dimer center-of-mass","category":"page"},{"location":"diffusion/#Tutorial","page":"Interaction-Diffusion","title":"Tutorial","text":"","category":"section"},{"location":"diffusion/#Running-the-Simulation","page":"Interaction-Diffusion","title":"Running the Simulation","text":"","category":"section"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"The main interface to the interaction-diffusion simulator is smoluchowski(). All arguments are keyword args and can be changed from default.  ","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"Simply run:","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"state_history, args = SMLMSim.smoluchowski()","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"The output is a MoleculeHistory structure and args, which is a structure that contains all the keyword input arguments that are used in the simulation.  ","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"The MoleculeHistory contains the position and monomer/dimer state  of each particle at each time frame.","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"The full list of options are in the docstring. ","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"SMLMSim.smoluchowski()","category":"page"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.smoluchowski-Tuple{}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.smoluchowski","text":"smoluchowski(; kwargs...)\n\nSimulate a system of molecules using the Smoluchowski model.\n\nKeyword Arguments\n\ndensity::Float64: the number density of the system (default: 1.0)\nbox_size::Float64: the size of the simulation box (default: 10.0)\ndiff_monomer::Float64: the diffusion coefficient of a monomer (default: 0.1)\ndiff_dimer::Float64: the diffusion coefficient of a dimer (default: 0.05)\ndiff_dimer_rot::Union{Nothing,Float64}: the rotational diffusion coefficient of a dimer (default: diff_dimer/d_dimer^2)\nk_off::Float64: the unbinding rate of a dimer (default: 0.2)\nr_react::Float64: the reaction radius for dimerization (default: 0.01)\ndt::Float64: the time step used in the simulation (default: 0.01)\nt_max::Float64: the maximum simulation time (default: 10.0)\nndims::Int64: the number of dimensions of the system (default: 2)\nd_dimer::Float64: Monomer seperation in the dimer (default: 0.05)\nboundary::String: the type of boundary conditions to apply (default: \"periodic\",  or \"reflecting\")\n\nReturns\n\nstate_history::MoleculeHistory: a MoleculeHistory object containing the states of the system at each time step\nargs::ArgsSmol: a struct containing the simulation parameters\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#Visualizing-the-Simulation","page":"Interaction-Diffusion","title":"Visualizing the Simulation","text":"","category":"section"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"To see positions of the particles in a single frame:","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"using SMLMSim\nstate_history, args = SMLMSim.smoluchowski()\nframenum = 100\nSMLMSim.show_frame(state_history,framenum,args)\nSMLMSim.show_frame(state_history,framenum,args,\"IDFrame.png\") # hide","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"<img src=\"IDFrame.png\" alt=\"Simulation Frame\" width=\"400\"/>","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"Monomers appear as blue dots. Dimers appear as red dots. ","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"To generate an mp4 movie of all frames:","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"using SMLMSim\nstate_history, args = SMLMSim.smoluchowski()\nSMLMSim.gen_movie(state_history,args; filename=\"defaultsim.mp4\")","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"<video width=\"320\" height=\"240\" controls>\n  <source src=\"defaultsim.mp4\" type=\"video/mp4\">\n  Your browser does not support the video tag.\n</video>","category":"page"},{"location":"diffusion/#Simulating-Microscope-Data","page":"Interaction-Diffusion","title":"Simulating Microscope Data","text":"","category":"section"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"The MoleculeHistory structure returned by smoluchowski() can  be use to simulate data as if it was collected in microscope.  This is implemented by the function gen_image_stack and requires a MicroscopePSFs.PSF and a SMLMSim.Camera.  ","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"A camera has a finite integration time that may be long compared to the dynamics of the diffusion simulation.  The  The example below shows the simulation of a system and the generation of data stack where each camera image is an integration over sub_sampling simulation frames. ","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"using SMLMSim\nusing MicroscopePSFs\nusing Images\n\ncamera_exposure_time = 0.05 # s\nbox_size = 10 # microns\nsub_sampling = 10\ndt = camera_exposure_time/sub_sampling # s\n\nstate_history, args = SMLMSim.InteractionDiffusion.smoluchowski(; \n    dt=dt, box_size=box_size);\n\n# Setup a Camera\npixelsize = 0.1\npixels = Int64(round(box_size/pixelsize))\ncamera = SMLMSim.IdealCamera(; xpixels=pixels, ypixels=pixels, pixelsize=pixelsize)\n\n# Setup a PSF\nna = 1.3\nwavelength = 0.6 # micron\npsf = MicroscopePSFs.Airy2D(na,wavelength,pixelsize)\n\nimage_stack = SMLMSim.gen_image_stack(psf, state_history, camera; \n    photons=1000.0, \n    bg=5.0, \n    poissonnoise=true, \n    frame_integration=sub_sampling);\n\n# Look at one frame\nim = Gray.(image_stack[:,:,1]./maximum(image_stack[:,:,1]))\n\nsave(\"simimageframe.png\",im) # hide","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"<img src=\"simimageframe.png\" alt=\"Image Frame\" width=\"400\"/>","category":"page"},{"location":"diffusion/#Extracting-Dimers","page":"Interaction-Diffusion","title":"Extracting Dimers","text":"","category":"section"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"We can extract the dimers into a MoleculeHistory structure, which can then be used with the above visualization and image generation tools. ","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"using SMLMSim\nusing MicroscopePSFs\nusing Images\n\ncamera_exposure_time = 0.05 # s\nbox_size = 10 # microns\nsub_sampling = 10\ndt = camera_exposure_time/sub_sampling # s\n\nstate_history, args = SMLMSim.InteractionDiffusion.smoluchowski(; \n    dt=dt, box_size=box_size);\n\n# Setup a Camera\npixelsize = 0.1\npixels = Int64(round(box_size/pixelsize))\ncamera = SMLMSim.IdealCamera(; xpixels=pixels, ypixels=pixels, pixelsize=pixelsize)\n\n# Setup a PSF\nna = 1.3\nwavelength = 0.6 # micron\npsf = MicroscopePSFs.Airy2D(na,wavelength,pixelsize)\n\nimage_stack = SMLMSim.gen_image_stack(psf, state_history, camera; \n    photons=1000.0, \n    bg=5.0, \n    poissonnoise=true, \n    frame_integration=sub_sampling);","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"dimer_history = SMLMSim.get_dimers(state_history)\n\ndimer_stack = SMLMSim.gen_image_stack(psf, dimer_history, camera; \n    photons=1000.0, bg=5.0, poissonnoise=true, frame_integration=10);\n\nframenum = length(state_history.frames)\nSMLMSim.show_frame(state_history,framenum,args)\nSMLMSim.show_frame(state_history,framenum,args,\"DimerFrame.png\") # hide\n\nall_image = Gray.(image_stack[:,:,end]./maximum(image_stack[:,:,end]))\n\ndimer_image = Gray.(dimer_stack[:,:,end]./maximum(dimer_stack[:,:,end]))\n\nsave(\"all_image.png\",all_image) # hide\nsave(\"dimer_image.png\",dimer_image) # hide","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"<div style=\"display: flex; flex-direction: row;\">\n    <img src=\"DimerFrame.png\" alt=\"DimerFrame\" width=\"300\" style=\"margin-right: 10px;\"/>\n    <img src=\"all_image.png\" alt=\"all_image\" width=\"300\" style=\"margin-right: 10px;\"/>\n    <img src=\"dimer_image.png\" alt=\"dimer_image\" width=\"300\"/>\n</div>","category":"page"},{"location":"diffusion/#API","page":"Interaction-Diffusion","title":"API","text":"","category":"section"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"SMLMSim.InteractionDiffusion","category":"page"},{"location":"diffusion/#SMLMSim.InteractionDiffusion","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion","text":"This module provides simulation tools for diffusion     and interaction between particles.\n\n\n\n\n\n","category":"module"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"Modules = [InteractionDiffusion]","category":"page"},{"location":"diffusion/","page":"Interaction-Diffusion","title":"Interaction-Diffusion","text":"Modules = [InteractionDiffusion]","category":"page"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.MoleculeFrame","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.MoleculeFrame","text":"MoleculeFrame(frame::Int64, molecules::Vector{<:AbstractOligomer})\n\nA struct representing a time frame of a simulation of a system of molecules.\n\nFields\n\nframe::Int64: The frame number of the simulation.\n'molecules::Vector{<:AbstractOligomer}': A vector of AbstractOligomer objects representing the state of the system at the current time step.\n\n\n\n\n\n","category":"type"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.MoleculeFrame-Tuple{Int64, Int64}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.MoleculeFrame","text":"MoleculeFrame(framenum::Int64, nmolecules::Int64)\n\nCreate a MoleculeFrame object with nmolecules molecules.\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.MoleculeHistory","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.MoleculeHistory","text":"MoleculeHistory(dt::Float64, frames::Vector{MoleculeFrame})\n\nA history of the state of a system of molecules over time.\n\nFields\n\ndt::Float64: The time step used in the simulation.\nframes::Vector{MoleculeFrame}: A vector of MoleculeFrame objects representing the state of the system at each time step.\n\n\n\n\n\n","category":"type"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.MoleculeHistory-Tuple{Float64, Int64, Int64}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.MoleculeHistory","text":"MoleculeHistory(dt::Float64, nframes::Int64, nmolecules::Int64)\n\nCreate a MoleculeHistory object with nframes frames, each containing nmolecules molecules.\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.Monomer","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.Monomer","text":"struct Monomer <: AbstractOligomer\n\nA struct representing a single monomer in an oligomer.\n\nFields\n\nx::Float64: the x-coordinate of the monomer's position\ny::Float64: the y-coordinate of the monomer's position\nz::Float64: the z-coordinate of the monomer's position\nstate::Int64: the state of the monomer (1 for monomer, 2 for dimer)\nlink::Union{Monomer,Nothing}: the monomer that this monomer is linked to, if any\nupdated::Bool: a flag indicating whether this monomer has been updated in the current iteration\n\n\n\n\n\n","category":"type"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.apply_boundary!-Tuple{SMLMSim.InteractionDiffusion.Monomer, SMLMSim.InteractionDiffusion.ArgsSmol}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.apply_boundary!","text":"apply_boundary!(mol::Monomer, args::ArgsSmol)\n\nApply boundary conditions to a single monomer in the system.\n\nArguments\n\nmol::Monomer: the monomer to apply boundary conditions to\nargs::ArgsSmol: a struct containing the simulation parameters\n\nReturns\n\nnothing\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.apply_boundary!-Tuple{SMLMSim.InteractionDiffusion.Monomer, SMLMSim.InteractionDiffusion.Monomer, SMLMSim.InteractionDiffusion.ArgsSmol}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.apply_boundary!","text":"apply_boundary!(mol1::Monomer, mol2::Monomer, args::ArgsSmol)\n\nApply boundary conditions to a dimer using center of mass.\n\nNote\n\nThis may result in molecules being placed outside of the box.\n\nArguments\n\nmol1::Monomer: the first monomer in the dimer to apply boundary conditions to\nmol2::Monomer: the second monomer in the dimer to apply boundary conditions to\nargs::ArgsSmol: a struct containing the simulation parameters\n\nReturns\n\nnothing\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.apply_boundary!-Tuple{Vector{<:SMLMSim.InteractionDiffusion.AbstractOligomer}, SMLMSim.InteractionDiffusion.ArgsSmol}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.apply_boundary!","text":"apply_boundary!(molecules::Vector{<:AbstractOligomer}, args::ArgsSmol)\n\nApply boundary conditions to all molecules in the system.\n\nArguments\n\nmolecules::Vector{<:AbstractOligomer}: a vector of AbstractOligomer objects representing the molecules in the system\nargs::ArgsSmol: a struct containing the simulation parameters\n\nReturns\n\nnothing\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.apply_boundary-Tuple{Float64, Float64, Float64, SMLMSim.InteractionDiffusion.ArgsSmol}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.apply_boundary","text":"apply_boundary(x::Float64, y::Float64, z::Float64, args::ArgsSmol)\n\nApply boundary conditions to a set of coordinates.\n\nArguments\n\nx::Float64: the x-coordinate of the point\ny::Float64: the y-coordinate of the point\nz::Float64: the z-coordinate of the point\nargs::ArgsSmol: a struct containing the simulation parameters\n\nReturns\n\nx::Float64: the x-coordinate of the point after applying boundary conditions\ny::Float64: the y-coordinate of the point after applying boundary conditions\nz::Float64: the z-coordinate of the point after applying boundary conditions\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.calc_r-Tuple{SMLMSim.InteractionDiffusion.Monomer, SMLMSim.InteractionDiffusion.Monomer}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.calc_r","text":"calc_r(mol1::Monomer, mol2::Monomer)\n\nCalculate the Euclidean distance between two monomers.\n\nArguments\n\nmol1::Monomer: the first monomer\nmol2::Monomer: the second monomer\n\nReturns\n\nr::Float64: the Euclidean distance between the two monomers\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.calc_θ-Tuple{SMLMSim.InteractionDiffusion.Monomer, SMLMSim.InteractionDiffusion.Monomer}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.calc_θ","text":"calc_θ(mol1::Monomer, mol2::Monomer)\n\nCalculate the angle between the z-axis and the vector connecting two monomers.\n\nArguments\n\nmol1::Monomer: the first monomer\nmol2::Monomer: the second monomer\n\nReturns\n\nθ::Float64: the angle between the z-axis and the vector connecting the two monomers, in radians\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.calc_ϕ-Tuple{SMLMSim.InteractionDiffusion.Monomer, SMLMSim.InteractionDiffusion.Monomer}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.calc_ϕ","text":"calc_ϕ(mol1::Monomer, mol2::Monomer)\n\nCalculate the azimuthal angle between two monomers.\n\nArguments\n\nmol1::Monomer: the first monomer\nmol2::Monomer: the second monomer\n\nReturns\n\nϕ::Float64: the azimuthal angle between the two monomers, in radians\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.dimerize!-Tuple{SMLMSim.InteractionDiffusion.Monomer, SMLMSim.InteractionDiffusion.Monomer, Float64}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.dimerize!","text":"dimerize!(mol1::Monomer, mol2::Monomer, distance::Float64)\n\nUpdate the state and position of two monomers to form a dimer.\n\nArguments\n\nmol1::Monomer: the first monomer to dimerize\nmol2::Monomer: the second monomer to dimerize\ndistance::Float64: the distance between the monomers in the dimer\n\nReturns\n\nnothing\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.gen_image-Tuple{MicroscopePSFs.PSF, SMLMSim.InteractionDiffusion.MoleculeHistory, SMLMSim.Camera, Int64}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.gen_image","text":"gen_image(psf::MicroscopePSFs.PSF, states::MoleculeHistory, camera::SMLMSim.Camera, framenum::Int64;\n    photons::Float64=1000.0, bg::Float64=5.0, poissonnoise::Bool=true)\n\nGenerate an image of a simulated system of molecules using a microscope PSF and camera.\n\nArguments\n\npsf::MicroscopePSFs.PSF: A MicroscopePSFs.PSF object representing the point spread function of the microscope.\nstates::MoleculeHistory: A MoleculeHistory object representing the history of the simulated system of molecules.\ncamera::SMLMSim.Camera: A SMLMSim.Camera object representing the camera used to capture the image.\nframenum::Int64: An integer representing the frame number of the simulation to generate an image for.\n\nOptional Arguments\n\nphotons::Float64=1000.0: The number of photons to simulate for each molecule.\nbg::Float64=5.0: Background photons per pixel.\npoissonnoise::Bool=true: A boolean indicating whether to corrupt with Poisson noise.\n\nReturns\n\nAn image of the simulated system of molecules as captured by the microscope.\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.gen_image_stack-Tuple{MicroscopePSFs.PSF, SMLMSim.InteractionDiffusion.MoleculeHistory, SMLMSim.Camera}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.gen_image_stack","text":"gen_image_stack(psf::MicroscopePSFs.PSF, states::MoleculeHistory, camera::SMLMSim.Camera;\n    photons::Float64=1000.0, bg::Float64=5.0, poissonnoise::Bool=true, frame_integration::Int64=1)\n\nGenerate a stack of images of a simulated system of molecules using a microscope PSF and camera.\n\nArguments\n\npsf::MicroscopePSFs.PSF: A MicroscopePSFs.PSF object representing the point spread function of the microscope.\nstates::MoleculeHistory: A MoleculeHistory object representing the history of the simulated system of molecules.\ncamera::SMLMSim.Camera: A SMLMSim.Camera object representing the camera used to capture the images.\n\nOptional Arguments\n\nphotons::Float64=1000.0: The number of photons emitted from each particle over the integration period.\nbg::Float64=5.0: Background photons per pixel in output images.\npoissonnoise::Bool=true: A boolean indicating whether to corrupt with Poisson noise.\nframe_integration::Int64=1: The number of simulation frames to integrate into each output image.\n\nReturns\n\nA stack of images of the simulated system of molecules as captured by the microscope.\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.gen_movie-Tuple{SMLMSim.InteractionDiffusion.MoleculeHistory, Any}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.gen_movie","text":"gen_movie(states::MoleculeStates, args::ArgsSmol; filename::String=\"smoluchowski.mp4\")\n\nGenerate an animation of the positions of all molecules in the system over time.\n\nArguments\n\nstates::MoleculeHistory: a MoleculeHistory object containing the states of the system at each time step\nargs::ArgsSmol: a struct containing the simulation parameters\nfilename::String: the name of the output file (default: \"smoluchowski.mp4\")\n\nReturns\n\nnothing\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.monomerize!-Tuple{SMLMSim.InteractionDiffusion.Monomer}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.monomerize!","text":"monomerize!(mol::Monomer)\n\nUpdate the state and links of a monomer and its linked molecule.\n\nArguments\n\nmol::Monomer: the monomer to monomerize\n\nReturns\n\nnothing\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.record_positions!-Tuple{Vector{<:SMLMSim.InteractionDiffusion.AbstractOligomer}, SMLMSim.InteractionDiffusion.MoleculeHistory, Int64}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.record_positions!","text":"record_positions!(molecules::Vector{<:AbstractOligomer}, state_history::MoleculeHistory, t::Int64)\n\nRecord the positions of all molecules in the system at a given time step.\n\nArguments\n\nmolecules::Vector{<:AbstractOligomer}: a vector of AbstractOligomer objects representing the molecules in the system\nstate_history::MoleculeHistory: a MoleculeStates object containing the states of the system at each time step\nt::Int64: the current time step\n\nReturns\n\nnothing\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.show_frame-Tuple{SMLMSim.InteractionDiffusion.MoleculeHistory, Any, Any, String}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.show_frame","text":"show_frame(states::MoleculeStates, framenum::Int64, args::ArgsSmol, filename::String)\n\nDisplay a scatter plot of the positions of all molecules in the system at a given time step and save it to a file.\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.show_frame-Tuple{SMLMSim.InteractionDiffusion.MoleculeHistory, Any, Any}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.show_frame","text":"show_frame(states::MoleculeStates, framenum::Int64, args::ArgsSmol)\n\nDisplay a scatter plot of the positions of all molecules in the system at a given time step.\n\nArguments\n\nstates::MoleculeHistory: a MoleculeHistory object containing the states of the system at each time step\nframenum::Int64: the time step to display\nargs::ArgsSmol: a struct containing the simulation parameters\n\nReturns\n\nf::Figure\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.update_position!-Tuple{SMLMSim.InteractionDiffusion.Monomer, SMLMSim.InteractionDiffusion.ArgsSmol}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.update_position!","text":"update_position!(mol::Monomer, args::ArgsSmol)\n\nUpdate the position of a single monomer using isotropic Brownian motion.\n\nArguments\n\nmol::Monomer: the monomer to update\nargs::ArgsSmol: a struct containing the simulation parameters\n\nReturns\n\nnothing\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.update_position!-Tuple{SMLMSim.InteractionDiffusion.Monomer, SMLMSim.InteractionDiffusion.Monomer, SMLMSim.InteractionDiffusion.ArgsSmol}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.update_position!","text":"update_position!(mol1::Monomer, mol2::Monomer, args::ArgsSmol)\n\nUpdate the position of a dimer using isotropic Brownian motion.\n\nThe center of mass of the dimer is updated using isotropic Brownian motion, and the orientation of the dimer is updated using rotational Brownian motion.\n\nArguments\n\nmol1::Monomer: the first monomer in the dimer to update\nmol2::Monomer: the second monomer in the dimer to update\nargs::ArgsSmol: a struct containing the simulation parameters\n\nReturns\n\nnothing\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.update_positions!-Tuple{Vector{<:SMLMSim.InteractionDiffusion.AbstractOligomer}, SMLMSim.InteractionDiffusion.ArgsSmol}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.update_positions!","text":"update_positions!(molecules::Vector{<:AbstractOligomer}, args::ArgsSmol)\n\nUpdate the positions of all molecules in the system according to the Smoluchowski model.\n\nArguments\n\nmolecules::Vector{<:AbstractOligomer}: a vector of AbstractOligomer objects representing the molecules in the system\nargs::ArgsSmol: a struct containing the simulation parameters\n\nReturns\n\nnothing\n\n\n\n\n\n","category":"method"},{"location":"diffusion/#SMLMSim.InteractionDiffusion.update_species!-Tuple{Vector{<:SMLMSim.InteractionDiffusion.AbstractOligomer}, SMLMSim.InteractionDiffusion.ArgsSmol}","page":"Interaction-Diffusion","title":"SMLMSim.InteractionDiffusion.update_species!","text":"update_species!(molecules::Vector{<:AbstractOligomer}, args::ArgsSmol)\n\nUpdate the state of a system of molecules according to the Smoluchowski model.\n\nArguments\n\nmolecules::Vector{<:AbstractOligomer}: a vector of AbstractOligomer objects representing the molecules in the system\nargs::ArgsSmol: a struct containing the simulation parameters\n\nReturns\n\nnothing\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = SMLMSim","category":"page"},{"location":"#SMLMSim","page":"Home","title":"SMLMSim","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SMLMSim.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SMLMSim]","category":"page"},{"location":"#SMLMData.SMLD2D-Union{Tuple{T}, Tuple{SMLMSim.Camera, Vector{T}, Vector{T}, Vector{T}, Vector{Int64}, Vector{Int64}, Vector{Int64}}} where T<:Real","page":"Home","title":"SMLMData.SMLD2D","text":"function SMLD2D(cam::Camera, y_microns::Vector{T}, x_microns::Vector{T}, photons::Vector{T},\nframenum::Vector{Int}, datasetnum::Vector{Int}, connectID::Vector{Int} ) where {T<:Real}\n\nGenerate SMLD2D using camera structure and physical units.\n\n\n\n\n\n","category":"method"},{"location":"#SMLMData.SMLD3D-Union{Tuple{T}, Tuple{SMLMSim.Camera, Vector{T}, Vector{T}, Vector{T}, Vector{T}, Vector{Int64}, Vector{Int64}, Vector{Int64}}} where T<:Real","page":"Home","title":"SMLMData.SMLD3D","text":"function SMLD3D(cam::Camera, y_microns::Vector{T}, x_microns::Vector{T}, z_microns::Vector{T}, photons::Vector{T},\nframenum::Vector{Int}, datasetnum::Vector{Int}, connectID::Vector{Int} ) where {T<:Real}\n\nGenerate SMLD3D using camera structure and physical units.\n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.CTMC","page":"Home","title":"SMLMSim.CTMC","text":"CTMC\n\nContinous Time Markov Chain    \n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Camera","page":"Home","title":"SMLMSim.Camera","text":"Camera\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.GenericFluor","page":"Home","title":"SMLMSim.GenericFluor","text":"GenericFluor\n\nDefines a fluorophore\n\nFields\n\nγ: photon emission rate in Hz, Default: 1e3\nq: state transision matrix. Default: q=[1.0]\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.IdealCamera","page":"Home","title":"SMLMSim.IdealCamera","text":"IdealCamera <: Camera\n\nA camera with no added noise. \n\nFields\n\npixelsize \nxpixels\nypixels\ngain\noffset\n\nIdealCamera(;\npixelsize=0.1,\nxpixels::Int=256,\nypixels::Int=256,\ngain=1.0,\noffset=0.0,\n)\n\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Images","page":"Home","title":"SMLMSim.Images","text":"Images\n\nOutput data type\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Line2D","page":"Home","title":"SMLMSim.Line2D","text":"Line2D <: Pattern2D\n\nPoints with uniform random distribution between 2 endpoints.    \n\nLine2D(;λ::AbstractFloat=10.0, endpoints=[(-1.0,0.0),(1.0,0.0)])\n\nFields\n\nλ: linear molecule density\n'endpoints': Vector of Tuple \n'n': Numbor of Points = 1\n'x': X position\n'y': Y position\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Molecule","page":"Home","title":"SMLMSim.Molecule","text":"Molecule\n\nPhotophysical properties of a molecule. \n\nThis is the most general type of luminecent or scattering single molecule.   Inherited types will defines the properties of a class of molecules. \n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Nmer2D","page":"Home","title":"SMLMSim.Nmer2D","text":"Nmer2D <: Pattern2D\n\nN molecules symmetricaly organized around a circle with diameter d    \n\nNmer2D(;n::Int=8, d::AbstractFloat=.1)\n\nFields\n\n'n': Numbor of Points = 1\n'd': Diameter\n'x': X position\n'y': Y position\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Nmer3D","page":"Home","title":"SMLMSim.Nmer3D","text":"Nmer3D <: Pattern3D\n\nN molecules symmetricaly organized around a circle with diameter d    \n\nNmer3D(;n::Int=8, d::AbstractFloat=.1)\n\nFields\n\n'n': Numbor of Points = 1\n'd': Diameter\n'x': X position\n'y': Y position\n'z': Z position\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Pattern","page":"Home","title":"SMLMSim.Pattern","text":"Pattern\n\nAbstract type for structured patterns of molecules    \n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Pattern2D","page":"Home","title":"SMLMSim.Pattern2D","text":"Pattern2D\n\nAbstract type for structured patterns of molecules    \n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Pattern3D","page":"Home","title":"SMLMSim.Pattern3D","text":"Pattern3D\n\nAbstract type for structured patterns of molecules    \n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Point2D","page":"Home","title":"SMLMSim.Point2D","text":"Point2D <: Pattern2D\n\nA single 2D point.\n\nPoint2D() = new(1, [0.0], [0.0])\n\nFields\n\n'n': Numbor of Points = 1\n'x': X position\n'y': Y position\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Point3D","page":"Home","title":"SMLMSim.Point3D","text":"Point3D <: Pattern3D\n\nA single 3D point.\n\nPoint3D() = new(1, [0.0], [0.0],[0.0])\n\nFields\n\n'n': Numbor of Points = 1\n'x': X position\n'y': Y position\n'z': Z position\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.getnext-Tuple{SMLMSim.CTMC, AbstractFloat}","page":"Home","title":"SMLMSim.getnext","text":"getnext(ctmc::CTMC,t::AbstractFloat)\n\nreturn the time and state of next transision\n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.getstate-Tuple{SMLMSim.CTMC, AbstractFloat}","page":"Home","title":"SMLMSim.getstate","text":"getstate(ctmc::CTMC,t::AbstractFloat)\n\nreturn the state at time t   \n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.intensitytrace-Tuple{SMLMSim.GenericFluor, Int64, Real}","page":"Home","title":"SMLMSim.intensitytrace","text":"intensitytrace(f::GenericFluor, nframes::Int, framerate::AbstractFloat;state1=1)\n\nCalculate an intensity trace.     \n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.kineticmodel","page":"Home","title":"SMLMSim.kineticmodel","text":"function kineticmodel(smd_true::SMLMData.SMLD,f::Molecule,nframes::Int,framerate::AbstractFloat;ndatasets::Int=1,minphotons=50.0)\n\ngenerate noise-free blinking model from smd_true\n\n\n\n\n\n","category":"function"},{"location":"#SMLMSim.noise","page":"Home","title":"SMLMSim.noise","text":"noise(smd_model::SMLMData.SMLD2D,σ_psf::AbstractFloat)\n\nAdd zero mean Gaussian noise to coordinates with σ = σ_pdf/sqrt(photons) \n\n\n\n\n\n","category":"function"},{"location":"#SMLMSim.noise-Tuple{SMLMData.SMLD3D, Vector{<:AbstractFloat}}","page":"Home","title":"SMLMSim.noise","text":"noise(smd_model::SMLMData.SMLD3D,σ_psf::Vector{<:AbstractFloat})\n\n3D data requries σ_psf = [σ_x,σ_y,σ_z]\n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.place2D-Tuple{Vector{<:SMLMSim.Pattern2D}, Real, Real}","page":"Home","title":"SMLMSim.place2D","text":"place2D(p::Vector{Pattern}, xsize::Real, ysize::Real)\n\nPlace the input patterns and return SMLD2D.\n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.rotate!-Tuple{SMLMSim.Pattern, Array{Real}}","page":"Home","title":"SMLMSim.rotate!","text":"Rotate a Pattern in 3D by premultiplying with the rotation matrix r.\n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.rotate!-Tuple{SMLMSim.Pattern, Real, Real, Real}","page":"Home","title":"SMLMSim.rotate!","text":"Rotate a Pattern in 3D by the improper Euler angles \\alpha \\beta \\gamma.\n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.rotate!-Tuple{SMLMSim.Pattern, Real}","page":"Home","title":"SMLMSim.rotate!","text":"rotate!(p::Pattern,θ::Real)\n\nRotate a Pattern in 2D by \\theta radians.\n\nBoth molecule positions and reference positions are rotated (e.g. endpoints of a line)\n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.sim-Tuple{}","page":"Home","title":"SMLMSim.sim","text":"sim(;\n    ρ=1.0,\n    σ_PSF=.13,\n    minphotons=50,\n    ndatasets=10,\n    nframes=1000,\n    framerate=50.0, \n    pattern::Pattern=Nmer2D(),\n    molecule::Molecule=GenericFluor(;q=[0 50; 1e-2 0]),\n    camera::Camera=IdealCamera(),\n    zrange::Vector{<:Real}=[-1.0, 1.0]\n)\n\nGenerate SMLD using simulation parmeters.      \n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.uniform2D-Tuple{Any, SMLMSim.Pattern2D, Real, Real}","page":"Home","title":"SMLMSim.uniform2D","text":"uniform2D(ρ, p::Pattern2D, xsize::AbstractFloat,ysize::AbstractFloat)\n\nCreate positions of molecules from uniformly randomly placed and rotated patterns.\n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.uniform2D-Tuple{Vector{<:SMLMSim.Pattern2D}, Real, Real, Real}","page":"Home","title":"SMLMSim.uniform2D","text":"uniform2D(p::Vector{Pattern}, xsize::Real, ysize::Real,θ::Real)\n\nRandomly place the input patterns with fixed rotation \\theta.\n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.uniform2D-Tuple{Vector{SMLMSim.Pattern2D}, Real, Real}","page":"Home","title":"SMLMSim.uniform2D","text":"uniform2D(p::Vector{Pattern2D}, xsize::Real, ysize::Real)\n\nRandomly place and rotate the input patterns.\n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.uniform3D-Tuple{Any, SMLMSim.Pattern3D, Real, Real}","page":"Home","title":"SMLMSim.uniform3D","text":"function uniform3D(ρ,p::Pattern, xsize::Real,ysize::Real; zrange::Vector{<:Real}=[-1.0,1.0])\n\nCreate positions of molecules from uniformly randomly placed and rotated patterns.\n\n!!! Note ρ is 2D density. 3D density is ρ/(zrange[2]-zrange[1]). \n\n\n\n\n\n","category":"method"}]
}
