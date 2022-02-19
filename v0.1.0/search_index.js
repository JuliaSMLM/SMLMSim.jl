var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = SMLMSim","category":"page"},{"location":"#SMLMSim","page":"Home","title":"SMLMSim","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for SMLMSim.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [SMLMSim]","category":"page"},{"location":"#SMLMSim.CTMC","page":"Home","title":"SMLMSim.CTMC","text":"CTMC\n\nContinous Time Markov Chain    \n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Camera","page":"Home","title":"SMLMSim.Camera","text":"Camera\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.GenericFluor","page":"Home","title":"SMLMSim.GenericFluor","text":"GenericFluor\n\nDefines a fluorophore\n\nFields\n\nγ: photon emission rate in Hz, Default: 1e3\nq: state transision matrix. Default: q=[1.0]\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.IdealCamera","page":"Home","title":"SMLMSim.IdealCamera","text":"IdealCamera <: Camera\n\nA camera with no added noise. \n\n#Fields\n\npixelsize \nxpixels\nypixels\ngain\noffset\n\nIdealCamera(;\npixelsize=0.1,\nxpixels::Int=256,\nypixels::Int=256,\ngain=1.0,\noffset=0.0,\n)\n\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Images","page":"Home","title":"SMLMSim.Images","text":"Images\n\nOutput data type\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Molecule","page":"Home","title":"SMLMSim.Molecule","text":"Molecule\n\nPhotophysical properties of a molecule. \n\nThis is the most general type of luminecent or scattering single molecule.   Inherited types will defines the properties of a class of molecules. \n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Nmer2D","page":"Home","title":"SMLMSim.Nmer2D","text":"Nmer2D <: Pattern\n\nN molecules symmetricaly organized with diameter d    \n\nNmer2D(n::Int,d::AbstractFloat)\n\nNote d is in physical units (e.g. microns)\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Pattern","page":"Home","title":"SMLMSim.Pattern","text":"Pattern\n\nAbstract type for structured patterns of molecules    \n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.Point2D","page":"Home","title":"SMLMSim.Point2D","text":"Uniform <: Pattern\n\nGenerate data with uniform random distribution.    \n\nFields\n\nρ: molecule density\n\n\n\n\n\n","category":"type"},{"location":"#SMLMSim.getnext-Tuple{SMLMSim.CTMC, AbstractFloat}","page":"Home","title":"SMLMSim.getnext","text":"getnext(ctmc::CTMC,t::AbstractFloat)\n\nreturn the time and state of next transision\n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.getstate-Tuple{SMLMSim.CTMC, AbstractFloat}","page":"Home","title":"SMLMSim.getstate","text":"getstate(ctmc::CTMC,t::AbstractFloat)\n\nreturn the state at time t   \n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.intensitytrace-Tuple{SMLMSim.GenericFluor, Int64, Real}","page":"Home","title":"SMLMSim.intensitytrace","text":"intensitytrace(f::GenericFluor, nframes::Int, framerate::AbstractFloat;state1=1)\n\nCalculate an intensity trace.     \n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.kineticmodel-Tuple{SMLMData.SMLD2D, SMLMSim.Molecule, Int64, Real}","page":"Home","title":"SMLMSim.kineticmodel","text":"function kineticmodel(smd_true::SMLMData.SMLD2D,f::Molecule,nframes::Int,framerate::AbstractFloat;ndatasets::Int=1,minphotons=50.0)\n\ngenerate noise-free blinking model from smd_true\n\n\n\n\n\n","category":"method"},{"location":"#SMLMSim.uniform2D-Tuple{Any, SMLMSim.Pattern, Real, Real}","page":"Home","title":"SMLMSim.uniform2D","text":"function uniformPattern2D(ρ,xsize::AbstractFloat,ysize::AbstractFloat, p::Pattern)\n\ncreate true positions of molecules from uniformly randomly placed patterns\n\n\n\n\n\n","category":"method"}]
}
