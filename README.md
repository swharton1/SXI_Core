# SXI_Core
This contains core functions used by a range of CMEM and SXI projects. 

It is not intended you would call this code directly. These are supporting functions for other CMEM projects. 

I've made it because I got fed up of copying the same scripts, and having to fix bugs in lots of copies, or just having keep lots of versions of things the same. 

You would need this module to run any of the following modules: 
CMEM/
CMEM_Image/
CMEM2/
CMEM_Ops/
Emissivity_Cubes/ 

THINGS YOU NEED TO DO: 

1. I recommend you put the top-level folder in your Pythonpath, say in your .bashrc file if you're linux. Mine looks like this: 
PYTHONPATH=$PYTHONPATH:~/Code/SXI_Core/

2. Adapt paths to emissivity cubes. read_fits_cube uses paths stored as environment variables to find the emissivity cubes. You need to set these in the init file. Change whichever of these variables you need too: PPMLR_PATH, OPENGGCM_PATH or BATSRUS_PATH. 
