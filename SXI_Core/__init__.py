#This is the init file for SXI_Core. 
#It will list all the modules in SXI_Core. 

#This is a list of general functions. 
from . import get_earth #Makes either a 2D or 3D representation of the earth 
from . import calc_pressures #Calculate dynamic and magnetic pressure 
from . import get_meridians #Calculates Y=0 and Z=0 planes in an emissivity cube
from . import colourbars #Makes colourbars
from . import coord_conv #Converts between spherical and cartesian coordinates
from . import sig_figs #Rounds to certain number of significant figures 

#This is a list of core functions for CMEM and Emissivity Cubes. 
from . import read_ppmlr #This reads in the original PPMLR ASCII cubes
from . import read_fits_cube #This reads in the new FITS format 
from . import set_initial_params #This sets the initial parameters for CMEM
from . import get_names_and_units #This gets the names and units of each CMEM parameter 
