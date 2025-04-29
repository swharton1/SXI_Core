#This is the init file for SXI_Core. 
#It will list all the modules in SXI_Core. 

#This is a list of general functions. 
from . import add_earth 
from . import calc_pressures 
from . import get_meridians
from . import colourbars 
from . import coord_conv 
from . import sig_figs

#This is a list of core functions for CMEM and Emissivity Cubes. 
from . import read_ppmlr #This reads in the original PPMLR ASCII cubes
from . import read_fits_cube #This reads in the new FITS format 
from . import set_initial_params #This sets the initial parameters for CMEM
from . import get_names_and_units #This gets the names and units of each CMEM parameter 
