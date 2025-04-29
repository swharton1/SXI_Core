#THIS WILL CONTAINS FUNCTIONS TO CONVERT BETWEEN COORDINATE SYSTEMS. 
import numpy as np
#Cartesian-SHUE coordinates. 

def convert_xyz_to_shue_coords(x, y, z):
    '''This will convert the x,y,z coordinates to those used in the Shue model 
     of the magnetopause and bowshock. 

    Parameters
    ---------
    x, y, z - now 3D. Must be entered as numpy arrays.  

    Returns
    -------
    r, theta (rad) and phi (rad)
    '''

   
    # r 
    r = (x**2 + y**2 + z**2)**0.5
       
    # theta - only calc. where coordinate singularities won't occur. 
    theta = np.zeros(r.shape)
    i = np.where(r != 0)
    theta[i] =  np.arccos(x[i]/r[i])

    # phi - only calc. where coordinate singularities won't occur. 
    phi = np.zeros(r.shape)
    phi[i] = np.arctan2(z[i], y[i]) 
     
    return r, theta, phi
    
def convert_shue_to_xyz_coords(r, theta, phi):
    '''This will convert the Shue coordinates back to xyz coordinates. 
        
    Parameters
    ----------
    r, theta (rad), phi (rad)
        
    Returns
    -------
    x,y,z
    '''

    x = r*np.cos(theta)
    y = r*np.sin(theta)*np.cos(phi)
    z = r*np.sin(theta)*np.sin(phi)

    return x,y,z 

#Cartesian-Aberrated coordinates. 

def convert_xyz_to_aberrated_shue_coords(x, y, z, dgamma, ddelta): 
    '''This will convert the x,y,z coordinates to those used in the Shue model 
     of the magnetopause and bowshock. 

    Parameters
    ---------
    x, y, z - now 3D. Must be entered as numpy arrays.  
    dgamma - rotation in the horizontal plane (rad)
    ddelta - rotation in the vertical plane (rad) 
    
    Returns
    -------
    r, theta (rad) and phi (rad)
    '''

    #Horizontal rotation. 
    x2 = x*np.cos(-dgamma) - y*np.sin(-dgamma)
    y2 = x*np.sin(-dgamma) + y*np.cos(-dgamma) 
    z2 = z 
    
    #Vertical rotation. 
    x3 = x2*np.cos(-ddelta) - z2*np.sin(-ddelta)
    y3 = y2 
    z3 = x2*np.sin(-ddelta) + z2*np.cos(-ddelta) 


    # r 
    r = (x3**2 + y3**2 + z3**2)**0.5
       
    # theta - only calc. where coordinate singularities won't occur. 
    theta = np.zeros(r.shape)
    i = np.where(r != 0)
    theta[i] =  np.arccos(x3[i]/r[i])

    # phi - only calc. where coordinate singularities won't occur. 
    phi = np.zeros(r.shape)
    phi[i] = np.arctan2(z3[i], y3[i]) 
     
    return r, theta, phi


