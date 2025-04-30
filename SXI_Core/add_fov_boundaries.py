#This is used in CMEM_Image to add the FOV boundaries to 3D viewing plots. 

def add_fov_boundaries(ax, xpos, ypos, zpos, color='k', lw=2):
    '''This will add the FOV boundaries in black/white to a 3D plot. 
    
    Parameters
    ----------
    ax - 3D axis object. 
    xpos - x positions.
    ypos - y positions.
    zpos - z positions. 
    color - linecolor.
    lw - linewidth. 
    
    '''
        
    #For corner pixels only. 
    ax.plot(xpos[0][0], ypos[0][0], zpos[0][0], color, lw=lw)
    ax.plot(xpos[0][-1], ypos[0][-1], zpos[0][-1], color, lw=lw)
    ax.plot(xpos[-1][0], ypos[-1][0], zpos[-1][0], color, lw=lw)
    ax.plot(xpos[-1][-1], ypos[-1][-1], zpos[-1][-1], color, lw=lw)
        
    #Join corners together. 
    ax.plot([xpos[0][0][-1],xpos[0][-1][-1]], [ypos[0][0][-1],ypos[0][-1][-1]], [zpos[0][0][-1],zpos[0][-1][-1]], color, lw=lw)
    ax.plot([xpos[0][-1][-1],xpos[-1][-1][-1]], [ypos[0][-1][-1],ypos[-1][-1][-1]], [zpos[0][-1][-1],zpos[-1][-1][-1]], color, lw=lw)
    ax.plot([xpos[-1][-1][-1],xpos[-1][0][-1]], [ypos[-1][-1][-1],ypos[-1][0][-1]], [zpos[-1][-1][-1],zpos[-1][0][-1]], color, lw=lw)
    ax.plot([xpos[-1][0][-1],xpos[0][0][-1]], [ypos[-1][0][-1],ypos[0][0][-1]], [zpos[-1][0][-1],zpos[0][0][-1]], color, lw=lw)
