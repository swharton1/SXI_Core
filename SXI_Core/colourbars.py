#Samuel Wharton 
#November 2017

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def make_colourbar(ax, limits, title='', cmap='seismic', levels=1000,
	fontsize=10, log = False, vertical=True, vpos='right', hpos='top'):
	'''Creates a colourbar in the axes. 
	
	Parameters
	---------------
	ax - Axis to turn into a colour bar. 
	limits - limits on colour bar. 
	title - Title on colour bar, including any units. Def is blank
	cmap - colour map. Default is 'seismic'
	levels - number of divisions on the colour scale. Def = 1000
	fontsize - def = 10 
	log - whether to plot number on axes as powers if it is a log scale. Def = False 
	
	Returns
	---------------
	ax
	
	'''
	
	levels = np.linspace(limits[0], limits[1], levels)
	
	cbar_x = list()
	cbar_y = list()
	
	for i in range(10):
		cbar_x.append(np.arange(10)+1)
		cbar_y.append(np.linspace(limits[0],limits[1],10))
		
	cbar_y = np.swapaxes(cbar_y,0,1)
	
	if vertical:
		ax.contourf(cbar_x, cbar_y, cbar_y, levels, cmap = mpl.cm.get_cmap(cmap))
		if vpos == 'right':
			ax.yaxis.tick_right()
			ax.yaxis.set_label_position("right")
		else:
			ax.yaxis.tick_left()
			ax.yaxis.set_label_position('left')
		ax.xaxis.set_major_formatter(plt.NullFormatter()) #No labels on x axis. 
		ax.xaxis.set_ticks_position('none') #No ticks on x axis. 
		ax.set_ylabel(title, fontsize=fontsize+2)
		for label in (ax.get_xticklabels() + ax.get_yticklabels()): 
			label.set_fontsize(fontsize)
	
		if log: 
			ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
			ticks = np.array(ax.get_yticks())#.astype(int)
			ax.set_yticklabels([r'$10^{'+str(int(i))+'}$' for i in ticks])
			
	else:
		ax.contourf(cbar_y, cbar_x, cbar_y, levels, cmap = mpl.cm.get_cmap(cmap))
		if hpos == 'top':
			ax.xaxis.tick_top()
			ax.xaxis.set_label_position("top")
		else:
			ax.xaxis.tick_bottom()
			ax.xaxis.set_label_position('bottom')
		ax.yaxis.set_major_formatter(plt.NullFormatter()) #No labels on x axis. 
		ax.yaxis.set_ticks_position('none') #No ticks on x axis. 
		ax.set_xlabel(title, fontsize=fontsize+2)
		for label in (ax.get_xticklabels() + ax.get_yticklabels()): 
			label.set_fontsize(fontsize)
		if log: 
			ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
			ticks = np.array(ax.get_xticks())#.astype(int)
			ax.set_xticklabels([r'$10^{'+str(i)+'}$' for i in ticks])
		
	
	
	
			
		#ax.set_yscale('log')
	#Matt's code for making a colourbar.
	#sm = plt.cm.ScalarMappable(cmap=mpl.cm.get_cmap(comap), 
		#norm=plt.Normalize(vmin=ymin, vmax=ymax))
		
	#sm._A = []
	#ax3 = plt.subplot2grid((1,10), (0,9))
	#ax3.colorbar(sm, shrink=0.9, pad = 0.01)
	#ax3.set_label(cbar_title, fontsize=fontsize+2)
	
	return ax  
	
def log_colourbar(ax, limits, title='',comap='seismic', levels=100,
	fontsize=10, log = False, vertical=True):
	'''Creates a colourbar in the axes. 
	
	Parameters
	---------------
	ax - Axis to turn into a colour bar. 
	limits - limits on colour bar. 
	title - Title on colour bar, including any units. Def is blank
	comap - colour map. Default is 'seismic'
	levels - number of divisions on the colour scale. Def = 1000
	fontsize - def = 10 
	log - whether to plot number on axes as powers if it is a log scale. Def = False 
	
	Returns
	---------------
	ax
	
	'''
	#levels = np.linspace(limits[0], limits[1], levels)
	#levels = np.logspace(0, limits[1], levels)
	
	cbar_x = list()
	cbar_y = list()
	
	for i in range(100):
		cbar_x.append(np.arange(100)+1)
		cbar_y.append(np.linspace(limits[0],limits[1],100))	
		
	cbar_y = np.swapaxes(cbar_y,0,1)
	cbar_z = np.log10(cbar_y) 
	
	#Need to remove -inf from array 
	for j,jval in enumerate(cbar_z[0]):
		if np.isneginf(jval):
			print ('Removing Negative Infinities')
			cbar_z[0][j] = cbar_z[1][j]
	
	ax.contourf(cbar_y, cbar_x, cbar_z, levels, cmap = mpl.cm.get_cmap(comap))
	ax.xaxis.tick_top()
	ax.xaxis.set_label_position("top")
	ax.yaxis.set_major_formatter(plt.NullFormatter()) #No labels on x axis. 
	ax.yaxis.set_ticks_position('none') #No ticks on x axis. 
	ax.set_xlabel(title, fontsize=fontsize+2)
	for label in (ax.get_xticklabels() + ax.get_yticklabels()): 
		label.set_fontsize(fontsize)
	#ax.set_xscale('log')
	#ticks = np.array(ax.get_xticks())#.astype(int)
	#ax.set_xticklabels([r'$10^{'+str(i)+'}$' for i in ticks])
		
	
	return ax 




