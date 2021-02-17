## This script sets plot parameters for the figures
# use exec(open("figure_layout.py").read()) in the scipt you want it in

from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib as mpl
 
sns.set_style("ticks")

#set font sizes

plt.rcParams['font.family'] = 'sans-serif' #'serif' (e.g., Times), 'sans-serif' (e.g., Helvetica), 'cursive' (e.g., Zapf-Chancery), 'fantasy' (e.g., Western), and 'monospace' (e.g., Courier).
plt.rcParams['font.serif'] = ['Helvetica']
plt.rcParams['mathtext.fontset'] = 'cm'	#this ensures you use computer modern as mathtext
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 11
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize']= 6

### Legend
#plt.rcParams['legend.frameon']='False'
plt.rcParams['legend.handletextpad']=0.5	#pad between legend handle and text
#plt.rcParams['legend.handlelength']=1.5 	#length of handles
#plt.rcParams['legend.fancybox']=True 	  	#round box
plt.rcParams['legend.shadow']=False 	 	#cast shadow

### Padding
#plt.rcParams['axes.titlepad'] = 6

##Ticks
plt.rcParams['xtick.direction']= 'in'
plt.rcParams['ytick.direction']= 'in'
plt.rcParams['xtick.major.width']= 1.2		#tick width
plt.rcParams['ytick.major.width']= 1.2
plt.rcParams['xtick.major.size']= 4		#tick length
plt.rcParams['ytick.major.size']= 4
plt.rcParams['xtick.major.pad']= 4	 	  	# distance to major tick label in points
plt.rcParams['ytick.major.pad']= 4
#plt.rcParams['xtick.minor.visible']= True 	#set visible
#plt.rcParams['ytick.minor.visible']= True 	#set visible
plt.rcParams['axes.formatter.useoffset'] = False 	#don't label ticks with respect to offset

##axes
plt.rcParams['axes.facecolor']= 'white'
plt.rcParams['figure.facecolor']= 'white'

##lines and markers
plt.rcParams['lines.linewidth']= 2#old = 2.5
plt.rcParams['axes.linewidth']= 1.2
plt.rcParams['lines.markersize']=3
#plt.rcParams['lines.markeredgewidth']=0.001     # the line width around the marker symbol

## colormap for figures
plt.rcParams['image.cmap']="viridis"
# plt.rcParams['image.cmap']="gnuplot"
# I'm overriding this cmap by a custom one in the script

###saving
plt.rcParams['savefig.transparent']= True
plt.rcParams['savefig.bbox']= 'tight'	 	# 'tight' or 'standard'.
plt.rcParams['savefig.pad_inches']= 0.04	# Padding to be used when bbox is set to 'tight'
plt.rcParams['savefig.dpi']= 300	# Padding to be used when bbox is set to 'tight'
#plt.rcParams['savefig.format']= 'pdf'	 	#standard format


## kwargs for contour and mode profile plots
contourkwargs = { 'linewidths':1, 'linestyles':'solid', 'colors':'white'}
fieldkwargs = {'vmin':0, 'vmax':1, 'cmap':'viridis'}
logfieldkwargs = {'vmin':0.05, 'vmax':1, 'cmap':'viridis', 'norm': LogNorm()}
subfiglabelkwargs = {'fontsize':plt.rcParams['axes.labelsize'], 'bbox':dict(facecolor='white', pad=2)}
 
temp = mpl.cm.viridis(np.arange(256))
mycols = temp[np.array([0, 120, 180, 220, 255])]
mycols[1] = np.array([80,178,208, 255])/255
mycols[3] = np.array([52,203,64, 255])/255

sns.set_palette(mycols)



#sns.palplot(mycols)
#sns.set_palette(mycols)












