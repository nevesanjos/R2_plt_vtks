import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import glob
from scipy.interpolate import griddata
import os
from shutil import copyfile
import re




# Define the axes
fig, ax = plt.subplots(1, 1, figsize=(14, 6))
ax2 = fig.add_axes([0.93, 0.1, 0.02, 0.8])
ax.set_xlim(0,60)
ax.set_ylim(-15,15)
ax.set_xlabel("Distance [$m$]",size=12)
ax.set_ylabel("Elevation [$m$]",size=12)
ax2.set_ylabel('Resistivity [$\Omega m$]',size=12)
ax2.set_title('$\log_{10}$',size=12)


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Define the colormap
number_of_color = 2**3
cmap = plt.cm.gist_earth  # define the colormap
#cmaplist = [cmap(i) for i in np.linspace(0,cmap.N+1,number_of_color)]
cmaplist = [cmap(i) for i in range(cmap.N)]
#cmaplist[0] = (.5, .5, .5, 1.0)
cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)

print(cmap.N/number_of_color)

# Working with directores and files
path_d = "../Example_files/"
ref = path_d + "invdir/ref/f001_res.vtk" 
copyfile(ref,path_d + "invdir/" + "f000_res.vtk") 
files = [f for f in glob.glob(path_d + "invdir/f???_res.vtk")]
files.sort()
data = [d for d in glob.glob(path_d + "data/*.tx0")]
data.sort()

if not os.path.exists(path_d + "pngs"):
        os.makedirs(path_d + "pngs")
#os.mkdir(path_+"pngs")

pngs = [p for p in glob.glob(path_d + "pngs/*.png")]

# define the global min and max values 
mm = np.zeros(2)
mm[0] = 10
for l in files:
    print(l)
    with open(l,'r') as gh:
        a = ' '
        while 'Resistivity(log10)' not in a:
            a = gh.readline()
        trash = gh.readline()
        d = gh.readline().split()
        if(mm[0] > float(min(d))):
            mm[0] = float(min(d))
        if(mm[1] < float(max(d))):
            mm[1] = float(max(d))            
gh.close()

# palete color
bounds = np.linspace(mm[0], mm[1], cmap.N)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=[np.linspace(np.around(mm[0]),np.around(mm[1]),5)], boundaries=bounds)

for w,f in enumerate(files):
    print(data[w])
    if(str(f) not in str(pngs)):
        print("Working on", f)
       
        with open(f,'r') as fh:
            a = ' '
            while 'POINTS' not in a:
                a = (fh.readline().strip())
            string = re.sub('double',' ',re.sub('POINTS',' ', a))
            nnodal = int(string)
            Nodal = []
            Elem = []
            rho_truelog = []
            rho = []
            S = []
            dump = fh.readlines()
            for i, line in enumerate(dump):
                if(i < nnodal):
                    Nodal.append(line.split())
                elif (i == nnodal):
                    n = line.split()
                    print(n[1])
                elif (i < int(n[1])+nnodal+1):
                    Elem.append(line.split())
                elif(i < int(n[1])+nnodal+7):
                    f = line.split()
                elif(i < int(n[1])+nnodal+8):
                    rho = line.split()
                elif(i < int(n[1])+nnodal+10):
                    f = line.split()
                elif(i < int(n[1])+nnodal+11):
                    rho_truelog = line.split()
                elif(i < int(n[1])+nnodal+13):
                    f = line.split()
                elif(i < int(n[1])+nnodal+14):
                    dif = line.split()
                elif(i < int(n[1])+nnodal+16):
                    f = line.split()
                elif(i < int(n[1])+nnodal+17):
                    S = line.split()
        Elem = np.array(Elem)
        Nodal = np.array(Nodal)
        S = np.array(S)
        rho_truelog = np.array(rho_truelog)
       
       
        for i in range(0,len(Nodal)):
            Nodal[i,-1]=i+1
       
        mshl = len(Elem)
       
        for i,e in enumerate(Elem):
            line1 = int(e[1])+1
            line2 = int(e[2])+1
            line3 = int(e[3])+1
           
            index1 = np.where(Nodal[:,-1] == str(line1))
            index2 = np.where(Nodal[:,-1] == str(line2))
            index3 = np.where(Nodal[:,-1] == str(line3))
            indexC = rho_truelog[i]
       
            plotX = [Nodal[index1,1], Nodal[index2,1], Nodal[index3,1],Nodal[index1,1]]
            plotY = [Nodal[index1,2], Nodal[index2,2], Nodal[index3,2],Nodal[index1,2]]
       
            plotC = cmaplist[int((float(indexC) - mm[0])*len(cmaplist)/(mm[1] - mm[0])) - 1] 
       
            polygon = ((Nodal[index1,0][0][0], Nodal[index1,1][0][0]), (Nodal[index2,0][0][0],Nodal[index2,1][0][0]), (Nodal[index3,0][0][0], Nodal[index3,1][0][0]),(Nodal[index1,0][0][0], Nodal[index1,1][0][0]))
            codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY]#patches.append(polygon)
            path = Path(polygon,codes)
       
            patch = patches.PathPatch(path, facecolor=plotC, edgecolor=None, lw=0.2)
            ax.add_patch(patch)
       
        fname = path_d + 'pngs/' + re.sub(path_d,' ',data[w])[6:-3] + 'png'
        ax.set_title(re.sub(path_d,' ',fname)[6:-4],fontsize=16, color='gray')
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None, metadata=None)

plt.close(fig)
