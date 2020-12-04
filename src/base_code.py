#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 09:24:51 2020

@author: felipe
"""

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
from scipy.optimize import minimize
import meshio

def line_(*argv):
    for i, arg in enumerate(argv):
        s = 0
        for j, a in enumerate(arg):
            s = s + a * depth_ ** j
    return s

def err(x):
    return np.sum((line_(x) - ff)**2)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

path_d = "../Example_files/"

if os.path.exists(path_d + "invdir/ref"):
    ref = path_d + "invdir/ref/f001_res.vtk" 
    copyfile(ref,path_d + "invdir/" + "f000_res.vtk") 

files = [f for f in glob.glob(path_d + "invdir/f???_res.vtk")]
files.sort()

## Please, name the files in the format yyyymmdd_hhss
data = [d for d in glob.glob(path_d + "data/*.Data")]
data.sort()

if not os.path.exists(path_d + "pngs/dif"):
        os.makedirs(path_d + "pngs/dif")

pngs = [p for p in glob.glob(path_d + "pngs/*.png")]

if not os.path.exists(path_d + "pngs/dif"):
        os.makedirs(path_d + "pngs/dif")

if os.path.exists(path_d + "nn_bh_temp.txt"):
    temperature = True
    sensor_depth = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    with open(path_d + "nn_bh_temp.txt",'r') as fh:
        dump = fh.readlines()
        tt = np.zeros((len(dump), 8))
        for i, line in enumerate(dump):
            d = line.split()[0:2]
            dd=(d[0][0:4]+d[0][5:7]+d[0][8:10]+d[1][0:2])
            tt[i,0] = dd
            tt[i,1:8] = line.split()[2:]        
    temp = np.zeros(len(data))
    for i, d in enumerate(data):
        for t in dump:
            if(t.split()[0]+"_"+t.split()[1][0:2] == d[11:-11]):
                temp[i] = float(t.split()[3])
    
    if not os.path.exists(path_d + "pngs/t1D"):
        os.makedirs(path_d + "pngs/t1D")
    
    if os.path.exists(path_d + "noname.txt"):
        el = np.loadtxt(path_d + "noname.txt", delimiter=',')
                

if not os.path.exists(path_d + "out"):
        os.makedirs(path_d + "out")


# Define the colormap_
cmap = plt.cm.seismic
cmaplist2 = [cmap(i) for i in range(cmap.N)]

factor = 10 # increase this to drop the namber of colors
cmaplist = []
for i, c in enumerate(cmaplist2):
    if (i== 0 or i%factor == 0):
        cmaplist.append(cmaplist2[i]) 
cmapN = len(cmaplist)
cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmapN)
mm =[-100, 100] #you can choose diffferent thresholds
bounds = np.linspace(mm[0], mm[1], cmapN)
norm = mpl.colors.BoundaryNorm(bounds, cmapN)

# Define the axes
fig_size = (14, 6)
axes2 = [0.93, 0.1, 0.02, 0.8] # palete
text_size = 12
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

ref={}
ttt_ref = {}
for w,f in enumerate(files):
    
    fname_t = path_d + 'pngs/t1D/' + re.sub(path_d,' ',data[w])[6:-4] + 'png'
    title = re.sub(path_d,' ',fname_t)[10:-4]
    print("Working on " + f + ' (' + title[4:6] + '/' + title[6:8] + '/' + title[:4] + ' ' + title[9:11] + 'h)')
    
    dd = [0, 0]
    mesh = meshio.read(f)
    rho_truelog = mesh.cell_data['Resistivity(log10)'][0]
    rho = mesh.cell_data['Resistivity(ohm.m)'][0]
    Elem = np.array(mesh.cells).all()
    Nodal = mesh.points
    
    fname_out = path_d + 'out/' + re.sub(path_d,' ',data[w])[6:-4] + 'dat'
    
       
    if(temperature):
        a = 0
        for t in tt:
            if(title[:8] + title[9:11] in str(t[0]) and a == 0):
                a = 1
                # Polynomium of 3rd degree. skip the first measure
                depth_ = sensor_depth[1:]
                ff = t[2:]
                x0 = [t[1], 1, 0, 0]
                bnds = ((None, None),(None, None),(None, None), (None, None),)
                s = minimize(err, x0, method='SLSQP', bounds=bnds)
                z1 = np.linspace(0,3,30)
                temp1 = s.x[0] + s.x[1]*z1 + s.x[2]*z1**2 + s.x[3]*z1**3
                plt.plot(temp1,z1,'b',t[1:],sensor_depth,'*r')
            
                # line
                depth_ = sensor_depth[-3:-1]
                ff = t[-3:-1]
                x0 = [ff[0], 1]
                bnds = ((None, None),(None, None))
                s2 = minimize(err, x0, method='SLSQP', bounds=bnds)           
                z2 = np.linspace(3,10,70)
                temp2 = s2.x[0] + s2.x[1]*z2
                plt.plot(temp2,z2,'b')
                plt.gca().invert_yaxis()
                plt.show()
                plt.xlim((-5, 20))
                plt.xlabel("Temperature",size=text_size)
                plt.ylabel("depth",size=text_size)
                plt.title(title[4:6] + '/' + title[6:8] + '/' + title[:4] + ' ' + title[9:11] + 'h',fontsize=16, color='gray')
                plt.savefig(fname_t)
                plt.close()                
    
    
    
    if(w == 0):
        for i,e in enumerate(Elem):
            
            index1 = int(e[0]) 
            index2 = int(e[1])
            index3 = int(e[2])

            polygon = ((Nodal[index1,0], Nodal[index1,1]),
                       (Nodal[index2,0],Nodal[index2,1]),
                       (Nodal[index3,0], Nodal[index3,1]),
                       (Nodal[index1,0], Nodal[index1,1]))
            
            x = (Nodal[index1,0] + Nodal[index2,0] + Nodal[index3,0]) / 3.0
            z = (Nodal[index1,1] + Nodal[index2,1] + Nodal[index3,1]) / 3.0
                
            elevation = el[find_nearest(el[:,0],x),1]
            z_ = elevation - (Nodal[index1,1] + Nodal[index2,1] + Nodal[index3,1]) / 3.0
            if(z_ < 0):
                z_ = 0.0
                #print(el[:,0])
                #print(elevation)
            if( z_ <= 3):
                ttt = s.x[0] + s.x[1]*z_ + s.x[2]*z_**2 + s.x[3]*z_**3
            else:
                ttt = s2.x[0] + s2.x[1]*z_
            
            
            ref.update({polygon:rho[i]})
            ttt_ref.update({polygon:ttt})
    else:
        
        xmin = min(Nodal[:,0])
        xmax = max(Nodal[:,0]) 
        ymin = min(Nodal[:,1])
        ymax = max(Nodal[:,1]) + 2
        
        fig, ax = plt.subplots(1, 1, figsize = fig_size)
        
        ax.set_xlim(xmin,xmax)
        ax.set_ylim(ymin,ymax)
        ax.set_xlabel("Distance [$m$]",size=text_size)
        ax.set_ylabel("Elevation [$m$]",size=text_size)
        
        with open(fname_out,'w') as gg:
            gg.write("       x_coord        y_coord           Rho           Rho_Ref      %     Depth    Temp    Temp_Ref    Rho_cT     Rho_ref_cT  %_cT\n")
        
        
            for i, e in enumerate(Elem):
                index1 = int(e[0]) 
                index2 = int(e[1])
                index3 = int(e[2])
                
                polygon = ((Nodal[index1,0], Nodal[index1,1]),
                            (Nodal[index2,0],Nodal[index2,1]),
                            (Nodal[index3,0], Nodal[index3,1]),
                            (Nodal[index1,0], Nodal[index1,1]))
                
                x = (Nodal[index1,0] + Nodal[index2,0] + Nodal[index3,0]) / 3.0
                z = (Nodal[index1,1] + Nodal[index2,1] + Nodal[index3,1]) / 3.0
                    
                elevation = el[find_nearest(el[:,0],x),1]
                z_ = elevation - (Nodal[index1,1] + Nodal[index2,1] + Nodal[index3,1]) / 3.0
                if(z_ < 0):
                    z_ = 0.0
                    #print(el[:,0])
                    #print(elevation)
                if( z_ <= 3):
                    ttt = s.x[0] + s.x[1]*z_ + s.x[2]*z_**2 + s.x[3]*z_**3
                else:
                    ttt = s2.x[0] + s2.x[1]*z_
                
                rhocT = rho[i]*(1 + 0.025*(ttt - 18))
                rhoRefcT = ref[polygon]*(1 + 0.025*(ttt_ref[polygon] - 18))
    
        
                indexC = (100.0* ( float(rho[i]) - float(ref[polygon])) / float(ref[polygon]))
                indexCT = (100.0* ( float(rhocT) - float(rhoRefcT)) / float(rhoRefcT)  )
                indice = int( (indexCT - mm[0])*(len(cmaplist) - 1) / (mm[1] - mm[0]) )
                
                if(indexCT < dd[0]):
                    dd[0] = indexC
                if(indexCT > dd[1]):
                    dd[1] = indexC
                
                gg.write('{:15.6f} {:15.6f} {:15.6f} {:15.6f} {:5.1f} {:5.3f} {:5.3f} {:5.3f} {:15.6f} {:15.6} {:5.1f}\n'.format(x, z, rho[i], ref[polygon], indexC, z_, ttt, ttt_ref[polygon], rhocT, rhoRefcT, indexCT))
                
                if (indice >= len(cmaplist)):
                    indice = len(cmaplist) - 1;
                elif(indice < 0):
                    indice = 0
                
                plotC = (cmaplist[indice])  
                            
                codes = [Path.MOVETO,Path.LINETO,Path.LINETO,Path.CLOSEPOLY]
                path = Path(polygon,codes)
           
                patch = patches.PathPatch(path, facecolor=plotC, edgecolor=None, lw=0.2)
                ax.add_patch(patch)

        ax2 = fig.add_axes(axes2)
        ax2.set_title('$ (\\rho - \\rho_0) /\\rho_0 \, [\%]$',size=text_size)
        cb = mpl.colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', ticks=[np.linspace(np.around(mm[0]),np.around(mm[1]),5)], boundaries=bounds)
       
        fname = path_d + 'pngs/dif/' + re.sub(path_d,' ',data[w])[6:-4] + 'png'
        title = re.sub(path_d,' ',fname)[10:-4]
        ax.set_title(title[4:6] + '/' + title[6:8] + '/' + title[:4] + ' ' + title[9:11] + 'h',fontsize=16, color='gray')
        ax.text(xmin+2,ymax-1.5,"[min max] = " + str(np.round(dd)))
        
        plt.savefig(fname)
        #plt.show()
        plt.close(fig)
        plt.clf()
        plt.cla()