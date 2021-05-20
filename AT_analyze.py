import numpy as np
import glob
import matplotlib.pyplot as plt
import math

# set fontsize
plt.rcParams.update({'font.size': 20})

# choose what to plot, 1: dist 2: Rzz, 3: dist vs Rzz
PLOT = 3

summary_file = 'summary_dist_Rzz.dat' 
with open(summary_file,'w') as f:
            f.write('%-10s %-10s %-10s %10s %10s %10s %14s %10s %10s %10s %14s \n' % ('folder','lip_name','rep','dist_mean','dist_start','dist_max','move_dist_max','Rzz_mean','Rzz_start','Rzz_min','move_Rzz_max'))

dist_list = []
Rzz_list = []
colors = ['red','darkorange','green']
markers = ['s','o','X','^','*','d','h']
folders = [22,27,28,18,23,25,15]
protein_names = ['KIBRA','PI3K (front mode)','PI3K (back mode)','RIM2','PTEN','SHIP2','SMURF2']

if PLOT == 3:
    plt.figure(figsize=(10,7))
    list = range(len(protein_names))
else:
    list = range(4,5)

for protein in list:
    folder = folders[protein]
    marker = markers[protein]
    protein_name = '%s-C2' % protein_names[protein]
    if PLOT < 3:
        plt.figure(figsize=(5,6))
    for lip in range(3):
        
        # get lipid name
        if lip == 0:
            lip_name = 'PC'
            lip_label = 'PC'
        elif lip == 1:
            lip_name = 'PCPS'
            lip_label = 'PC:PS'
        else:
            lip_name = 'PCPSPIP2'
            lip_label = 'PC:PS:PIP2'
        
        # get reps
        folder_name = folder
        if folder == 22:
            reps = [17,15,11]
        elif folder == 27:
            reps = [2,10,7]
        elif folder == 28:
            reps = [22,6,19]
            folder_name = 27
        elif folder == 18:
            reps = [10,23,11]
        elif folder == 23:
            reps = [3,19,1]
        elif folder == 25:
            reps = [15,22,23]
        elif folder == 15:
            reps = [23,1,24]
        rep = reps[lip]

        # import dist
        dist_file = '%d/%s_%d/CG*/FINAL/dist_com_cent.xvg' % (folder_name,lip_name,rep)
        for name in glob.glob(dist_file):
            time,dist = np.genfromtxt(name,skip_header=17,usecols=[0,1],unpack=True)
        ndx = np.where(dist>3.6)
        dist_start = dist[0]
        dist_mean = np.mean(dist[ndx])
        dist_max = np.amax(dist)
        move_dist_max = dist_max - dist_start
        D_dist= dist - dist_start

        # import Rzz
        Rzz_file = '%d/%s_%d/CG*/FINAL/Rzz.xvg' % (folder_name,lip_name,rep) 
        for name in glob.glob(Rzz_file):
            time,Rzz = np.genfromtxt(name,skip_header=32,usecols=[0,9],unpack=True)
        Rzz_start = Rzz[0]
        Rzz_mean = np.mean(Rzz)
        Rzz_min = np.amin(Rzz)
        move_Rzz_max = Rzz_start - Rzz_min
        D_Rzz = Rzz_start - Rzz
        
        # change time units
        time /= 1000

        # print results to output and file
        print('%d/%s_%d: %1.2f, %1.2f, %1.2f, %1.2f' % (folder,lip_name,rep,dist_mean,dist_start,dist_max,move_dist_max))
        print('%d/%s_%d: %1.2f, %1.2f, %1.2f, %1.2f' % (folder,lip_name,rep,Rzz_mean,Rzz_start,Rzz_min,move_Rzz_max))
        with open(summary_file,'a') as f:
            f.write('%-10d %-10s %-10d %10.2f %10.2f %10.2f %14.2f %10.2f %10.2f %10.2f %14.2f \n' % (folder,lip_name,rep,dist_mean,dist_start,dist_max,move_dist_max,Rzz_mean,Rzz_start,Rzz_min,move_Rzz_max))
        
        # append to varibles
        dist_list.append(move_dist_max)
        Rzz_list.append(move_Rzz_max)

        # plot results
        if PLOT == 1:
            plt.plot(time[ndx],D_dist[ndx],label=lip_label,linestyle='none',marker='.',color=colors[lip])
            plt.title(r'$\Delta$dist, %s' % protein_name)
        elif PLOT == 2:
            plt.plot(time,D_Rzz,label=lip_label,linestyle='none',marker='.',color=colors[lip])
            plt.title(r'$\Delta R_{zz}$, %s' % protein_name)
        else:
            plt.plot(move_Rzz_max,move_dist_max,color=colors[lip],marker=marker,markersize=16)

if PLOT == 1 or PLOT == 2:
    if PLOT == 1:
        plt.ylabel(r'$\Delta$dist [nm]')
    else:
        plt.ylabel(r'$\Delta R_{zz}$')
    plt.xlabel('time [ns]')
    #plt.legend()
else:
    plt.xlim([0.0,2.1])
    ymax = np.amax(dist_list)
    plt.ylim([0.0,math.ceil(ymax)])
    plt.ylabel(r'max($\Delta$dist) [nm]')
    plt.xlabel(r'max($\Delta R_{zz}$)')
plt.tight_layout()
if PLOT == 1:
    plt.savefig('../../Seafile/C2_Manus/AT/AT_Dist_vs_time.svg',format='svg') # libreoffice cannot work with eps
elif PLOT == 2:
    plt.savefig('../../Seafile/C2_Manus/AT/AT_Rzz_vs_time.svg',format='svg')
elif PLOT == 3:
    plt.savefig('../../Seafile/C2_Manus/AT/AT_Dist_vs_Rzz.svg',format='svg')
plt.show()
            #plt.plot(time[ndx],dist[ndx],label=label)
