#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot Rzz against distance
Andreas Haahr Larsen
"""

# import libraries
import numpy as np
import pylab
import matplotlib as mpl
from matplotlib import pyplot as plt

#fontsize
plt.rcParams.update({'font.size': 22  })
    
# folder
Folder = 28

# plot a single or selected trajectories
SELECTED_TRAJS = 1
trajs = [4]

# bins
bins = 100

# find frame corresponding to point on Dist vs. Rzz plot
aTol = 0 # increase tolerance area by a factor 
RzzInt = 0.99
DistInt = 4.03
RzzTol = 0.002*aTol
DistTol = 0.002*aTol

# include only last frames 286 frames (out of 5734), i.e. 100 ns out of 2 us
LAST = 0
last_frames = 500

# choose to save figures or not
SAVE = 0
# save high resolution (and without title and grid)
HIGH_RES = 0
# save high res with color bar as Rzz_color_bar.eps
HIGH_RES_COLOR_BAR = 0 

# plot data from all repeats
PLOT_ALL = 0

# lipid compositon, color, number of repetitions
Reps = [25,0,0] # PCPSPIP2, PCPS, PC

if Folder == 8:
    Protein = "PI3K"
elif Folder == 9:
    Protein = "FVa"
    #traj = [18,19,20,21,23] mode 1a Rzz 1.0
    #traj = [1,2,4,5,8,10,15] mode 1b Rzz 0.75
    #traj = [0,14] mode 2 Rzz 0.35    
    #traj = [6,7] mode 4 Rzz 0.0
    #traj = [9,24] mode 5 Rzz -0.93

    #traj = [3,11] mode 1b and XX (Rzz -0.3)
    #traj = [12,16] mode 5 and XX (Rzz -0.5)
    #traj = [13,22] mode 1 and mode 3 (Rzz 0.3)
elif Folder == 13:
    Protein = "MFGE8"
elif Folder == 15:
    Protein = "Smurf2"
    #trajs = [0,5,6,8,11,13,16,18,20,21] # mode 1 Rzz 0.99 (0 good examples)
    #trajs = [1,3,10,12,14,19,22,23,24] # mode 2 Rzz -0.3
    #trajs = [2] # mode 3 Rzz -1.0 (ends in mode 2)
    
    #trajs = [7] # mode 1 and ??
    #trajs = [4] # mode ?? and ??
    #trajs = [9] # mode 2 and 3
    #trajs = [15] # mode ?? and ?? and ???
elif Folder == 17:
    Protein = "FVIII"
    #trajs = [1,6,7,12,19,22] # mode 1 Rzz 1.0 (19 good example)
    #trajs = [8,14] # mode 2 Rzz 0.6 (14 good example)
    #trajs = [9] # mode 3 Rzz 0.2
    #trajs = [10,21,24] # mode 4 Rzz -0.3 (10 and 21 mix)
    #trajs = [0,2,3,13,17,18,20] # mode 5 from -0.5 to -0.9
    #trajs = [15,16] # mode 6 Rzz -1.0
elif Folder == 18:
    Protein = "RIM2" #t
    # trajs = [0,2,5,8,11,12,14,16,17,23,24] # mode 1
    # trajs = [3,6,9] # mode 2
    # trajs = [1,4,19,20,21] # mode 1 and 2
    # trajs = [7,10,15,18,22]  # mode 1 and 3
    # trajs = [13] # mode 2 and 3
elif Folder == 21:
    Protein = "DLL4" #t
    #trajs = [1,2,8,12,13] # mode 1 Rzz 0.99 Dist 3.42 (12 good example)
    #trajs = [4,9,11,19,21,22] # mode 2 Rzz -0.04 Dist 3.53 (11 good example)
    #trajs = [0,5,6,7,14,15,16,18,20,23,24] # mode 3 Rzz -0.67 Dist 3.95 (14 good example). not very different from mode 2 
    #trajs = [10] # mode 1 and 2 
    #trajs = [3] # mode 1 and 3 
elif Folder == 22:
    Protein = "KIBRA" #t
    #trajs = [1,12,13,17,19,21,23] # mode 1
elif Folder == 23:
    Protein = "PTEN" #t
    #trajs = [1,2,5,11,18,20] mode 1, top, productive mode
    #trajs = [0,3,4,6,7,8,9,10,13,14,15,16,17,19,21,22,23] mode 2, back, unproductive mode
elif Folder == 24:
    Protein = "SHIP2" #t
elif Folder == 25:
    Protein = "SHIP2" #tt
    #trajs = [6,7,16,17,23] # mode 1 Rzz 1.0
    #trajs = [1,5,8,10,11,12,13,18,19,21] # mode 2 Rzz -0.15
    #trajs = [2,14,20,22] # mode 3 Rzz -0.9
    
    #trajs = [9] # mode 1 and 2
    #trajs = [0] # mode 1 and 3 
    #trajs = [3,4,24] # mode 2 and 3 
    #trajs = [15] # mode 1, 2 and 3 (end in mode 2)
elif Folder == 26:
    Protein = "DLL1" #t
elif Folder == 27:
    Protein = "PI3K"
    #trajs = [0,1,6,7,12,13,15,20] # mode 1 Rzz 1.0 (7 good example)
    #trajs = [3,5,9,10,14,17,18,23] # mode 2 Rzz 0.2 to 0.3
    #trajs = [2,19] # mode 3 Rzz -0.9 Dist 3.2
    #trajs = [4,8,16,21,24] # mode 4 Rzz -0.9 Dist 3.6
    
    #trajs = [22] # mode 2 and 3
    #trajs = [11] # mode 3 and 4
elif Folder == 28:
    Protein = "SHIP2_FL"
    # trajs = [5,13,14,20,23,27] # mode 1, top
    # trajs = [0,3,11,12] # mode 4 --> mode 1, Ptase --> top
    # trajs = [1,2,7,9,10,16,17,18,19,21,22,26] # mode 2, back
    # trajs = [6,8,24,25] # mode 3, front/bottom
    # trajs = [15] # mode 4 --> mode 3, Ptase --> front/bottom
    # trajs = [4] # mode 4, Ptase

    # mode 1 and mode 4>1 : 10/28 = 36 %
    # mode 2              : 12/28 = 43 %
    # mode 3 and mode 4>3 :  5/28 = 18 % 
    # mode 4              :  1/28 =  4 % 
elif Folder == 29:
    Protein = "PTEN_FL"
    #trajs = [3,7,13,15,18,19,21,23] # mode 1, C2 top [8/11]
    #trajs = [6,17,22] # mode 2, C2 side, N-term down, intermediate mode towards mode 1 [3/25]
    #trajs = [0,1,2,3,4,9,10,11,14,20,24] # mode 3, C2 side, N-term up [11/25]
    #trajs = [5,8,12,16] # move from mode 2  (side) to mode 1 (top) [4/25]


Elements = [\
            ["PCPSPIP2","red",Reps[0],"PC:PS:PIP$_2$"],\
            ["PCPS","blue",Reps[1],"PC:PS"],\
            ["PC","green",Reps[2],"PC"],\
            ] 
for Element in Elements:
    Lip = Element[0]
    Color = Element[1]
    Repeats = Element[2]
    Lip_name = Element[3]
    if Repeats > 0:
        RzzAll = []
        DistAll = []
        print('\n##################################################\nLipid: %s' % Lip)
        if SELECTED_TRAJS == 0:
            trajs = np.arange(Repeats)
        for i in trajs:
            # print('i = %d' % i)
            # Rzz
            File = str(Folder) + "/" + Lip + "_" + str(i) + "/Rzz.xvg"
            Time,Rxx,Rxy,Rxz,Ryx,Ryy,Ryz,Rzx,Rzy,Rzz = np.genfromtxt(File,skip_header=32,usecols=[0,1,2,3,4,5,6,7,8,9],unpack=True)
            Time /= 1000 # convert from ps to ns
            
            if LAST:
                Time = Time[-last_frames:]
                Rzz = Rzz[-last_frames:]
                
            if PLOT_ALL:
                Figure1=pylab.figure(1)
                plot1=pylab.subplot(111)
                plot1.plot(Time,Rzz,label="Rzz",color="r")
                
                plot1.legend(loc=1,fontsize=9)
                plot1.set_xlabel("Time [ns]")
                plot1.set_ylabel("Rzz [nm]")
                pylab.suptitle("Rzz")
                pylab.show()
            
            RzzAll = np.append(RzzAll,Rzz)
        
            # Dist
            File = str(Folder) + "/" + Lip + "_" + str(i) + "/dist_com_cent.xvg"
            Time,Dist = np.genfromtxt(File,skip_header=17,usecols=[0,1],unpack=True)
            
            if LAST:
                Time = Time[-last_frames:]
                Dist = Dist[-last_frames:]
            
            # find structure of interest
            Rzz_filter_index = np.where((Rzz > RzzInt-RzzTol) & (Rzz < RzzInt + RzzTol))
            Dist_filter_index = np.where((Dist > DistInt-DistTol) & (Dist < DistInt + DistTol))
            Common_index = np.intersect1d(Rzz_filter_index,Dist_filter_index)
            l_Common = len(Common_index)
            
            if l_Common:
                for j in range(l_Common):
                    print('Structure of interest: lip %s, repeat %d, frame %d' % (Lip,i,Common_index[j]))
            # for debubbing:
            #else:    
                #len_Rzz = len(Rzz_filter_index[0])
                #len_Dist = len(Dist_filter_index[0])
                #print('left Rzz values = %d, left Dist values = %d, left common values = %d' % (len_Rzz,len_Dist,l_Common) )
            
            if PLOT_ALL:  
                Figure1=pylab.figure(1)
                plot1=pylab.subplot(111)
                plot1.plot(Time,Dist,label="Dist",color="g")
                
                plot1.legend(loc=1,fontsize=9)
                plot1.set_xlabel("Time [ns]")
                plot1.set_ylabel("Dist [nm]")
                pylab.suptitle("Dist")
                pylab.show()
                
            DistAll = np.append(DistAll,Dist)
            
            # Together
            if PLOT_ALL:
                Figure1=pylab.figure(1)
                plot1=pylab.subplot(111)
                plot1.plot(Rzz,Dist,label="Dist vs Rzz",color="k",marker=".",linestyle="none")
                
                plot1.legend(loc=1,fontsize=9)
                plot1.set_xlabel("Rzz")
                plot1.set_ylabel("Dist")
                pylab.suptitle("Dist vs Rzz")
                pylab.show()
            
        # Histogram Rzz vs Dist
        hist=plt.hist2d(RzzAll,DistAll,bins=bins,norm=mpl.colors.LogNorm())
        
        
        plt.clim(1,1000)
        
            
        plt.ylim((2.5,7.5))
        plt.xlim((-1.0,1.0))
        plt.tight_layout()
        if HIGH_RES == 0:
            plt.title("C2 from " + Protein + " and " + Lip + " membrane")
            plt.minorticks_on()
            plt.grid()
            plt.grid(which='minor')
            plt.xlabel('$R_{zz}$')
            plt.ylabel('Distance [nm]')
            plt.colorbar()
        else:
            if HIGH_RES_COLOR_BAR:
                plt.colorbar()
            else:
                plt.title(r'%s' % Lip_name)
                if Folder == 15:
                    plt.xlabel('$R_{zz}$')    
                if Lip == 'PC':
                    plt.ylabel('Distance [nm]')
#            if Lip == 'PCPSPIP2':
#                plt.colorbar()    
        plt.tight_layout()     
        if SAVE and SELECTED_TRAJS == 0:
            if LAST:
                plt.savefig('../../Seafile/C2_Manus/Rzz_%s_%s_last.png' % (Protein,Lip))
            else:
                if HIGH_RES:
                    if HIGH_RES_COLOR_BAR:
                        plt.savefig('../../Seafile/C2_Manus/Rzz_color_bar.eps', format='eps')
                    else:
                        plt.savefig('../../Seafile/C2_Manus/Rzz_%s_%s.eps' % (Protein,Lip), format='eps')
                    
                else:
                    plt.savefig('../../Seafile/C2_Manus/Rzz_%s_%s.png' % (Protein,Lip))
        plt.show()

        # find bin with highest occupancy
        maxbin=np.amax(hist[0])
        indc = np.where(maxbin==hist[0])
        if len(indc) > 1:
            indxRzz = indc[0][0]
            indxDist = indc[1][0]
        else:
            indxRzz = indc[0]
            indxDist = indc[1]
        Rzz_indx_min = hist[1][indxRzz]
        Rzz_indx_max = hist[1][indxRzz+1]
        Dist_indx_min = hist[2][indxDist]
        Dist_indx_max = hist[2][indxDist+1]
        Rzz_indx_mean = (Rzz_indx_min + Rzz_indx_max)/2
        Dist_indx_mean = (Dist_indx_min + Dist_indx_max)/2
        
        print('Bin with highest density: Rzz = %1.2f, Dist = %1.2f' % (Rzz_indx_mean,Dist_indx_mean))
        
