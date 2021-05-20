#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot contacts per residue
Finds frame with closest bound residue
Andreas Haahr Larsen
"""

# import libraries
import numpy as np
from matplotlib import pyplot as plt

# folder
Folder = 28

# lips, color, repeats, PS, PIP2, PIP3
if   Folder == 8 : Elements = [["PCPSPIP2","red",25,0,0,0,0.3,"PI3K",1]] # PI3K, 6btz
elif Folder == 15: Elements = [["PCPSPIP2","red",25,0,0,0,0.3,"SMURF2",1]] # 2jqz
elif Folder == 18: Elements = [["PCPSPIP2","red",25,0,0,0,0.3,"RIM2",1]] # RIM2, 2bwq_t
#elif Folder == 21: Elements = [["PCPSPIP2","red",25,0,0,0,0.2,"DLL4",1]] # DLL4, 5mvx_t
elif Folder == 22: Elements = [["PCPSPIP2","red",25,0,0,0,0.15,"KIBRA",1]] # KIBRA, 6fjd_t
elif Folder == 23: Elements = [["PCPSPIP2","red",25,0,0,0,0.45,"PTEN",1]] # PTEN, 1d5r_t
elif Folder == 25: Elements = [["PCPSPIP2","red",25,0,0,0,0.3,"SHIP2",1]] # SHIP2, 5okm_tt
elif Folder == 26: Elements = [["PCPSPIP2","red",25,0,0,0,0.3,"DLL1",1]] # DLL1, 4xbm_t
#elif Folder == 27: Elements = [["PCPSPIP2","red",5,0,0,0,0.3,"PI3K",1]] # PI3K, 6bu0
elif Folder == 28: Elements = [["PCPSPIP2","darkgreen",[5,13,14,20,23,27],0,1,0,0.3,"SHIP2, mode 1",418],["PCPSPIP2","darkred",[1,2,7,9,10,16,17,18,19,21,22,26],0,1,0,0.3,"SHIP2-FL, mode 2",418],["PCPSPIP2","blue",[6,8,24,25],0,1,0,0.3,"SHIP2-FL, mode 3",418]] # SHIP2, 5okm_FL
elif Folder == 29: Elements = [["PCPSPIP2","darkgreen",[3,7,13,15,18,19,21,23],0,1,0,0.3,"PTEN, mode 1",7],["PCPSPIP2","darkred",[0,1,2,3,4,9,10,11,14,20,24],0,1,0,0.3,"PTEN-FL, mode 2",7]] # PTEN, 155r_FL

# distance below CUTOFF is considered a contact [nm]
cutoff = 0.45 

ymax_previous = 0.0

# figure settings
plt.rcParams.update({'font.size': 14})
plt.figure(figsize=(8,6))

# for loop over lipids
for Element in Elements:  
    Lip = Element[0]
    Color = Element[1]
    Repeats = Element[2]
    PS = Element[3]
    PIP2 = Element[4]
    PIP3 = Element[5]
    threshold = Element[6]
    name = Element[7]
    first_res = Element[8]
    #Mode = Element[9]

    #active_site_min=20
    #active_iste_max=30
    bound_res_total=0
    Nres_total=0  
 
    # naming depending on PIP choice
    if PS:
        PIP = '_PS'
        extra_label = ' (binding to PS only)'
    elif PIP2:
        PIP = '_PIP2'
        extra_label = ' (binding to PIP2 only)'
    elif PIP3:
        PIP = '_PIP3'
        extra_label = ' (binding to PIP3 only)'
    else:
        PIP = ''
        extra_label = ' (binding to any lipid)'
    
    FIRST = 1
    # for loop over repeats
    # for rep in range(Repeats):
    #for rep in [3,7,13,15,18,19,21,23]: 
    for rep in Repeats:
        filename = str(Folder) + '/' + Lip + '_' + str(rep) + '/dist_res_min' + PIP + '.xvg'
        data = np.genfromtxt(filename,skip_header=24,unpack=True)
        time = data[0]
        Nres = len(data)-1
        Nframes = len(time)
        if FIRST:
            bound_res = np.zeros(Nres)
            FIRST = 0
        for res in range(1,Nres+1):
            # indices of timeslots for that residued with min dist below the cutoff
            ndx=np.where(data[res]<cutoff) 
            bound_res[res-1] = len(ndx[0])
        bound_res_total += bound_res    
        Nres_total =+ Nframes
        
    # fraction of measured residues with contact
    fraction = bound_res_total/Nres_total
    
    # make residues and frames arrays
    residue = np.arange(1,Nres+1) 
    frame = np.arange(1,Nframes+1)
    
    # plot histogram
    label = Lip + extra_label
    #plt.bar(residue+first_res-1,fraction,color=Color,alpha=0.5,label=label)
    plt.bar(residue+first_res-1,fraction,linewidth=0.8,color=Color,alpha=0.5,label=name)
    
    # find closest bound residue
    key_residue_index = np.where(fraction==np.amax(fraction))[0]
    key_residue = residue[key_residue_index]
    
    # find frame with closest bound residue (of all repetitions)
    min_dist_global = 99999 # very large number
    for rep in Repeats:
        filename = str(Folder) + '/' + Lip + '_' + str(rep) + '/dist_res_min' + PIP + '.xvg'
        data = np.genfromtxt(filename,skip_header=24,unpack=True)
        min_dist_local = np.amin(data[key_residue_index])
        if min_dist_local < min_dist_global:
            min_dist_global = min_dist_local
            key_frame_index = np.where(data[key_residue_index][0]==min_dist_local)
            rep_min = rep
    key_frame = frame[key_frame_index]
    
    # print
    print(label + '\t : tightest bound residue, no. ' + str(key_residue[0]) + ', is in frame ' + str(key_frame[0]) + ' (frame ' + str(key_frame_index[0][0]) + ' if index 0), in repetition ' + str(rep_min) + ' (index 0), distance: ' + str(min_dist_global) + ' A.' )
    
    # find several key binding residues
    key_residues_index = np.where(fraction>threshold)[0]
    key_residues = key_residues_index + 1
    print('     key residues    \t \t : ' + str(key_residues))
    
    ymax = np.amax([np.amax(fraction),ymax_previous])
    ymax_previous = ymax
# plot threshold line for selection of key residues
#plt.plot(residue,np.ones(len(residue))*threshold,color='grey',label='threshold',linestyle='dashed')

# plot area for the active site

y=ymax/2.
y1=0.9*y
y2=1.1*y
y3=1.2*y
lw=2
alf=0.7

if Folder == 28:

    plt.plot([424,431],[y,y],linewidth=lw,color='lightgray')
    plt.plot([441,447],[y,y],linewidth=lw,color='lightgray')
    plt.plot([467,472],[y,y],linewidth=lw,color='lightgray')
    plt.plot([477,494],[y,y],linewidth=lw,color='lightgray')
    #plt.plot([479,480],[y,y],linewidth=lw,color='red',alpha=alf)
    plt.plot([499,505],[y,y],linewidth=lw,color='lightgray')
    plt.plot([509,514],[y,y],linewidth=lw,color='lightgray')
    plt.plot([523,531],[y,y],linewidth=lw,color='lightgray')
    plt.plot([533,541],[y,y],linewidth=lw,color='red',alpha=alf) # active part of L??
    plt.plot([543,550],[y,y],linewidth=lw,color='lightgray')
    plt.plot([553,561],[y,y],linewidth=lw,color='lightgray')
    plt.plot([562,568],[y,y],linewidth=lw,color='red',alpha=alf) # active part of L??
    plt.plot([569,583],[y,y],linewidth=lw,color='lightgray')
    plt.plot([595,599],[y,y],linewidth=lw,color='lightgray')
    plt.plot([602,607],[y,y],linewidth=lw,color='lightgray')
    plt.plot([615,625],[y,y],linewidth=lw,color='lightgray')
    plt.plot([629,643],[y,y],linewidth=lw,color='lightgray')
    plt.plot([675,684],[y,y],linewidth=lw,color='red',alpha=alf) # active part of L4
    plt.plot([690,696],[y,y],linewidth=lw,color='lightgray')
    plt.plot([703,709],[y,y],linewidth=lw,color='lightgray')
    plt.plot([721,727],[y,y],linewidth=lw,color='lightgray')
    plt.plot([748,758],[y,y],linewidth=lw,color='lightgray')
    plt.plot([765,770],[y,y],linewidth=lw,color='lightgray')
    plt.plot([776,780],[y,y],linewidth=lw,color='lightgray')
    plt.plot([793,798],[y,y],linewidth=lw,color='lightgray')
    plt.plot([821,827],[y,y],linewidth=lw,color='lightgray')
    plt.plot([833,840],[y,y],linewidth=lw,color='lightgray')
    plt.plot([850,858],[y,y],linewidth=lw,color='lightgray')
    plt.plot([861,873],[y,y],linewidth=lw,color='lightgray')
    
    
    #plt.plot([first_res,731],[y,y],linewidth=lw,color='black',label='Ptase')

    #plt.plot([433,435],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='active site, motif 1 (LeCoq2017, fig 6B, Mg site)')
    #plt.plot([452,463],[y3,y3],linewidth=lw,color='skyblue',alpha=alf,label='loop 1, res 452-463')
    #plt.plot([474,475],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='active site, motif 2 (LeCoq2017, fig 6B, Mg site)')
    #plt.plot([474,475],[y2,y2],linewidth=lw,color='pink',alpha=alf,label='active site, motif 2')
    #plt.plot([498,505],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='beta strand 3')
    #plt.plot([498,505],[y2,y2],linewidth=lw,color='pink',alpha=alf,label='beta strand 3')
    #plt.plot([532,541],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='active site, motif 3, loop 2 (LeCoq2017, fig 6') 
    #plt.plot([587,594],[y2,y2],linewidth=lw,color='pink',alpha=alf,label='loop 3 (LeCoq2017)')
    #plt.plot([607,608],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='res 607, active site (LeCoq2017, fig 6 and fig 7)')
    #plt.plot([620,640],[y3,y3],linewidth=lw,color='skyblue',alpha=alf,label='alpha-helices, res 620-640')
    #plt.plot([661,662],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='res 661, active site (LeCoq2017, fig 6)')
    #plt.plot([674,684],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='active site, motif 3, loop 4 (LeCoq2017)')

    #plt.plot([732,744],[y2,y2],linewidth=lw,color='pink',label='linker between C2 and Ptase')
    #plt.plot([732,744],[y3,y3],linewidth=lw,color='skyblue',label='linker between C2 and Ptase')
    
    #plt.plot([746,878],[y,y],linewidth=lw,color='lightgray',label='C2')
    #plt.plot([758,763],[y1,y1],linewidth=lw,color='lightgreen',label='CBL1')
    #plt.plot([775,783],[y2,y2],linewidth=lw,color='pink',label='B3')
    #plt.plot([800,812],[y2,y2],linewidth=lw,color='pink',label='loop 45')
    #plt.plot([827,832],[y2,y2],linewidth=lw,color='pink',label='CBL2')
    #plt.plot([827,832],[y1,y1],linewidth=lw,color='lightgreen',label='CBL2')
    
elif Folder == 29:
    
    plt.plot([12,16],[y,y],linewidth=lw,color='red',alpha=alf)
    plt.plot([25,29],[y,y],linewidth=lw,color='lightgray',alpha=alf)
    plt.plot([32,35],[y,y],linewidth=lw,color='lightgray',alpha=alf)
    plt.plot([49,61],[y,y],linewidth=lw,color='lightgreen',alpha=alf)
    plt.plot([65,71],[y,y],linewidth=lw,color='lightgray',alpha=alf)
    plt.plot([85,90],[y,y],linewidth=lw,color='lightgray',alpha=alf)
    plt.plot([97,113],[y,y],linewidth=lw,color='lightgreen',alpha=alf)
    plt.plot([119,123],[y,y],linewidth=lw,color='lightgray',alpha=alf)
    plt.plot([128,142],[y,y],linewidth=lw,color='lightgreen',alpha=alf)
    plt.plot([147,160],[y,y],linewidth=lw,color='lightgreen',alpha=alf)
    plt.plot([168,185],[y,y],linewidth=lw,color='lightgreen',alpha=alf)

    plt.plot([192,200],[y,y],linewidth=lw,color='lightgray',alpha=alf)
    plt.plot([213,219],[y,y],linewidth=lw,color='lightgray',alpha=alf)
    plt.plot([222,226],[y,y],linewidth=lw,color='lightgray',alpha=alf)
    plt.plot([238,250],[y,y],linewidth=lw,color='lightgray',alpha=alf)
    plt.plot([252,259],[y,y],linewidth=lw,color='lightgray',alpha=alf)
    plt.plot([268,276],[y,y],linewidth=lw,color='lightgray',alpha=alf)
    plt.plot([316,321],[y,y],linewidth=lw,color='lightgray',alpha=alf)
    plt.plot([342,349],[y,y],linewidth=lw,color='lightgray',alpha=alf)
    
    
    """
    plt.plot([first_res,186],[y,y],linewidth=lw,color='black',label='Ptase')
    plt.plot([7,23],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='N-terminal loop res 7-23')
    plt.plot([40,47],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='loop res 40-47')
    #plt.plot([72,83],[y2,y2],linewidth=lw,color='darkred',alpha=alf,label='loop A72 to C83')
    #plt.plot([84,90],[y2,y2],linewidth=lw,color='brown',alpha=alf,label='beta strand R84 to F90')
    plt.plot([72,90],[y2,y2],linewidth=lw,color='pink',alpha=alf,label='loop A72 to C83 and beta R84 to F90')
    plt.plot([93,94],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='WPD loop (Lee 1999)')
    plt.plot([122,131],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='P loop (Lee2015)')
    plt.plot([166,171],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='TI loop (Lee1999)')
   
    #plt.plot([187,191],[y2,y2],linewidth=lw,color='pink',label='linker between Ptase and C2')


    plt.plot([192,353],[y,y],linewidth=lw,color='lightgray',label='C2')
    #plt.plot([192,200],[y2,y2],linewidth=lw,color='blue',alpha=alf,label='B1')
    plt.plot([202,212],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='CBL1, loop 12')
    plt.plot([213,219],[y2,y2],linewidth=lw,color='pink',alpha=alf,label='B2')
    plt.plot([222,225],[y2,y2],linewidth=lw,color='pink',alpha=alf,label='B3')
    plt.plot([238,250],[y2,y2],linewidth=lw,color='pink',alpha=alf,label='B4')
    plt.plot([259,267],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='CBL2, loop 56')
    #plt.plot([252,259],[y2,y2],linewidth=lw,color='darkgreen',alpha=alf,label='B5')
    #plt.plot([268,276],[y2,y2],linewidth=lw,color='gold',alpha=alf,label='B6')
    plt.plot([307,314],[y1,y1],linewidth=lw,color='lightgreen',alpha=alf,label='end of loop 67')
    #plt.plot([316,321],[y2,y2],linewidth=lw,color='darkgreen',alpha=alf,label='B7')
    #plt.plot([342,349],[y2,y2],linewidth=lw,color='red',alpha=alf,label='B8')
    #plt.plot([322,341],[y2,y2],linewidth=lw,color='pink',alpha=alf,label='loop 78, K322 to F341')
    plt.plot([322,331],[y2,y2],linewidth=lw,color='pink',alpha=alf,label='first half of loop 78, K322 to D331')

    print(len(residue))
    print(first_res)
    print(len(residue)+first_res-1)
    """

# plot settings
plt.xticks(fontsize=16)
plt.xlabel('residue', fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel('contact frequency', fontsize=16) 
plt.legend(loc=2,fontsize=16)
#plt.title('Probability of contacts between protein and membrane lipids')
#plt.title(name)
plt.tight_layout()
plt.show()
