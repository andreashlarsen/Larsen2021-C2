import numpy as np
import matplotlib.pyplot as plt
from itertools import chain

HIGH_RES = 1
SHOW_PLOT = 1

f_min_x = open('file_min_x.txt','w')
f_min_x.write('# %-20s %-20s %-10s %20s\n' % ('protein_folder', 'protein_name', 'lipid', 'x pos for min_y'))

# set fontsize
plt.rcParams.update({'font.size': 18})

# user input
#for protein_folder in [22,27,28,18,23,25,15]:
for protein_folder in [25]:
#    for lip in ['PC','PCPS','PCPSPIP2']:
    for lip in ['PCPS','PCPSPIP2']:
        # default values
        skip_hist = 0
        umb_rep = 0
        protein_folder_actual = protein_folder

        # set params, depeding on input
        if lip == 'PCPSPIP2':
            lip_name = 'PC:PS:PIP$_2$'
            ylimits = [-26,-23.5]
        elif lip == 'PCPS':
            lip_name = 'PC:PS'
            ylimits = [-12,-8.5]
        elif lip == 'PC':
            lip_name = 'PC'

        if protein_folder == 22:
            protein_name = 'KIBRA'
            if lip == 'PCPSPIP2':
                choice = [10,8,6.0,11]
            elif lip == 'PCPS':
                choice = [11,11,5.9,15]
            elif lip == 'PC':
                choice = [25,14,5.8,17]
                skip_hist = 1
        elif protein_folder == 27:
            protein_name = 'PI3K'
            if lip == 'PCPSPIP2':
                choice = [4,8,5.8,7]
            elif lip == 'PCPS':
                choice = [14,8,5.8,10]
            elif lip == 'PC':
                choice = [19,10,5.9,2]
        elif protein_folder == 28:
            protein_name = 'PI3K_back'
            if lip == 'PCPSPIP2':
                choice = [5,8,5.9,19]
            elif lip == 'PCPS':
                choice = [13,6,6.1,6]
                skip_hist = 3
            elif lip == 'PC':
                choice = [12,11,5.8,22]
            protein_folder_actual = 27
        elif protein_folder == 18:
            protein_name = 'RIM2'
            if lip == 'PCPSPIP2':
                #choice = [32,9,5.5,11]
                choice = [8,9,6.0,11]
            elif lip == 'PCPS':
                #choice = [22,13,5.5,23]
                choice = [7,10,6.1,23]
            elif lip == 'PC':
                #choice = [38,12,6.0,10]
                choice = [9,12,6.0,10]
        elif protein_folder == 23:
            umb_rep = 1
            protein_name = 'PTEN'
            if lip == 'PCPSPIP2':
                choice = [36,9,6.0,1]
                skip_hist = 6
            elif lip == 'PCPS':
                choice = [39,8,6.2,19]
                skip_hist = 6
            elif lip == 'PC':
                choice = [49,10,5.9,3]
                skip_hist = 6
        elif protein_folder == 25:
            umb_rep = 1
            protein_name = 'SHIP2'
            if lip == 'PCPSPIP2':
                choice = [16,11,6.0,23]
                skip_hist = 1
            elif lip == 'PCPS':
                choice = [25,11,5.8,22]
                skip_hist = 3
            elif lip == 'PC':
                choice = [27,12,5.7,15]
        elif protein_folder == 15:
            umb_rep = 1
            protein_name = 'SMURF2'
            if lip == 'PCPSPIP2':
                choice = [34,15,6.0,24]
                skip_hist = 4
            elif lip == 'PCPS':
                choice = [37,15,5.9,1]
                skip_hist = 4
            elif lip == 'PC':
                choice = [51,10,6.0,23]
                skip_hist = 10

        skip_first = choice[0]
        skip_last = choice[1]
        plateau_start = choice[2]
        rep = choice[3]

        # import data
        file = '%s/%s_%s/umbrella_t1000_l5_k2000_n%d/bsResult.xvg' % (protein_folder_actual,lip,rep,umb_rep)
        print(file)
        skip_header = 18+skip_first
        skip_footer = skip_last
        x,y,dy = np.genfromtxt(file,skip_header=skip_header,skip_footer=skip_footer,usecols=[0,1,2],unpack=True)

        # import hist
        file = '%s/%s_%s/umbrella_t1000_l5_k2000_n%d/histo.xvg' % (protein_folder_actual,lip,rep,umb_rep)
        print(file)
        dist = np.genfromtxt(file,skip_header=17,usecols=[0],unpack=True)

        file_tpr = '%s/%s_%s/umbrella_t1000_l5_k2000_n%d/tpr_files_umbrella_trunc.dat' % (protein_folder_actual,lip,rep,umb_rep)
        print(file_tpr)
        x_hist = np.genfromtxt(file,skip_header=17,usecols=[0],unpack=True)

        Nfiles = sum(1 for line in open(file_tpr))
        hist = np.zeros((Nfiles,len(dist)))
        w_mean = np.zeros((Nfiles,1))
        for i in range(Nfiles):
            h = np.genfromtxt(file,skip_header=17,usecols=[i+1],unpack=True)
            h /= np.amax(h)
            hist[i,:] = h
            w_mean[i] = np.sum(h*x_hist)/np.sum(h)
        ndx_w_mean = int(np.where(w_mean==np.amin(w_mean))[0])

        # find minimum
        min_y = np.amin(y)
        indx_min_y = np.where(y==min_y)
        d_min_y = dy[indx_min_y]
        f_min_x.write('  %-20d %-20s %-10s %20.1f\n' % (protein_folder, protein_name, lip, x[indx_min_y]))
        print('PMF_min = %1.2f +/- %1.2f' % (min_y,d_min_y))

        # find max by user defined interval
        indx_plateau = np.where(x>plateau_start)
        plateau_x = x[indx_plateau]
        plateau_y = y[indx_plateau]
        plateau_y_mean = np.mean(plateau_y)

        plateau_dy = dy[indx_plateau]

        d_max_y = np.sqrt(np.sum(plateau_dy**2))/len(plateau_dy)
        d_max_y_2 = np.mean(plateau_dy)

        print('PMF_max = %1.2f +/- %1.2f' % (plateau_y_mean,d_max_y_2))

        # determine PMF and d_PMF
        PMF = plateau_y_mean - min_y
        d_PMF = np.sqrt(d_min_y**2 + d_max_y_2**2)
        print('PMF = %1.1f +/- %1.1f' % (PMF,d_PMF))

        # set plateau to zero for plotting
        y -= plateau_y_mean
        min_y -= plateau_y_mean

        # plot data
        fig,ax = plt.subplots(2,1,gridspec_kw={'height_ratios': [9,2]})
        xlimits = [3.3,3.8]

        ax[0].plot(x,y,color='#1B2ACC')
        ax[0].fill_between(x,y-dy,y+dy,alpha=1,color='skyblue') # .eps format does not support transparancy, so alpha should be 1
        ax[0].set_xlim(xlimits)
        ax[0].plot(xlimits,np.array([1.0,1.0])*min_y,linestyle='--',color='black')
        ax[0].plot(xlimits,[0.0,0.0],linestyle='--',color='black')
        #ax[0].axvline(x=plateau_start,color='red')
        plateau_end = np.amax(x)
        ax[0].set_xticks([])
        ax[0].set_title(r'%s' % lip_name)
        
        ax[0].plot([x[indx_min_y],x[indx_min_y]],ylimits,linestyle='-',color='#1B2ACC')
        if protein_folder == 25:
            if lip == 'PC':
                # plot vertical line for minimum
                ax[0].plot([x[indx_min_y],x[indx_min_y]],ylimits,linestyle='-',color='#1B2ACC')    
            
            else:
                if lip == 'PCPS':
                    productive_mode_x = 3.672
                else:
                    #productive_mode_x = 4.256
                    productive_mode_x = 3.622
                f_min_x.write('  %-7d (productive) %-20s %-10s %20.1f\n' % (protein_folder, protein_name, lip, productive_mode_x))

                # plot vertical lines for two minima
                ax[0].plot([x[indx_min_y],x[indx_min_y]],ylimits,linestyle='-',color='darkred')
                ax[0].plot([productive_mode_x,productive_mode_x],ylimits,linestyle='-',color='#1B2ACC')
                indx_min_y_prod = np.argmax(x>productive_mode_x) # find index of first element greater than productive_mode_x 
                #print(indx_min_y_prod)
                min_y_prod = y[indx_min_y_prod]
                PMF_prod = -min_y_prod
                d_min_y_prod = dy[indx_min_y_prod]
                d_PMF_prod = np.sqrt(d_min_y_prod**2 + d_max_y_2**2)
                print('PMF (productive) = %1.1f +/- %1.1f' % (PMF_prod,d_PMF_prod))
        else:
            # plot vertical line for minimum
            ax[0].plot([x[indx_min_y],x[indx_min_y]],ylimits,linestyle='-',color='#1B2ACC')
        ax[0].set_ylim(ylimits)

        list = chain(range(0,ndx_w_mean+1-skip_hist),range(ndx_w_mean+1,Nfiles))

        for i in list:
            ax[1].plot(dist,hist[i,:])
        ax[1].set_xlim(xlimits)
        ax[1].set_ylim([-0.1,1.2])
        if protein_folder == 25:
            ax[1].set_xlabel('Distance [nm]')
        #if lip == 'PC':
        ax[0].set_ylabel('E [kJ/mol]')
        ax[1].set_ylabel('Count')

        fig.tight_layout()

        figure_name = '../../Seafile/C2_Manus/PMF/ZOOM_PMF_%s_%s' % (protein_name, lip)
        print('save figure: %s' % figure_name)

        if HIGH_RES:
            print('high_res figure ON')
            plt.savefig('%s.eps' % figure_name,format='eps')
        else:
            plt.savefig(figure_name)
        
        if SHOW_PLOT:
            plt.show()

f_min_x.close()
