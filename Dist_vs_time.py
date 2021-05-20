import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

## user input
reps      = 25

HIGH_RES  = 1
SHOW_PLOT = 0
SAVE_PLOT = 1

## smooth funciton
def smooth(x,n):

    # ensure n is odd
    if (n % 2) == 0:
        n += 1

    # calculate wing size
    wing_size = int((n-1)/2)

    # initiate x_smooth (no need to change first and last value)
    x_smooth = x.copy() # brackets ensures list x is not edited
    #x_smooth = x_smooth.astype(float)
    j = 1
    rest = len(x)-j
    while j < len(x)-1:
        if (j < wing_size) or (rest-1 < wing_size) > j:
            i_max = j
        elif rest-1 < wing_size:
            i_max = rest-1
        else:
            i_max = wing_size
        sum = x[j]
        for i in range(1,i_max+1):
            sum += x[j-i] + x[j+i]
        x_smooth[j] = sum/(2*i_max+1)
        j += 1
        rest = len(x) - j

    return x_smooth

def exp_decay(t,A,B,D):

    return A*np.exp(-D*t)+B

# set fontsize
plt.rcParams.update({'font.size': 20})

# summary file
summary_filename = 'Dist_vs_time.txt'
f = open(summary_filename,'w')
f.write('%12s %12s %12s %12s %12s %12s %12s %12s\n' % ('Folder','Protein','D, PC','B, PC','D, PCPS', 'B, PCPS','D, PCPSPI','B, PCPSPI'))

# select system
for protein,protein_name in zip([15,18,22,23,25,27],['Smurf2','RIM2','KIBRA','PTEN','SHIP2','PI3KCa']):
#for protein in [15,18,22,25,27]:
    f.write('%12d %12s ' % (protein,protein_name))
    for lip in ['PC','PCPS','PCPSPIP2']:
#    for lip in ['PC']:
        plt.figure()
        file = '/sansom/s157/bioc1642/Desktop/prj_C2/%d/%s_%d/dist_min.xvg' % (protein,lip,0)
        time = np.genfromtxt(file,skip_header=24,usecols=[0,],unpack=True)
        dist_matrix = np.zeros((reps,len(time)))
#        for rep in range(25):
        for rep in range(reps):
            filename = '../../Seafile/C2_Manus/Dist/Dist_vs_time_%d_%s' % (protein,lip)
            file = '/sansom/s157/bioc1642/Desktop/prj_C2/%d/%s_%d/dist_min.xvg' % (protein,lip,rep)
            print('Import file:   %s' % file)
            time,dist = np.genfromtxt(file,skip_header=24,usecols=[0,1],unpack=True)
            plt.plot(time,dist,zorder=rep)
            dist_matrix [rep,:] = dist 
        avg_dist = np.mean(dist_matrix,axis=0)
        std_dist = np.std(dist_matrix,axis=0)
        std_mean_dist = std_dist/np.sqrt(reps)
#        plt.fill_between(time,avg_dist-std_mean_dist,avg_dist+std_mean_dist,alpha=1.0,edgecolor='black',facecolor='black',linewidth=1,zorder=rep+1) #label='average distance'
#        plt.plot(time,avg_dist,            linestyle='-',color='black',label='avg. dist.',zorder=rep+2)
#        plt.plot(time,smooth(avg_dist,121),linestyle='-',color='black',label='avg dist')
        
        # fit
        A0 = avg_dist[0]
        B0 = 0.5 
        D0 = 0.01
        print('A0 = %1.2f, B0 = %1.2f, D0 = %1.4f' % (A0,B0,D0))
        dist_0 = exp_decay(time,A0,B0,D0)
#        plt.plot(time,dist_0,color='green',label='exp, initial',zorder=rep+3)
        
        
        popt,pcov = curve_fit(exp_decay,time,avg_dist,p0=[A0,B0,D0])
        A,B,D = popt[0],popt[1],popt[2]
        print('A  = %1.2f, B  = %1.2f, D  = %1.4f' % (A,B,D))
        dist_fit = exp_decay(time,*popt)
#        plt.plot(time,dist_fit,color='red',label='decay rate = %1.3f' % D,zorder=rep+4)

        # print fit results to file
        f.write('%12.5f %12.3f ' % (D,B))

        plt.xlabel('Time [ns]')
        if lip == 'PC':
            plt.ylabel('Min prot-lip dist [nm]')
            lipname = lip
        elif lip == 'PCPS':
            lipname = 'PC:PS'
        else:
            lipname = 'PC:PS:PIP$_2$'
        plt.ylim([0,4.7])
        plt.title('%s' % lipname)
        plt.tight_layout()
        #plt.legend(prop={'size':16})
        
        if SAVE_PLOT:
            if reps == 25:
               filename = '../../Seafile/C2_Manus/Dist/Dist_vs_time_%d_%s' % (protein,lip)
            else:
               filename = '../../Seafile/C2_Manus/Dist/Dist_vs_time_%d_%s_%d' % (protein,lip,reps)
            print('Writing file:   %s' % filename)
            if HIGH_RES:
                plt.savefig('%s.png' % filename,format='png',dpi=600)
                print('HIGH_RES ON')
            else:
                plt.savefig(filename)

        if SHOW_PLOT:
            plt.show()
    f.write('\n')
f.close()

DC,BC,DS,BS,DI,BI = np.genfromtxt(summary_filename,skip_header=1,usecols=[2,3,4,5,6,7],unpack=True)

with open(summary_filename,'a') as f:
    f.write('%12s %12s %12.5f %12.3f %12.5f %12.3f %12.5f %12.3f\n' % ('mean',' ',np.mean(DC),np.mean(BC),np.mean(DS),np.mean(BS),np.mean(DI),np.mean(BI)))
    f.write('%12s %12s %12.5f %12.3f %12.5f %12.3f %12.5f %12.3f\n' % ('std dev',' ',np.std(DC),np.std(BC),np.std(DS),np.std(BS),np.std(DI),np.std(BI)))
