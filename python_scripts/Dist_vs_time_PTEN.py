import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from boxsize import get_area

## user input
reps      = 24

HIGH_RES  = 0
SHOW_PLOT = 1
SAVE_PLOT = 1

plt.rcParams["figure.figsize"] = (9,12)

# set fontsize
plt.rcParams.update({'font.size': 10})

## calculate surface charge densities
charge_CH = 0
charge_PC = 0
charge_PS = -1
charge_P1 = -3
charge_P2 = -5
charge_P3 = -7

"""
Area = 7*7 # surface area in nm
Area_PC = 54.7
Area_PCPS = 54.1
Area_CHPCPS = 50.6
Area_PCPSP1 = 51.2
Area_PCPSP1 = 51.2
Area_PCPSP3 = 51.2
Area_CHPCPSP2 = 51.0
Area_PCPS40 = 53.7
"""


list_A = get_area() #A_PC,A_PCPS,A_PCPS40,A_CHPCPS,A_PCPSPIP1,A_PCPSPIP2,A_PCPSPIP3,A_CHPCPSPIP2 = get_area()

print(list_A)

PC = 81*charge_PC/list_A[0]
PCPS = (64*charge_PC+16*charge_PS)/list_A[1]
PCPS40 = (48*charge_PC+32*charge_PS)/list_A[2]
PCPSP1 = (64*charge_PC+12*charge_PS+4*charge_P1)/list_A[4]
PCPSP2 = (64*charge_PC+12*charge_PS+4*charge_P2)/list_A[5]
PCPSP3 = (64*charge_PC+12*charge_PS+4*charge_P3)/list_A[6]
CHPCPS = (20*charge_CH+44*charge_PC+16*charge_PS)/list_A[3]
CHPCPSP2 = (20*charge_CH+44*charge_PC+12*charge_PS+4*charge_P2)/list_A[7]
CHPC = (20*charge_CH+60*charge_PC)/list_A[3]

print('PC = %s/nm^2' % PC)
print('PCPS = %s/nm^2' % PCPS)
print('PCPS40 = %s/nm^2' % PCPS40)
print('PCPSP1 = %s/nm^2' % PCPSP1)
print('PCPSP2 = %s/nm^2' % PCPSP2)
print('PCPSP3 = %s/nm^2' % PCPSP3)
print('CHPCPS = %s/nm^2' % CHPCPS)
print('CHPCPSP2 = %s/nm^2' % CHPCPSP2)
print('CHPC = %s/nm^2' % CHPC)

# summary file
#summary_filename = 'Dist_vs_time_PTEN.txt'
#f = open(summary_filename,'w')
#f.write('%12s %12s %12s %12s %12s %12s %12s %12s\n' % ('Folder','Protein','D, PC','B, PC','D, PCPS', 'B, PCPS','D, PCPSPI','B, PCPSPI'))

# select system
for protein,protein_name in zip([23],['PTEN']):
#for protein in [15,18,22,25,27]:
#    f.write('%12d %12s ' % (protein,protein_name))
    surface_charge_density = [PC,CHPC,PCPS,CHPCPS,PCPSP2,CHPCPSP2,PCPS40,PCPSP1,PCPSP3]
    lips = ['PC','CHPC','PCPS','CHPCPS','PCPSPIP2','CHPCPSPIP2','PCPS40','PCPSPIP1','PCPSPIP3']
#    surface_charge_density = [PC,PCPS,PCPSP2,CHPCPSP2,PCPS40,PCPSP1,PCPSP3]
#    lips = ['PC','PCPS','PCPSPIP2','CHPCPSPIP2','PCPS40','PCPSPIP1','PCPSPIP3']
#    surface_charge_density = [PCPS]
#    lips = ['PCPS']
#    lips = ['PC','PCPS','PCPSPIP2','CHPCPSPIP2','PCPS40']
    Nlips = len(lips)
    count=1
#    f,(p1,p2,p3,p4,p5,p6,p7) = plt.subplot(7,1,figsize=[10,20])
    for lip,scd in zip(lips,surface_charge_density):
#    for lip in ['CHPCPSPIP2']:
        p = plt.subplot(Nlips,1,count)
        count += 1
        file = '/sansom/s157/bioc1642/Desktop/prj_C2/%d/%s_%d/dist_min.xvg' % (protein,lip,0)
        time = np.genfromtxt(file,skip_header=24,usecols=[0,],unpack=True)
        dist_matrix = np.zeros((reps,len(time)))
        #if lip != 'PCPSPIP3':
        if lip != 'CHPCPS':
            p.set_xticklabels([])

#       for rep in range(25):
        if lip == 'PC':
            lipname = lip
        if lip == 'CHPC':
            lipname = 'CH:PC (25:75)'
        elif lip == 'PCPS':
            lipname = 'PC:PS (80:20)'
        elif lip == 'PCPSPIP2':
            lipname = 'PC:PS:PIP$_2$ (80:15:5)'
        elif lip == 'CHPCPSPIP2':
            lipname = 'CHOL:PC:PS:PIP$_2$ (25:55:15:5)'
            p.set_ylabel('Min prot-lip dist [nm]')
        elif lip == 'PCPS40':
            lipname = 'PC:PS (60:40)'
        elif lip == 'PCPSPIP1':
            lipname = 'PC:PS:PIP$_1$ (80:15:5)'
        elif lip == 'PCPSPIP3':
            lipname = 'PC:PS:PIP$_3$ (80:15:5)'
        elif lip == 'CHPCPS':
            lipname = 'CHOL:PC:PS (25:55:20)'
            p.set_xlabel('Time [ns]')
        for rep in range(0,reps):
            filename = '../../Seafile/C2_Manus/Dist/Dist_vs_time_PTEN_%d_%s' % (protein,lip)
            file = '/sansom/s157/bioc1642/Desktop/prj_C2/%d/%s_%d/dist_min.xvg' % (protein,lip,rep)
            print('Import file:   %s' % file)
            time,dist = np.genfromtxt(file,skip_header=24,usecols=[0,1],unpack=True)
            if rep==reps-1:
                p.plot(time,dist,zorder=rep,color='black',label='%s' % lipname)
            else:
                p.plot(time,dist,zorder=rep)
            dist_matrix [rep,:] = dist 
        avg_dist = np.mean(dist_matrix,axis=0)
        std_dist = np.std(dist_matrix,axis=0)
        std_mean_dist = std_dist/np.sqrt(reps)
        
        #print surface charge density
        p.text(1300,3.5,r'%s' % lipname)
        p.text(1300,2.4,r'surface charge density = %1.2f/nm$^2$' % scd)
        p.set_xlim(-100,2100)
        p.set_ylim(0,4.5)
        # fit
#        A0 = avg_dist[0]
#        B0 = 0.5 
#        D0 = 0.01
#        print('A0 = %1.2f, B0 = %1.2f, D0 = %1.4f' % (A0,B0,D0))
#        dist_0 = exp_decay(time,A0,B0,D0)
#        plt.plot(time,dist_0,color='green',label='exp, initial',zorder=rep+3)
        
#        popt,pcov = curve_fit(exp_decay,time,avg_dist,p0=[A0,B0,D0])
#        A,B,D = popt[0],popt[1],popt[2]
#        print('A  = %1.2f, B  = %1.2f, D  = %1.4f' % (A,B,D))
#        dist_fit = exp_decay(time,*popt)
#        plt.plot(time,dist_fit,color='red',label='decay rate = %1.3f' % D,zorder=rep+4)

        # print fit results to file
#        f.write('%12.5f %12.3f ' % (D,B))
        
        #p.legend()

#        p.set_ylim([0,4.7])
#        p.set_title('%s' % lipname)
        #plt.legend(prop={'size':16})

#        if SAVE_PLOT:
#            if reps == 25:
#               filename = '../../Seafile/C2_Manus/Dist/Dist_vs_time_%d_%s' % (protein,lip)
#            else:
#               filename = '../../Seafile/C2_Manus/Dist/Dist_vs_time_%d_%s_%d' % (protein,lip,reps)
#            print('Writing file:   %s' % filename)
#            if HIGH_RES:
#                plt.savefig('%s.png' % filename,format='png',dpi=600)
#                print('HIGH_RES ON')
#            else:
#                plt.savefig(filename)
#
    if SHOW_PLOT:
#        plt.tight_layout()
        plt.show()
#    f.write('\n')
#f.close()

if SAVE_PLOT:
    plt.savefig('../../Seafile/C2_Manus/Dist/Dist_vs_time_PTEN.pdf',format='pdf')

#DC,BC,DS,BS,DI,BI = np.genfromtxt(summary_filename,skip_header=1,usecols=[2,3,4,5,6,7],unpack=True)

#with open(summary_filename,'a') as f:
#    f.write('%12s %12s %12.5f %12.3f %12.5f %12.3f %12.5f %12.3f\n' % ('mean',' ',np.mean(DC),np.mean(BC),np.mean(DS),np.mean(BS),np.mean(DI),np.mean(BI)))
#    f.write('%12s %12s %12.5f %12.3f %12.5f %12.3f %12.5f %12.3f\n' % ('std dev',' ',np.std(DC),np.std(BC),np.std(DS),np.std(BS),np.std(DI),np.std(BI)))
