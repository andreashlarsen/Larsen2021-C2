import numpy as np

number_of_rep = 5
number_of_lip = 4
FEP_Bilayer_mean   = 254.55

f = open('FEP_energies.dat','w')
f.write('# FEP energies calculated from %d repeats\n' % number_of_rep)
f.write('# \n')
f.write('# %-12s  %-13s  %-13s  %-13s  %-13s  %-13s\n' % ('Protein','lip1','lip2','lip3','lip4','mean'))

for prot in ['KIBRA','PI3K','PI3K_back','RIM2','PTEN','SHIP2','SMURF2']:
    print('protein  : %s' % prot)
    FEP_mean            = np.zeros(number_of_lip)
    FEP_std             = np.zeros(number_of_lip)
    FEP_std_mean_indep  = np.zeros(number_of_lip)
    FEP_sub             = np.zeros(number_of_lip)
    f.write('  %-12s' % prot)
    for lip in range(number_of_lip):
        print('  lip    : %d' % lip)
        FEP = np.zeros(number_of_rep)
        dFEP = np.zeros(number_of_rep)
        for rep in range(number_of_rep):
            print('    rep  : %d' % rep)
            file = '%s_1000ns/ANALYSIS/lip%d_rep%d/results.txt' % (prot,lip+1,rep)
            #print(file)

            FEP[rep],dFEP[rep] = np.genfromtxt(file,skip_header=28,usecols=[16,18],unpack=True)
            #print('prot = %s; lip = %d; rep = %d;FEP = %1.2f +/- %1.2f' % (prot,lip,rep,FEP[rep],dFEP[rep]))
        FEP_std[lip] = np.std(FEP)
        FEP_std_mean_indep[lip] = FEP_std[lip]/np.sqrt(number_of_rep)
        FEP_sub[lip] = np.mean(FEP) - FEP_Bilayer_mean
        #FEP_std_mean_indep = FEP_std/np.sqrt(number_of_rep)
        #print('prot = %s; lip = %d:   FEP_mean = %1.1f, FEP_std = %1.1f, FEP_std_mean = %1.1f' % (prot,lip,FEP_mean[lip],FEP_std[lip],FEP_std_mean_indep[lip]))
        #f.write('  %-12s %-5d %8.1f +/- %-8.1f\n' % (prot,lip[lip],FEP_mean[lip],FEP_std[lip]))
        f.write('  %4.1f +/- %-4.1f' % (FEP_sub[lip],FEP_std[lip]))
    FEP_sum = np.sum(FEP_sub)
    dFEP_sum = np.sum(FEP_std)
    dFEP_sum_indep = np.sqrt(np.sum(FEP_std**2))
    f.write('  %4.1f +/- %-4.1f ( %-4.1f)\n' % (FEP_sum,dFEP_sum,dFEP_sum_indep))
f.close()
