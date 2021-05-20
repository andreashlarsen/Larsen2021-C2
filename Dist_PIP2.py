import numpy as np
import matplotlib.pyplot as plt

Threshold = 0.6 # contact below this value, nm
contacts = []
Nreps=1
for i in range(Nreps):
    print('repetition %d' % i)
    filename = '15/PCPSPIP2_%d/dist_res_min_PIP2.xvg' % i
    time = np.genfromtxt(filename,skip_header=24,usecols=[0],unpack=True)
    Ntime = len(time)
    j = 1
    STOP = 0
    while STOP == 0:
        try:
            dist = np.genfromtxt(filename,skip_header=24,usecols=[j],unpack=True)
            idx = np.where(dist<Threshold)
            Ncontacts = len(dist[idx])
            if i == 0:
                contacts.append(Ncontacts)
            else:
                contacts[j-1] = Ncontacts + contacts[j-1]
            #print('res %d: %d/%d contacts' % (j,Ncontacts,len(dist)))
            j += 1
            #print(j)
        except:
            #print('no more residues')
            STOP = 1
#print(j)    
Nres = j-1
res = np.linspace(1,Nres,Nres)
total = Nres * Ntime
#print(res)
#print(contacts)
plt.bar(res,np.asarray(contacts)/total,color='black')
    
plt.xlabel(r'residue')
plt.ylabel(r'PIP$_2$ contact frequency')

plt.show()
