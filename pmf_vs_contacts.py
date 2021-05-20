import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

names = [r'KIBRA',r'PI3KC2$\alpha$ (front)',r'PI3KC2$\alpha$ (back)',r'RIM2',r'PTEN',r'SHIP2',r'Smurf2']
PMF = -np.asarray([48.0,64.1,115.2,120.4,77.7,23.1,81.5])
dPMF = [1.4,1.2,2.6,2.4,4.3,1.0,2.0]

Ncontacts = np.asarray([4,6,5,8,6,3,6.5])

def func(Ncontacts,a,b):
    y = a*Ncontacts+b
    return y

popt,pcov = curve_fit(func,Ncontacts,PMF,sigma=dPMF)
xmin,xmax=0,10
x = np.linspace(xmin,xmax,5)
PMF_fit = func(x,*popt)

fig,ax = plt.subplots()
ax.errorbar(Ncontacts,PMF,yerr=dPMF,linestyle='none',marker='.',markersize=10,color='blue')
ax.plot(x,PMF_fit,color='darkgrey')

ax.set_xlabel(r'Number of PIP$_2$/(Arg,Lys) contacts')
ax.set_ylabel('PMF well depth [kJ/mol]')
ax.set_xlim(2,9)
ax.set_ylim(-140,0)
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top') 

for i in range(len(PMF)):
    if names[i] == 'PTEN':
        displace = -0.61
    elif names[i] == 'KIBRA':
        displace = -0.68
    else:
        displace = 0.13
    plt.text(Ncontacts[i]+displace,PMF[i]-1.9,names[i])

plt.savefig('pmf_vs_contacts.pdf',format='pdf')
plt.show()
