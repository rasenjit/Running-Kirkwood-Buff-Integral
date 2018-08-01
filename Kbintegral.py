#!/usr/bin/python/
# Python code for calculation of Kirkwood-Buff Integral using Ganguly corrected RDF and Kruger corrected KBI # 
# Ref: J. Phys. Chem. B 2018, 122, 5515-5526 # 
# Rasenjit Ghosh Dastidar # Lab-202, Dept of Chemistry, IISER Bhopal # Last modified Date: 06-07-18 # 

import numpy as np 
import math 
import sys 
from scipy.integrate import *
from scipy.interpolate import *
import matplotlib.pyplot as plt 

# reading the RDF file in XY format #
file = np.loadtxt(sys.argv[1])
r=[];gr=[]
for i in range(len(file)): 
    r.append(float(file[i,0]))
    gr.append(float(file[i,1]))

r=np.array(r); gr=np.array(gr)
Rmax = np.max(r)            # Maximum length of the solvation shell 

Nc = float(sys.argv[2])     # Number of cosolvent molecules 
boxl = float(sys.argv[3])   # length of sides of the cubic box 
outfile = open(sys.argv[4],'w') # output file to save the corrected KBIs 
vol = boxl**3               # Volume of the box 
Rho_c = Nc/vol              # Number Density of the cosolvent molecules

# Correction of RDF using Ganguly method (J. Chem. Theo. Comput 2013, 9, 1347) # 

dr = r[1]-r[0]  # width of radius 
DNc = Rho_c*4*math.pi*np.cumsum((gr-1.0)*dr*(r**2.0)) # calculation of number fluctuation of cosolvents 
Vol_frac = 1.0-(4.0*math.pi*(r**3.0)/3.0/vol) 
gr_cor = gr*Nc*Vol_frac/(Nc*Vol_frac-DNc-1.0) # corrected RDF 

gr_cor = np.array(gr_cor) 

# Computation of KBIntegral with correction by Kruger method (J. Phys. Chem. Lett 2013, 4, 235-238) # 

KBI = 4*math.pi*np.cumsum((1.0-(r/Rmax)**3.0)*(gr_cor-1.0)*dr*(r**2.0))  # KBI multiplied by (1.0-(r/Rmax)**3.0) : Kruger correction # 
KBI2 = 4*math.pi*np.cumsum((gr_cor-1.0)*dr*(r**2.0)) # corrected KBI using ganguly corrected RDF # 
KBI3 = 4*math.pi*np.cumsum((gr-1.0)*dr*(r**2.0)) # uncorrected  

# Plotting of KBIntegral of RDF # 

plt.plot(r,KBI,'r-',label='corrected: Ganguly+Kruger method')
plt.plot(r,KBI2,'b-',label='corrected: Ganguly method')
plt.plot(r,KBI3,'g-',label='uncorrected: with raw RDF')
plt.xlabel('r (Angstrom)')
plt.ylabel('RKBI (A**3)')
legend = plt.legend(loc='upper right')
frame = legend.get_frame()
frame.set_facecolor('0.90')
plt.show()

# Saving the corrected KBIntegral using both methods # 

k=0
while k<len(r):
    outfile.write(str(r[k])+'  '+str(KBI[k])+'\n')
    k=k+1 
outfile.close()

