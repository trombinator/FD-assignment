# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 14:01:24 2019

@author: Boris Englebert
"""

#------------------BOOM AREA----------------
##-----------------

from math import*
from matplotlib import pyplot as plt
from cs import crosssec

#-------------INPUT----------
##NOTE: USE CROSS SECTION FUNCTION FROM CROSS SECTION SCRIPT TO FIND BOOM COORDS

c = 0.515  # aileron chord (m)
h = 0.248  # aileron height (m)
r = h / 2  # leading edge section radius (m)
le = sqrt((c - r) ** 2 + r ** 2)  # length of linear section (m)
circ = pi * r + 2 * le  # circumference of cross section (m)
phi = circ /11  # stiffener spacing (m)
tsk = 0.0011 #Skin thickness (m)
tsp = 0.0022 #spar thickness (m)
Astiff = 5.4*10**(-5) #Stiffener point area (m^2)

B = []  #BOOM AREAS (m^2)



#--------FUNCTIONS-------

#CROSS SECTION COORDS (YPOS AND ZPOS)

nstiff, zpos, ypos = crosssec(11)


i = 0

for i in range(len(zpos)):
    if i==4 or i==8:
        Bnew = (tsp*h/6)+(tsk*sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2)/6)*(2+ypos[4]/ypos[3])+(tsk*(phi-sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2))/6)*(2+ypos[5]/ypos[4])
        B.append(Bnew)
    elif i==6:
        B.append(Astiff)
    elif i==5 or i==7:
        Bnew = Astiff+(tsk*sqrt((ypos[5]-ypos[4])**2+(zpos[5]-zpos[4])**2)/6)*(2+ypos[5]/ypos[4])+(tsk*2*phi/6)
        B.append(Bnew)
    elif i==3 or i==9:
        Bnew = Astiff+(tsk*sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2)/6)*(2+ypos[4]/ypos[3])+(tsk*phi/6)*(2+ypos[3]/ypos[2])
        B.append(Bnew)
    elif i==0 or i==12:
        Bnew = Astiff + (tsk*phi/6)*(2+ypos[0]/ypos[12]) + (tsk*phi/6)*(2+ypos[1]/ypos[0])
        B.append(Bnew)
    elif i==1 or i==11:
        Bnew = Astiff + (tsk*phi/6)*(2+ypos[1]/ypos[0]) + (tsk*phi/6)*(2+ypos[2]/ypos[1])
        B.append(Bnew)
    elif i==2 or i==10:
        Bnew = Astiff + (tsk*phi/6)*(2+ypos[2]/ypos[1]) + (tsk*phi/6)*(2+ypos[3]/ypos[2])
        B.append(Bnew)
    print (i)      

print(B, len(zpos))

plt.subplot(121)
plt.plot(zpos,ypos, "-ro")
plt.gca().invert_xaxis()
plt.axis("equal")
plt.grid(linestyle='-', linewidth='0.5')

plt.subplot(122)
plt.plot(nstiff,B)

plt.show()