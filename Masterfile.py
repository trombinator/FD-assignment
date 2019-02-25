#--------------MASTER FILE NUMERICAL TOOL-------

#This file will provide the LE and TE deflection and the maximum shear stress in the ribs.
#Input functions are the twistrate, the shear and moment distribution, the beam deflection
#and the shear flow in the ribs. All units in meters(Standard SI).

from math import*
from shearflows import twistrate
import numpy as np
from intshearmom import deflections, sheardiagram, momentdiagram
from ribabd import ribABD
from ribc import ribC


n = 50 #number of sections (Even number!)

c = 0.515  # aileron chord
h = 0.248  # aileron height
r = h / 2  # leading edge section radius
le = sqrt((c - r) ** 2 + r ** 2)  # length of linear section
circ = pi * r + 2 * le  # circumference of cross section
phi = circ / 11  # stiffener spacing
Astiff = 5.4 * 10 ** (-5)  # Stiffener point area (m^2)
tsk = 0.0011  # skin thickness
tsp = 0.0022  # spar thickness
G = 28 * 10 ** 9  # Shear modulus Al-2024 T3
A1 = 0.5 * h * (c - r)  # Area of section 1
A2 = 0.5 * pi * r ** 2  # Area of section 2
theta = radians(25)  # Aileron deflection angle (rad)
la = 2.691 #Aileron span
lsec = la/n



P =  20.6 * 10**3

x = np.linspace(0,la,n)
v = []
for i in x:
    v.append(deflections(i)[0])



Sy = []
Sz = []
for j in x:
    Sy.append(sheardiagram(j)[0])
    Sz.append(sheardiagram(j)[1])

My = []
Mz = []
for k in x:
    My.append(momentdiagram(k)[0])
    Mz.append(momentdiagram(k)[1])

twist1 = []
tw1 = 0
for f in range(0,20):
    ind = 19-f
    rot = twistrate(Sz[ind],Sy[ind],Mz[ind],My[ind],P)
    if f==0:
        twist1.append(0)
    else:
        tw = rot*lsec
        tw1 += tw
        twist1.append(tw1)

twist2 = []
tw2 = 0
for f in range(20,len(x)):
    rot = twistrate(Sz[f],Sy[f],Mz[f],My[f],P)
    tw = rot*lsec
    tw2 += tw
    twist2.append(tw2)

twist = twist1[::-1]+twist2

dLE = []
dTE = []

for a in range(len(twist)):
    dle = -r*sin(twist[a])*sin(theta)
    dte = (c-r)*sin(twist[a])*sin(theta)
    dLE.append(dle)
    dTE.append(dte)



defle = []
defte = []
for b in range(len(v)):
    def1 = dLE[b]+v[b]
    def2 = dTE[b]+v[b]
    defle.append(def1)
    defte.append(def2)


print("Max Leading edge deflection =", max(defle), "Max trailing edge deflection =",max(defte))




###-------SHEAR FLOW IN RIBS----------

qmax = []

x1 = [0.174]
Sy1 = []
Sz1 = []
for j in x1:
    Sy1.append(sheardiagram(j)[0])
    Sz1.append(sheardiagram(j)[1])

qmaxA = ribABD(Sz1[0],Sy1[0])
print("qmax in rib A =", qmaxA)

x2 = [1.051-0.015]
Sy2 = []
Sz2 = []
for j in x2:
    Sy2.append(sheardiagram(j)[0])
    Sz2.append(sheardiagram(j)[1])

qmaxB = ribABD(Sz2[0],Sy2[0])
print("qmax in rib B =", qmaxB)

x3 = [1.051+0.015]
Sy3 = []
Sz3 = []
for j in x3:
    Sy3.append(sheardiagram(j)[0])
    Sz3.append(sheardiagram(j)[1])

qmaxC = ribC(Sz3[0],Sy3[0])
print("qmax in rib C =", qmaxC)

x4 = [2.512]
Sy4 = []
Sz4 = []
for j in x4:
    Sy4.append(sheardiagram(j)[0])
    Sz4.append(sheardiagram(j)[1])

qmaxD = ribC(Sz4[0],Sy4[0])
print("qmax in rib D =", qmaxD)

















