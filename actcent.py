#---------------ACTUAL CENTROID----------



##This script finds the actual centroid of the cross section based on the
##stiffener point areas and thin walled sections. This centroid is later used
##in the computation of the sigma ratio to find the boom areas. After this, the
##centroid is updated based on the newly found boom areas. All units in meters.

#-----FUNCTION INPUTS------
##Requires the cross section geometry as input
def centactual(m):
    from math import sqrt, pi
    from cs import crosssec
    c = 0.515  # aileron chord
    h = 0.248  # aileron height
    r = h / 2  # leading edge section radius
    tsk = 0.0011##skin thickness
    tsp = 0.0022 #spar thickness
    le = sqrt((c - r) ** 2 + r ** 2)  # length of linear section
    circ = pi * r + 2 * le  # circumference of cross section
    phi = circ /m # stiffener spacing
    Astiff = 5.4*10**(-5) #Stiffener point area (m^2)

    Alin = le*tsk
    zlin = le*0.5*(c-r)

    Acirc = pi*r*tsk
    zcirc = (c-r)+2*r/pi

    Asp = h*tsp
    zsp = (c-r)

    nstiff, zpos, ypos = crosssec(11)

    zA = []
    A = []

    for n in range(len(zpos)):
        if n==4 or n==8:
            break
        else:
            za = zpos[n]*Astiff
            zA.append(za)
            A.append(Astiff)


    zcent = (sum(zA)+2*Alin*zlin+Acirc*zcirc+Asp*zsp)/(sum(A)+2*Alin+Acirc+Asp)
    return zcent
