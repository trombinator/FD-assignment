##------------------ACTUAL MOI----------


##This script finds the actual moments of inertia Iz'z' and Iy'y' of the cross section, based on the
##stiffener point areas (that only contribute due to the Steiner terms) and the angled thin wall
##rectangular sections and the thin-walled semicircle. These moments of inertia can then be used
##to find the sigma ratios and as such the simplified section. All units in meters.


#-------FUNCTION INPUTS-------
##This module requires for both the actual centroid and the cross section geometry to be imported
def actualmoi(m):
    from math import sqrt, pi
    from cs import crosssec
    from actcent import centactual

    c = 0.515  # aileron chord
    h = 0.248  # aileron height
    r = h / 2  # leading edge section radius
    le = sqrt((c - r) ** 2 + r ** 2)  # length of linear section
    circ = pi * r + 2 * le  # circumference of cross section
    phi = circ /m  # stiffener spacing
    Astiff = 5.4*10**(-5) #Stiffener point area (m^2)
    tsk = 0.0011 #skin thickness
    tsp = 0.0022 #spar thickness

    zc = centactual(11)
    nstiff, zpos, ypos = crosssec(11)

    #-----Thin walled sections Izz-----
    Izz1 = le*tsk*r**2/12+le*tsk*(0.5*(c-r)-zc)**2 #Izz of linear sections
    Izz2 = pi*r**3*tsk/2+pi*r*tsk*((c-r)+2*r/pi-zc)**2 #Izz of semicircle
    Izzsp = tsp*h**3/12 #Izz of spar

    #-----Thin walled section Iyy-----
    Iyy1 = le*tsk*(c-r)**2/12+(0.5*r)**2 #Iyy of linear section
    Iyy2 = pi*r**3*tsk/2+pi*r*tsk*((c-r)) #Iyy of semicircle
    Iyysp = tsp*h*((c-r)-zc)**2

    Izzst = []
    Iyyst = []

    for n in range(len(zpos)):
        if n==4 or n==8:
            break
        else:
            izz = Astiff*(zpos[n]-zc)**2
            Izzst.append(izz)

    for k in range(len(zpos)):
        if k==4 or k==8:
            break
        else:
            iyy = Astiff*(ypos[k])**2
            Iyyst.append(iyy)

    Izztot = sum(Izzst)+2*Izz1+Izz2+Izzsp
    Iyytot = sum(Iyyst)+2*Iyy1+Iyy2+Iyysp
    return Izztot, Iyytot

