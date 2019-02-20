##---------SIMPLIFIED CROSS SECTION---------

##This script finds the section properties of the simplified boom cross section, based on
##the actual found moments of inertia and centroid of the thin-walled cross section, and
##the applied moments on the cross section. The outputs are the boom areas at a given
##position on the aileron (according to stiffener numbering). All units are in meters.

#-----FUNCTION INPUTS-----
##This script requires as function inputs the aileron cross section geometry,
##the thin-walled moments of inertia and centroid, the sigma ratio function 
##and the local internal moments in z'and y' direction.
def aboom(Mz, My):
    from math import sqrt, pi
    from matplotlib import pyplot as plt
    from actmoi import actualmoi
    from sigmarat import sigratio
    from actcent import centactual
    from cs import crosssec

    c = 0.515  # aileron chord
    h = 0.248  # aileron height
    r = h / 2  # leading edge section radius
    le = sqrt((c - r) ** 2 + r ** 2)  # length of linear section
    circ = pi * r + 2 * le  # circumference of cross section
    phi = circ /11  # stiffener spacing
    Astiff = 5.4*10**(-5) #Stiffener point area (m^2)
    tsk = 0.0011 #skin thickness
    tsp = 0.0022 #spar thickness

    B = []

    nstiff, zpos, ypos = crosssec(11)
    Izz, Iyy = actualmoi(11)

    for i in range(len(nstiff)):
        if i != 12:
            sr1 = sigratio(Mz, My,zpos[i+1],zpos[i],ypos[i+1],ypos[i])
        sr2 = sigratio(Mz, My,zpos[i-1],zpos[i],ypos[i-1],ypos[i])
        if i==3:
            Bnew = Astiff+(tsk*sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2)/6)*(2+sr1)+(tsk*phi/6)*(2+sr2)
            B.append(Bnew)
        elif i==4:
            sr3 = sigratio(Mz, My,zpos[8],zpos[4],ypos[8],ypos[4])
            Bnew = (tsp*h/6)*(2+sr3)+(tsk*sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2)/6)*(2+sr2)+(tsk*(phi-sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2))/6)*(2+sr1)
            B.append(Bnew)
        elif i==5:
            Bnew = Astiff+(tsk*(phi-sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2)/6))*(2+sr2)+(tsk*2*phi/6)*(2+sr1)
            B.append(Bnew)
        elif i==6:
            Bnew = Astiff+(tsk*phi/6)*(2+sr2)+(tsk*phi/6)*(2+sr1)
            B.append(Bnew)
        elif i==7:
            Bnew = Astiff+(tsk*(phi-sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2)/6))*(2+sr1)+(tsk*phi/6)*(sr2)
            B.append(Bnew)
        elif i==8:
            sr3 = sigratio(Mz, My,zpos[4],zpos[8],ypos[4],ypos[8])
            Bnew = (tsp*h/6)*(2+sr3)+(tsk*sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2)/6)*(2+sr1)+(tsk*(phi-sqrt((ypos[4]-ypos[3])**2+(zpos[4]-zpos[3])**2))/6)*(2+sr2)
            B.append(Bnew)
        elif i==12:
            sr1 = sigratio(Mz, My,zpos[0],zpos[12],ypos[0],ypos[12])
            Bnew = Astiff + (tsk*phi/6)*(2+sr2) + (tsk*phi/6)*(2+sr1)
            B.append(Bnew)
        else:
            Bnew = Astiff + (tsk*phi/6)*(2+sr2) + (tsk*phi/6)*(2+sr1)
            B.append(Bnew)

    return B


