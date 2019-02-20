##---------SIGMA RATIO-------

##This script finds the sigma ratio that can be used to find the simplified
##boom cross section. It requires as inputs both the Iz'z' and Iy'y' of the
##thin-walled section, the thin-walled centroid, and the moments about z'
##and y'. All units in meters.

#------FUNCTION INPUTS------
##Inputs are the functions that find the thin-walled moments of inertia
##and cross section

def sigratio(Mz,My,z2,z1,y2,y1):
    from cs import crosssec
    from actmoi import actualmoi
    from actcent import centactual

    Mz = 10000 #moment in z' direction (Nm)
    My = 10000 #moment in y' direction (Nm)


    Izz, Iyy = actualmoi(11)
    zc = centactual(11)

    sigrat = (Izz*My*(z2-zc)+Iyy*Mz*y2)/(Izz*My*(z1-zc)+Iyy*Mz*y1)
    return sigrat

