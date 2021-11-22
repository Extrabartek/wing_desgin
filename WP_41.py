import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import interpolate

# ---------------- Values & Constants ------------------
# Values known from previous work packages, and constants needed
cr = 6.46
ct = 1.83464
b = 31.1
rho = 1.225
V = 256
q = 0.5 * rho * V ** 2

# Reading the data files and making an array with columns y (wing span [m]), cl, cd, cm
data1 = "MainWing_a=0.00_v=10.00ms.txt"
data2 = "MainWing_a=10.00_v=10.00ms.txt"
data3 = "MainWing_a=-10.00_v=10.00ms.txt"
list1 = np.genfromtxt(data1, dtype='float', skip_header=40, skip_footer=1029, usecols=(0, 3, 5, 7))
list2 = np.genfromtxt(data2, dtype='float', skip_header=40, skip_footer=1029, usecols=(0, 3, 5, 7))
list3 = np.genfromtxt(data3, dtype='float', skip_header=40, skip_footer=1029, usecols=(0, 3, 5, 7))

# Making lists for values of y (wing span [m]), cl, cd, cm at each angle of attack (AoA)
y = []
cl1 = []
cd1 = []
cm1 = []
cl2 = []
cd2 = []
cm2 = []
cl3 = []
cd3 = []
cm3 = []
for i in range(19):
    y.append(list1[i][0])
    cl1.append(list1[i][1])
    cd1.append(list1[i][2])
    cm1.append(list1[i][3])
    cl2.append(list2[i][1])
    cd2.append(list2[i][2])
    cm2.append(list2[i][3])
    cl3.append(list3[i][1])
    cd3.append(list3[i][2])
    cm3.append(list3[i][3])

# Interpolating cl, cd, cm to find a function describing each cl, cd, cm, at each AoA
cly1 = sp.interpolate.interp1d(y, cl1, kind='cubic', fill_value="extrapolate")
cdy1 = sp.interpolate.interp1d(y, cd1, kind='cubic', fill_value="extrapolate")
cmy1 = sp.interpolate.interp1d(y, cm1, kind='cubic', fill_value="extrapolate")

cly2 = sp.interpolate.interp1d(y, cl2, kind='cubic', fill_value="extrapolate")
cdy2 = sp.interpolate.interp1d(y, cd2, kind='cubic', fill_value="extrapolate")
cmy2 = sp.interpolate.interp1d(y, cm2, kind='cubic', fill_value="extrapolate")

cly3 = sp.interpolate.interp1d(y, cl3, kind='cubic', fill_value="extrapolate")
cdy3 = sp.interpolate.interp1d(y, cd3, kind='cubic', fill_value="extrapolate")
cmy3 = sp.interpolate.interp1d(y, cm3, kind='cubic', fill_value="extrapolate")


# -----------------Functions-------------------
def c(y):
    '''
    This function describes the variation of the chord along the wingspan.
    :param y: Wing span [m]
    :type y: float
    :return: Chord length [m]
    :rtype: float
    '''
    b = 31.1
    c = cr - (cr - ct) * 2 * y / b
    return c


def Ldis(y, cly):
    """
    This function describes the lift distribution
    :param y: Wing span [m]
    :type y: float
    :param cly: Interpolated function of the lift coefficient versus wing span
    :type cly: array
    :return: Lift distribution at a certain point along the wing span
    :rtype: array
    """
    Ldis = c(y) * q * cly(y)
    return Ldis


def Ddis(y, cdy):
    """
    This function describes the drag distribution
    :param y: Wing span [m]
    :type y: float
    :param cdy: Interpolated function of the drag coefficient versus wing span
    :type cdy: array
    :return: Drag distribution at a certain point along the wing span
    :rtype: array
    """
    Ddis = c(y) * q * cdy(y)
    return Ddis


def Mdis(y, AoA):
    """
    This function describes the moment distribution
    :param y: Wing span [m]
    :type y: float
    :param AoA: Angle of attack [deg]
    :type AoA: float
    :return: Moment distribution at a certain point along the wing span, based upon the angle of attack
    :rtype: array
    """
    Mdis = 0
    if AoA == 0:
        Mdis = (c(y)) **2 * q * cmy1(y)
    if AoA == 10:
        Mdis = (c(y)) **2 * q * cmy2(y)
    if AoA == -10:
        Mdis = (c(y) **2 * q * cmy3(y))
    return Mdis


def Ndis0(y, AoA):
    """
    This function describes the normal force distribution (contribution of lift and drag at an angle)
    :param y: Wing span [m]
    :type y: float
    :param AoA: Angle of attack [deg]
    :type AoA: float
    :return: Normal force distribution at a certain point along the wing span, based upon the angle of attack
    :rtype: array
    """
    Ndis = 0
    if AoA == 0:
        Ndis = Ldis(y, cly1) * np.cos(np.radians(0)) + Ddis(y, cdy1) * np.sin(np.radians(0))
    if AoA == 10:
        Ndis = Ldis(y, cly2) * np.cos(np.radians(10)) + Ddis(y, cdy2) * np.sin(np.radians(10))
    if AoA == -10:
        Ndis = Ndis = Ldis(y, cly3) * np.cos(np.radians(-10)) + Ddis(y, cdy3) * np.sin(np.radians(-10))
    return Ndis


# ------------------- Distribution calculations ------------------
range = np.linspace(0, 15.55, 500)  # Wing span [m]

L1 = Ldis(range, cly1)              # Lift distribution at AoA of 0 [deg]
D1 = Ddis(range, cdy1)              # Drag distribution at AoA of 0 [deg]
M1 = Mdis(range, cmy1)              # Moment distribution at AoA of 0 [deg]

L2 = Ldis(range, cly2)              # Lift distribution at AoA of 10 [deg]
D2 = Ddis(range, cdy2)              # Drag distribution at AoA of 10 [deg]
M2 = Mdis(range, cmy2)              # Moment distribution at AoA of 10 [deg]

# --------- Effect of velocity and angle of attack----------
# Following values result from XFLR5
CL0 = 0.176815
CL10 = 0.970715
CD0 = 0.001319
CD10 = 0.038922
CM0 = -0.267948
CM10 = -1.16705

# Common values used for CL, CD, CM
CLd = 1
CDd = 0.035
CMd = 0.00000005


def CLdy(y):
    '''
    :param y: Wing span [m]
    :type y: float
    :return: Lift coefficient at a certain point along the wing span
    :rtype: float
    '''
    CLdyvalue = cly1(y) + (CLd - CL0) / (CL10 - CL0) * (cly2(y) - cly1(y))
    return CLdyvalue


def CDdy(y):
    '''
    :param y: Wing span [m]
    :type y: float
    :return: Drag coefficient at a certain point along the wing span
    :rtype: float
    '''
    CDdyvalue = cdy1(y) + (CDd - CD0) / (CD10 - CD0) * (cdy2(y) - cdy1(y))
    return CDdyvalue


def CMdy(y):
    '''
    :param y: Wing span [m]
    :type y: float
    :return: Moment coefficient at a certain point along the wing span
    :rtype: float
    '''
    CMdyvalue = cmy1(y) + (CMd - CM0) / (CM10 - CM0) * (cmy2(y) - cmy1(y))
    return CMdyvalue


def alphad():
    '''
    Function that gives the design angle of attack for a design lift coefficient.
    :return: Design angle of attack
    :rtype: float
    '''
    alphadvalue = np.degrees(np.arcsin((CLd - CL0) / (CL10 - CL0) * np.sin(np.radians(10))))
    return alphadvalue

'''
# --------------------- Shear moment and bending diagrams -----------------------
def Ny(AoA):
    
    :param AoA: Angle of attack [deg]
    :type AoA: float
    :return: Normal force distribution
    :rtype: array
    
    if AoA == 0:
        L = L1
        D = D1
    else:
        L = L2
        D = D2
    return np.cos(np.radians(AoA)) * L + np.sin(np.radians(AoA)) * D


def Ty(AoA):
    
    :param AoA: Angle of attack [deg]
    :type AoA: float
    :return: Tangential force distribution
    :rtype: array
    
    if AoA == 0:
        L = L1
        D = D1
    else:
        L = L2
        D = D2
    return np.sin(np.radians(AoA)) * L + np.cos(np.radians(AoA)) * D
'''
