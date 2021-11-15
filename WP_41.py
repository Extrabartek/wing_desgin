import matplotlib.pyplot as plt
import scipy as sp
from scipy import interpolate
from scipy import integrate
import numpy as np
from io import StringIO

# ----------------Constants------------------
cr = 6.46
ct = 1.83464
b = 31.1
rho = 1.225
V = 256
q = 0.5 * rho * V ** 2

# Reading file and making an array with columns y, cl, cd, cm
data1 = "MainWing_a=0.00_v=10.00ms.txt"
data2 = "MainWing_a=10.00_v=10.00ms.txt"
list1 = np.genfromtxt(data1, dtype='float', skip_header=40, skip_footer=1029, usecols=(0, 3, 5, 7))
list2 = np.genfromtxt(data2, dtype='float', skip_header=40, skip_footer=1029, usecols=(0, 3, 5, 7))

# Making individual lists for y, cl, cd, cm
y = []
cl1 = []
cd1 = []
cm1 = []
cl2 = []
cd2 = []
cm2 = []
for i in range(19):
    y.append(list1[i][0])
    cl1.append(list1[i][1])
    cd1.append(list1[i][2])
    cm1.append(list1[i][3])
    cl2.append(list2[i][1])
    cd2.append(list2[i][2])
    cm2.append(list2[i][3])

cltot1 = sum(cl1)
cltot2 = sum(cl2)

# Interpolating cl, cd, cm
cly1 = sp.interpolate.interp1d(y, cl1, kind='cubic', fill_value="extrapolate")
cdy1 = sp.interpolate.interp1d(y, cd1, kind='cubic', fill_value="extrapolate")
cmy1 = sp.interpolate.interp1d(y, cm1, kind='cubic', fill_value="extrapolate")

cly2 = sp.interpolate.interp1d(y, cl2, kind='cubic', fill_value="extrapolate")
cdy2 = sp.interpolate.interp1d(y, cd2, kind='cubic', fill_value="extrapolate")
cmy2 = sp.interpolate.interp1d(y, cm2, kind='cubic', fill_value="extrapolate")


# -----------------Functions-------------------
def c(y):
    '''

    :param y: y distance along the wingspan starting from root
    :return: chordlength
    '''
    b = 31.1
    c = cr - (cr - ct) * 2 * y / b
    return c


def Ldis(y, cly):
    '''

    :param y:
    :param cly:
    :return:
    '''
    Ldis = c(y) * q * cly(y)
    return Ldis


print(Ldis(5, cly1))  # Del later


def Ddis(y, cdy):
    '''

    :param y:
    :param cdy:
    :return:
    '''
    Ddis = c(y) * q * cdy(y)
    return Ddis


def Mdis(y, AoA):
    '''

    :param y:
    :param cmy:
    :return:
    '''
    Mdis = 0
    if AoA == 0:
        Mdis = (c(y)) **2 * q * cmy1(y)
    if AoA == 10:
        Mdis = (c(y)) **2 * q * cmy2(y)
    return Mdis


def Ndis0(y):
    '''

    :param y:
    :return:
    '''
    Ndis = Ldis(y, cly1) * np.cos(np.radians(0)) + Ddis(y, cdy1) * np.sin(np.radians(0))
    return Ndis


# -------------------Distribution calculations------------------
range = np.linspace(0, 15.55, 500)
L1 = Ldis(range, cly1)
D1 = Ddis(range, cdy1)
M1 = Mdis(range, cmy1)

L2 = Ldis(range, cly2)
D2 = Ddis(range, cdy2)
M2 = Mdis(range, cmy2)

# ---------------------Plots------------------
'''
plt.plot(range, L1)
plt.plot(range, L2)
plt.xlabel("Wingspan")
plt.ylabel("CL")
plt.show()

plt.plot(range, D1)
plt.plot(range, D2)
plt.xlabel("Wingspan")
plt.ylabel("Cd")
plt.show()

plt.plot(range, M1)
plt.plot(range, M2)
plt.xlabel("Wingspan")
plt.ylabel("CM")
plt.show()
'''

# --------- Effect of velocity and angle of attack----------
# From MainWing_a=0.00_v=10.00ms.txt, MainWing_a=10.00_v=10.00ms.txt
CL0 = 0.176815
CL10 = 0.970715
CD0 = 0.001319
CD10 = 0.038922
CM0 = -0.267948
CM10 = -1.16705

CLd = 1
CDd = 0.035
CMd = 0.00000005


def CLdy(y):  # y here would be the wingspan, "range"
    '''

    :param y:
    :return:
    '''
    CLdyvalue = cly1(y) + (CLd - CL0) / (CL10 - CL0) * (cly2(y) - cly1(y))
    return CLdyvalue


def CDdy(y):  # y here would be the wingspan, "range"
    '''

    :param y:
    :return:
    '''
    CDdyvalue = cdy1(y) + (CDd - CD0) / (CD10 - CD0) * (cdy2(y) - cdy1(y))
    return CDdyvalue


def CMdy(y):  # y here would be the wingspan, "range"
    '''

    :param y:
    :return:
    '''
    CMdyvalue = cmy1(y) + (CMd - CM0) / (CM10 - CM0) * (cmy2(y) - cmy1(y))
    return CMdyvalue


def alphad():  # Design angle of attack
    '''

    :return:
    '''
    alphadvalue = np.degrees(np.arcsin((CLd - CL0) / (CL10 - CL0) * np.sin(np.radians(10))))
    return alphadvalue


print("Design angle of attack = ", alphad())

# -----------------Plotting lift distribution----------------------
'''
plt.plot(range, CLdy(range))
plt.title("Wingspan vs. CLd(y)")
plt.show()

plt.plot(range, CDdy(range))
plt.title("Wingspan vs. CDd(y)")
plt.show()

plt.plot(range, CMdy(range))
plt.title("Wingspan vs. CMd(y)")
plt.show()
'''

# ---------------------Shear moment and bending diagrams-----------------------
def Ny(AoA):
    '''

    :param AoA:
    :return:
    '''
    if AoA == 0:
        L = L1
        D = D1
    else:
        L = L2
        D = D2
    return np.cos(np.radians(AoA)) * L + np.sin(np.radians(AoA)) * D


def Ty(AoA):
    """ # Example:
    This function returns the sweep of the planform at any given percentage of the cord
    :param c_percent: The percent of the cord (range from 0 to 1)
    :type c_percent: float
    :return: The sweep of the planform at a percentage of the cord [rad]
    :rtype: float
    """
    if AoA == 0:
        L = L1
        D = D1
    else:
        L = L2
        D = D2
    return np.sin(np.radians(AoA)) * L + np.cos(np.radians(AoA)) * D

'''
plt.plot(range, Ny(0))
plt.title("Normal force distribution at 0 deg")
# plt.show()

plt.plot(range, Ty(0))
plt.title("Tangential force distribution at 0 deg")
# plt.show()

plt.plot(range, Ny(10))
plt.title("Normal force distribution at 10 deg")
# plt.show()

plt.plot(range, Ty(10))
plt.title("Tangential force distribution at 10 deg")
# plt.show()
'''

'''
# ----------- Shear force diagram ------------
def Ay(y):
    t = 0.0005
    return c(y) * 2.0580 * t


def weight(y):
    

    :param y:
    :return:
    
    [Vtot, Verr] = sp.integrate.quad(Ay, 0, 15.55)
    Vfuel = 29174 / (800 * 2)  # Fuel density is 800 kg/m^3 as per WP 1
    mwing = 1980 / 2

    # rho = a * rho_fuel + b * rho_wing, a + b = 1
    a = Vfuel / Vtot
    b = 1 - a
    rhofuel = 800  # kg/m^3
    rhowing = 2700  # kg/m^3
    rhotot = a * rhofuel + b * rhowing

    Wy = rhotot * Ay(range) * 9.81
    return Wy  # Return a single value from a function in the shape of "y=5*x+10"


# -------------------- Total weight --------------------
def wy0(y):
    return Ny(0) - weight(y)


def wy10(y):
    return Ny(10) - weight(y)


wy0inter = sp.interpolate.interp1d(list(range), list(wy0(range)), kind='cubic', fill_value="extrapolate")
wy10inter = sp.interpolate.interp1d(list(range), list(wy10(range)), kind='cubic', fill_value="extrapolate")


def wtotfunc(y):
    return wy10inter(y)



plt.plot(range, wy0inter(range))
plt.title('Interpolated Wy0')
plt.show()

plt.plot(range, Ay(range))
plt.title('Airfoil Area')
plt.show()

plt.plot(range, wy10inter(range))
plt.title('Interpolated Wy10')
plt.show()



# -------------------------------CHANGE THIS CODE A LITTLE BIT, COPIED FROM INTERNET ---------------------------------
def V0(y):
    res = np.zeros_like(y)
    for i, val in enumerate(y):
        [V0, err] = sp.integrate.quad(wtotfunc, 0, val)
        res[i] = V0
    return res
# -------------------------------CHANGE THIS CODE A LITTLE BIT, COPIED FROM INTERNET ---------------------------------
'''

'''
plt.plot(range, V0(range))
plt.title('Shear force')
plt.show()

plt.plot(range, weight(range))
plt.title("Weight diagrams")
plt.show()
'''

'''
Assumptions in this model:
A(y) = P*t

Problems:
Shear force reaction force at root not included --> starts at 0 goes to max
Solve by subtracting by total shear? gives negative shear, but it should be positive.
'''
