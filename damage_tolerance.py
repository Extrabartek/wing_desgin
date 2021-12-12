import numpy as np
import matplotlib.pyplot as plt
import main as fn
import testing_numbers as tn
'''
# ------------------- Parameters ---------------------------------
# Wing (box) properties (copied from main.py for ease of access)
aluminum = fn.Material(2700, 276 * (10 ** 6), 310 * (10 ** 6), 68.9 * (10 ** 9), 26 * (10 ** 9))
planform = fn.Planform(31.11, 6.46, 1.84, 0.63, 0.6, 0.15)
stringer_full = fn.Stringer(0.007, 0.15, 0.15, planform.b / 2, aluminum)
list_stringers = []
for x in range(4):
    list_stringers.append(stringer_full)
wingbox = fn.WingBox(0.005, list_stringers, aluminum)
'''
# Material Properties AL6061-T6
sigma_ult = 310 * 10 ** 6  # Ultimate tensile strength [Pa]
sigma_yield = 275 * 10 ** 6  # Yield strength [Pa]
tau_ult = 205 * 10 ** 6  # Ultimate shear strength [Pa]

E = 68.9 * 10 ** 9  # Young's modulus [Pa]
K_IC = 29 * 10 ** 6  # Fracture toughness [Pa-m½] [Ashby # We take the safest case, lowest K_IC]

# ------------------- Tensile strength failure -------------------
...

# ------------------- Crack effect -------------------------------
c = 0.005  # Crack size [m]
a = c/2  # Half-crack size [m]

# Following makes use of K_IC = K_I = sigma_a * np.sqrt(np.pi * a) where sigma_a is tension stress applied


sigma_a = np.sqrt(K_IC ** 2 / (np.pi * a))  # Maximal tensile loading [Pa]
print("Crit. tension stress =", sigma_a / (10 ** 6), "[MPa]")
c_max = 2 * (((K_IC / sigma_ult) ** 2) / np.pi)
print("Maximum size of crack such that wing skin fails at limit load = ", c_max, "[m]")

sigmapos = tn.normal_list[0]
sigmaneg = tn.normal_list[1]
x = 0
rib_placement = []
i = 0
while x <= 15.5:
    rib_placement.append(x + fn.rib_spacing(abs(fn.normal_stress(x, tn.wingbox, tn.planform, tn.wingbox.height(tn.planform, x)/2))))
    x = rib_placement[i]
    i = i + 1


print('Web buckling is', fn.web_buckling(13.5, tn.wingbox, tn.planform)/10**6)
print('Shear stress is', fn.shear_stress(13.5,tn.wingbox.height(tn.planform, 13.5), tn.wingbox, tn.planform)/10**6)

print('Rib placement =', rib_placement)


'''
del sigmapos[-1]  # MOS goes to infinity (very high in graph), del last point for clearer graph
del sigmaneg[-1]  # Same as the line above
b = np.arange(0, 15.55, 15.55/(len(sigmapos)-1))

MOSpos = np.abs(sigma_a / sigmapos)  # Calculation of safety factor for positive loading
MOSneg = np.abs(sigma_a / sigmaneg)  # Calculation of safety factor for negative loading
print('len b = ', len(b), 'len sigma = ', len(sigmapos))

print("b =", b)
print("sigmapos =", sigmapos)

MOS_testpos = np.abs(sigma_a/sigmapos[0])  # Calculate lowest safety factor (positive loading)
MOS_testneg = np.abs(sigma_a/sigmaneg[1])  # Calculate lowest safety factor (negative loading)
print('The margin of safety at the root for positive loading is', MOS_testpos)
print('The margin of safety at the root for negative loading is', MOS_testneg)

# ------------------ Plots ----------------------
plt.plot(b, MOSpos, label="Positive loading")
plt.plot(b, MOSneg, label='Negative loading')
plt.yscale("log")
plt.xlabel("Wing span [m]")
plt.ylabel("Safety Margin [-]")
plt.legend()
plt.show()
'''