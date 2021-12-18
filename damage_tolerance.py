import numpy as np
import matplotlib.pyplot as plt
import main as fn
import testing_numbers as tn
import WP_41
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
'''
x = 0
rib_placement = []
vertstringer_placement = []
i = 0
while x <= 15.5:
    dummy1 = fn.rib_spacing_column(abs(fn.normal_stress(x, tn.wingbox, tn.planform, tn.wingbox.height(tn.planform, x)/2)))
    rib_placement.append(x + dummy1)
    x = rib_placement[i]
    i = i + 1
i = 0
x = 0
while x <= 15.5:
    vertstringer_placement.append(x + fn.vertstringer_spacing_web(tn.aluminum, tn.wingbox, tn.planform, x))
    x = vertstringer_placement[i]
    i = i + 1

'''
vertstringer_placement = []
rib_placement = []
fn.load_factor = 2.5
fn.AoA = 10
WP_41.q = fn.dynamic_pressure(tn.wingbox, tn.planform)
x = 0
y = 0
i = 0
j = 0
while x <= tn.planform.b/2:
    rib_placement.append(x + fn.rib_spacing_column(x, tn.wingbox, tn.planform))
    x = rib_placement[i]
    if x > tn.planform.b / 2:
        print("I broke the infinite loop")
        break
    while y + fn.vertstringer_spacing_web(tn.aluminum, tn.wingbox,tn.planform, y) <= x:
        vertstringer_placement.append(y + fn.vertstringer_spacing_web(tn.aluminum, tn.wingbox, tn.planform, y))
        y = vertstringer_placement[j]
        j = j + 1
    y = rib_placement[i]
    i = i + 1



locvertstringers = []
locribs = []
for i in range(len(vertstringer_placement)-1):
    if vertstringer_placement[i] <= tn.planform.b/2:
        locvertstringers.append(vertstringer_placement[i])

for i in range(len(rib_placement)-1):
    if rib_placement[i] <= tn.planform.b/2:
        locribs.append(rib_placement[i])


stiffening_elements = locvertstringers
'''
for element in locvertstringers:
    stiffening_elements.append(element)
'''

for i in range(len(locribs)):
    stiffening_elements.append(locribs[i])
stiffening_elements.sort()
print(stiffening_elements)

x = 0.00
critstressweb = []
spacing = []
while x <= tn.planform.b/2:
    number = fn.web_buckling(x, tn.wingbox, tn.planform, tn.aluminum, stiffening_elements) / fn.shear_stress(x, tn.wingbox.height(tn.planform, x), tn.wingbox, tn.planform)
    if number > 4:
        critstressweb.append(4)
    else:
        critstressweb.append(number)
    x += 0.01
#del critstressweb[-1]
critstressweb = np.array(critstressweb)
spacing = np.array(spacing)

critstresscolumn = []
while x <= 11:
    number = fn.column_buckling(x, tn.wingbox.stringers_top[0]) / fn.normal_stress(x, tn.wingbox, tn.planform, tn.wingbox.height(tn.planform, x)*2)
    if number > 4:
        critstresscolumn.append(4)
    else:
        critstresscolumn.append(number)
    x += 0.01

del critstresscolumn[-1]
critstresscolumn = np.array(critstresscolumn)
b = np.arange(0, 11, 0.01)

plt.plot(b, critstressweb, label = "MOS web buckling")
plt.plot(b, critstresscolumn, label = "MOS Column Buckling")
plt.axis([0, tn.planform.b / 2, 1, 2])
plt.legend()
plt.title('Margin of safety for web buckling')
plt.xlabel("Distance from root [m]")
plt.ylabel("Margin of safety [-]")
plt.grid()
plt.show()

print('Web buckling critical stress is', fn.web_buckling(0, tn.wingbox, tn.planform, tn.aluminum, locvertstringers)/10**6)
print('Shear stress is', fn.shear_stress(0, tn.wingbox.height(tn.planform, 0), tn.wingbox, tn.planform)/10**6)
print('Rib placement =', locribs)
print('Number of ribs is', (len(locribs)))
print('Vertical stringer placement is', locvertstringers)
print('Number of vertical stringers is', (len(locvertstringers) - len(locribs)))
print('Vertical stringer placement is', fn.vertstringer_spacing_web(tn.aluminum, tn.wingbox, tn.planform, 0))


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