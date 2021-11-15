import math
import time

import matplotlib.pyplot as plt
import scipy.integrate as integrate
import WP_41 as WP41
import main as fn

aluminum = fn.Material(2700, 276 * 10 ** 6, 310 * 10 ** 6, 68.9 * 10 ** 9, 26 * 10 ** 9)
planform = fn.Planform(31.11, 6.46, 1.84, 0.63, 0.6, 0.15)
stringer = fn.Stringer(0.005, 0.08, 0.08, planform.b / 2, aluminum)
list_stringers = []
for x in range(4):
    list_stringers.append(stringer)
wingbox = fn.WingBox(0.005, list_stringers, aluminum)
# print(WP41.Ndis0(0))

# bending_list = []
# for x in range(0, int(100*planform.b / 2)):
#   start = time.time()
#  number = fn.shear_force(x/100, wingbox, planform)
# end = time.time()
# print("The force and ... ", number, x/100, (end-start))
# bending_list.append(number)
# plt.plot(range(0, int(100*planform.b / 2)), bending_list)
# plt.show()

print(fn.twist_angle(planform.b / 2, wingbox, aluminum, planform) * (180 / math.pi))
twist_list = []
g_list = []
torsion_list = []
accuracy_thing = 1
rangy_range = range(0, accuracy_thing * int(planform.b / 2))
lift_list = []
bending_list = []

print(fn.lateral_deflection(planform.b/2, aluminum, wingbox, planform))

for x in rangy_range:
    start = time.time()
    number = fn.twist_angle(x / accuracy_thing, wingbox, aluminum, planform)
    number_1 = wingbox.torsional_constant(x / accuracy_thing, planform)
    number_2 = fn.torsion(x, planform)
    number_3 = WP41.Ndis0(x)
    twist_list.append(number)
    g_list.append(number_1)
    torsion_list.append(number_2)
    lift_list.append(number_3)
    #number_4 = fn.bending_moment(x, wingbox, planform)
    #bending_list.append(number_4)
    end = time.time()
    #print("Value .. ", number_4, "Time... ", end-start)
    # lift_list.append(WP41.Ny(10, x))

plt.plot(rangy_range, torsion_list)
plt.plot(rangy_range, lift_list)
plt.show()

plt.plot(rangy_range, g_list)
plt.show()

plt.plot(rangy_range, twist_list)
plt.show()


#plt.plot(rangy_range, bending_list)
#plt.show()

#print(integrate.quad(lambda x: WP41.Ndis0(x), 0, planform.b / 2))
