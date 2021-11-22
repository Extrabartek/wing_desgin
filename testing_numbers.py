import math
import time

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import WP_41 as WP41
import main as fn

aluminum = fn.Material(2700, 276 * 10 ** 6, 310 * 10 ** 6, 68.9 * 10 ** 9, 26 * 10 ** 9)
planform = fn.Planform(31.11, 6.46, 1.84, 0.63, 0.6, 0.15)
stringer = fn.Stringer(0.005, 0.08, 0.08, planform.b / 2, aluminum)
weight_final = 26299
list_stringers = []
for x in range(4):
    list_stringers.append(stringer)
wingbox = fn.WingBox(0.00075, list_stringers, aluminum)
WP41.b = planform.b
WP41.cr = planform.cr
WP41.ct = planform.ct
WP41.q = fn.load_factor * fn.g * (weight_final/2 + integrate.quad(lambda a: wingbox.mass_distribution(planform, a), 0, planform.b/2)[0])/(integrate.quad(lambda a: WP41.c(a) * ((fn.AoA == 0) * WP41.cly1(a) + (fn.AoA == 10) * WP41.cly2(a)) + (fn.AoA == -10) * WP41.cly3(a), 0, planform.b/2)[0])
#print(WP41.q)

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

#print(fn.twist_angle(planform.b / 2, wingbox, aluminum, planform) * (180 / math.pi))
twist_list = []
g_list = []
torsion_list = []
accuracy_thing = 5
rangy_range = range(0, accuracy_thing*int(planform.b/2))
lift_list = []
bending_list = []
shear_list = []
deflection_list = []
normal_stress_list = []
shear_stress_list = []
number_of_strigers_list = []
'''
start_test = time.time()
for x in rangy_range:
    test_shear_list.append(fn.test_shear_force(x / accuracy_thing, wingbox, planform))
end_test = time.time()
'''
'''
for x in rangy_range:
    #bending_list.append(fn.bending_moment(x / accuracy_thing, wingbox, planform))
    start = time.time()
    number = fn.vertical_deflection(x / accuracy_thing, aluminum, wingbox, planform)
    deflection_list.append(number)
    end = time.time()
    print("Value: ", number, "Time: ", end-start)
'''

#number = fn.vertical_deflection(planform.b/2, aluminum, wingbox, planform)

#print(number, "Time to execute:", end-start, "s")


for x in rangy_range:
    while abs(fn.normal_stress(x / accuracy_thing, wingbox, planform)) > aluminum.sig_yld:
        for a in range(4):
            list_stringers.append(stringer)
        wingbox = fn.WingBox(0.001, list_stringers, aluminum)
    number_of_strigers_list.append(len(list_stringers))
    list_stringers = []
    for a in range(4):
        list_stringers.append(stringer)
    wingbox = fn.WingBox(0.001, list_stringers, aluminum)
    #start = time.time()
    #number = fn.twist_angle(x / accuracy_thing, wingbox, aluminum, planform)
    #number_1 = wingbox.torsional_constant(x / accuracy_thing, planform)
    #number_2 = fn.torsion(x / accuracy_thing, planform)
    #number_3 = fn.shear_force(x / accuracy_thing, wingbox, planform)
    #number_4 = fn.normal_stress(x /accuracy_thing, wingbox, planform)
    #number_5 = fn.shear_stress(x /accuracy_thing, wingbox, planform)
    #normal_stress_list.append(number_4)
    #shear_stress_list.append(number_5)
    #twist_list.append(np.degrees(number))
    #g_list.append(number_1)
    #torsion_list.append(number_2)
    #shear_list.append(number_3)
    #number_4 = fn.bending_moment(x, wingbox, planform)
    #bending_list.append(number_4)
    #end = time.time()
    #print("Value .. ", number_4, "Time... ", end-start)
    #lift_list.append(WP41.Ndis0(x / accuracy_thing))

#plt.plot(rangy_range, torsion_list)
#plt.plot(rangy_range, normal_stress_list)
#plt.show()
#plt.plot(rangy_range, shear_stress_list)
#plt.show()

plt.plot(rangy_range, number_of_strigers_list)
plt.show()

#plt.plot(rangy_range, shear_list)
#plt.show()
#plt.plot(rangy_range, twist_list)
#plt.show()

#print("test speed", end_test-start_test)

#print(shear_list[0])
#print(shear_list[-1])

#print("Actual tip:", fn.shear_force(planform.b/2, wingbox, planform))


#plt.plot(rangy_range, bending_list)
#plt.show()

#print(integrate.quad(lambda x: WP41.Ndis0(x), 0, planform.b / 2))




print(len(list_stringers))