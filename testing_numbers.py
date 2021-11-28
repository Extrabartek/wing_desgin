import logging
import math
import time

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import WP_41 as WP41
import main as fn

'''
Parameter input 
'''

aluminum = fn.Material(2700, 276 * 10 ** 6, 310 * 10 ** 6, 68.9 * 10 ** 9, 26 * 10 ** 9)
planform = fn.Planform(31.11, 6.46, 1.84, 0.63, 0.6, 0.15)
stringer_full = fn.Stringer(0.007, 0.15, 0.15, planform.b / 2, aluminum)
weight_final = 23528
list_stringers = []
for x in range(4):
    list_stringers.append(stringer_full)
wingbox = fn.WingBox(0.005, list_stringers, aluminum)
WP41.b = planform.b
WP41.cr = planform.cr
WP41.ct = planform.ct
WP41.q = fn.load_factor * fn.g * (
        weight_final / 2 + integrate.quad(lambda a: wingbox.mass_distribution(planform, a), 0, planform.b / 2)[
    0]) / (integrate.quad(
    lambda a: WP41.c(a) * ((fn.AoA == 0) * WP41.cly1(a) + (fn.AoA == 10) * WP41.cly2(a)) + (fn.AoA == -10) * WP41.cly3(
        a), 0, planform.b / 2)[0])

twist_list = []
g_list = []
torsion_list = []
step_size = 0.1
rangy_range = np.arange(0, planform.b / 2, step_size)
lift_list = []
bending_list = []
shear_list = []
deflection_list = []
normal_stress_list = []
shear_stress_list = []
mass_list = []
MMOI_list = []
torsional_list = []
number_of_strigers_list = []
stringer_len = []
list_of_box_thickness = np.arange(0.001, 0.006, 0.001)
list_of_stringer_thickness = np.arange(0.001, 0.006, 0.001)
list_of_base_len = np.arange(0.05, 0.15, 0.01)
list_of_flange_len = np.arange(0.05, 0.15, 0.01)
list_of_combinations = []
i = 0


#Optimisation
#"""
for t_b in list_of_box_thickness:
    print(f"I'm doing wingbox thickness {1000 * t_b: 0.3f} mm. We at i = {i}")
    for t_s in list_of_stringer_thickness:
        for s_b in list_of_base_len:
            for s_f in list_of_flange_len:
                number_of_strigers_list = []
                i += 1
                test_stringer = fn.Stringer(t_s, s_b, s_f, planform.b / 2, aluminum)
                test_list_of_stringer = []
                stringer_len = []
                for u in range(4):
                    test_list_of_stringer.append(test_stringer)
                test_wing = fn.WingBox(t_b, test_list_of_stringer, aluminum)

                number = [0, True]

                for x in rangy_range:
                    number = fn.optimize_stringers(x, test_wing, planform)
                    if not number[1]:
                        break
                    else:
                        number_of_strigers_list.append(number[0])
                if not number[1]:
                    break
                stringer_len = fn.stringer_length_conversion(number_of_strigers_list, test_stringer, step_size, rangy_range)
                wingbox.stringers = stringer_len
                wingbox.thickness = t_b
                len_list = []
                for stringer in wingbox.stringers:
                    len_list.append(stringer.x_stop)

                try:
                    if wingbox.total_weight(planform) < list_of_combinations[0][1]:
                        print(
                            f"Wingbox number {i} is completed with weight {wingbox.total_weight(planform):.3f} kg, wingbox thickness {1000 * t_b:.3f} mm, stringer "
                            f"thickness {1000 * t_s:.3f} mm, base length {100 * s_b:.3f} cm, flange height {100 * s_f:.3f} cm, number of stri"
                            f"ngers {len(wingbox.stringers):.3f}, distribution {len_list}")
                        plt.plot(rangy_range, number_of_strigers_list)
                        plt.show()
                        print(f"The maximum vertical deflection is {fn.vertical_deflection(planform.b / 2, aluminum, wingbox, planform)} m, allowed is 4.7 m")
                except LookupError:
                    print(f"First wingbox")
                    print(
                        f"Wingbox number {i} is completed with weight {wingbox.total_weight(planform):.3f} kg, wingbox thickness {1000 * t_b:.3f} mm, stringer "
                        f"thickness {1000 * t_s:.3f} mm, base length {100 * s_b:.3f} cm, flange height {100 * s_f:.3f} cm, number of stri"
                        f"ngers {len(wingbox.stringers):.3f}, distribution {len_list}")
                    plt.plot(rangy_range, number_of_strigers_list)
                    plt.show()
                    print(f"The maximum vertical deflection is {fn.vertical_deflection(planform.b / 2, aluminum, wingbox, planform)} m, allowed is 4.7 m")
                list_of_combinations.append([wingbox, wingbox.total_weight(planform)])
                list_of_combinations = sorted(list_of_combinations, key=lambda u: u[1])
                print(
                    f"Wingbox number {i} is completed with weight {wingbox.total_weight(planform):.3f} kg, wingbox thickness {1000 * t_b:.3f} mm, stringer "
                    f"thickness {1000 * t_s:.3f} mm, base length {100 * s_b:.3f} cm, flange height {100 * s_f:.3f} cm, number of stri"
                    f"ngers {len(wingbox.stringers):.3f}, distribution {len_list}")

"""
#analysis
t_w = 0.013
t_s = 0.003
b = 0.05
f = 0.05

#Make wingbox
print(WP41.q)
print(WP41.q * (integrate.quad(
    lambda a: WP41.c(a) * ((fn.AoA == 0) * WP41.cly1(a) + (fn.AoA == 10) * WP41.cly2(a)) + (fn.AoA == -10) * WP41.cly3(
        a), 0, planform.b / 2)[0]))
number_of_stringers_list = []
test_stringer = fn.Stringer(t_s, b, f, planform.b / 2, aluminum)
test_list_of_stringer = []
stringer_len = []
for u in range(4):
    test_list_of_stringer.append(test_stringer)
test_wing = fn.WingBox(t_w, test_list_of_stringer, aluminum)
number = [0, True]

for x in rangy_range:
    number = fn.optimize_stringers(x, test_wing, planform)
    if not number[1]:
        break
    else:
        number_of_strigers_list.append(number[0])
stringer_len = fn.stringer_length_conversion(number_of_strigers_list, test_stringer, step_size, rangy_range)
wingbox.stringers = stringer_len
len_list = []
for stringer in wingbox.stringers:
    len_list.append(stringer.x_stop)

print(
    f"Wingbox is completed with weight {wingbox.total_weight(planform):.3f} kg, wingbox thickness {1000 * t_w:.3f} mm, stringer "
    f"thickness {1000 * t_s:.3f} mm, base length {100 * b:.3f} cm, flange height {100 * f:.3f} cm, number of stri"
    f"ngers {len(wingbox.stringers):.3f}, distribution {len_list}")

# Analysis plots
plt.plot(rangy_range, number_of_strigers_list)
plt.show()

fn.load_factor = 1
fn.AoA = 0
fn.fuel = 0
print(f"The maximum vertical deflection is {fn.vertical_deflection(planform.b / 2, aluminum, wingbox, planform)} m, allowed is 4.7 m")
print(f"The twist angle at the tip is {np.degrees(fn.twist_angle(planform.b/2, wingbox, aluminum, planform))} degrees")
for x in rangy_range:
    lift_list.append(WP41.Ndis0(x, fn.AoA))
    torsion_list.append(fn.torsion(x, planform))
    bending_list.append(fn.bending_moment(x, wingbox, planform))
    shear_list.append(fn.shear_force(x, wingbox, planform))
    shear_stress_list.append(fn.shear_stress(x, wingbox, planform))
    mass_list.append(wingbox.mass_distribution(planform, x))
    MMOI_list.append(wingbox.moment_of_inertia(planform, x))
    torsional_list.append(wingbox.torsional_constant(x, planform))

plt.plot(rangy_range, lift_list)
plt.show()
plt.plot(rangy_range, torsion_list)
plt.show()
plt.plot(rangy_range, bending_list)
plt.show()
plt.plot(rangy_range, shear_list)
plt.show()
plt.plot(rangy_range, shear_stress_list)
plt.show()
plt.plot(rangy_range, mass_list)
plt.show()
plt.plot(rangy_range, MMOI_list)
plt.show()
plt.plot(rangy_range, torsional_list)
plt.show()
"""