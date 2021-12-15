import logging
import math
import time

import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
import WP_41 as WP41
import main as fn

# Standard variables
aluminum = fn.Material(2700, 276 * (10 ** 6), 310 * (10 ** 6), 68.9 * (10 ** 9), 26 * (10 ** 9), 0.33)
planform = fn.Planform(31.11, 6.46, 1.84, 0.63, 0.6, 0.15)
stringer_full = fn.Stringer(0.003, 0.08, planform.b / 2, aluminum)
list_stringers = []
rib_list =[]
for x in range(4):
    list_stringers.append(stringer_full)
wingbox = fn.WingBox([0.005, 0.005, 0.005], list_stringers, list_stringers, rib_list, aluminum)
WP41.b = planform.b
WP41.cr = planform.cr
WP41.ct = planform.ct
WP41.q = fn.dynamic_pressure(wingbox, planform)

# List statements
step_size = 0.1
rangy_range = np.arange(0, planform.b / 2, step_size)
number_of_strigers_list = []
stringer_len = []
list_of_box_thickness = np.arange(0.003, 0.010, 0.001)
list_of_stringer_thickness = np.arange(0.001, 0.006, 0.001)
list_of_base_len = np.arange(0.15, 0.31, 0.05)
list_of_flange_len = np.arange(0.15, 0.31, 0.05)
list_of_combinations = []

# Optimisation
i = 0
for t_b in list_of_box_thickness:
    print(f"I'm doing wingbox thickness {1000 * t_b: 0.3f} mm. We at i = {i}")
    for t_s in list_of_stringer_thickness:
        for s_b in list_of_base_len:
            for s_f in list_of_flange_len:
                print(f"We are at iteration = {i}")
                number_of_strigers_list = []
                number_of_strigers_list1 = []
                number_of_strigers_list2 = []
                i += 1
                test_stringer = fn.Stringer(t_s, s_b, s_f, planform.b / 2, aluminum)
                test_list_of_stringer = []
                stringer_len = []
                for u in range(4):
                    test_list_of_stringer.append(test_stringer)
                test_wing = fn.WingBox(t_b, test_list_of_stringer, aluminum)

                number = [0, True]
                fn.load_factor = 2.5
                fn.AoA = 10
                fn.fuel = 1
                WP41.q = fn.dynamic_pressure(wingbox, planform)
                for x in rangy_range:
                    number = fn.optimize_stringers(x, test_wing, planform)
                    if not number[1]:
                        break
                    else:
                        number_of_strigers_list1.append(number[0])
                if not number[1]:
                    break

                number = [0, True]
                fn.load_factor = -1
                fn.AoA = -10
                fn.fuel = 1
                WP41.q = fn.dynamic_pressure(wingbox, planform)
                for x in rangy_range:
                    number = fn.optimize_stringers(x, test_wing, planform)
                    if not number[1]:
                        break
                    else:
                        number_of_strigers_list2.append(number[0])
                if not number[1]:
                    break
                number_of_strigers_list = list(map(max, number_of_strigers_list1,
                                                   number_of_strigers_list2))
                stringer_len = fn.stringer_length_conversion(number_of_strigers_list,
                                                             test_stringer, step_size,
                                                             rangy_range)
                wingbox.stringers = stringer_len
                wingbox.thickness = t_b
                len_list = []
                for stringer in wingbox.stringers:
                    len_list.append(stringer.x_stop)
                fn.fuel = 0
                try:
                    if wingbox.total_weight(planform) < list_of_combinations[0][1]:
                        print(
                            f"Wingbox number {i} is completed with weight {wingbox.total_weight(planform):.3f} kg, wingbox thickness {1000 * t_b:.3f} mm, stringer "
                            f"thickness {1000 * t_s:.3f} mm, base length {100 * s_b:.3f} cm, flange height {100 * s_f:.3f} cm, number of stri"
                            f"ngers {len(wingbox.stringers):.3f}, distribution {len_list}")
                except LookupError:
                    print(f"First wingbox")
                    print(
                        f"Wingbox number {i} is completed with weight {wingbox.total_weight(planform):.3f} kg, wingbox thickness {1000 * t_b:.3f} mm, stringer "
                        f"thickness {1000 * t_s:.3f} mm, base length {100 * s_b:.3f} cm, flange height {100 * s_f:.3f} cm, number of stri"
                        f"ngers {len(wingbox.stringers):.3f}, distribution {len_list}")
                list_of_combinations.append([wingbox, wingbox.total_weight(planform)])
                list_of_combinations = sorted(list_of_combinations, key=lambda u: u[1])

# t_w = 0.004, t_s = 0.003, b = 0.05, f = 0.05 stringers 34 (max 60) done
# t_w = 0.004, t_s = 0.003, b = 0.15, f = 0.05 max 20 done
# t_W = 0.005, t_s = 0.001, b = 0.05, f = 0.05 max 8
# t_w = 0.004, t_s = 0.003, b = 0.15, f = 0.15 max 20