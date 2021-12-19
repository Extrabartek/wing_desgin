import main as fn
import numpy as np
import WP_41 as WP41
import matplotlib.pyplot as plt
import math as math
# this is pure agony
'''
def does_it_work(wingbox, planform):
    i = 0
    step_size = 0.25
    rangy_range = np.arange(0, planform.b / 2, step_size)
    sig_strengh_max = 99999999999999
    does_it_buckle = True # change it to false if the whole thing is good

    while sig_strengh_max > wingbox.material.sig_yld and does_it_buckle and i < 10:
        for y in rangy_range:
            n_top = len(wingbox.stringers_top)
            n_bottom = len(wingbox.stringers_bottom)
            j = 0
            while fn.tau_max(y, wingbox, planform) > wingbox.material.sig_yld:
                if j < n_top:
                    wingbox.stringers_top[j].x_stop = y
                    j += 1
                else:
                    wingbox.stringers_top.append(fn.Stringer(wingbox.stringers_top[0].t, wingbox.stringers_top[0].a, y, wingbox.stringers_top[0].material))
                if j < n_top:
                    wingbox.stringers_top[j].x_stop = y
                    j += 1
                else:
                    wingbox.stringers_top.append(fn.Stringer(wingbox.stringers_top[0].t, wingbox.stringers_top[0].a, y, wingbox.stringers_top[0].material))


    return "Yess"
'''
planform = fn.Planform(31.11, 6.46, 1.84, 0.63, 0.6, 0.15)

t = 0.001  # meter
a = 0.10  # meter
t_top = 0.005
t_bottom = 0.005
t_spar = 0.005
step_size = 0.01
rangy_range_2 = np.arange(0, planform.b / 2, step_size)


aluminum = fn.Material(2700, 276 * (10 ** 6), 310 * (10 ** 6), 68.9 * (10 ** 9), 26 * (10 ** 9), 0.33)
# list_of_stringer_lengths = [0.7, 1.2, 1.6, 1.9, 2.2, 2.5, 2.7, 3.0, 3.2, 3.5, 3.7, 3.9, 4.1, 4.3, 4.4, 4.7, 4.9, 5.2, 15.5, 15.5]
list_of_stringer_top = []
list_of_stringer_bottom = []
list_of_stringer_lengths1 = [15.56, 15.56,14.73, 14.73, 14.09, 13.31, 13.31, 11.19, 11.19, 11.19, 11.19, 11.19, 11.19, 6.6899999999999995, 6.6899999999999995, 6.6899999999999995, 6.6899999999999995, 6.6899999999999995, 6.6899999999999995, 6.6899999999999995, 6.6899999999999995, 6.6899999999999995, 6.6899999999999995, 6.6899999999999995]
list_of_stringer_lengths2 = [15.56, 15.56,14.3, 14.3, 13.75, 12.25, 12.25, 9.13, 9.13, 9.13, 9.13, 9.13, 9.13, 2.9299999999999997, 2.9299999999999997, 2.9299999999999997, 2.9299999999999997, 2.9299999999999997, 2.9299999999999997, 2.9299999999999997, 2.9299999999999997, 2.9299999999999997, 2.9299999999999997, 2.9299999999999997]

for x in range(len(list_of_stringer_lengths1)):
    list_of_stringer_top.append(fn.Stringer(t, a, list_of_stringer_lengths1[x], aluminum))
for x in range(len(list_of_stringer_lengths2)):
    list_of_stringer_bottom.append(fn.Stringer(t, a, list_of_stringer_lengths2[x], aluminum))

wingbox = fn.WingBox([t_top, t_bottom, t_spar], list_of_stringer_top, list_of_stringer_bottom, [0.5], aluminum)


tau_max_list = [[], []]
buckling_stringer_count = [[], []]
actual_stringer_count = [[], []]
spacing_list = [[], []]

fn.load_factor = 2.5
fn.AoA = 10
WP41.q = fn.dynamic_pressure(wingbox, planform)

for x in rangy_range_2:
    # tau_max_list[0].append(fn.tau_max(x, wingbox, planform)[0])
    buckling_stringer_count[0].append(fn.skin_buckling_stringer_count(x, wingbox, planform)[0])
    number = 0
    for stringer in wingbox.stringers_top:
        number += (stringer.x_stop > x)
    actual_stringer_count[0].append(number)
    stringer_count = 0
    for stringer in wingbox.stringers_top:
        stringer_count += (stringer.x_stop >= x) * 1
    spacing_list[0].append((wingbox.width(planform, x) - stringer_count * wingbox.stringers_top[0].a / 2)/(stringer_count-1))

fn.load_factor = -1
fn.AoA = -10
WP41.q = fn.dynamic_pressure(wingbox, planform)

for x in rangy_range_2:
    # tau_max_list[1].append(fn.tau_max(x, wingbox, planform)[0])
    buckling_stringer_count[1].append(fn.skin_buckling_stringer_count(x, wingbox, planform)[0])
    number = 0
    for stringer in wingbox.stringers_bottom:
         number += (stringer.x_stop > x)
    actual_stringer_count[1].append(number)
    stringer_count = 0
    for stringer in wingbox.stringers_bottom:
        stringer_count += (stringer.x_stop >= x) * 1
    spacing_list[1].append(
        (wingbox.width(planform, x) - stringer_count * wingbox.stringers_top[0].a / 2) / (stringer_count - 1))


'''
yid_list = []
proc80_yid_list = []
for x in rangy_range_2:
    yid_list.append(aluminum.sig_yld / 2)
    proc80_yid_list.append((aluminum.sig_yld / 2) * 0.8)

plt.plot(rangy_range_2, tau_max_list[0], label="Load factor: 2.5")
# plt.plot(rangy_range_2, tau_max_list[1], label="Load factor: -1.0")
# plt.axis([0, planform.b / 2, int(min(min(tau_max_list[0]), min(tau_max_list[1])) * 1.1), int(yid_list[0] * 1.1)])
plt.plot(rangy_range_2, yid_list, label="Half yield stress")
plt.plot(rangy_range_2, proc80_yid_list, label="80% of half yield stress")
plt.title("Design Option 2: Maximum Shear Stress Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Shear stress [Pa]")
plt.tight_layout()
plt.grid()
plt.legend()
plt.show()
'''


plt.plot(rangy_range_2, actual_stringer_count[0], label="Real stringer count top")
plt.plot(rangy_range_2, actual_stringer_count[1], label="Real stringer count bottom")
plt.plot(rangy_range_2, buckling_stringer_count[0], label="Needed stringers for buckling top")
plt.plot(rangy_range_2, buckling_stringer_count[1], label="Needed stringers for buckling bottom")
# plt.axis([0, planform.b / 2, int(min(min(tau_max_list[0]), min(tau_max_list[1])) * 1.1), int(yid_list[0] * 1.1)])
# plt.plot(rangy_range_2, yid_list, label="Half yield stress")
# plt.plot(rangy_range_2, proc80_yid_list, label="80% of half yield stress")
plt.title("Design Option 2: Maximum Shear Stress Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Stringer count [#]")
plt.tight_layout()
plt.grid()
plt.legend()
plt.show()

list_of_length_req1 = []
list_of_length_req2 = []
len_list = []

for stringer in fn.stringer_length_conversion(buckling_stringer_count[0], wingbox.stringers_top[0], step_size, rangy_range_2):
    list_of_length_req1.append(stringer.x_stop)

for stringer in fn.stringer_length_conversion(buckling_stringer_count[1], wingbox.stringers_top[0], step_size, rangy_range_2):
    list_of_length_req2.append(stringer.x_stop)

list_of_length_req1[-1] = planform.b / 2
list_of_length_req2[-1] = planform.b / 2

print(f"{list_of_length_req1} top side")
print(f"{list_of_length_req2} bottom side")
'''
'''
plt.plot(rangy_range_2, spacing_list[0], label="Top spacing")
plt.plot(rangy_range_2, spacing_list[1], label="Bottom spacing")
# plt.axis([0, planform.b / 2, int(min(min(tau_max_list[0]), min(tau_max_list[1])) * 1.1), int(yid_list[0] * 1.1)])
plt.title("Design Option 2: Maximum Shear Stress Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Spacing [m]")
plt.tight_layout()
plt.grid()
plt.legend()
plt.show()

fn.load_factor = 2.5
fn.AoA = 10
WP41.q = fn.dynamic_pressure(wingbox, planform)

# print(f"Twist is {math.degrees(fn.twist_angle(planform.b / 2, wingbox, aluminum, planform))} degree")
# print(f"Deflection is {fn.vertical_deflection(planform.b / 2, aluminum, wingbox, planform)} meters")


# [15.56, 14.71, 14.71, 13.37, 13.37, 13.37, 11.29, 11.29, 11.29, 11.29, 11.29, 6.89, 6.89, 6.89, 6.89, 6.89, 6.89, 6.89, 6.89, 6.89, 6.89, 6.89] top side
# [15.56, 14.28, 14.28, 12.95, 12.34, 12.34, 9.4, 9.4, 9.4, 9.4, 9.4, 3.2399999999999998, 3.2399999999999998, 3.2399999999999998, 3.2399999999999998, 3.2399999999999998, 3.2399999999999998, 3.2399999999999998, 3.2399999999999998, 3.2399999999999998, 3.2399999999999998, 3.2399999999999998] bottom side

