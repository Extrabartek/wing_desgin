import main as fn
import numpy as np
import WP_41 as WP41
import matplotlib.pyplot as plt
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
step_size = 0.5
rangy_range_2 = np.arange(0, planform.b/2, step_size)


aluminum = fn.Material(2700, 276 * (10 ** 6), 310 * (10 ** 6), 68.9 * (10 ** 9), 26 * (10 ** 9), 0.33)
# list_of_stringer_lengths = [0.7, 1.2, 1.6, 1.9, 2.2, 2.5, 2.7, 3.0, 3.2, 3.5, 3.7, 3.9, 4.1, 4.3, 4.4, 4.7, 4.9, 5.2, 15.5, 15.5]
list_of_stringer_top = []
list_of_stringer_bottom = []
list_of_stringer_lengths1 = [15.0, 15.0, 15.0, 14.0, 14.0, 12.0, 12.0, 12.0, 12.0, 12.0, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5, 8.5]
list_of_stringer_lengths2 = [14.5, 14.5, 14.5, 13.0, 13.0, 10.5, 10.5, 10.5, 10.5, 10.5, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0]

for x in range(len(list_of_stringer_lengths1)):
    list_of_stringer_top.append(fn.Stringer(t, a, list_of_stringer_lengths1[x], aluminum))
for x in range(len(list_of_stringer_lengths2)):
    list_of_stringer_bottom.append(fn.Stringer(t, a, list_of_stringer_lengths2[x], aluminum))

wingbox = fn.WingBox([t_top, t_bottom, t_spar], list_of_stringer_top, list_of_stringer_bottom, [0.5], aluminum)


tau_max_list = [[], []]
buckling_stringer_count = [[], []]
actual_stringer_count = [[], []]

fn.load_factor = 2.5
fn.AoA = 10
WP41.q = fn.dynamic_pressure(wingbox, planform)

for x in rangy_range_2:
    tau_max_list[0].append(fn.tau_max(x, wingbox, planform)[0])
    buckling_stringer_count[0].append(fn.skin_buckling_stringer_count(x, wingbox, planform)[0])
    number = 0
    for stringer in wingbox.stringers_top:
        number += (stringer.x_stop > x)
    actual_stringer_count[0].append(number)

fn.load_factor = -1
fn.AoA = -10
WP41.q = fn.dynamic_pressure(wingbox, planform)

for x in rangy_range_2:
    tau_max_list[1].append(fn.tau_max(x, wingbox, planform)[0])
    buckling_stringer_count[1].append(fn.skin_buckling_stringer_count(x, wingbox, planform)[0])
    number = 0
    for stringer in wingbox.stringers_bottom:
        number += (stringer.x_stop > x)
    actual_stringer_count[1].append(number)



yid_list = []
proc80_yid_list = []
for x in rangy_range_2:
    yid_list.append(aluminum.sig_yld / 2)
    proc80_yid_list.append((aluminum.sig_yld / 2) * 0.8)

plt.plot(rangy_range_2, tau_max_list[0], label="Load factor: 2.5")
plt.plot(rangy_range_2, tau_max_list[1], label="Load factor: -1.0")
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
plt.plot(rangy_range_2, actual_stringer_count[0], label="Load factor real: 2.5")
plt.plot(rangy_range_2, actual_stringer_count[1], label="Load factor real: -1.0")
plt.plot(rangy_range_2, buckling_stringer_count[0], label="Buckling Load factor needed: 2.5")
plt.plot(rangy_range_2, buckling_stringer_count[1], label="Buckling Load factor needed: -1.0")
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

print(f"{list_of_length_req1} top side")
print(f"{list_of_length_req2} bottom side")
'''