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
for x in range(20):
    list_stringers.append(stringer_full)
wingbox = fn.WingBox([0.005, 0.005, 0.005], list_stringers, list_stringers, rib_list, aluminum)
WP41.b = planform.b
WP41.cr = planform.cr
WP41.ct = planform.ct
WP41.q = fn.dynamic_pressure(wingbox, planform)

# List statements
twist_list = [[], []]
torsion_list = [[], []]
lift_list = [[], []]
bending_list = [[], []]
shear_list = [[], []]
deflection_list = [[], []]
tau_max_list = [[], []]
mass_list = []
normal_list = [[], []]
critical_column_stress = [[], []]
MMOI_list = []
torsional_list = []
stringer_count_skin_buck = [[], []]
stringer_count_skin_buck1 = [[], []]


# analysis
step_size = 0.1
rangy_range = np.arange(0, planform.b / 2 * 0.8, step_size)

'''
t_w = 0.004
t_w_list = [0.004, 0.004, 0.004, 0.004]
t_s = 0.003
b = 0.15
f = 0.15

# Make wingbox
number_of_strigers_list = []
number_of_strigers_list1 = []
number_of_strigers_list2 = []
test_stringer = fn.Stringer(t_s, b, f, planform.b / 2, aluminum)
test_list_of_stringer = []
stringer_len = []
for u in range(4):
    test_list_of_stringer.append(test_stringer)
test_wing = fn.WingBox(t_w_list, test_list_of_stringer, aluminum)

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

number_of_strigers_list = list(map(max, number_of_strigers_list1, number_of_strigers_list2))
stringer_len = fn.stringer_length_conversion(number_of_strigers_list, test_stringer, step_size, rangy_range)
wingbox.stringers = stringer_len
wingbox.thickness = t_w
len_list = []
for stringer in wingbox.stringers:
    len_list.append(stringer.x_stop)
fn.fuel = 0
print(
    f"Wingbox is completed with weight {wingbox.total_weight(planform):.3f} kg, wingbox thickness {1000 * t_w:.3f} mm, stringer "
    f"thickness {1000 * t_s:.3f} mm, base length {100 * b:.3f} cm, flange height {100 * f:.3f} cm, number of stri"
    f"ngers {len(wingbox.stringers):.3f}, distribution {len_list}")

# Analysis plots

plt.plot(rangy_range, number_of_strigers_list)
plt.axis([0, planform.b / 2, 0, int(max(number_of_strigers_list) * 1.1)])
plt.grid()
plt.title("Design Option 2: Stringer Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Number of stringers")
plt.tight_layout()
plt.show()


step_size = 0.5
rangy_range = np.arange(0, planform.b / 2, step_size)
'''
fn.load_factor = 2.5
fn.AoA = 10
fn.fuel = 1
WP41.q = fn.dynamic_pressure(wingbox, planform)
for x in rangy_range:
    # lift_list[0].append(WP41.Ndis0(x, fn.AoA))
    # torsion_list[0].append(fn.torsion(x, planform))
    # bending_list[0].append(fn.bending_moment(x, wingbox, planform))
    # shear_list[0].append(fn.shear_force(x, wingbox, planform))
    # tau_max_list[0].append(fn.tau_max(x, wingbox, planform))
    # mass_list.append(wingbox.mass_distribution(planform, x))
    # MMOI_list.append(wingbox.moment_of_inertia(planform, x))
    # torsional_list.append(wingbox.torsional_constant(x, planform))
    # twist_list[0].append(np.degrees(fn.twist_angle(x, wingbox, aluminum, planform)))
    #deflection_list[0].append(fn.vertical_deflection(x, aluminum, wingbox, planform))
    # normal_list[0].append(abs(fn.normal_stress(x, wingbox, planform, wingbox.height(planform, x)/2)))
    #critical_column_stress[0].append(fn.column_buckling(x, test_stringer))
    stringer_count_skin_buck[0].append(fn.skin_buckling_stringer_count(x, wingbox, planform)[0])
    # stringer_count_skin_buck1[0].append(fn.skin_buckling_stringer_count(x, wingbox, planform)[1])


fn.load_factor = -1
fn.AoA = -10
fn.fuel = 1
WP41.q = fn.dynamic_pressure(wingbox, planform)
for x in rangy_range:
    # lift_list[1].append(WP41.Ndis0(x, fn.AoA))
    # torsion_list[1].append(fn.torsion(x, planform))
    # bending_list[1].append(fn.bending_moment(x, wingbox, planform))
    # shear_list[1].append(fn.shear_force(x, wingbox, planform))
    # tau_max_list[1].append(fn.tau_max(x, wingbox, planform))
    # twist_list[1].append(np.degrees(fn.twist_angle(x, wingbox, aluminum, planform)))
    #deflection_list[1].append(fn.vertical_deflection(x, aluminum, wingbox, planform))
    # normal_list[1].append(fn.normal_stress(x, wingbox, planform, wingbox.height(planform, x)/2))
    #critical_column_stress[1].append(fn.column_buckling(x, test_stringer))
    stringer_count_skin_buck[1].append(fn.skin_buckling_stringer_count(x, wingbox, planform)[0])
    # stringer_count_skin_buck1[1].append(fn.skin_buckling_stringer_count(x, wingbox, planform)[1])
'''
# print(f"The maximum vertical deflection is {max(abs(deflection_list[0][-1]), abs(deflection_list[1][-1]))} m, allowed is 4.7 m")
# print(f"The twist angle at the tip is {max(abs(twist_list[0][-1]), abs(twist_list[1][-1]))} degrees")

#creating lines for the tau_max graph

plt.plot(rangy_range, normal_list[0], label="Load factor: 2.5")
plt.plot(rangy_range, normal_list[1], label="Load factor: -1.0")
plt.plot(rangy_range, critical_column_stress[0], label="Load factor: 2.5 Buckle")
plt.plot(rangy_range, critical_column_stress[1], label="Load factor: -1 Buckle")
plt.axis([0, planform.b / 2, int(min(min(normal_list[0]), min(normal_list[1])) * 1.1), int(max(max(normal_list[0]), max(normal_list[1])) * 1.1)])
plt.title("Design Option 2: Normal Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Lift per unit length[N/m]")
plt.tight_layout()
plt.grid()
plt.legend()
plt.show()


yid_list = []
proc80_yid_list = []
for x in rangy_range:
    yid_list.append(aluminum.sig_yld / 2)
    proc80_yid_list.append((aluminum.sig_yld / 2) * 0.8)

plt.plot(rangy_range, lift_list[0], label="Load factor: 2.5")
plt.plot(rangy_range, lift_list[1], label="Load factor: -1.0")
plt.axis([0, planform.b / 2, int(min(min(lift_list[0]), min(lift_list[1])) * 1.1), int(max(max(lift_list[0]), max(lift_list[1])) * 1.1)])
plt.title("Design Option 2: Lift Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Lift per unit length[N/m]")
plt.tight_layout()
plt.grid()
plt.legend()
plt.show()

plt.plot(rangy_range, torsion_list[0], label="Load factor: 2.5")
plt.plot(rangy_range, torsion_list[1], label="Load factor: -1.0")
plt.axis([0, planform.b / 2, int(min(min(torsion_list[0]), min(torsion_list[1])) * 1.1), int(max(max(torsion_list[0]), max(torsion_list[1])) * 1.1)])
plt.title("Design Option 2: Internal Torque Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Torque per unit length [Nm/m]")
plt.tight_layout()
plt.grid()
plt.legend()
plt.show()

plt.plot(rangy_range, bending_list[0], label="Load factor: 2.5")
plt.plot(rangy_range, bending_list[1], label="Load factor: -1.0")
plt.axis([0, planform.b / 2, int(min(min(bending_list[0]), min(bending_list[1])) * 1.1), int(max(max(bending_list[0]), max(bending_list[1])) * 1.1)])
plt.title("Design Option 2: Internal Bending Moment Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Internal bending moment [Nm]")
plt.tight_layout()
plt.grid()
plt.legend()
plt.show()

plt.plot(rangy_range, shear_list[0], label="Load factor: 2.5")
plt.plot(rangy_range, shear_list[1], label="Load factor: -1.0")
plt.axis([0, planform.b / 2, int(min(min(shear_list[0]), min(shear_list[1])) * 1.1), int(max(max(shear_list[0]), max(shear_list[1])) * 1.1)])
plt.title("Design Option 2: Internal Shear Force Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Internal shear force [N]")
plt.tight_layout()
plt.grid()
plt.legend()
plt.show()

plt.plot(rangy_range, tau_max_list[0], label="Load factor: 2.5")
plt.plot(rangy_range, tau_max_list[1], label="Load factor: -1.0")
plt.axis([0, planform.b / 2, int(min(min(tau_max_list[0]), min(tau_max_list[1])) * 1.1), int(yid_list[0] * 1.1)])
plt.plot(rangy_range, yid_list, label="Half yield stress")
plt.plot(rangy_range, proc80_yid_list, label="80% of half yield stress")
plt.title("Design Option 2: Maximum Shear Stress Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Shear stress [Pa]")
plt.tight_layout()
plt.grid()
plt.legend()
plt.show()

plt.plot(rangy_range, mass_list)
plt.axis([0, planform.b / 2, 0, int(max(mass_list) * 1.1)])
plt.title("Design Option 2: Mass Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Mass per unit area [kg/m]")
plt.tight_layout()
plt.grid()
plt.show()

plt.plot(rangy_range, MMOI_list)
plt.axis([0, planform.b / 2, 0, float(max(MMOI_list)) * 1.1])
plt.title("Design Option 2: Moment of Inertia Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Moment of inertia [$m^4$]")
plt.tight_layout()
plt.grid()
plt.show()

plt.plot(rangy_range, torsional_list)
plt.axis([0, planform.b / 2, 0, float(max(torsional_list)) * 1.1])
plt.title("Design Option 2: Torsional Constant Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Torsional Constant [$m^4$]")
plt.tight_layout()
plt.grid()
plt.show()

# plt.plot(rangy_range, deflection_list[0], label="Load factor: 2.5")
# plt.plot(rangy_range, deflection_list[1])
# plt.axis([0, planform.b / 2, min(min(deflection_list[0]), min(deflection_list[1])) * 1.1, max(max(deflaction_list[0]), max(deflaction_list[1])) * 1.1)])
# plt.legend()
# plt.show()

plt.plot(rangy_range, twist_list[0], label="Load factor: 2.5")
plt.plot(rangy_range, twist_list[1], label="Load factor: -1")
plt.axis([0, planform.b / 2, min(min(twist_list[0]), min(twist_list[1])) * 1.1, max(max(twist_list[0]), max(twist_list[1])) * 1.1])
plt.title("Design Option 2: Wing twist Angle Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Twist angle [deg]")
plt.tight_layout()
plt.grid()
plt.legend()
plt.show()


plt.plot(rangy_range, deflection_list[0], label="Load factor: 2.5")
plt.plot(rangy_range, deflection_list[1], label="Load factor: -1")
plt.axis([0, planform.b / 2, min(min(deflection_list[0]), min(deflection_list[1])) * 1.1, max(max(deflection_list[0]), max(deflection_list[1])) * 1.1])
plt.title("Design Option 3: Vertical deflection Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Vertical deflection [m]")
plt.legend()
plt.grid()
plt.show()
'''
print(wingbox.width(planform, 8))

plt.plot(rangy_range, stringer_count_skin_buck[0], label="Load factor: 2.5")
plt.plot(rangy_range, stringer_count_skin_buck[1], label="Load factor: -1")
# plt.plot(rangy_range, stringer_count_skin_buck1[0], label="Load factor: 2.5")
# plt.plot(rangy_range, stringer_count_skin_buck1[1], label="Load factor: -1")
# plt.axis([0, planform.b / 2, min(min(stringer_count_skin_buck[0]), min(stringer_count_skin_buck[1])) * 1.1, max(max(stringer_count_skin_buck[0]), max(stringer_count_skin_buck[1])) * 1.1])
plt.title("Design Option 3: Vertical deflection Distribution")
plt.xlabel("Distance from root [m]")
plt.ylabel("Vertical deflection [m]")
plt.legend()
plt.grid()
plt.show()