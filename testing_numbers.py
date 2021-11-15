import math
import time

import WP_41
import WP_41 as WP41
import main as fn
import matplotlib.pyplot as plt


aluminum = fn.Material(2700, 276 * 10 ** 6, 310 * 10 ** 6, 68.9 * 10 ** 9, 26 * 10 ** 9)
planform = fn.Planform(31.11, 6.46, 1.84, 0.63, 0.6, 0.15)
stringer = fn.Stringer(0.005, 0.08, 0.08, planform.b / 2, aluminum)
list_stringers = []
for x in range(4):
    list_stringers.append(stringer)
wingbox = fn.WingBox(0.005, list_stringers, aluminum)
print(WP41.Ndis0(0))

#bending_list = []
#for x in range(0, int(100*planform.b / 2)):
 #   start = time.time()
  #  number = fn.shear_force(x/100, wingbox, planform)
   # end = time.time()
    # print("The force and ... ", number, x/100, (end-start))
    #bending_list.append(number)
#plt.plot(range(0, int(100*planform.b / 2)), bending_list)
#plt.show()

print(fn.twist_angle(planform.b/2, wingbox, aluminum, planform)*(180/math.pi))
torque_list = []
#lift_list = []
for x in range(0, int(10*planform.b/2)):
    start = time.time()
    number = fn.bending_moment(x/10, wingbox, planform)
    torque_list.append(fn.bending_moment(x/10, wingbox, planform))
    end = time.time()
    print("Value .. ", number, "Time... ", end-start)
    #lift_list.append(WP41.Ny(10, x))

plt.plot(range(0, int(10*planform.b / 2)), torque_list)
plt.show()

