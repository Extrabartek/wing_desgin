import main as fn
import matplotlib.pyplot as plt


aluminum = fn.Material(2700, 276 * 10 ** 6, 310 * 10 ** 6, 68.9 * 10 ** 9, 26 * 10 ** 9)
planform = fn.Planform(31.11, 6.46, 1.84, 0.63, 0.6, 0.15)
stringer = fn.Stringer(0.005, 0.08, 0.08, planform.b / 2, aluminum)
list_stringers = []
for x in range(4):
    list_stringers.append(stringer)
wingbox = fn.WingBox(0.005, list_stringers, aluminum)

bending_list = []
for x in range(0, int(1000*planform.b/2)):
    bending_list.append(fn.bending_moment(x/1000, wingbox, planform))

plt.plot(range(0, int(1000*planform.b/2)), bending_list)
plt.show()