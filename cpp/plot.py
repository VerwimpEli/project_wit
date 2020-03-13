import matplotlib.pyplot as plt
import numpy as np

VW = 10
VH = 10
file_name = "test.out"

with open(file_name, "r") as f:
    v = [float(x) for x in f.readline().rstrip().split(",")]
    t = [float(x) for x in f.readline().rstrip().split(",")]

v = np.array(v).reshape((VW, VH))
t = np.array(t).reshape((VW, VH))

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.imshow(v, extent=[0, 1, 0, 1])
ax2.imshow(t, extent=[0, 1, 0, 1])

plt.show()
