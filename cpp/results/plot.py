import matplotlib.pyplot as plt
import numpy as np
import argparse

def main():

	parser = argparse.ArgumentParser()

	parser.add_argument('-f', type=str, default='test.out')
	parser.add_argument('-u', type=int, default=25)
	parser.add_argument('-w', type=int, default=-1)

	args = parser.parse_args()

	args.w = args.u if args.w == -1 else args.w

	VW = args.w
	VH = args.u
	file_name = args.f
	
	with open(file_name, "r") as f:
		v = [float(x) for x in f.readline().rstrip().split(",")]
		t = [float(x) for x in f.readline().rstrip().split(",")]

	t_min, t_max = np.min(t), np.max(t)

	v = np.array(v).reshape((VW, VH))
	t = np.array(t).reshape((VW, VH))

	fig, (ax1, ax2) = plt.subplots(1, 2)
	vs = ax1.imshow(v, extent=[0, 1, 0, 1], vmin=0, vmax=1)
	ts = ax2.imshow(t, extent=[0, 1, 0, 1], vmin=t_min, vmax=t_max)

	fig.colorbar(vs, ax=ax1)
	fig.colorbar(ts, ax=ax2)

	plt.show()


if __name__ == '__main__':
	main()
