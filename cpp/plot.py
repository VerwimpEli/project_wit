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
	#    t = [float(x) for x in f.readline().rstrip().split(",")]

	v = np.array(v).reshape((VW, VH))
	#t = np.array(t).reshape((VW, VH))

	fig, ax1 = plt.subplots(1, 1)
	vs = ax1.imshow(v, extent=[0, 1, 0, 1], vmin=0, vmax=1)
	fig.colorbar(vs, ax=ax1)

	plt.show()


if __name__ == '__main__':
	main()
