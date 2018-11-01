import sys
sys.path.append(sys.path[0] + "/..")

import numpy as np
from python_tools.plotting_routines import plot_with_slider


def main():
    t = np.loadtxt(sys.path[0] + "/t.dat")
    x = np.loadtxt(sys.path[0] + "/x.dat")
    u = np.loadtxt(sys.path[0] + "/u.dat")

    plot_with_slider(t, x, u)


if __name__ == "__main__":
    main()
