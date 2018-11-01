import sys
sys.path.append(sys.path[0] + "/..")

import numpy as np
from python_tools.plotting_routines import plot_2_with_slider


def main():
    t = np.loadtxt(sys.path[0] + "/t.dat")
    x = np.loadtxt(sys.path[0] + "/x.dat")
    u = np.loadtxt(sys.path[0] + "/u.dat")
    x_proj = np.loadtxt(sys.path[0] + "/x_proj.dat")
    u_proj = np.loadtxt(sys.path[0] + "/u_proj.dat")

    plot_2_with_slider(t, x, u, x_proj, u_proj, "u", "P(u)")


if __name__ == "__main__":
    main()
