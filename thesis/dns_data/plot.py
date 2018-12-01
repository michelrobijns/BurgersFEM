import sys
sys.path.append(sys.path[0] + "/../..")

import numpy as np
from python_tools.plotting_routines import plot_2_with_slider, plot_3_with_slider


def main():
    t_dns = np.loadtxt(sys.path[0] + "/t_dns.dat")
    x_dns = np.loadtxt(sys.path[0] + "/x_dns.dat")
    u_dns = np.loadtxt(sys.path[0] + "/u_dns.dat")
    x_proj = np.loadtxt(sys.path[0] + "/x_proj.dat")
    u_proj = np.loadtxt(sys.path[0] + "/u_proj.dat")

    plot_2_with_slider(t_dns, x_dns, u_dns, x_proj, u_proj, "u", "P(u)")


if __name__ == "__main__":
    main()
