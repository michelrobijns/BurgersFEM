import numpy as np
from plotting_routines import plot_2_with_slider


def main():
    t = np.loadtxt('t.dat')
    x = np.loadtxt('x.dat')
    u = np.loadtxt('u.dat')
    u_e = np.loadtxt('u_e.dat')

    plot_2_with_slider(t, x, u, u_e)


if __name__ == '__main__':
    main()
