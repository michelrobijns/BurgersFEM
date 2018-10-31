import numpy as np
from plotting_routines import plot_with_slider


def main():
    t = np.loadtxt('t.dat')
    x = np.loadtxt('x.dat')
    u = np.loadtxt('u.dat')

    plot_with_slider(t, x, u)


if __name__ == '__main__':
    main()
