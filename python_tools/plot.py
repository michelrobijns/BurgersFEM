import numpy as np
from plotting_routines import plot_with_slider


def main():
    time = np.loadtxt('time.dat')
    x = np.loadtxt('x.dat')
    u = np.loadtxt('u.dat')

    plot_with_slider(time, x, u)


if __name__ == '__main__':
    main()
