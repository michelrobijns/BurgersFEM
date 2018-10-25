import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider


def plot(x, u):
    plt.figure()

    if x.size < 50:
        fmt = "-o"
    else:
        fmt = "-"

    plt.plot(x, u, fmt, markerfacecolor='none')

    plt.xlabel("x")
    plt.ylabel("u")

    plt.show()


def plot_with_slider(time, x, u, pad_factor=0.05):
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.2)
    fig.suptitle("t = {0:.3f}".format(time[0]))

    if x.size < 50:
        fmt = "-o"
    else:
        fmt = "-"

    # Plot the first entry of `u'
    line, = plt.plot(x, u[0, :], fmt, markerfacecolor='none')

    plt.xlabel("x")
    plt.ylabel("u")

    # Setup limits
    x_min = np.min(x)
    x_max = np.max(x)
    y_min = np.min(np.min(u))
    y_max = np.max(np.max(u))

    x_pad = (x_max - x_min) * pad_factor
    y_pad = (y_max - y_min) * pad_factor

    x_min -= x_pad
    x_max += x_pad
    y_min -= y_pad
    y_max += y_pad

    ax = plt.axis([x_min, x_max, y_min, y_max])

    # Setup slider
    slider_axes = plt.axes([0.25, 0.025, 0.65, 0.03],
                           facecolor="lightgoldenrodyellow")
    slider = Slider(slider_axes, 'timestep', 0, time.size-1, valinit=0,
                    valfmt='%0.0f')

    def update_plot_on_slider_changed(val):
        # Get the integer value of the slider
        timestep = int(round(slider.val))

        # Update lines
        line.set_ydata(u[timestep, :])

        # Update the timestamp in the title of the plot
        fig.suptitle("t = {0:.3f}".format(time[timestep]))

        # Redraw canvas while idle
        fig.canvas.draw_idle()

    # Call update function on slider value change
    slider.on_changed(update_plot_on_slider_changed)

    plt.show()
