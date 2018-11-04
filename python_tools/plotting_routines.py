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


def plot_2_with_slider(time, x0, u0, x1, u1, label0="u0", label1="u1", pad_factor=0.05):
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.2)
    fig.suptitle("t = {0:.3f}".format(time[0]))

    if x0.size < 50:
        fmt0 = "-o"
    else:
        fmt0 = "--"

    if x1.size < 50:
        fmt1 = "-o"
    else:
        fmt1 = "-"

    # Plot the first entry of `u'
    line0, = plt.plot(x0, u0[0, :], fmt0, markerfacecolor='none', label=label0)
    line1, = plt.plot(x1, u1[0, :], fmt1, markerfacecolor='none', label=label1)

    plt.xlabel("x")
    plt.legend(loc="best")

    # Setup limits
    x_min = np.min([np.min(x0), np.min(x1)])
    x_max = np.max([np.max(x0), np.max(x1)])
    y_min = np.min([np.min(u0), np.min(u1)])
    y_max = np.max([np.max(u0), np.max(u1)])

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
        line0.set_ydata(u0[timestep, :])
        line1.set_ydata(u1[timestep, :])

        # Update the timestamp in the title of the plot
        fig.suptitle("t = {0:.3f}".format(time[timestep]))

        # Redraw canvas while idle
        fig.canvas.draw_idle()

    # Call update function on slider value change
    slider.on_changed(update_plot_on_slider_changed)

    plt.show()


def plot_3_with_slider(time, x0, u0, x1, u1, x2, u2, label0="u0", label1="u1", label2="u2", pad_factor=0.05):
    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.2)
    fig.suptitle("t = {0:.3f}".format(time[0]))

    if x0.size < 50:
        fmt0 = "-o"
    else:
        fmt0 = "--"

    if x1.size < 50:
        fmt1 = "-o"
    else:
        fmt1 = "-"

    if x2.size < 50:
        fmt2 = "-o"
    else:
        fmt2 = "-"

    # Plot the first entry of `u'
    line0, = plt.plot(x0, u0[0, :], fmt0, markerfacecolor='none', label=label0)
    line1, = plt.plot(x1, u1[0, :], fmt1, markerfacecolor='none', label=label1)
    line2, = plt.plot(x2, u2[0, :], fmt2, markerfacecolor='none', label=label2)

    plt.xlabel("x")
    plt.legend(loc="best")
    plt.grid(True)

    # Setup limits
    x_min = np.min([np.min(x0), np.min(x1), np.min(x2)])
    x_max = np.max([np.max(x0), np.max(x1), np.max(x2)])
    y_min = np.min([np.min(u0), np.min(u1), np.min(u2)])
    y_max = np.max([np.max(u0), np.max(u1), np.max(u2)])

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
        line0.set_ydata(u0[timestep, :])
        line1.set_ydata(u1[timestep, :])
        line2.set_ydata(u2[timestep, :])

        # Update the timestamp in the title of the plot
        fig.suptitle("t = {0:.3f}".format(time[timestep]))

        # Redraw canvas while idle
        fig.canvas.draw_idle()

    # Call update function on slider value change
    slider.on_changed(update_plot_on_slider_changed)

    plt.show()
