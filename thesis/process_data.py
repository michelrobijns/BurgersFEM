import numpy as np
import matplotlib.pyplot as plt


def compute_st1(time, x_dns, u_dns, x_proj, u_proj):
    num_timesteps = time.size
    num_elements = x_proj.size - 1  # Assuming piecewise-linear elements

    element_average_stabilization = np.zeros((num_timesteps, num_elements))

    for i in range(num_timesteps):
        for j in range(num_elements):
            # Edge coordinates and element spacing
            x_l = x_proj[j]
            x_r = x_proj[j+1]
            h = x_r - x_l

            # Sampling locations within element
            x = np.linspace(x_l, x_r, 1001)

            u = np.interp(x, x_dns, u_dns[i, :])

            ubar = np.interp(x, x_proj, u_proj[i, :])

            uprime = u - ubar

            element_average_stabilization[i, j] = \
                np.mean(ubar * uprime + 0.5 * np.power(uprime, 2))# - 0.002 * (uprime[-1] - uprime[0]) / h

    return element_average_stabilization


def main():
    # Load the DNS data and its projection
    print()
    print("Loading raw data...")

    input_directory = "./raw_data"
    output_directory = "./training_data"

    t = np.load('{}/t.npy'.format(input_directory))
    x_dns = np.load('{}/x.npy'.format(input_directory))
    u_dns = np.load('{}/u.npy'.format(input_directory))
    x_proj = np.load('{}/x_proj.npy'.format(input_directory))
    u_proj = np.load('{}/u_proj.npy'.format(input_directory))
    f_proj = np.load('{}/f_proj.npy'.format(input_directory))

    # Print a summary

    print()
    print("    Found {} timesteps".format(t.shape[0]))
    print("    Found {} elements in the fine mesh".format(x_dns.shape[0] - 1))
    print("    Found {} elements in the coarse mesh".format(x_proj.shape[0] - 1))

    # Build the training datasets

    print()
    print("Building dataset...")
    data, targets = build_dataset(t, x_dns, u_dns, x_proj, u_proj, f_proj)

    # Randomly shuffle the dataset
    indices = np.arange(data.shape[0])
    np.random.shuffle(indices)

    data = data[indices]
    targets = targets[indices]

    # Print a summary

    print()
    print("    Created a dataset of {} examples".format(data.shape[0]))
    print("    Every example has {} features and 1 target".format(data.shape[1]))
    print()
    print("        Feature 1  = ubar_l")
    print("        Feature 2  = ubar_r")
    print("        Feature 3  = dubar_dt_l")
    print("        Feature 4  = dubar_dt_r")
    print("        Feature 5  = f_l")
    print("        Feature 6  = f_r")
    print()
    print("        Target  1  = element-average of ubar * uprime + 0.5 * uprime * uprime")

    # Save datasets to files

    print()
    print("Saving dataset in '{}'...".format(output_directory))

    np.save("{}/data.npy".format(output_directory), data)
    np.save("{}/targets.npy".format(output_directory), targets)

    print()
    print("    Done")
    print()


if __name__ == '__main__':
    np.set_printoptions(
        precision=4,
        linewidth=200,
        suppress=True)

    main()
