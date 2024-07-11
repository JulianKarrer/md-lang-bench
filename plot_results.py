# type: ignore
from itertools import cycle
import matplotlib.pyplot as plt

NUMBER_OF_TIMESTEPS = 1000


def data_from_csv(filename: str, name_x: str, name_y: str) -> tuple[list[float], list[float]]:
    """From a given `filename`, parse a csv file containing data 
    for different runs of a measurement by inspecting the column names in the first row of the file:

    - `name_x` is the name of the column used for x-values
    - `name_y` is the name of the column used for y-values
    - returns a tuple of lists of x and y values
    """
    with open(filename) as run_rs:
        lines: list[str] = run_rs.readlines()
        fl: str = lines[0].replace("\n", "")
        x_index = [i for (i, x) in enumerate(
            fl.split(",")) if x == name_x][0]
        y_index = [i for (i, y) in enumerate(
            fl.split(",")) if y == name_y][0]
        xs: list[float] = []
        ys: list[float] = []
        for line in lines[1:]:
            split = line.split(",")
            xs += [float(split[x_index])]
            ys += [float(split[y_index])]
        return xs, ys


linestyles = cycle(["-", "--", "-.", ":"])
markerstyles = cycle(["o", "v", "s", "D", "X", "^", "*"])

if __name__ == "__main__":
    print("Plotting the results of both benchmarks")

    # load data from csvs
    xs_cpp, ys_cpp = data_from_csv(
        "./cpp/runtimes.csv", "nb_atoms", "runtime_micros")
    xs_rs, ys_rs = data_from_csv(
        "./rust/runtimes.csv", "nb_atoms", "runtime_micros")

    # prepare plot and set its parameters
    plt.rcParams['text.usetex'] = True
    plt.rcParams.update({'font.size': 15})
    fig, ax = plt.subplots()
    fig.tight_layout()
    fig.set_size_inches(16., 9.)
    plt.yscale("log")
    plt.xscale("log")

    # describe the plot, label axes
    ax.set_title(
        r"Runtime per Simulation Step as a Function of the Number of Atoms")
    ax.set_xlabel(r"Number of Atoms $N$")
    ax.set_ylabel(r"Runtime $t$ per Simulation Step $(\mu s)$")
    fig.text(0.99, 0.01, r"Average across 100 subsequent timesteps each, Lennard-Jones direct summation, Regular grid of atoms",
             horizontalalignment='right', fontsize="10")

    # actually plot the curves
    ax.plot(xs_cpp, ys_cpp, next(linestyles),
            label="C++", marker=next(markerstyles))
    ax.plot(xs_rs, ys_rs, next(linestyles),
            label="Rust", marker=next(markerstyles))

    # save to file
    plt.legend(fontsize="12")
    plt.savefig("runtimes.png", dpi=400)
    plt.yscale("linear")
    plt.xscale("linear")
    plt.savefig("runtimes_lin.png", dpi=400)
