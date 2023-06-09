import argparse
import numpy as np
import os
import matplotlib.pyplot as plt
from sklearn.datasets import make_blobs

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--n-samples",
        default=100,
        type=int,
        help="the number of data points to generate"
    )
    parser.add_argument(
        "--dim",
        default=2,
        type=int,
        help="the dimension of the points"
    )
    parser.add_argument(
        "--spread",
        default=1.0,
        type=float,
        help="spread of each cluster. (higher number will make the clusters more diffuse)"
    )
    parser.add_argument(
        "--output-dir",
        default="tsp_examples",
        type=str,
        help="the path where the points should be written"
    )
    parser.add_argument(
        "--plot",
        action="store_true",
        default=False,
        help="whether or not to plot the points (only works if dim is 2)"
    )
    parser.add_argument(
        "--plot-path",
        default="plot.png",
        type=str,
        help="the path to save the plot to. (only used if the --plot argument is provided)"
    )

    args = parser.parse_args()

    points, _ = make_blobs(
        n_samples=args.n_samples,
        n_features=args.dim,
        cluster_std=args.spread,
    )

    full_output_dir = os.path.join(args.output_dir, f"spread={args.spread}")
    if not os.path.exists(full_output_dir):
        os.mkdir(full_output_dir)

    # write the output file
    output_path = os.path.join(full_output_dir, f"{args.n_samples}.txt")
    with open(output_path, "w") as out:
        out.write(f"{args.n_samples} {args.dim}\n")
        for point in points:
            out.write(" ".join([str(x) for x in np.round(point, 5)]) + "\n")

    if args.dim == 2 and args.plot:
        # plot the points and save
        plt.scatter(points[:, 0], points[:, 1])
        plt.title("TSP points")
        plt.xlabel("x")
        plt.ylabel("y")
        plt.savefig(args.plot_path)
