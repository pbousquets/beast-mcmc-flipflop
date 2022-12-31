"""
This function will run the simulator for all the combinations we considered for the parameters.

IMPORTANT: this script was a one-use script. In order to run, it needs that the click decorators are removed from the main function.
Otherwise it won't work. If the function is going to be imported in multiple occasions, it should better be a standalone function, allowing
the import, and the main function would just be a wrapper for the cli calling this function.
"""

from flipflopphylogenysimulate import main
from os import path
import itertools
import tqdm

N = 2000
delta = 0.05
eta = 0.93

eff_pop_size = 5
time = 35
root_limits = "5,30"
seed = None
random_tree = False
no_tree = True
directory = "simulation_datasets/"

cells = [3, 10, 20]
Lambdas = [0.1, 1.0, 10.0]
mus = [0.001, 0.01, 0.1]
gammas = [0.001, 0.01, 0.1]
kappas = [20, 50, 100]
samples = [4, 9]


runs = list(itertools.product(*[cells, Lambdas, mus, gammas, kappas, samples]))
for S, lam, mu, gamma, kappa, sample in tqdm.tqdm(runs):
    output = f"sim_{S}_{lam}_{mu}_{gamma}_{kappa}_{sample}"

    if path.exists(f"{directory}/{output}"):
        continue

    tree_file = f"{directory}/test_{sample}samples/test_{sample}samples.tree"
    main(N, S, lam, mu, gamma, delta, eta, kappa, random_tree, tree_file, seed, samples, eff_pop_size, time, root_limits, directory, output, no_tree)


