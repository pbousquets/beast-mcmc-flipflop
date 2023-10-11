"""
This function will run the simulator for all the combinations we considered for the parameters.
"""

from flipflopphylogenysimulate import simulate_fcpgs
from os import path, urandom
import tqdm
import numpy as np


def get_random_seed():
    r = urandom(4)
    return int.from_bytes(r, byteorder="big")


N = 2000  # Sites
delta = 0.04
eta = 0.92
kappa = 100
sample = 5  # leaves

eff_pop_size = 10000
time = 45
root_limits = "5,40"
seed = None
random_tree = True

S = 3
lambdas = [0.5, 0.75, 1.0, 1.25, 1.5]
multipliers = 10
mus = [round(l * multiplier, 5) for multiplier in 1 / np.linspace(10, 500, multipliers) for l in lambdas]  # multipliers: [0.1, 0.01552, 0.00841, 0.00577, 0.00439, 0.00354, 0.00297, 0.00256, 0.00224, 0.002]
gammas = [round(l * multiplier, 5) for multiplier in 1 / np.linspace(8, 450, multipliers) for l in lambdas]  # multipliers: [0.125, 0.01751, 0.00941, 0.00644, 0.00489, 0.00394, 0.0033, 0.00284, 0.00249, 0.00222]

replicates = list(range(3))

directory = "new_simulations/"
no_tree = False

x = 0
for lam, mu, gamma in tqdm.tqdm(zip(lambdas * multipliers, mus, gammas)):
    for replicate in replicates:
        output = f"sim_{S}_{lam}_{mu}_{gamma}_{kappa}_{sample}_{replicate+1}"
        print(output)
        if path.exists(f"{directory}/{output}"):
            continue

        tree_file = f"{directory}/test_{sample}samples/test_{sample}samples.tree"
        simulate_fcpgs(N, S, lam, mu, gamma, delta, eta, kappa, random_tree, tree_file, get_random_seed(), sample, eff_pop_size, time, root_limits, directory, output, no_tree)
