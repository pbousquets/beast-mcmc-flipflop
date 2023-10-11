"""
This function will run the simulator for all the combinations we considered for the parameters.
"""

from flipflopphylogenysimulate import simulate_fcpgs
from os import path, urandom
import tqdm
import numpy as np
import itertools

np.random.seed(42)


def get_random_seed():
    r = urandom(4)
    return int.from_bytes(r, byteorder="big")


N = 2000  # Sites
delta = 0.04
eta = 0.92
kappa = 100
samples = [5, 9]  # leaves

eff_pop_size = 10000
time = 45
root_limits = "5,40"
seed = None
random_tree = True

cells = [4, 7, 10]
mus = [0.005, 0.01]

# We select the number of fixed mutations we want in our trees. Lambda will depend on this. We use: tm = (s-1)**2/(stemCell*lambda) and fixations = tree_height/tm
mean, std_dev, fixations = 10, 5, []  # We sample from a normal distribution with mean 10 and std_dev 5 for the number of fixed mutations we want in our trees.

while len(fixations) < 3: # Ensure that only positive values are used
    x = np.random.normal(mean, std_dev)
    if x > 0:
        fixations.append(x)

fixations = np.array(fixations)

directory = "newest_simulations/"
no_tree = False

for sample, S, mu in tqdm.tqdm(list(itertools.product(*[samples, cells, mus]))):
    lambdas = int((S - 1) ** 2) * fixations / (time * S)
    for lam in lambdas:
        lam = round(lam, 3)
        gamma = round(mu * np.random.normal(1, 0.04), 5)
        output = f"sim_{S}_{lam}_{mu}_{gamma}_{kappa}_{sample}"

        if path.exists(f"{directory}/{output}"):
            continue

        tree_file = f"{directory}/test_{sample}samples/test_{sample}samples.tree"
        simulate_fcpgs(N, S, lam, mu, gamma, delta, eta, kappa, random_tree, tree_file, get_random_seed(), sample, eff_pop_size, time, root_limits, directory, output, no_tree)
