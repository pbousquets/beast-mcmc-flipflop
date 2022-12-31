import click
from collections import OrderedDict
import numpy as np
import pandas as pd
from scipy import linalg, stats
from Bio import Phylo
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
import os
import json

import msprime
import re

## Custom click to allow grouping arguments
class CustomOption(click.Option):
    def __init__(self, *args, **kwargs):
        kwargs.update({'show_default': True})
        self.help_group = kwargs.pop('help_group', None)
        super(CustomOption, self).__init__(*args, **kwargs)


class CustomCommand(click.Command):
    def format_options(self, ctx, formatter):
        """Writes all the options into the formatter if they exist."""
        opts = OrderedDict()
        for param in self.get_params(ctx):
            rv = param.get_help_record(ctx)
            if rv is not None:
                if hasattr(param, 'help_group') and param.help_group:
                    opts.setdefault(str(param.help_group), []).append(rv)
                else:
                    opts.setdefault('Options', []).append(rv)

        for name, opts_group in opts.items():
            with formatter.section(name):
                formatter.write_dl(opts_group)

# Random tree
def generate_random_tree(k, N, t, root_height, seed):
    # Simulate the coalescent tree using msprime
    ts = msprime.simulate(sample_size=k, Ne=N, mutation_rate=0, random_seed=seed)
    tree = next(ts.trees())

    # Get the tree in Newick format
    newick = tree.newick()

    # Calculate the generation time in years
    total_length_generations = tree.time(tree.root)
    generation_time = t / total_length_generations

    # Convert the branch lengths from number of generations to years
    branch_lengths = re.findall(r":(\d+(?:\.\d+)?)", newick)
    for i, length in enumerate(branch_lengths):
        length = float(length)
        length *= generation_time
        newick = newick.replace(f":{branch_lengths[i]}", f":{length:.6f}", 1)

    # Add the root height
    newick = newick[:-1] + f":{root_height:.6f};"

    # Print the modified Newick string
    return newick

def generate_state_var(S):
    return [(k, m) for m in range(S+1) for k in range(S+1) if k+m <=S]

def generate_rate_matrix(S, lam, mu, gamma, stateVar=None):
    if stateVar is None:
        stateVar = generate_state_var(S)

    RateMatrix = np.zeros((len(stateVar), len(stateVar)))

    for down, (k_down, m_down) in enumerate(stateVar):
        for across, (k, m) in enumerate(stateVar):
            if k == k_down-1 and m == m_down :
                RateMatrix[down, across] = (S-m-k) * (k*lam/(S-1) + 2*mu)
            elif k == k_down and m == m_down-1:
                RateMatrix[down, across] = m * (S-m-k) * lam / (S-1)
            elif k == k_down+1 and m == m_down-1:
                RateMatrix[down, across] = k * (m*lam/(S-1)+mu)
            elif k == k_down+1 and m == m_down:
                RateMatrix[down, across] = k * ((S-m-k)*lam/(S-1)+gamma)
            elif k == k_down and m == m_down+1:
                RateMatrix[down, across] = m * (S-m-k) * lam / (S-1)
            elif k == k_down-1 and m == m_down+1:
                RateMatrix[down, across] = m * (k*lam/(S-1)+2*gamma)
            elif k == k_down and m == m_down:
                RateMatrix[down, across] = -(2*((k+m)*(S-m-k)+k*m)*lam/(S-1) +
                                            (k+2*m)*gamma + (2*S-(k+2*m))*mu)

    return RateMatrix

def multinomial_rvs(counts, p, rng=None):
    if rng is None:
        rng = np.random.default_rng()
    if not isinstance(counts, (np.ndarray)):
        counts = np.full(p[0, ...].shape, counts)

    out = np.zeros(np.shape(p), dtype=int)
    ps = np.cumsum(p[::-1, ...], axis=0)[::-1, ...]
    # Conditional probabilities
    with np.errstate(divide='ignore', invalid='ignore'):
        condp = p / ps
    condp[np.isnan(condp)] = 0.0

    for i in range(p.shape[0]-1):
        binsample = rng.binomial(counts, condp[i, ...])
        out[i, ...] = binsample
        counts -= binsample

    out[-1, ...] = counts

    return out

def evolve_state(RateMatrix, state, dt, rng = None):
    if rng is None:
        rng = np.random.default_rng()

    probStates = linalg.expm(RateMatrix * dt) @ state
    stateOut = multinomial_rvs(1, probStates, rng)

    return stateOut

def convert_state_beta(state, S, stateVar=None):
    if stateVar is None:
        stateVar = generate_state_var(S)

    methState = np.zeros((2*S+1, state.shape[1]))
    for index, (k, m) in enumerate(stateVar):
        methState[k+2*m, :] += state[index, :]

    betaOut = np.linspace(0, 1, 2*S+1) @ methState

    return betaOut

def beta_convert_params(mu, kappa):
    """
    Convert mean/dispersion parameterization of a beta distribution to the ones scipy supports

    """

    if np.any(kappa <= 0):
        raise Exception("kappa must be greater than 0")
    elif np.any(mu <= 0) or np.any(mu >= 1):
        raise Exception("mu must be between 0 and 1")

    alpha = kappa * mu
    beta = kappa * (1- mu)

    return alpha, beta

def rescale_beta(beta, delta, eta):
    # Linear transform of beta values from between
    # 0 and 1 to between delta and eta
    return (eta - delta) * beta + delta

def add_noise(beta, delta, eta, kappa):
    beta_rescale = rescale_beta(beta, delta, eta)
    a, b = beta_convert_params(beta_rescale, kappa)

    return stats.beta.rvs(a, b)

def all_parents(tree):
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade:
            parents[child] = clade
    return parents

def simulate_beta(tree, S, lam, mu, gamma, delta, eta, kappa, N):
    rng = np.random.default_rng()

    stateVar = generate_state_var(S)
    RateMatrix = generate_rate_matrix(S, lam, mu, gamma, stateVar)

    initialConditions = np.zeros((int(0.5 * (S+1) * (S+2)), N))
    initialConditions[[0, -1], :] = 0.5

    stateInit = multinomial_rvs(1, initialConditions, rng)
    tree.root.fcpgs = evolve_state(RateMatrix, stateInit, tree.root.branch_length, rng)

    n = 0
    for clade in tree.find_clades():
        if clade.name is None:
            clade.name = f'{n}'
            n += 1

    parent_dict = all_parents(tree)

    for clade in tree.find_clades():
        if clade.name != '0':
            parent_clade = parent_dict[clade]
            clade.fcpgs = evolve_state(RateMatrix, parent_clade.fcpgs, clade.branch_length, rng)
            clade.beta = convert_state_beta(clade.fcpgs, S, stateVar)
            clade.betaNoisy = add_noise(clade.beta, delta, eta, kappa)

    return tree

def file_exists_callback(ctx, param, value):
    if value and not os.path.isfile(value):
        raise click.BadParameter(f"The provided tree file ({value}) does not exist")
    return value

@click.command(cls=CustomCommand)
@click.option("--nSites", "N", default=1000, help="Number of fCpGs to simulate", type=int, cls=CustomOption, help_group="Model parameters")
@click.option("--nCells", "S", default=5, help="Number of stem cells to consider", type=int, cls=CustomOption, help_group="Model parameters")
@click.option("--lam", default=1.3, help="Fixed lambda parameter value", type=float, cls=CustomOption, help_group="Model parameters")
@click.option("--mu", default=0.03, help="Fixed mu parameter value", type=float, cls=CustomOption, help_group="Model parameters")
@click.option("--gamma", default=0.03, help="Fixed gamma parameter value", type=float, cls=CustomOption, help_group="Model parameters")
@click.option("--delta", default=0.05, help="Fixed delta parameter value", type=float, cls=CustomOption, help_group="Model parameters")
@click.option("--eta", default=0.93, help="Fixed eta parameter value", type=float, cls=CustomOption, help_group="Model parameters")
@click.option("--kappa", default=100, help="Fixed kappa parameter value", type=float, cls=CustomOption, help_group="Model parameters")
@click.option("--random_tree", default=False, help="Use a random tree", is_flag=True, cls=CustomOption, help_group="Phylogenetic tree")
@click.option("--tree_file", help="Tree file in Newick format to use", type=str, callback=file_exists_callback, cls=CustomOption, help_group="Phylogenetic tree")
@click.option("--samples", default=10, help="Number of samples in the simulated tree", type=int, cls=CustomOption, help_group="Phylogenetic tree")
@click.option("--eff_pop_size", default=10_000, help="Effective population size for the simulated tree", type=float, cls=CustomOption, help_group="Phylogenetic tree")
@click.option("--time", default=30, help="Number of years for the simulated tree", type=int, cls=CustomOption, help_group="Phylogenetic tree")
@click.option("--root_limits", default="5,25", help="Low and high limit the size of the root height. Comma-separated total years.", type=str, cls=CustomOption, help_group="Phylogenetic tree")
@click.option("--seed", default=123, help="Set a seed", type=int, cls=CustomOption, help_group="Other")
@click.option("--directory", default=".", help="Alternative root directory for the output", type=click.Path(exists=True), cls=CustomOption, help_group="Other")
@click.option("--output", required = True, help="Prefix that output files will carry", type=str, cls=CustomOption, help_group="Other")
@click.option("--no_tree", default=False, help="Whether to avoid plotting the phylogenetic tree", is_flag=True, cls=CustomOption, help_group="Other")
def main(N, S, lam, mu, gamma, delta, eta, kappa, random_tree, tree_file, seed, samples, eff_pop_size, time, root_limits, directory, output, no_tree):
    np.random.seed(seed)

    arguments = locals() # Get current variables as a dictionary

    directory = f"{directory}/{output}"
    os.makedirs(directory, exist_ok=True)

    if random_tree:
        tree_file = f"{directory}/{output}.tree"
        low_limit, high_limit = root_limits.split(",")
        root_height = np.random.uniform(float(low_limit), float(high_limit))
        tree_time = time - root_height

        arguments["tree_time"] = tree_time
        arguments["root_height"] = root_height

        tree = generate_random_tree(samples, eff_pop_size, tree_time, root_height, seed)
        with open(tree_file, "w") as f:
            f.write(tree)
    else:
        assert tree_file, "Error. By default it's expected a tree. Otherwise, please activate the --random_tree flag."

    tree = Phylo.read(tree_file, "newick")
    assert tree.root.branch_length, "The tree root's branch length is null"

    tree.ladderize()  # Flip branches so deeper clades are displayed at top
    tree.rooted = True

    if not no_tree:
        with PdfPages(f'{directory}/{output}Tree.pdf') as export_pdf:
            fig, ax = plt.subplots()
            Phylo.draw(tree, axes=ax, do_show=False)
            plt.tick_params(
                axis='y',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                left=False,      # ticks along the bottom edge are off
                right=False,         # ticks along the top edge are off
                labelleft=False) # labels along the bottom edge are off
            plt.ylabel('')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['left'].set_visible(False)
            export_pdf.savefig()
            plt.close()

    tree = simulate_beta(tree, S, lam, mu, gamma, delta, eta, kappa, N)

    df = pd.DataFrame({})

    for clade in tree.get_terminals():
        df[clade.name] = clade.betaNoisy

    df.to_csv(f'{directory}/{output}Beta.csv')

    with PdfPages(f'{directory}/{output}SimulatedClustermap.pdf') as export_pdf:
        sns.clustermap(df, cmap="vlag")
        export_pdf.savefig()
        plt.close()

    with PdfPages(f'{directory}/{output}BetaHist.pdf') as export_pdf:
        for col in df:
            fig, ax = plt.subplots()
            plt.hist(df[col], np.linspace(0, 1, 51), density=True, alpha=0.4)
            plt.xlabel("Fraction methylated")
            plt.ylabel("Probability density")
            sns.despine()
            plt.tight_layout()
            export_pdf.savefig()
            plt.close()

    with open(f'{directory}/{output}Parameters.json', 'w') as fp:
        json.dump(arguments, fp, indent=1)

if __name__ == "__main__":
    main()