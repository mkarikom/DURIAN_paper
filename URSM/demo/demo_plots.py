##
## Check the results for demo data
##

from __future__ import with_statement
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import logging, sys

def plot_scatter(ax, prng, nb_samples=100):
    """Scatter plot.
    """
    for mu, sigma, marker in [(-.5, 0.75, 'o'), (0.75, 1., 's')]:
        x, y = prng.normal(loc=mu, scale=sigma, size=(2, nb_samples))
        ax.plot(x, y, ls='none', marker=marker)
    ax.set_xlabel('X-label')
    return ax


## simulated data
data_dir = "demo_data/"
true_A = np.loadtxt(data_dir + "demo_profile_matrix.csv", dtype=float, delimiter=",")
true_W = np.loadtxt(data_dir + "demo_bulk_mix_proportions.csv", dtype=float, delimiter=",")

## resulting data
res_dir = "demo_out/"
est_A = np.loadtxt(res_dir + "gemout_est_A.csv", dtype=float, delimiter=",")
est_W = np.loadtxt(res_dir + "gemout_exp_W.csv", dtype=float, delimiter=",")

## scatter plot for eatimated A
plt.scatter(true_A.flatten(), est_A.flatten(), alpha=0.5)
plt.plot([-0.05, np.max(est_A)+0.05], [-0.05, np.max(est_A)+0.05], 'k-')
plt.xlabel('true A')
plt.ylabel('estimated A')
plt.title('profile matrix')
plt.savefig('estimation_A.png')
plt.close()

## scatter plot for eatimated W
plt.scatter(true_W.flatten(), est_W.transpose().flatten(), alpha=0.5)
plt.plot([-0.05, np.max(est_W)+0.05], [-0.05, np.max(est_W)+0.05], 'k-')
plt.xlabel('true W')
plt.ylabel('estimated W')
plt.title('mixing proportions')
plt.savefig('estimation_W.png')
plt.close()

