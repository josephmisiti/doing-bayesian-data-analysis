from IPython.core.pylabtools import figsize
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import scipy.stats as stats

plt.close("all")

figsize(11, 9)

x = np.linspace(0,1,10000)

# n_trials = [0, 1, 2, 3, 4, 5, 8, 15, 50, 500]
# data = stats.bernoulli.rvs(0.5, size=n_trials[-1])
# x = np.linspace(0, 1, 100)
#
# # For the already prepared, I'm using Binomial's conj. prior.
# for k, N in enumerate(n_trials):
#     sx = plt.subplot(len(n_trials) / 2, 2, k + 1)
#     plt.xlabel("$p$, probability of heads") \
#         if k in [0, len(n_trials) - 1] else None
#     plt.setp(sx.get_yticklabels(), visible=False)
#     heads = data[:N].sum()
#     y = dist.pdf(x, 1 + heads, 1 + N - heads)
#     plt.plot(x, y, label="observe %d tosses,\n %d heads" % (N, heads))
#     plt.fill_between(x, 0, y, color="#348ABD", alpha=0.4)
#     plt.vlines(0.5, 0, 4, color="k", linestyles="--", lw=1)
#
#     leg = plt.legend()
#     leg.get_frame().set_alpha(0.4)
#     plt.autoscale(tight=True)
#
#
# plt.suptitle("Bayesian updating of posterior probabilities",
#              y=1.02,
#              fontsize=14)





index = 0
col = 0
for a in [.5,1,2,3,4]:
	row = 0
	for b in [.5,1,2,3,4]:
		y = stats.beta(a, b).pdf(x)
		sx = plt.subplot(5,5,index)
		row += 1
		index+=1
		plt.setp(sx.get_yticklabels(), visible=False)
		plt.setp(sx.get_xticklabels(), visible=False)
		plt.xlabel("$theta$")
		plt.ylabel("$p(theta)$")
		plt.plot(x, y, label="a=%s, b=%s" % (a,b))
		leg = plt.legend()
		leg.get_frame().set_alpha(0.1)
		plt.autoscale(tight=True)
		
	col += 1

plt.suptitle("Beta Distribution For Different a and b Values",
             y=1.02,
             fontsize=14)
plt.tight_layout()
plt.show()

	