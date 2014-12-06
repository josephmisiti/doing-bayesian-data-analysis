##################
# Plots the beta distributions, page 82
####################

from IPython.core.pylabtools import figsize
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import scipy.stats as stats

plt.close("all")

figsize(11, 9)

a,b = 4,4
a1,b1 = 1,6

x = np.linspace(0,1,10000)
y = stats.beta(a, b).pdf(x)
plt.plot(x, y, label="a=%s, b=%s" % (a,b),lw=3)

y = stats.beta(a1, b1).pdf(x)
plt.plot(x, y, label="a=%s, b=%s" % (a1,b1),lw=3)

leg = plt.legend()
leg.get_frame().set_alpha(0.1)
plt.autoscale(tight=True)
plt.title("The Beta Distribution",
             y=1.02,
             fontsize=14)
plt.grid(True)``
plt.tight_layout()
plt.show()