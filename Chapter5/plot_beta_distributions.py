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

x = np.linspace(0,1,10000)

index = 0
col = 0

As = [.5,1,2,3,4]
Bs = [.5,1,2,3,4]
Bs.reverse()

for a in As:
	row = 0
	for b in Bs:
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

	