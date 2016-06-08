import os
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------
# Open and read in data file, define parameters
# ----------------------------------------------------------------------------------

data = np.loadtxt('outputs/simMagsR.txt')
N = len(data) # number of data points
# generate time axis
time = np.linspace (0, 1, N) # in s

# ----------------------------------------------------------------------------------
# Plot results
# ----------------------------------------------------------------------------------

fig = plt.figure()
fig.subplots_adjust(hspace = .5)

plotMags = fig.add_subplot(1, 1, 1)
plt.plot(time, data, 'k')
# hide every 2nd tick label on y axis to make it look nicer:
for label in plotMags.yaxis.get_ticklabels()[::2]:
	label.set_visible(False)
plt.xlabel('Time [s]')
plt.ylabel('Amplitude [arb]')

plt.show()
#fig.savefig("outputs/mags.png")
#fig.savefig("outputs/mags.pdf")