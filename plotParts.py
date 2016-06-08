import os
import numpy as np
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------------
# Open and read in data file, define parameters
# ----------------------------------------------------------------------------------

dataReal = np.loadtxt('outputs/simRealR.txt')
dataImag = np.loadtxt('outputs/simImagR.txt')
N = len(dataReal) # number of data points
# generate time axis
time = np.linspace (0, 1, N) # in s

# ----------------------------------------------------------------------------------
# Plot results
# ----------------------------------------------------------------------------------

fig = plt.figure()
fig.subplots_adjust(hspace = .5)

plotMags = fig.add_subplot(1, 1, 1)
plt.plot(time, dataReal, 'k', label='Real part')
plt.plot(time, dataImag, 'k--', label='Imaginary part')
# hide every 2nd tick label on y axis to make it look nicer:
for label in plotMags.yaxis.get_ticklabels()[::2]:
	label.set_visible(False)
# get legend
lg = plt.legend(loc='upper right', prop={'size':12})
lg.draw_frame(False)
plt.xlabel('Time [s]')
plt.ylabel('Amplitude [arb]')

plt.show()
#fig.savefig("outputs/parts.png")
#fig.savefig("outputs/parts.pdf")