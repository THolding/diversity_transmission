import numpy as np
import matplotlib.pyplot as plt


###Import all data
time = np.loadtxt("default_timesteps.csv", delimiter="\n")

hostPrev = np.loadtxt("default_host_prevalences.csv", delimiter="\n")
hostMOI = np.loadtxt("default_host_moi.csv", delimiter="\n")
meanImmunity = np.loadtxt("default_host_immunity_mean.csv", delimiter="\n")
immunityVariance = np.loadtxt("default_host_immunity_variance.csv", delimiter="\n")
#relativeImmunity = np.loadtxt("default_host_immunity_relative.csv", delimiter="\n")

hostDiversity = np.loadtxt("default_host_diversity.csv", delimiter="\n")
antigenHostDiversity = np.loadtxt("default_antigenic_host_diversity.csv", delimiter="\n")
#antigenHostFrequency = np.loadtxt("default_antigen_frequencies.csv", delimiter=", ")

mosPrev = np.loadtxt("default_mosquito_prevalences.csv", delimiter="\n")

eir = np.loadtxt("default_eir.csv", delimiter="\n")

immuneProportion = meanImmunity / antigenHostDiversity;

### Plots
fig, ax1 = plt.subplots()
plt.title("Host prevalence with M:H ratio ?unknown?")
ax1.plot(time, hostPrev, 'r', linewidth=2, label="host prevalence")
ax1.plot(time, antigenHostDiversity, 'k', linewidth=2, label="antigen diversity in host infections")
plt.legend(loc=0)
plt.ylabel("prevalence or proportion of max antigens in circulation");
plt.xlabel("time (days)")
ax2 = ax1.twinx()
ax2.plot(time, eir*365.0, 'b', linewidth=2, label="EIR (projected annual)")
plt.ylabel("Projected annual EIR")
plt.legend(loc=0)
plt.savefig("ratio_m-h=0.5.pdf")

#
#fig, ax1 = plt.subplots()
#plt.title(titleStr)
#ax1.plot(x, ylast[numStrains:numStrains*2], 'ro', linewidth=3)
#ax1.set_xlabel('Cumulative exposure')
#ax1.set_ylabel('Infected proportion', color='r')
#ax2 = ax1.twinx()
#ax2.plot(x, ylast[0:numStrains], 'go', linewidth=3)
#ax2.set_ylabel('Susceptible proportion', color='g')
#plt.show()
#plt.savefig(filePath+runName+".pdf")
