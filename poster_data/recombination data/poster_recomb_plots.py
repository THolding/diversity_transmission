import numpy as np
import matplotlib.pyplot as plt

lowPath = "recomb_low/"
highPath = "recomb_high/"

###Import all data

#hostMOI = np.loadtxt("default_host_moi.csv", delimiter="\n")
#meanImmunity = np.loadtxt("default_host_immunity_mean.csv", delimiter="\n")
#immunityVariance = np.loadtxt("default_host_immunity_variance.csv", delimiter="\n")

#hostDiversity = np.loadtxt("default_host_diversity.csv", delimiter="\n")

#mosPrev = np.loadtxt("default_mosquito_prevalences.csv", delimiter="\n")

time = np.loadtxt(lowPath+"default_timesteps.csv", delimiter="\n")

antigenHostDiversityLow = np.loadtxt(lowPath+"default_antigenic_host_diversity.csv", delimiter="\n")
hostPrevLow = np.loadtxt(lowPath+"default_host_prevalences.csv", delimiter="\n")
eirLow = np.loadtxt(lowPath+"default_eir.csv", delimiter="\n")

antigenHostDiversityHigh = np.loadtxt(highPath+"default_antigenic_host_diversity.csv", delimiter="\n")
hostPrevHigh = np.loadtxt(highPath+"default_host_prevalences.csv", delimiter="\n")
eirHigh = np.loadtxt(highPath+"default_eir.csv", delimiter="\n")


meanImmunityLow = np.loadtxt(lowPath+"default_host_immunity_mean.csv", delimiter="\n")
meanImmunityHigh = np.loadtxt(highPath+"default_host_immunity_mean.csv", delimiter="\n")

immuneProportionLow = meanImmunityLow / antigenHostDiversityLow;
immuneProportionHigh = meanImmunityHigh / antigenHostDiversityHigh;

###new plots
plt.figure()
plt.title("host immunity")
plt.xlabel("time (years)")
plt.ylabel("immunity")
plt.plot(time/365.0, meanImmunityHigh, 'r', linewidth=2, label="high intragenic recombination")
plt.plot(time/365.0, meanImmunityLow, 'b', linewidth=2, label="low intragenic recombination")
plt.legend(loc=0)
#plt.ylim(0,1)
plt.savefig("recomb_immunity.pdf")


plt.figure()
plt.title("diversity")
plt.xlabel("time (days)")
plt.ylabel("proportion of max diversity used")
plt.plot(time, antigenHostDiversityHigh, 'r', linewidth=2, label="high intragenic recombination")
plt.plot(time, antigenHostDiversityLow, 'b', linewidth=2, label="low intragenic recombination")
plt.legend(loc=0)
plt.ylim(0,1)
plt.savefig("recomb_diversity.pdf")

plt.figure()
plt.title("prevalence")
plt.xlabel("time (days)")
plt.ylabel("prevalence")
plt.plot(time, hostPrevHigh, 'r', linewidth=2, label="high intragenic recombination")
plt.plot(time, hostPrevLow, 'b', linewidth=2, label="low intragenic recombination")
plt.legend(loc=0)
plt.savefig("recomb_prevalence.pdf")

plt.figure()
plt.title("EIR")
plt.xlabel("time (days)")
plt.ylabel("EIR")
plt.plot(time, eirHigh*365, 'r', linewidth=2, label="high intragenic recombination")
plt.plot(time, eirLow*365, 'b', linewidth=2, label="low intragenic recombination")
plt.legend(loc=0)
#plt.ylim(0,1)
plt.savefig("recomb_eir.pdf")


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