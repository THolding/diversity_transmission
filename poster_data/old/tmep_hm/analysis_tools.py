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



#plot some timeseries
plt.figure();
plt.plot(time, mosPrev, 'y--', linewidth=2, label="mosquito prevalence")
plt.plot(time, meanImmunity, 'k', linewidth=2, label="host immunity (mean)")
plt.plot(time, immunityVariance, 'k--', linewidth=2, label="host immunity (variance)")
plt.plot(time, hostPrev, 'r', linewidth=2, label="host prevalence")
plt.title("time series data (prevalence and immunity)")
plt.xlabel("time (days)")
#plt.ylim(0,0.01)
#plt.xlim(200,1000)
plt.legend(loc=0)
plt.savefig("prevalence_immunity.pdf")

plt.figure()
plt.plot(time, hostMOI, 'b', linewidth=2, label="host MOI")
plt.title("time series (multiplicity of infection)")
plt.xlabel("time (days)")
plt.legend(loc=0)
plt.savefig("MOI.pdf")

plt.figure()
plt.plot(time, hostDiversity)
plt.title("diversity in host infections")
plt.xlabel("time (days)")
plt.ylabel("number of antigenically distinct strains")
plt.savefig("diversity.pdf")

plt.figure()
plt.plot(time, antigenHostDiversity)
plt.title("proportion of antigenic space represented in host infections")
plt.xlabel("time (days)")
plt.ylabel("Proportion of antigenic space")
plt.savefig("antigen_representation.pdf")

plt.figure()
plt.plot(time, immuneProportion)
plt.title("mean immunity / antigen representation in population")

plt.figure();
plt.plot(time, eir);
plt.title("daily eir")



#plt.figure()
#plt.pcolor(antigenHostFrequency)
#plt.colorbar()
#plt.clim(0,5)