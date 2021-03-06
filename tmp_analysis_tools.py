# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 07:57:09 2016

@author: rr
"""

import numpy as np
import matplotlib.pyplot as plt
import os

def read_dataset(_filePath, runName=None):
    filePath = _filePath;    
    if filePath != "" and filePath[-1] != '/':
        filePath+='/';
        
    #if no runName is given, intuit it.
    if runName == None:
        filenames = os.walk(filePath).next()[2];
        for filename in filenames:
            if "_seed.txt" in filename:
                runName = filename[0:-9];
                break;
    if runName == None:
        runName = "";
        
    
    data = {}
    ###Import all data
    data["t_days"] = np.array(np.loadtxt(filePath+runName+"_timesteps.csv", delimiter="\n"))
    data["t_years"] = np.array(data["t_days"]/365.0)
    data["mosPrev"] = np.array(np.loadtxt(filePath+runName+"_mosquito_prevalences.csv", delimiter="\n"))
    data["hostPrev"] = np.array(np.loadtxt(filePath+runName+"_host_prevalences.csv", delimiter="\n"))
    data["hostMOI"] = np.array(np.loadtxt(filePath+runName+"_host_moi.csv", delimiter="\n"))
    data["immunityVariance"] = np.array(np.loadtxt(filePath+runName+"_host_immunity_variance.csv", delimiter="\n"))
    data["meanImmunity"] = np.array(np.loadtxt(filePath+runName+"_host_immunity_mean.csv", delimiter="\n"))
    data["hostDiversity"] = np.array(np.loadtxt(filePath+runName+"_host_diversity.csv", delimiter="\n"))
    data["antigenHostDiversity"] = np.array(np.loadtxt(filePath+runName+"_antigenic_host_diversity.csv", delimiter="\n"))
    data["eir"] = np.array(np.loadtxt(filePath+runName+"_eir.csv", delimiter="\n"))
    
    ###Calculate other metrics
    data["immuneProportion"] = np.array(data["meanImmunity"] / data["antigenHostDiversity"]);
    
    ###File path and run name
    data["filePath"] = filePath
    data["runName"] = runName
    
    
    #Some output files are only there under certain circumstances!
    #Timeseries of number of mosquitoes is only outputted if mosquito population is dynamic:
    if os.path.isfile(filePath+runName+"_num_mosquitoes.csv"):
        data["num_mosquitoes"] = np.array(np.loadtxt(filePath+runName+"_num_mosquitoes.csv"));
    #Timeseries of bite rate is only outputted if bite rate is dynamic
    if os.path.isfile(filePath+runName+"_bite_rate.csv"):
        data["bite_rate"] = np.array(np.loadtxt(filePath+runName+"_bite_rate.csv"));
    #Antigen frequencies are computensively intensive to calculate so are only outputted if specifically requested:
    if os.path.isfile(filePath+runName+"_antigen_frequency.csv"):
        data["antigen_frequency"] = np.array(np.genfromtxt(filePath+runName+"_antigen_frequency.csv", delimiter=","))
    
    return data

def summary_plot(data, saveFig=False, figSavePath=None, burnInPeriod=-1):
    if burnInPeriod > 0:
        burnY = [0, 1]
        burnX = [burnInPeriod/365.0, burnInPeriod/365.0]
    
    plt.figure(figsize=(12,10))
    plt.subplot(2,2,1);
    plt.title("prevalence")
    plt.plot(data["t_years"], data["hostPrev"], 'k', linewidth=3, label="host prevalence");
    plt.plot(data["t_years"], data["mosPrev"], 'r--', linewidth=3, label="mosquito prevalence");
    plt.legend(loc=0, fontsize=9)
    plt.ylim(0,1)
    plt.ylabel("prevalence")
    plt.xlabel("time (years)")
    if burnInPeriod > 0:
        plt.plot(burnX, burnY, 'y:', linewidth = 2, label="burn in period")
    
    plt.subplot(2,2,2);
    plt.title("transmission")
    plt.plot(data["t_years"], data["eir"], 'k', linewidth=3)
    plt.ylabel("EIR (daily)")
    plt.xlabel("time (years)")
    plt.ylim(0, 0.2)#np.max(data["eir"]))
    if burnInPeriod > 0:
        plt.plot(burnX, burnY, 'y:', linewidth = 2, label="burn in period")
    
    plt.subplot(2,2,3);
    plt.title("diversity")
    plt.plot(data["t_years"], data["antigenHostDiversity"], 'k', linewidth=3, label="diversity")
    plt.ylim(0,1)
    plt.ylabel(r"diversity (proportion of $a_{max}$)")
    plt.xlabel("time (years)")
    if burnInPeriod > 0:
        plt.plot(burnX, burnY, 'y:', linewidth = 2, label="burn in period")
    
    ax1 = plt.subplot(2,2,4);
    plt.title("immunity")    
    plt.xlabel("time (years)")
    
    immuneProportion = np.array(data["immuneProportion"])
    immuneProportion[np.isnan(immuneProportion)] = 0.0
    plt.plot(data["t_years"], immuneProportion, 'k', linewidth=3)
    plt.ylabel("immunity (proportion of circulating diversity)")
    plt.ylim(0, np.max(immuneProportion));
    ax1.twinx();
    plt.plot(data["t_years"], data["meanImmunity"], 'r--', linewidth=3)
    print immuneProportion
    plt.ylabel(r"immunity (proportion of $a_{max}$)", color="red")
    plt.ylim(0,1);
    if burnInPeriod > 0:
        plt.plot(burnX, burnY, 'y:', linewidth = 2, label="burn in period")    
        #plt.plot(burnX, [0, max(1, np.max(data["immuneProportion"]))], 'y:', linewidth = 2, label="burn in period")    
    
    plt.tight_layout(pad=2.5)#, w_pad=0.5, h_pad=1.0)
    
    if figSavePath == None:
        figSavePath = data["filePath"]
    
    if figSavePath == "":
        figSavePath = "/";
    elif figSavePath[-1] != '/':
        figSavePath+='/'
    
    if saveFig:
        print "***SAVING FIGURE***"
        plt.savefig(figSavePath+data["runName"]+"_summaryplot.pdf")


#Plot_script to plot summary_plot and be used with run_tools for remote plotting.
#Plot_scripts must read data, plot data, save figure in correct place and exit plot by itself.
def summary_plot_script(args):
    if "file_path" in args and "run_name in args":
        filePath = args[args.index("file_path")+1];
        runName = args[args.index("run_name")+1];
        data = read_dataset(filePath, runName);
        
        if "burn_in_period" in args:
            summary_plot(data, burnInPeriod=int(args[args.index("burn_in_period")+1]), saveFig=True);
        else:
            summary_plot(data, saveFig=True);
        
        plt.close(); #Close fig after call to summary_plot


def phase_plot(_time, _xParam, _yParam, title="default title", xlabel="paramX", ylabel="paramY", saveFig=False, figSavePath="/", figFilename="phasePlot.pdf", burnInPeriod=None):
    time = np.array(_time);
    xParam = np.array(_xParam);
    yParam = np.array(_yParam);
    
    if figSavePath == "":
        figSavePath = "/";
    elif figSavePath[-1] != '/':
        figSavePath+='/';
    
    if figFilename[-4:] != ".pdf":
        figFilename+=".pdf";
    
    #discard burn in period
    if burnInPeriod != None:
        mask = time >= (int(burnInPeriod)/365.0);
        time = time[mask]
        xParam = xParam[mask]
        yParam = yParam[mask]
    
    segmentsT = []
    segmentsX = []
    segmentsY = []
    
    for i in range(1,len(time)):
        segmentsT.append(time[i-1:i+1])
        segmentsX.append(xParam[i-1:i+1])
        segmentsY.append(yParam[i-1:i+1])
    
    plt.figure();
    for i in range(0,len(segmentsT)):
        red = i/float(len(time));
        blue = 1.0-red;
        green = 0.25;
        plt.plot(segmentsX[i], segmentsY[i], color=(red, green, blue), linewidth=3)
    
    plt.title(title);
    plt.xlabel(xlabel);
    plt.ylabel(ylabel);
    
    if saveFig == True:
        plt.savefig(figSavePath+figFilename);
        plt.close();
        







#filePath = ""
#run_sim_reps(["run_name", "default", "num_hosts", "10", "run_time", "1000"], 100, 3)
#d = read_dataset(filePath+"rep0/", "default")
#summary_plot(d)




####Could create a namedtuple type to store data in... but easier to use dictionary!
#import collections
#DataSet = collections.namedtuple("DataSet", "run_name file_path time hostPrev mosPrev hostMOI immunityVariance meanImmunity hostDiversity antigenHostDiversity eir immuneProportion")
#
#myData = DataSet(run_name = runName, file_path = filePath,
#    time = np.loadtxt(filePath+runName+"_timesteps.csv", delimiter="\n"),
#    mosPrev = np.loadtxt(filePath+runName+"_mosquito_prevalences.csv", delimiter="\n"),
#    hostPrev = np.loadtxt(filePath+runName+"_host_prevalences.csv", delimiter="\n"),
#    hostMOI = np.loadtxt(filePath+runName+"_host_moi.csv", delimiter="\n"),
#    immunityVariance = np.loadtxt(filePath+runName+"_host_immunity_variance.csv", delimiter="\n"),
#    meanImmunity = np.loadtxt(filePath+runName+"_host_immunity_mean.csv", delimiter="\n"),
#    hostDiversity = np.loadtxt(filePath+runName+"_host_diversity.csv", delimiter="\n"),
#    antigenHostDiversity = np.loadtxt(filePath+runName+"_antigenic_host_diversity.csv", delimiter="\n"),
#    eir = np.loadtxt(filePath+runName+"_eir.csv", delimiter="\n"),
#    immuneProportion = meanImmunity / antigenHostDiversity)
#plt.plot(myData.time, myData.hostPrev)
