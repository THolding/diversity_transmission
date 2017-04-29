# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 09:50:49 2017

@author: rr
"""

import numpy as np;
import matplotlib.pyplot as plt;

immunityMask = np.loadtxt("immunity_mask.csv", delimiter=",");
immuneState = np.loadtxt("test_immune_state.csv", delimiter=",");

plt.figure()
plt.plot(immunityMask)

plt.figure();
plt.plot(immuneState);
#plt.xlim(0,200);



x = np.array(range(0,200));
mu = 100.0;
var = 1.0;
amp = 1.0;

y = amp * np.exp(-(((x-mu)**2)/(2*var*var)))

plt.figure()
plt.plot(x,y)
plt.xlim(0,200)