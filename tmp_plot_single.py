import numpy as np
import matplotlib.pyplot as plt
import tmp_analysis_tools as antools

data = antools.read_dataset("", runName="default")

antools.summary_plot(data, burnInPeriod=2000)
