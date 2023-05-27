import csv
import re
import numpy as np
import matplotlib.pyplot as plt


# list of csv files 
files = [".\Water_3_Assigned.csv",
         ".\Water_3_Unssigned.csv",
         ".\Water_3_Unmatched.csv"]


all_mzs = []
all_file_names = []
all_int = []

# load data from csv files
for file in files:
    mz = []
    Intensity = []
    with open(file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            mz.append(float(row["mz"]))
            Intensity.append(float(row["Intensity"]))
    
    
    x = np.array(mz)
    file_name = re.sub(".csv","",file.rsplit("/",1)[-1])
    all_mzs.append(x)
    all_int.append(Intensity)
    all_file_names.append(file_name)
    

# Define parameters of matplotlib 
plt.style.use("seaborn-bright")
plt.style.use("seaborn")
plt.rcParams['axes.facecolor']='white'
plt.rcParams['axes.edgecolor']="#000000"
plt.rcParams["axes.linewidth"]=.5
plt.rcParams["xtick.major.size"]=3
plt.rcParams["ytick.major.size"]=3


# Geneate plots for each file
for i in range(len(all_mzs)):
    x = all_mzs[i]
    Intensity=all_int[i]
    file_name = all_file_names[i]
    plt.subplots()
    plt.scatter(x,x%1,s=.1)
    plt.plot([0,2000],[0,2000*0.0005],color="r") # Peptide slope
    plt.ylabel("m/z mod 1")
    plt.ylim(0,1)
    plt.xlim(0,3000)
    plt.xlabel("m/z")
    plt.title(re.sub("_"," ",file_name))
    plt.show()