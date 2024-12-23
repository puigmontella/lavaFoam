import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# File paths
path = "postProcessing/sampleSets/"  # Replace with your file name



inclinations=[0,30]
listTime=np.arange(50,50000,50)
totalL=80
L_max=[]
L_deg=np.zeros(len(inclinations),len(listTime))
T_deg=[]

for theta in range(len(inclinations)):

    for i in range(len(listTime)):
        input_file=path+str(listTime[i])+"/gauge_1_alpha.lava.xy"

        # Load the data
        data = pd.read_csv(input_file, delim_whitespace=True, header=None)
        for k in range(len(data[1])):
            # print ("L: ", data[0][k], ", alpha: ", data[1][k])
            if data[1][k]>0.5:
                L_0deg.append(totalL-data[0][k])
                T_0deg.append(listTime[i]/3600)
                break





# Plot
plt.figure(figsize=(8, 6))
plt.plot(T_0deg, L_0deg, marker='o', label='$\\theta=0^o$')
plt.xlabel('$T$ [h]}')
plt.ylabel('$L_{front}$ [m]')
# plt.title('Alpha vs. Distance')
plt.grid()
plt.legend()
plt.show()
