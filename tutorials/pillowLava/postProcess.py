import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# File paths

path = "postProcessing/sampleSets/"


gradient_colors = ['#0d47a1', '#1565c0', '#1976d2', '#1e88e5', '#42a5f5', '#90caf9']

inclinations=[0,10,20,30,40,50]
listTime=np.arange(50,20000,50)
totalL=80
L_max=[]
L_deg = [[] for _ in inclinations]  # Initialize empty lists for each inclination
T_deg = [[] for _ in inclinations]

for theta_idx, theta in enumerate(inclinations):
    for i in range(len(listTime)):
        input_file="deg"+str(theta)+"_PillowLava/"+path+str(listTime[i])+"/gauge_1_alpha.lava.xy"

        # Load the data
        data = pd.read_csv(input_file, delim_whitespace=True, header=None)
        for k in range(len(data[1])):
            # print ("L: ", data[0][k], ", alpha: ", data[1][k])
            if data[1][k]>0.5:
                    L_deg[theta_idx].append(totalL - data[0][k])
                    T_deg[theta_idx].append(listTime[i] / 3600)  # Convert time to hours
                    break
    L_max.append(L_deg[theta_idx][-1])



# Plot
# Plot
fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2,
                                figsize=(10, 6))

# Plot for each inclination
for i, theta in enumerate(inclinations):
    ax0.plot(np.array(T_deg[i]), np.array(L_deg[i]), color=gradient_colors[i], label=f'$\\theta={theta}^o$')

ax0.set_xlabel('$T$ [h]')
ax0.set_ylabel('$L_{front}$ [m]')
ax0.set_xlim([listTime[0] / 3600-0.2, listTime[-1] / 3600+0.5])
ax0.grid()
ax0.legend()

# Plot the final L_max values vs inclination angle
ax1.plot(inclinations, L_max, marker='o', color="k", label='$\\theta=0^o$')
ax1.set_xlabel('$\\theta$ [$^o$]')
ax1.set_ylabel('Final $L_{front}$ [m]')
ax1.set_xlim([inclinations[0] - 2, inclinations[-1] + 2])


ax1.grid()

plt.show()
