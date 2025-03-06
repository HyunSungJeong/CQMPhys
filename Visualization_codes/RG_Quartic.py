import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Define parameters
lambda_z_critical = 1 - 1/np.sqrt(2)
#C_values = np.append(np.linspace(-0.5*(1 + np.log(2)), 0.1, 7), [-0.8555, -0.9, -1, -0.7])  # Added new C values
C_values1 = [-0.9, -0.88, -0.86, -0.5*(1 + np.log(2))]
color1a = ['black', 'black', 'black', 'black']
color1b = ['black', 'black', 'black', 'black']
C_values2 = [-0.82, -0.79]
color2 = ['black', 'black', 'black']

# Define the function for flow lines
def lambda_x(lambda_z, C):
    
    result = 2.24 * np.sqrt((lambda_z - 1)**2 - np.log(np.abs(lambda_z - 1)) + C)

    result[np.isnan(result)] = 0  # Avoid NaNs for better plotting
    return result

# Create figure
fig, ax = plt.subplots(figsize=(8, 5))
ymax = 0.85
ax.set_xlim([0, 0.5])
ax.set_ylim([0, ymax])
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
ax.set_xlabel("$\\lambda_z$", fontsize=22)
ax.set_ylabel("$\\lambda_x$", fontsize=22)

# Generate curved flow lines
for it in range(0,len(C_values1),1):
    C = C_values1[it]
    lambda_z_vals = np.linspace(0, lambda_z_critical, 300)
    lambda_x_vals = lambda_x(lambda_z_vals, C)
    lambda_z_vals = lambda_z_vals[lambda_x_vals>0]
    lambda_x_vals = lambda_x_vals[lambda_x_vals>0]
    lambda_z_vals = lambda_z_vals[lambda_x_vals<ymax]
    lambda_z_vals = np.append(lambda_z_vals, lambda_z_vals[-1]+0.001)
    lambda_x_vals = lambda_x(lambda_z_vals, C)
    valid = ~np.isnan(lambda_x_vals)
    ax.plot(lambda_z_vals[valid], lambda_x_vals[valid], color=color1a[it], linewidth=1)
    
    # Adding flow direction arrows along the curves
    for i in range(len(lambda_z_vals[valid])//4, len(lambda_z_vals[valid])-1, len(lambda_z_vals[valid])//2-1):
        dx0 = lambda_z_vals[valid][i+1] - lambda_z_vals[valid][i-1]
        dy0 = lambda_x_vals[valid][i+1] - lambda_x_vals[valid][i-1]
        dx = 0.01*dx0/np.sqrt(dx0**2+dy0**2)
        dy = 0.01*dy0/np.sqrt(dx0**2+dy0**2)
        ax.quiver(lambda_z_vals[valid][i], lambda_x_vals[valid][i], dx, dy,
                  angles='xy', scale=0.5, width=0.01, headwidth=2, headlength=4, headaxislength=4, color=color1a[it])
        
for it in range(0,len(C_values1),1):
    C = C_values1[it]
    lambda_z_vals = np.linspace(lambda_z_critical, 0.5, 300)
    lambda_x_vals = lambda_x(lambda_z_vals, C)
    lambda_z_vals = lambda_z_vals[lambda_x_vals>0]
    lambda_x_vals = lambda_x_vals[lambda_x_vals>0]
    lambda_z_vals = lambda_z_vals[lambda_x_vals<ymax]
    lambda_z_vals = np.append(lambda_z_vals[0]-0.001, lambda_z_vals)
    lambda_x_vals = lambda_x(lambda_z_vals, C)
    valid = ~np.isnan(lambda_x_vals)
    ax.plot(lambda_z_vals[valid], lambda_x_vals[valid], color=color1b[it], linewidth=1)
    ax.plot(lambda_z_vals[valid], -lambda_x_vals[valid], color=color1b[it], linewidth=1)  # Reflect for symmetry
    
    # Adding flow direction arrows along the curves
    for i in range(len(lambda_z_vals[valid])//4, len(lambda_z_vals[valid])-1, len(lambda_z_vals[valid])//2-1):
        dx0 = lambda_z_vals[valid][i+1] - lambda_z_vals[valid][i-1]
        dy0 = lambda_x_vals[valid][i+1] - lambda_x_vals[valid][i-1]
        dx = 0.01*dx0/np.sqrt(dx0**2+dy0**2)
        dy = 0.01*dy0/np.sqrt(dx0**2+dy0**2)
        ax.quiver(lambda_z_vals[valid][i], lambda_x_vals[valid][i], dx, dy,
                  angles='xy', scale=0.5, width=0.01, headwidth=2, headlength=4, headaxislength=4, color=color1b[it])
        
for it in range(0, len(C_values2), 1):
    C = C_values2[it]
    lambda_z_vals = np.linspace(0, 0.5, 300)
    lambda_x_vals = lambda_x(lambda_z_vals, C)
    lambda_z_vals = lambda_z_vals[lambda_x_vals>0]
    lambda_x_vals = lambda_x_vals[lambda_x_vals>0]
    lambda_z_vals = lambda_z_vals[lambda_x_vals<ymax]
    lambda_x_vals = lambda_x(lambda_z_vals, C)
    valid = ~np.isnan(lambda_x_vals)
    ax.plot(lambda_z_vals[valid], lambda_x_vals[valid], color=color2[it], linewidth=1)
    
    # Adding flow direction arrows along the curves
    for i in range(len(lambda_z_vals[valid])//4, len(lambda_z_vals[valid])-1, len(lambda_z_vals[valid])//2-1):
        dx0 = lambda_z_vals[valid][i+1] - lambda_z_vals[valid][i-1]
        dy0 = lambda_x_vals[valid][i+1] - lambda_x_vals[valid][i-1]
        dx = 0.01*dx0/np.sqrt(dx0**2+dy0**2)
        dy = 0.01*dy0/np.sqrt(dx0**2+dy0**2)
        ax.quiver(lambda_z_vals[valid][i], lambda_x_vals[valid][i], dx, dy,
                  angles='xy', pivot='middle', scale=0.5, width=0.01, headwidth=2, headlength=4, headaxislength=4, color=color2[it])

# Annotate critical dashed line
ax.axvline(lambda_z_critical, linestyle='dashed', color='k')
ax.text(lambda_z_critical - 0.07, ymax - 0.05, "$\\lambda_z = 1 - \\frac{1}{\sqrt{2}}$", fontsize=18)

# Define custom colormap for gradation effect
cmap_brown = LinearSegmentedColormap.from_list("brown_grad", [(1, 0.6, 0.4), (1, 1, 1)], N=100)
cmap_blue = LinearSegmentedColormap.from_list("blue_grad", [(1, 1, 1), (0, 0, 1)], N=100)

# Shade regions with gradient effect
Z_vertical_SFL = np.linspace(0, 1, 100).reshape(-1, 1)  # Vertical gradation for SFL
Z_diagonal_FL = np.zeros((1,100))                       # Diagonal gradation for FL

for it in range(1,150):
    Z_col = np.zeros((1,100-2*it//3))
    fadelen = 50
    if 2*it//3 > fadelen:
        Z_col = np.concatenate((Z_col, np.linspace(0,100,fadelen).reshape(1,fadelen)), axis=1)
        Z_col = np.concatenate((Z_col, 100*np.ones((1,2*it//3-fadelen))), axis=1)
    else:
        Z_col = np.concatenate((Z_col, np.linspace(0,(100/fadelen)*2*it//3, 2*it//3).reshape(1,2*it//3)), axis=1)
    
    Z_diagonal_FL = np.concatenate((Z_diagonal_FL, Z_col), axis=0)

# Shade SFL region with vertical gradient
ax.imshow(Z_vertical_SFL, extent=[0, 0.2929, 0, 0.05], origin='lower', cmap=cmap_brown, aspect='auto', alpha=0.6)

# Shade FL region with horizontal gradient
ax.imshow(Z_diagonal_FL, extent=[0.4, 0.5, ymax-0.6, ymax], origin='lower', cmap=cmap_blue, aspect='auto', alpha=0.3)

# Labels
ax.text(0.05, 0.01, "SFL", fontsize=14, color='brown', fontweight='bold')
ax.text(0.475, ymax-0.085, "FL", fontsize=14, color='blue', fontweight='bold')

plt.show()
