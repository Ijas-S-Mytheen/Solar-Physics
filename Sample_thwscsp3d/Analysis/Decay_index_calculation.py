import numpy as np
import pandas as pd

data = np.loadtxt('ac64.dat')

print(data.shape)
data2 = data[:,3]
F = data2.reshape(160,288,511)
r = data[:511,1]
phi = data[:,2]
phi2 = phi.reshape(160,288,511)
p = phi2[0,:,0]
theta = data[:,0]
theta2 = theta.reshape(160,288,511)
t = theta2[:,0,0]

Bt64 = np.load('/home/ijas/Fortran_IIA/testing_code/res/Components2/64Bt.npy')
Bp64 = np.load('/home/ijas/Fortran_IIA/testing_code/res/Components2/64Bp.npy')
Br64 = np.load('/home/ijas/Fortran_IIA/testing_code/res/Components2/64Br.npy')
print(Bt64.shape)

import math
trail_panels = np.linspace(0,510,511).astype(int)
Bh64_values = []
panel_height = []
Bh64_by_height = []
side_length = 5
rows, cols = Bt64[:,:,0].shape
grid_points = side_length * side_length
start_row = (rows - side_length) // 2
start_col = (cols - side_length) // 2
for i in trail_panels:
    Btt64 = Bt64[:,:,i]
    Bpt64 = Bp64[:,:,i]
    Brt64 = Br64[:,:,i]
    square_region_Btt64 = Btt64[start_row:start_row + side_length, start_col:start_col + side_length]
    square_region_Bpt64 = Bpt64[start_row:start_row + side_length, start_col:start_col + side_length]
    square_region_Brt64 = Brt64[start_row:start_row + side_length, start_col:start_col + side_length]
    avg_Btt64 = np.mean(square_region_Btt64)
    avg_Bpt64 = np.mean(square_region_Bpt64)
    avg_Brt64 = np.mean(square_region_Brt64)
    Bh64 = math.sqrt((avg_Btt64)**2 + (avg_Bpt64)**2)
    grid_points = side_length * side_length
    Bh64_values.append(Bh64)
    panel_height.append(r[i])
# Bh_by_height.append(Bh/r[i])
# print("Bh @ radius",r[i],'is',Bh)
print('number of entries of Bh:',len(Bh64_values))
#print(panel_height)
#print(Bh_by_height)
panel_height2 = np.array(panel_height[1:510]) - 1
log_Bh64 = np.log(np.array(Bh64_values[1:510])) # I elemated the last entri as iit was zero
log_h = np.log(panel_height2)


d_log_Bh64 = np.gradient(log_Bh64)
d_log_h = np.gradient(log_h)
decay_index64 = -d_log_Bh64/d_log_h

D_I = pd.DataFrame({'VAR64' : decay_index64[:100]})
D_I.to_csv('Decay_index.csv')

import matplotlib.pyplot as plt

fig,ax = plt.subplots(figsize = (12,8))
ax.plot(panel_height[:100], decay_index64[:100],label='VAR64')
ax.axhline(y = 1.5, color = 'red', linestyle = '--', label = 'Decay Index = 1.5')

for label in (ax.get_xticklabels() + ax.get_yticklabels()):
    label.set_fontsize(10)
ax.spines['bottom'].set_linewidth(1.5)
ax.spines['top'].set_linewidth(1.5)
ax.spines['right'].set_linewidth(1.5)
ax.spines['left'].set_linewidth(1.5)
ax.tick_params(direction = 'inout', length = 10, width = 1.5, colors = 'k')

plt.xlabel('panel_height in solar radii', fontsize = 13)
plt.ylabel('Decay Index',fontsize = 13)
plt.title('Decay Index vs panel_height',fontsize = 13)
plt.legend(fontsize = 13)
plt.grid(False)
plt.show()
fig.savefig('decay_index_plot.png',dpi = 500) #bbox_inches = 'tight'
