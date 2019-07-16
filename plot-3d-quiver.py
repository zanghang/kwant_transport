import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#===========================================================================
#  below for reading the data, merge the data,  and copy the data for plot
df = pd.read_csv('current_density.dat',sep='\s+',header=None,engine='python')

df_sum = df.groupby([0,1,2,3,4,5],as_index=False).sum()

#for i in range(len(df_sum)):
#    print(df_sum[0][i],df_sum[1][i],df_sum[2][i],df_sum[3][i],df_sum[4][i],df_sum[5][i],df_sum[6][0]) 

data = np.zeros((len(df_sum),7))

for i in range(len(df_sum)):
    for j in range(7):
        data[i][j] = df_sum[j][i]

x = []
y = []
z = []
u = []
v = []
w = []
value = []

for i in range(len(data)):
    x.append(data[i][3])
    y.append(data[i][4])
    z.append(data[i][5])
    u.append(data[i][0])
    v.append(data[i][1])
    w.append(data[i][2])
    value.append(data[i][6])
#===========================================================================
#==below for modifying the data, make all the hopping values to be positive
u_tmp = np.zeros(len(data))
v_tmp = np.zeros(len(data))
w_tmp = np.zeros(len(data))

for i in range(len(data)):
    if value[i] < 0:
        u_tmp[i] = u[i]
        v_tmp[i] = v[i]
        w_tmp[i] = w[i]
        u[i] = x[i]
        v[i] = y[i]
        w[i] = z[i]
        x[i] = u_tmp[i]
        y[i] = v_tmp[i]
        z[i] = w_tmp[i]
        value[i] = -value[i]
#===========================================================================
#==below for modifying the vector u v w according to the hopping value
new_u = []
new_v = []
new_w = []
r_value = []

for i in range(len(data)):
    r_value.append(  np.sqrt((u[i]-x[i])**2 + (v[i]-y[i])**2 + (w[i]-z[i])**2 ) )
    new_u.append( np.sqrt(value[i]) * (u[i]-x[i])/r_value[i] )
    new_v.append( np.sqrt(value[i]) * (v[i]-y[i])/r_value[i] )
    new_w.append( np.sqrt(value[i]) * (w[i]-z[i])/r_value[i] )
#===========================================================================
# below is for the color map , and for extracting partial data

fig = plt.figure()
#fig = plt.figure(figsize=(28,38))

cmap = 'jet'
#cmap = 'gnuplot'

value_norm = value/np.max(value)


x_part = []
y_part = []
z_part = []
new_u_part = []
new_v_part = []
new_w_part = []
color = []

for i in range(len(x)):
    if value_norm[i] > 0.0:
        x_part.append(x[i])
        y_part.append(y[i])
        z_part.append(z[i])
        new_u_part.append(new_u[i])
        new_v_part.append(new_v[i])
        new_w_part.append(new_w[i])
        color.append(value_norm[i])

color = np.concatenate((color, np.repeat(color,2))) # Repeat for each body line and two head lines

color = getattr(plt.cm,cmap)(color)   # it cmap = 'jet' , means color = plt.cm.jet(color)

ax = fig.gca(projection='3d')

q=ax.quiver(x_part, y_part, z_part, new_u_part, new_v_part, new_w_part, length=1, normalize=False ,arrow_length_ratio=0.4, linewidths=0.5, cmap=cmap)
#q=ax.quiver(x_part, y_part, z_part, 0, 0, new_w_part, length=20, normalize=False ,arrow_length_ratio=0.4, linewidths=0.1, cmap=cmap)
#q=ax.quiver(x, y, z, new_u, new_v, new_w, length=1, normalize=False, arrow_length_ratio=0.6, linewidths=0.5, cmap=cmap)

# for the colorbar
q.set_array(np.linspace(0,np.max(color),10))
q.set_edgecolor(color)
q.set_facecolor(color)
fig.colorbar(q)



#===========================================================================
coord_x = []
coord_y = []
coord_z = []
atom_type = []

with open('../cp2k-device-2Au/device.xyz','r') as f:
    f.readline()
    f.readline()
    for line in f:
        a, x1, x2, x3 = line.split()
        coord_x.append(float(x1))
        coord_y.append(float(x2))
        coord_z.append(float(x3))
        atom_type.append(a)
f.close()

atom_type_color = ['r'  for i in range(len(atom_type))]

for i in range(len(atom_type_color)):
    if atom_type[i] == 'Au' :
        atom_type_color[i] = '#FF1493'
    elif atom_type[i] == 'C' :
        atom_type_color[i] = '#00FF00'
    elif atom_type[i] == 'H' :
        atom_type_color[i] = '#FFFFFF'
    elif atom_type[i] == 'N' :
        atom_type_color[i] = '#8F8FFF'
    elif atom_type[i] == 'O' :
        atom_type_color[i] = '#F00000'
    elif atom_type[i] == 'S' :
        atom_type_color[i] = '#FFC832'

ax.scatter(coord_x,coord_y,coord_z,c=atom_type_color,s=0.9,depthshade=False,edgecolor='')#,marker='o')


ax.set_xlim(-10,10)
ax.set_ylim(-10,10)
ax.set_zlim(-15,66)

ax.set_xlabel('X (Angstrom)')
ax.set_ylabel('Y (Angstrom)')
ax.set_zlabel('Z (Angstrom)')

ax.view_init(elev=5,azim=40)

#plt.show()
plt.savefig('figure.png',format='png',dpi=1400)
