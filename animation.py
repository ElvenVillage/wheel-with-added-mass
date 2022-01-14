import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

m_0 = 100  #Масса однородного колеса
m_ext = 30 #Масса точечной нагрузки
R = 20     #Радиус колеса
m = m_0 + m_ext
R_c = m_ext*R/(m_0+m_ext)   

data = pd.read_csv('data.csv')

phi = data.iloc[:,2]
dt = 0.025

fig = plt.figure()
fig.set_dpi(100)
fig.set_size_inches(5,5)

ax = plt.axes(xlim=(-50, 400), ylim=(-100, 100))
point = plt.Circle((R_c, R), radius=0.5, fc='y', color='red')
patch = plt.Circle((0, R), radius=R, fc='y',)

def init():
    patch.center = (0, R)
    point.center = (R_c, R)
    ax.add_patch(patch)
    ax.add_patch(point)
    return patch, point, 

def animate(i):
    x, y = patch.center
    x_c, y_c = point.center
    
    x = data.iloc[:,0][i]
    y = data.iloc[:,1][i]
    x_c = x + R_c * np.cos(phi[i])
    y_c = y + R_c * np.sin(phi[i])



    patch.center = (x, y)
    point.center = (x_c, y_c)

    return patch, point, 

anim = animation.FuncAnimation(fig, animate, 
                               init_func=init, 
                               frames=len(phi), 
                               repeat=False,
                               interval=dt*100*6,
                               blit=True)

plt.plot([x for x in range(-50, 400)], [0 for x in range(-50, 400)])
plt.show()