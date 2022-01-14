import numpy as np
from scipy import integrate

g = 9.81

m_0 = 100  #Масса однородного колеса
m_ext = 30 #Масса точечной нагрузки
R = 20     #Радиус колеса


m = m_0 + m_ext

Theta = (m_0 + m_ext*m_0/(m_0+m_ext))*R**2 #Момент инерции
R_c = m_ext*R/(m_0+m_ext)                  #Смещение цетра тяжести


dt = 0.025 #Шаг интегрирования

#Начальные условия
phi_0 = np.pi
dot_phi_0 = -2

def N_force(phi, ddot_phi):
	return m*(ddot_phi*R_c*np.cos(phi) - phi**2*R_c*np.sin(phi) + g)

#Режим качения без проскальзывания
def rollingSolver(Y, _):
    enumerator = R_c * np.cos(Y[0]) * (R * Y[1]**2 + g)
    denumerator = Theta/m + R_c**2 + 2*R*R_c *np.cos(Y[0]) + R**2

    return [Y[1], (-1)*enumerator / denumerator ]


def checkLiftoff(phi, dot_phi):
    N = N_force(phi, dot_phi)
    return N < 0

def checkLanding(y):
    return y < R 

solution = []

isFlying = False

def handleFlying(t_start, phi_0, dot_phi_0, x_0):
    t = np.arange(t_start, t_start + 30, dt)
    phi = phi_0
    v_b = dot_phi_0 * R * (-1)
    v_c_x = v_b - dot_phi_0*R_c*np.sin(phi_0)
    v_c_y = dot_phi_0*R_c*np.cos(phi_0)

    x = x_0; y = R #Геометрический центр
    x_c = x + R_c * np.cos(phi); y_c = y + R_c * np.sin(phi) #Центр масс

    for i in range(len(t)):
        cond = checkLanding(y)
        if cond:
            return t[i]
        
        phi += dot_phi_0 * dt
        x_c += v_c_x * dt; y_c += v_c_y * dt
        v_c_y -= g*dt
        x = x_c - R_c*np.cos(phi); y = y_c - R_c*np.sin(phi)

        solution.append([x, y, phi, dot_phi_0])

def handleRolling(t_start, phi_0, dot_phi_0, x_0):
    if (len(solution) > 99):
        pass
    t_0 = t_start
    iterations = 0
    while True:
        iterations += 1
        if (iterations > 10):
            return t_start
        
        t = np.arange(t_0, t_0 + 30, dt)
        sol = integrate.odeint(rollingSolver, [phi_0, dot_phi_0], t)

        phi = sol[:, 0]
        dot_phi = sol[:, 1]
        v_b = dot_phi * R * (-1)
        x = x_0

        for i in range(len(phi)):
            cond = checkLiftoff(phi[i], dot_phi[i])
            if (cond):
                return t[i]
            x += v_b[i] * dt
            solution.append([x, R, phi[i], dot_phi[i]])
        t_0 += 30

def get_last_entry(t_start):
    return solution[int(t_start/dt)-1]


t_start = 0; x_start = 0; phi_0_start = phi_0
for i in range(13):
    t_start = handleRolling(t_start, phi_0_start, dot_phi_0, x_start)
    last_entry = get_last_entry(t_start)
    t_start = handleFlying(t_start, last_entry[2], last_entry[3], last_entry[0])
    last_entry = get_last_entry(t_start)
    x_start = last_entry[0]; phi_0_start = last_entry[2]

np.savetxt('data.csv', solution, delimiter=',', fmt='%.5f')