import numpy as np
import matplotlib.pyplot as plt

def f(t, y):
    return -2*y + t**2

def g(k ,y):
    return -k*y

def h(t, theta, omega, g, l):
    dtheta_dt = omega
    domega_dt = - (g / l) * np.sin(theta)
    return dtheta_dt, domega_dt

def euler_method(f, theta0, omega0, t0, tf, h, g, l):
    t_values = np.arange(t0, tf + h, h)
    theta_values = np.zeros(len(t_values))
    omega_values = np.zeros(len(t_values))

    theta_values[0] = theta0
    omega_values[0] = omega0
    t = t0
    theta = theta0
    omega = omega0
    for i in range(1, len(t_values)):
        dtheta_dt, domega_dt = f(t, theta, omega, g, l)
        theta = theta + h * dtheta_dt
        omega = omega + h * domega_dt
        t = t_values[i]
        theta_values[i] = theta
        omega_values[i] = omega

    return t_values, theta_values, omega_values


def runge_kutta_4(f, y0, t0, tf, h):
    t_values = np.arange(t0, tf + h, h)  
    y_values = np.zeros(len(t_values))  

    y_values[0] = y0
    t = t0
    y = y0

    #runge kutta ordre 4
    for i in range(1, len(t_values)):
        k1 = f(t, y)
        k2 = f(t + h/2, y + h/2 * k1)
        k3 = f(t + h/2, y + h/2 * k2)
        k4 = f(t + h, y + h * k3)
        y = y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
        t = t_values[i]
        y_values[i] = y

    return t_values, y_values

y0 = 1    
t0 = 0    
tf = 2    
h = 0.1   
k=1

#Euler

t_values, y_values = euler_method(f, y0, t0, tf, h)
plt.plot(t_values, y_values, label='euler')
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.title('resultat')
plt.grid(True)
plt.show()

t_values, y_values = runge_kutta_4(f, y0, t0, tf, h)

#Runge Kutta
plt.plot(t_values, y_values, label='Runge-Kutta ordre 4')
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.title('Résolution de l\'EDO par la méthode de Runge-Kutta d\'ordre 4')
plt.grid(True)
plt.show()