import numpy as np
import matplotlib.pyplot as plt



#-----------FONCTION------------

def Fonction_ex(y, t):
    return -2*y + t**2

#Q2

def Fonction_Q2(k ,y):
    return -k*y

def Fonction_Q2_Ana(k,t):
    return 

#Q3

def Fonction_Q3(t, theta, omega, g, l):
    dtheta_dt = omega
    domega_dt = - (g / l) * np.sin(theta)
    return dtheta_dt, domega_dt

#--------------------------------

#------------METHODE-------------
#EULER
def euler_method2(f, theta0, omega0, t0, tf, h, g, l):
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

def euler_method_simple(f,y0,t):
    y = np.zeros(len(t))
    y[0] = y0
    for i in range(0,len(t)-1):
        y[i+1] = y[i] + f(y[i],t[i])*(t[i+1]-t[i])
    return y

#RANGE KUTTA ORDRE 4
def RK_method(f,y0,t):
    y = np.zeros(len(t))
    y[0] = y0
    for i in range(0,len(t)-1):
        h = t[i+1]-t[i]
        F1 = h*f(y[i],t[i])
        F2 = h*f((y[i]+F1/2),(t[i]+h/2))
        F3 = h*f((y[i]+F2/2),(t[i]+h/2))
        F4 = h*f((y[i]+F3),(t[i]+h))
        y[i+1] = y[i] + 1/6*(F1 + 2*F2 + 2*F3 + F4)
    return y

#--------------------------------

#--------------INIT--------------

#Euler
    #Init
depart=0
arret=5
y0 = 1
points = 210
t = np.linspace(depart,arret,points)

    #Analytique
t_ana = np.linspace(depart,arret,points)
y_ana = np.exp(-y0*t_ana)

    #Plot
plt.plot(t,euler_method_simple(Fonction_Q2,y0,t),'b-', label = "Euler")
plt.plot(t_ana,y_ana,'ro',label = "Analytical Solution")
plt.plot(t,RK_method(Fonction_Q2,y0,t),'k-', label = "Runge Kutta")
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.title('resultat Euler simple')
plt.grid(True)
plt.show()

#RK4
y_RK = RK_method(Fonction_Q2,y0,t)






