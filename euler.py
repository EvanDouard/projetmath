import numpy as np
import matplotlib.pyplot as plt



#-----------FONCTION------------

def Fonction_ex(y, t):
    return -2*y + t**2

#Q2

def Fonction_Q2(y,t):
    return -k*y


#Q3

def Fonction_Q3(y, t, g, l):#y t b c
    theta=y[0]
    omega=y[1]
    dtheta_dt = omega
    domega_dt = - (g / l) * np.sin(theta)
    return np.array([dtheta_dt,domega_dt])

#--------------------------------

#------------METHODE-------------
#EULER
def euler_method2(f, y0, t, g, l):
    y = np.zeros((len(t), 2)) 

    y[0] = y0
    for i in range(0, len(t) - 1):
        dt = t[i + 1] - t[i]
        y[i + 1] = y[i] + f(y[i], t[i], g, l) * dt

    return y

def RK_method2(f, y0, t, g, l):
    y = np.zeros((len(t), 2))  # Tableau pour theta et omega
    y[0] = y0

    for i in range(0, len(t) - 1):
        h = t[i + 1] - t[i]

        F1 = h * f(y[i], t[i], g, l)
        F2 = h * f(y[i] + F1 / 2, t[i] + h / 2, g, l)
        F3 = h * f(y[i] + F2 / 2, t[i] + h / 2, g, l)
        F4 = h * f(y[i] + F3, t[i] + h, g, l)

        y[i + 1] = y[i] + (F1 + 2 * F2 + 2 * F3 + F4) / 6

    return y

def euler_method_simple(f,y0,t,args=()):
    if (isinstance(y0,int)):
       y = np.zeros(len(t), ) 
    else:
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
k=1
points = 21
t = np.linspace(depart,arret,points)

    #Analytique
t_ana = np.linspace(depart,arret,points)
y_ana = y0*np.exp(-k*t_ana)

    #Plot
plt.plot(t,euler_method_simple(Fonction_Q2,y0,t),'b-', label = "Euler")
plt.plot(t_ana,y_ana,'ro',label = "Solution analytique")
plt.plot(t,RK_method(Fonction_Q2,y0,t),'k-', label = "Runge Kutta ordre 4") #RK4
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.title('resultat Euler simple')
plt.grid(True)
plt.show()

# Init
g = 1
l = 1
points=100
t = np.linspace(0, 10, points)
y0 = [np.pi / 4, 0]  #theta = pi/4 omega = 0


y_euler = euler_method2(Fonction_Q3, y0, t, g, l)
y_rk = RK_method2(Fonction_Q3, y0, t, g, l)


# Plot
plt.plot(t, y_euler[:, 0], 'b-', label="Euler - θ")
plt.plot(t, y_euler[:, 1], 'b--', label="Euler - ω")
plt.plot(t, y_rk[:, 0], 'k-', label="Runge Kutta ordre 4 - θ")
plt.plot(t, y_rk[:, 1], 'k--', label="Runge Kutta ordre 4 - ω")
plt.xlabel('Temps')
plt.ylabel('Valeurs')
plt.legend()
plt.show()






