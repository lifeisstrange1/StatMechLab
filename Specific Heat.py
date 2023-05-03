import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

# Dulong and petit plot
Cv_3R = [1,1,1,1,1,1,1,1,1,1,1]
T = [0,10,20,30,40,50,60,70,80,90,100]
plt.plot(T,Cv_3R,label = "Cv = 3R")
plt.xlabel("T(K)")
plt.ylabel("Cv/3R")
plt.title("Dulong and petit plot")
plt.legend()
plt.grid()
plt.show()


#Einstein plot
def G(x):
    y = []
    for i in range(len(x)):
        if x[i]==0.5:
            y.append(1)
        else:
            y.append(0)
    return y
        
V_Ve = np.linspace(0,2,21)
plt.plot(V_Ve,G(V_Ve))
plt.xlabel("V/Ve")
plt.ylabel("G(V)/N_a*f")
plt.title("Einstein Density plot")
plt.legend()
plt.grid()
plt.show()

def Cv(x):
    y = (((1/x)**2)*(np.exp(1/x)))/(np.exp(1/x) - 1)**2
    return y

T = np.linspace(0,2,50)
plt.plot(T,Cv(T),marker = ".")
plt.xlabel("x = T/theta_e")
plt.ylabel("Cv(x)")
plt.title("Einstein Theory of specific heat")
plt.legend()
plt.grid()
plt.show()

#Debye plot
def D(x):
    y = x**2
    return y

x = np.linspace(0,2,20)
plt.plot(x,D(x))
plt.xlabel("V/Vd")
plt.ylabel("D(V)")
plt.title("Debye Density plot")
plt.legend()
plt.grid()
plt.show()

x = np.linspace(0,2,21)
plt.plot(x,D(x),label = "Debye Density")
plt.plot(x,G(x),label = "Einstein Density")
plt.title("Density plot")
plt.legend()
plt.grid()
plt.show()

x = np.linspace(0,2,20)

f = np.vectorize(lambda x : x**3/(np.exp(x)-1))
def Cv1(x):
    result = []
    for i in range(len(x)):
        y = (-3)*(x[i])/(np.exp(x[i])-1) + (12/(x[i])**3)*quad(f,0,x[i])[0]
        result.append(y)
    return result

x = np.linspace(0,2,100)
y = Cv1(1/x)
plt.plot(x,y,marker = ".")
plt.xlabel("theta/T")
plt.ylabel("Cv/3R")
plt.title("Debye model for specific heat")
plt.grid()
plt.show()


T = np.linspace(0,2,50)
Cv_3R = np.ones(50)
plt.plot(T,Cv_3R,label ="Dulong and petit plot" ,marker = ".")
plt.plot(T,Cv(T),label = "Einstein Theory of specific heat",marker = ".")
plt.plot(x,y,label = "Debye model for specific heat",marker = ".")
plt.legend(loc = "4")
plt.grid()
plt.show()








import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

def dulong_petit(x):
    return np.ones(len(x))

def einstein(x):
    return ((1/x)**2) * ((np.exp(1/x))/(np.exp(1/x) - 1)**2)

einstein = np.vectorize(einstein)

def debye(x):
    u1 = -3*x/(np.exp(x)-1)
    u2 = 12/x**3
    ine = quad(lambda x : (x**3/(np.exp(x)-1)),0,x)
    return u1 + u2*ine[0]

debye = np.vectorize(debye)

def dos_d(v):
    if v <= 1/r:
       return v**2
    else:
       return 0

def dos_e(v):
    if abs(v-1) <= 0.02:
       return 1
    else:
       return 0

dos_d = np.vectorize(dos_d)
dos_e = np.vectorize(dos_e)
    

if __name__ == "__main__":
    theta_e = 230
    theta_d = 343
    r = theta_e/theta_d

    # EINSTEIN LAW
    
    T_einstein = np.linspace(1e-8,2*theta_e,100)/theta_e

    Cv_einstein = einstein(T_einstein)
    
    # DEBYE LAW

    T_debye = np.linspace(1e-8,2*theta_d,100)/theta_d

    Cv_debye = debye(1/T_debye)
    
    # DULONG AND PETIT LAW
    
    T_dulong = np.linspace(T_debye[0],T_debye[-1],100) 

    Cv_dulong = dulong_petit(T_dulong)
    
    # DENSITY OF STATES 
    
    v = np.linspace(0,2,100)

    fig,ax = plt.subplots(figsize = (16,8))
    
    ax.plot(v,dos_e(v),label = 'Einstein theory')
    ax.plot(v,dos_d(v)/dos_d(1/r),label = 'Debye theory')
    ax.scatter(v,dos_d(v)/dos_d(1/r),marker = ".",c = "darkorange")
    ax.set_ylim([0,1.1])
    ax.set(xlabel = '$\\nu$/$\\nu$*',title = "Density of states",ylabel = "No. of states")
    ax.legend()
    ax.grid(ls = "--")
    fig.savefig(f"exp3_1.png")
    plt.show()
    
    # PLOTTING

    fig1,ax1 = plt.subplots(figsize = (16,8))

    ax1.plot(T_dulong,Cv_dulong,label = "Dulong & Petit's Law",ls = "--")
    ax1.plot(T_debye,Cv_debye,label = "Debye's Law")
    ax1.plot(T_einstein,Cv_einstein,label = "Einstein's Law",)
    ax1.scatter(T_debye,Cv_debye,marker = ".",c = "darkorange")
    ax1.scatter(T_einstein,Cv_einstein,marker = ".",c = "green")
    ax1.set(xlabel = "T/T*",ylabel = "$\dfrac{C_v}{3R}$",title = "Specific Heat of Solids")
    ax1.legend()
    ax1.grid(ls = "--")
    fig1.savefig(f"exp3_2.png")
    plt.show()
