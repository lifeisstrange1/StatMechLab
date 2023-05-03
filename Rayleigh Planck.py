import numpy as np
import matplotlib.pyplot as plt

x=np.linspace(0,12,1000)
K=1.380649*10**(-23)
T=1000
c=3*10**8
h=6.626*10**(-34)
l = 1 * 10**(-10)
T=[1000,1200,1800]

def f_rj(x):
    return (x**2)
def f_p(x):
    return (x**3)/(np.exp(x)-1)

plt.scatter(x,np.pi*f_rj(x))
plt.title("Density of States (Rayleigh-Jeans)")
plt.xlabel("x")
plt.ylabel("G(x)")
plt.legend()
plt.grid(ls='--')
plt.show()

plt.scatter(x,f_p(x))
plt.title('Density of States (Planck)')
plt.xlabel("x")
plt.ylabel("f_rj(x)")
plt.grid(ls='--')
plt.show()

for i in T:
    A = (8*np.pi)*(K*i)**4/(h*c)**3
    plt.scatter(x*(K*i)/h,A*np.pi*f_rj(x)*h/(K*i),label = 'T='+ str(i)+"K")
plt.legend()
plt.title("Energy Density (Rayleigh-Jeans)")
plt.xlabel('frequency')
plt.ylabel('Energy density (Rayliegh-Jeans)')
plt.grid(ls='--')
plt.show()

for i in T:
    A = (8 * np.pi) * (K * i) ** 4 / (h * c) ** 3
    plt.scatter(x*(K*i)/h,A*f_p(x)*h/(K*i), label = 'T='+ str(i)+"K")
plt.legend()
plt.xlabel('frequency')
plt.ylabel('Energy density (Planck)')
plt.title("Energy Density (Planck)")
plt.grid(ls='--')
plt.show()



def density(v):
    return (8*np.pi)/(c**3)*(v**2)*(l**3)

v_complete=np.logspace(10,30,3000)
v_visible=np.linspace(4e14 ,8e14,20)

plt.scatter(v_complete,density(v_complete))
plt.xscale("log")
plt.yscale("log")
plt.xlabel("Frequency(v)")
plt.ylabel("Density of states")
plt.title("Density of states with frequency")
plt.grid()
plt.show()

plt.scatter(v_visible,density(v_visible))
plt.xlabel("Frequency(v)")
plt.ylabel("Density of states")
plt.title("Density of states with frequency in visible region")
plt.grid()
plt.show()
