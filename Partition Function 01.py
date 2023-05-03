import numpy as np
import matplotlib.pyplot as plt

k = 8.617 * 10**(-5) # eV/K

def Z(g,E,T):
    n = len(E)
    part = []
    
    for j in (T):
        
        z = 0

        for i in range(n):

            z = z + g[i]*np.exp(-E[i]/(k*j))

        part.append(z)

    return part

def frac(g,E,T):
    n = len(E)
    FRAC = np.zeros(n*len(T)).reshape(n, len(T))
    z = Z(g,E,T)
    
    for i in range(len(T)):

        for j in range(n):
            
            f1 = (g[j]*np.exp(-E[j]/(k*T[i]))) / z[i]
            FRAC[j][i] = f1
            
        
        
    return FRAC


T1 = np.linspace(0.0001,5000,1000)
T2 = np.linspace(5000,10**6,1000)

# 2 level
E = [0,1,]
g=[1,1,]


plt.plot(T1,Z(g,E,T1))
plt.xlabel("T")
plt.ylabel("Z")
plt.title("Z vs T (low T)")
plt.grid()
plt.show()

plt.plot(T2,Z(g,E,T2))
plt.xlabel("T")
plt.ylabel("Z")
plt.grid()
plt.title("Z vs T (high T)")
plt.show()



plt.plot(T1,frac(g,E,T1)[0], label = "for energy = "+str(E[0]) + " eV")
plt.plot(T1,frac(g,E,T1)[1], label = "for energy = "+str(E[1]) + " eV")
#plt.plot(T1,frac(g,E,T1)[2], label = "for energy = "+str(E[2]) + " eV")
plt.xlabel("T")
plt.ylabel("N_j / N")
plt.legend(loc="best")
plt.grid()
plt.title("N_j/N vs T (low T)")
plt.show()


plt.plot(T2,frac(g,E,T2)[0], label = "for energy = "+str(E[0]) + " eV")
plt.plot(T2,frac(g,E,T2)[1], label = "for energy = "+str(E[1]) + " eV")
#plt.plot(T2,frac(g,E,T2)[2], label = "for energy = "+str(E[2]) + " eV")
E_mid = np.sum(E)/len(E)
E_mid = np.full(shape = len(T2) ,fill_value = E_mid)
plt.plot(T2,E_mid,'--')
plt.xlabel("T")
plt.ylabel("N_j / N")
plt.legend()
plt.grid()
plt.title("N_j/N vs T (high T)")
plt.show()


#INTERNAL ENERGY

population_1 = frac(g,E,T1)
population_2 = frac(g,E,T2)
U_1,U_2 = 0,0
for i in range(len(population_1)) :
    U_1 += population_1[i]*E[i]

for i in range(len(population_2)) :
    U_2 += population_2[i]*E[i]


plt.plot(T1,U_1)
plt.xlabel("T")
plt.ylabel("U")
plt.title("U vs T (low T)")
plt.grid()
plt.show()

plt.plot(T2,U_2)
plt.xlabel("T")
plt.ylabel("U")
plt.grid()
plt.title("U vs T (high T)")
plt.show()



#ENTROPY
z1 = Z(g,E,T1)
z2 = Z(g,E,T2)
N = 1

S1 = N*k*np.log(np.array(z1) / N) + U_1 / T1 + N*k
S2 = N*k*np.log(np.array(z2) / N) + U_2 / T2 + N*k

plt.plot(T1,S1)
plt.xlabel("T")
plt.ylabel("S")
plt.title("S vs T (low T)")
plt.grid()
plt.show()

plt.plot(T2,S2)
plt.xlabel("T")
plt.ylabel("S")
plt.grid()
plt.title("S vs T (high T)")
plt.show()


#HELMHOLTZ
F1 = -N*k*np.array(T1) * np.log(np.array(z1))
F2 = -N*k*np.array(T2) * np.log(np.array(z2))

plt.plot(T1,F1)
plt.xlabel("T")
plt.ylabel("F")
plt.title("F vs T (low T)")
plt.grid()
plt.show()

plt.plot(T2,F2)
plt.xlabel("T")
plt.ylabel("F")
plt.grid()
plt.title("F vs T (high T)")
plt.show()








import numpy as np
import matplotlib.pyplot as plt

k = 8.617 * 10**(-5) # eV/K

def Z(g,E,T):
    n = len(E)
    part = []
    
    for j in (T):
        
        z = 0

        for i in range(n):

            z = z + g[i]*np.exp(-E[i]/(k*j))

        part.append(z)

    return part

def frac(g,E,T):
    n = len(E)
    FRAC = np.zeros(n*len(T)).reshape(n, len(T))
    z = Z(g,E,T)
    
    for i in range(len(T)):

        for j in range(n):
            
            f1 = (g[j]*np.exp(-E[j]/(k*T[i]))) / z[i]
            FRAC[j][i] = f1
            
        
        
    return FRAC


T1 = np.linspace(0.0001,5000,1000)
T2 = np.linspace(5000,10**6,1000)

# 3 level
E = [0,1,2]
g=[1,1,1]


plt.plot(T1,Z(g,E,T1))
plt.xlabel("T")
plt.ylabel("Z")
plt.title("Z vs T (low T)")
plt.grid()
plt.show()

plt.plot(T2,Z(g,E,T2))
plt.xlabel("T")
plt.ylabel("Z")
plt.grid()
plt.title("Z vs T (high T)")
plt.show()



plt.plot(T1,frac(g,E,T1)[0], label = "for energy = "+str(E[0]) + " eV")
plt.plot(T1,frac(g,E,T1)[1], label = "for energy = "+str(E[1]) + " eV")
plt.plot(T1,frac(g,E,T1)[2], label = "for energy = "+str(E[2]) + " eV")
plt.xlabel("T")
plt.ylabel("N_j / N")
plt.legend(loc="best")
plt.grid()
plt.title("N_j/N vs T (low T)")
plt.show()

plt.plot(T2,frac(g,E,T2)[0], label = "for energy = "+str(E[0]) + " eV")
plt.plot(T2,frac(g,E,T2)[1], label = "for energy = "+str(E[1]) + " eV")
plt.plot(T2,frac(g,E,T2)[2], label = "for energy = "+str(E[2]) + " eV")
plt.xlabel("T")
plt.ylabel("N_j / N")
plt.legend()
plt.grid()
plt.title("N_j/N vs T (high T)")
plt.show()


#INTERNAL ENERGY

population_1 = frac(g,E,T1)
population_2 = frac(g,E,T2)
U_1,U_2 = 0,0
for i in range(len(population_1)) :
    U_1 += population_1[i]*E[i]

for i in range(len(population_2)) :
    U_2 += population_2[i]*E[i]


plt.plot(T1,U_1)
plt.xlabel("T")
plt.ylabel("U")
plt.title("U vs T (low T)")
plt.grid()
plt.show()

plt.plot(T2,U_2)
plt.xlabel("T")
plt.ylabel("U")
plt.grid()
plt.title("U vs T (high T)")
plt.show()



#ENTROPY
z1 = Z(g,E,T1)
z2 = Z(g,E,T2)
N = 1

S1 = N*k*np.log(np.array(z1) / N) + U_1 / T1 + N*k
S2 = N*k*np.log(np.array(z2) / N) + U_2 / T2 + N*k

plt.plot(T1,S1)
plt.xlabel("T")
plt.ylabel("S")
plt.title("S vs T (low T)")
plt.grid()
plt.show()

plt.plot(T2,S2)
plt.xlabel("T")
plt.ylabel("S")
plt.grid()
plt.title("S vs T (high T)")
plt.show()


#HELMHOLTZ
F1 = -N*k*np.array(T1) * np.log(np.array(z1))
F2 = -N*k*np.array(T2) * np.log(np.array(z2))

plt.plot(T1,F1)
plt.xlabel("T")
plt.ylabel("F")
plt.title("F vs T (low T)")
plt.grid()
plt.show()

plt.plot(T2,F2)
plt.xlabel("T")
plt.ylabel("F")
plt.grid()
plt.title("F vs T (high T)")
plt.show()






import numpy as np
import matplotlib.pyplot as plt

def partition_function(T,epsilon,g):
    Z_list = []
    for i in T:
        Z = 0
        for j,m in zip(epsilon,g):
            Z = Z + m * np.exp(-j/(k*i))
        Z_list.append(Z)
    return np.array(Z_list)

def fraction_population(g,T,epsilon,Z):
    frac_pop_list = []
    for j in range(len(epsilon)):                                                                                             
        frac_pop = (g[j] * np.exp(-epsilon[j]/(k*T)))/Z
        frac_pop_list.append(frac_pop)
    frac_pop_list = np.array(frac_pop_list)
    return frac_pop_list

def internal_energy(frac_pop,N,epsilon,T):
    N_j = N * frac_pop
    inte_energy = np.zeros([len(T)])
    for i in range(len(N_j)):
        inte_energy = inte_energy + N_j[i]*epsilon[i]
    return inte_energy

def entropy(Z,N,T,U):
    S = (N*k*np.log(Z/N)) + (U/T) + (N*k)
    return S

def free_energy(N,T,Z):
    F = -N*k*T*np.log(Z)
    return F

def graph(x1,x2,y1,y2,title,y_label,frac_pop_low,frac_pop_high,key):
    fig1,ax1 = plt.subplots(1,2)
    fig1.suptitle(title)
    if key == 0:
        ax1[0].scatter(x1,y1,label = "Low Temperature")
        ax1[0].set_xlabel("T")
        ax1[0].set_ylabel(y_label)
        ax1[0].grid(ls = "--")
        ax1[0].legend()
        ax1[1].scatter(x2,y2,label = "High Temperature")
        ax1[1].set_xlabel("T")
        ax1[1].set_ylabel(y_label)
        ax1[1].grid(ls = "--")
        ax1[1].legend()
        plt.show()
    elif key == 1:
        for i in range(len(frac_pop_low)):
            ax1[0].scatter(x1,frac_pop_low[i],label = "LowTemperature")
            ax1[1].scatter(x2,frac_pop_high[i],label = "High Temperature")
        ax1[0].set_xlabel("T")
        ax1[0].set_ylabel("$\\dfrac{N_i}{N}$")
        ax1[1].set_xlabel("T")
        ax1[1].set_ylabel("$\\dfrac{N_i}{N}$")
        ax1[0].grid(ls = "--")
        ax1[1].grid(ls = "--")
        ax1[0].legend()
        ax1[1].legend()
        plt.show()
                    
if __name__ == "__main__":
    k = 8.617 * 10**(-5)
    
    T_low = np.linspace(10**(-18),5000,50)
    T_high = np.linspace(5000,10**(5),50)

    g = [1,1,1] ; epsilon = [0,1,2]

    Z_1 = partition_function(T_low,epsilon,g)
    Z_2 = partition_function(T_high,epsilon,g)
    
    frac_pop_low = fraction_population(g,T_low,epsilon,Z_1)
    frac_pop_high = fraction_population(g,T_high,epsilon,Z_2)
    
    U_low = internal_energy(frac_pop_low,1,epsilon,T_low)
    U_high = internal_energy(frac_pop_high,1,epsilon,T_high)
    
    S_low = entropy(Z_1,1,T_low,U_low)
    S_high = entropy(Z_2,1,T_high,U_high)
    
    F_low = free_energy(1,T_low,Z_1)
    F_high = free_energy(1,T_high,Z_2)
    
    graph(T_low,T_high,Z_1,Z_2,"Partition Function","Z",None,None,0)
    graph(T_low,T_high,None,None,"Fraction Population",None,frac_pop_low,frac_pop_high,1)
    graph(T_low,T_high,U_low,U_high,"Internal Energy","U",None,None,0)
    graph(T_low,T_high,S_low,S_high,"Entropy","S",None,None,0)
    graph(T_low,T_high,F_low,F_high,"Helmholtz free energy","F",None,None,0)
