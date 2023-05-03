import matplotlib.pyplot as plt 
import numpy as np
import pandas as pd 
import scipy.constants as const
import scipy
from scipy.integrate import quad 

h = const.h
m = const.m_e
k = const.k

def integrand(x):
 return pow(x,2)*np.exp(-pow(h*x,2)/(8*6.6464731*pow(10,-27)*pow(V[i],2/3)*k*T[j]))

def Z(V,T):
 result, error = quad(integrand, 0, pow(10,11))
 return (np.pi/2)*result 

V = np.linspace(2*pow(10,-4),5*pow(10,-4),50)
T = np.linspace(150,450,50)

Z_matrix = []
row = []
for i in range(0,len(V)):
 for j in range(0,len(T)):
  row.append(round(Z(V[i],T[j]),4))
 Z_matrix.append(row)
 row = []

for i in range(len(T)):
 Z_matrix[i] = np.array(Z_matrix[i])
Z_matrix = np.array(Z_matrix) # V constant for each row of variable T

Z_matrix_transpose = Z_matrix.transpose() # T constant for each row of variable V

lnZ_matrix = []
lnZ_matrix_transpose = []
for i in range(0,len(V)):
 lnZ_matrix.append(np.log(Z_matrix[i]))
 lnZ_matrix_transpose.append(np.log(Z_matrix_transpose[i]))

lnZ_matrix = np.array(lnZ_matrix)
lnZ_matrix_transpose = np.array(lnZ_matrix_transpose)
f1 = plt.figure(1)
for i in range(0,len(T)):
 if i%10==0:
  plt.plot(np.log(T),lnZ_matrix[i],label=f'{round(V[i]*pow(10,4),2)}$*10^{-4} m^3$',marker='.')
plt.xlabel("ln(T)")
plt.legend()
plt.grid(True)
plt.title("lnZ v T (for different V)")
plt.ylabel("lnZ")
f1.show()

f2 = plt.figure(2)
for i in range(0,len(T)):
 if i%10==0:
  plt.plot(np.log(V),lnZ_matrix_transpose[i],label=f'{round(T[i],2)} K',marker='.')
plt.xlabel("ln(V)")
plt.ylabel("lnZ")
plt.legend()
plt.title("lnZ v V (for different T)")
plt.grid(True)
f2.show()

# lnV = np.log(V)
# lnT = np.log(T)

# f3 = plt.figure(3)
# plt.plot(T,lnT,marker='.')
# plt.xlabel("T")
# plt.ylabel("lnT")
# plt.grid(True)
# plt.title("lnT v T")
# f3.show()

# f4 = plt.figure(4)
# plt.plot(V,lnV,marker='.')
# plt.xlabel("V")
# plt.ylabel("lnV")
# plt.grid(True)
# plt.title("lnV v V")
# f4.show()

P_matrix = []

def for_der(x,y):
 del_x = np.array(x)
 del_y = np.array(y)
 for i in range(1,len(x)):
  np.append(del_x, x[i]-x[i-1])
  np.append(del_y, y[i]-y[i-1])
 arr = del_y/del_x
 return arr

P_matrix = []
for i in range(0,len(T)):
 row = k*T[i]*for_der(V,lnZ_matrix_transpose[i])
 P_matrix.append(row)

P_matrix = np.array(P_matrix)
P_matrix_transpose = np.transpose(P_matrix)
f5 = plt.figure(5)
for i in range(0,len(T)):
 if i%10==0:
  plt.plot(V,P_matrix[i],marker='.',label=f'{round(T[i],2)} K')
plt.xlabel("V")
plt.grid(True)
plt.title("pressure v V (for different T)")
plt.legend()
plt.ylabel("pressure")
f5.show()

f6 = plt.figure(6)
for i in range(0,len(T)):
 if i%10==0:
  plt.plot(T,P_matrix_transpose[i],marker='.',label=f'{round(V[i]*pow(10,4),2)}$*10^{-4} m^3$')
plt.xlabel("T")
plt.ylabel("pressure")
plt.title("pressure v T (for different V)")
plt.legend()
plt.grid(True)
f6.show()

E_matrix = []
row = []
for i in range(0,len(T)):
 row.append(k*pow(T,2)*for_der(T,lnZ_matrix[i]))
 E_matrix.append(row)
E_matrix = np.array(E_matrix)

f7 = plt.figure(7)
for i in range(0,len(V)):
 if i%10==0:
  plt.plot(T,E_matrix[i][0],marker='.',label=f'{round(V[i]*pow(10,4),2)}$*10^{-4} m^3$')
plt.xlabel("T")
plt.ylabel("<E>")
plt.title("<E> vs T")
plt.legend()
plt.grid(True)
f7.show()

# TAKING N = N_A 
N = const.N_A
f8 = plt.figure(8)
plt.plot(T,N*E_matrix[0][0],marker='.',label=f'{round(V[0]*pow(10,4),2)}$*10^{-4} m^3$')
plt.xlabel("T")
plt.ylabel("U")
plt.title("U vs T")
plt.legend()
plt.grid(True)
f8.show()

result = scipy.stats.linregress(T,E_matrix[0][0])
Cv = result[0]
print("Slope : C_v =", result[0])

S = []
row = []
for i in range(0,len(V)):
 S.append((E_matrix[i][0]/T) + N*k*(lnZ_matrix[i] - np.log(N) + 1))
S = np.array(S)
S_transpose = np.transpose(S)
f9 = plt.figure(9)
for i in range(0,len(V)):
 if i%10==0:
  plt.plot(T,S[i],marker='.',label=f'{round(V[0]*pow(10,4),2)}$*10^{-4} m^3$')
plt.xlabel("T")
plt.ylabel("S")
plt.grid(True)
plt.legend()
plt.title("S vs T from different V")
f9.show()

f10 = plt.figure(10)
for i in range(0,len(V)):
 if i%10==0:
  plt.plot(V,S_transpose[i],marker='.',label=f'{round(T[i],2)} K')
plt.xlabel("V")
plt.ylabel("S")
plt.grid(True)
plt.legend()
plt.title("S vs V from different T")
f10.show()

Energy_fluctuations = k*pow(T,2)*Cv
f11 = plt.figure(11)
plt.plot(T,Energy_fluctuations,marker='.')
plt.xlabel("T")
plt.ylabel(r"$<(\Delta E)^2>$")
plt.title(r"$<(\Delta E)^2>$ vs $T$")
plt.grid(True)

plt.show()






import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
from numpy import diff
from scipy.stats import linregress

def partition_function(V,T):

    def integral(n):
        const = h**2/(8*m*(V**(2/3))*k_b*T)
        return n**(2)* np.exp(-n**2 * const)

    return np.pi/2 * quad(integral,0,10**(11))[0]

def result_matrix(V_array,T_array):
    result = np.zeros((len(V_array),len(T_array)))
    result_log = np.zeros((len(V_array),len(T_array)))

    for i in range(len(V_array)):
        for j in range(len(T_array)):
            result[i][j] = partition_function(V_array[i],T_array[j])
            result_log[i][j] = np.log(partition_function(V_array[i],T_array[j]))
    return result.T,result_log.T

def pressure(V,T,Z):
    Z_diff = diff(Z)
    V_diff = diff(V)
    diff_term = Z_diff/V_diff
    p = N_a*k_b*T*diff_term
    return p

def energy(Z,T):
    Z_diff = diff(Z)
    T_diff = diff(T)
    diff_term = Z_diff/T_diff
    E = k_b * T[:-1]**2 * diff_term
    return E,E*N_a

def entropy(U,T,Z):
    term = N_a * k_b * (Z - np.log(N_a) + 1)
    en = U/T + term
    return en

def energy_fluctation(T,C_v):
    e_fluc = k_b * T**2 * C_v
    return e_fluc
    
if __name__ == "__main__":
    k_b = 1.38 * 10**(-23)
    N_a = 6.022 * 10**(23)
    m = 3.32 * 10**(-27)
    h = 6.63 * 10**(-34)

    V_array = np.linspace(20 * 10**(-3),50 * 10**(-3),10)
    T_array = np.linspace(150,450,10)
    p_array = [] ; e_array = [] ; u_array = [] ; en_array = []

    result,result_log = result_matrix(V_array,T_array)

    for i in range(len(V_array)):
        plt.plot(V_array,result_log[i,:],label = "At T = "+str(T_array[i]))
    plt.xlabel("V")
    plt.ylabel("log(Z)")
    plt.title("Partition Function (at constant temperature)")
    plt.grid(ls = "--")
    #plt.legend()
    plt.show()

    for i in range(len(V_array)):
        plt.plot(T_array,result_log[:,i],label = "At V = "+str(V_array[i]))
    plt.xlabel("T")
    plt.ylabel("log(Z)")
    plt.title("Partition Function (at constant volume)")
    plt.grid(ls = "--")
    #plt.legend()
    plt.show()

    # PRESSURE

    for i in range(len(T_array)):
        p = pressure(V_array,T_array[i],result_log[i,:])
        p_array.append(p)
        plt.plot(V_array[:-1],p,label = "At T = "+str(T_array[i]))
    plt.xlabel("Volume (meter^3)")
    plt.ylabel("Pressure (Pa)")
    plt.title("Pressure Vs Volume (At Constant Temperature)")
    plt.grid(ls = "--")
    #plt.legend()
    plt.show()
        
    p_array = np.array(p_array)
    
    for i in range(len(V_array[:-1])):
        plt.plot(T_array[:-1],p_array[:,i][:-1],label = "At V = "+str(V_array[i]))
    plt.xlabel("Temperature (K)")
    plt.ylabel("Pressure (Pa)")
    plt.title("Pressure Vs Tempertaure (At Constant Volume)")
    plt.grid(ls = "--")
    #plt.legend()
    plt.show()
        
    # ENERGY
  
    for i in range(len(T_array)):
        E,U = energy(result_log[:,i],T_array)
        e_array.append(E)
        u_array.append(U)
        #print(len(U))
        plt.plot(T_array[:-1],E)
    plt.xlabel("Temperature")
    plt.ylabel("Energy")
    plt.title("Energy Vs temperature (At constant Volume)")
    plt.grid(ls = "--")
    plt.show()

    slope = linregress(T_array[:-1],e_array[0])[0]
    print("\nCv :",slope)
        
    u_array = np.array(u_array)

    # ENTROPY

    for i in range(len(T_array[:-1])):
        en = entropy(u_array[i,:],T_array[:-1],result_log[:,i][:-1])
        en_array.append(en)
        plt.plot(T_array[:-1],en,label = "At V = "+str(V_array[i]))
    plt.xlabel("Temperature")
    plt.ylabel("Entropy")
    plt.title("Entropy Vs temperature (At constant Volume)")
    plt.grid(ls = "--")
    #plt.legend()
    plt.show()

    # ENERGY FLUCTATION

    e_fluc = energy_fluctation(T_array,slope)
    print("\nEnergy Fluctation :",e_fluc)

    plt.plot(T_array,e_fluc)
    plt.xlabel("Temperature")
    plt.ylabel("Variance")
    plt.title("Variance Vs temperature")
    plt.grid(ls = "--")
    plt.show()
