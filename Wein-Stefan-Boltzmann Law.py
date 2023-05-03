import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as int1
from scipy.integrate import simps
from scipy.stats import linregress
import pandas as pd
import scipy

def planck(x):
    return x**3/(np.exp(x)-1)


'''Wein's Displacement'''

if __name__ == "__main__":
    x = np.linspace(0.001,12,10000)
    h = 6.62 * (10 ** (-34))
    c= 3 * (10**(8))
    k = 1.380649 * (10 ** (-23))
    T = 5500

    y_p = planck(x)
    area = simps(y_p,x)
    print("Area = ", area)
    

    y_p_max =  np.max(y_p)
    x_p = 0

    for i in range(len(y_p)):
        if y_p[i] == y_p_max:
           x_p = x[i]
           print("x_p = ",x_p)

    b= (h*c)/(k*x_p)
    print("b for x_p = ",b)

    wavelength = b/T
    print("wavelength = ", wavelength)
   

    plt.scatter(x,y_p)
   # plt.plot(x, y_p, color = "green")
    plt.scatter(x_p,y_p_max)
    plt.grid()
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("F_p(x)")
    plt.show()


if __name__ == "__main__":
    x = np.linspace(0.000001,12,10000)
    h = 6.67 * (10 ** (-34))
    c= 3 * (10**(8))
    k = 1.380649 * (10 ** (-23))
    T = 5500
    
    y=planck(x)
    for i in range(0,len(x)):
        I1,error=int1.quad(planck,x[0],x[i])
        I2,error=int1.quad(planck,x[i],np.inf)
        if round(I1,2)==round(I2,2):
            xm=x[i]
            ym=y[i]

    print("x_m = ",xm)
    print("b for x_m = ", h*c/(k*xm))

    wavelength = (h*c/(k*xm))/T
    print("wavelength = ", wavelength)
    
    plt.scatter(x, planck(x))
  #  plt.plot(x, planck(x), color = "green") 
    plt.scatter(xm,planck(xm))
    plt.grid()
    plt.legend()
    plt.xlabel("x")
    plt.ylabel("F_p(x)")
    plt.show()


T_s = np.arange(1000, 20000, 1000)
dict = {"Temp (K)": T_s,
        "Lambda Peak (nm)": np.round(10**9*b/T_s, 1),
        "Lambda Median (nm)": np.round(10**9*h*c/(k*xm)/T_s, 1)}
df = pd.DataFrame(dict)
print(df)


'''Stefan's Boltzmann Law'''

T = np.arange(100, 10000, 500)
CT = 8*np.pi*(k*T)**4/(h*c)**3
FT = (c/4)*(np.pi**4/15)*CT

plt.scatter(T, FT)
plt.plot(T, FT, color = "green")
plt.xlabel("T (K)")
plt.ylabel("F(T)")
plt.grid()
plt.legend()
plt.show()


lnF = np.log(FT)
lnT = np.log(T)
fit_res = scipy.stats.linregress(lnT, lnF)

new_lnT = np.linspace(0, lnT[-1], 100)
calc_lnF = fit_res[0]*new_lnT + fit_res[1]


print("\nStefan-Boltzmann Constant\nSlope = {:.5}\nIntercept = {:.5}".format(fit_res[0], fit_res[1]))
sigma = np.exp(fit_res[1])
sigma_org = 2*np.pi**5*k**4/(15*h**3*c**2)
print("Calculated Stefan's Constant = {:.5}".format(sigma))
print("Stefan's Constant = {:.5}".format(sigma_org))

plt.scatter(lnT, lnF)
plt.plot(new_lnT, calc_lnF, color = "green")
plt.xlabel("ln(T)")
plt.ylabel("ln(F(T))")
plt.grid()
plt.legend()
plt.show()








import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.integrate import quad
from scipy.stats import linregress

def planck(x):
    return x**3/(np.exp(x) - 1)

def median_point(tol):
    area = simps(u,x)
    area_h = area/2
    for i in range(2,len(x)):
        area1 = simps((u[0:i]),(x[0:i]))
        
        if np.abs(area1 - area_h) <= tol:
            break
    
    return x[i]

def u_d(T):
    return (np.pi**4/15)*(8*np.pi*(k*T)**4/(h*c)**3)       

if __name__ == "__main__":
    k = 1.38 * 10**(-23)
    h = 6.624 * 10**(-34)
    c = 3 * 10**(8)
    
    x_in = 1e-2 ; x_fin = 12
    x = np.linspace(x_in,x_fin,5000)
    
    u = planck(x)
    
    median = median_point(1e-3)
    
    print('\nx_median:',median)

    # WEIN'S CONSTANT

    b = h*c/(k*median)
    print("\nWein's constant:",b)
    
    # STEFAN
    
    I_p = quad(planck,1e-15,20)
    
    print("\nValue of I_p is:",I_p[0],"and Standard value of I_p is:",np.pi**4/15)
    
    T = np.arange(100,10200,500)
    
    F = lambda T : (c/4)*u_d(T)
    
    res = linregress(np.log(T),np.log(F(T)))

    print('\nSlope:',res[0])
    print('\nIntercept:',res[1])

    print('\nStefan constant:',np.exp(res[1]),"\n")

    fig,ax = plt.subplots()
    ax.plot(x,u)
    ax.set_title('Planck law of radiation')
    ax.set_xlabel('x')
    ax.set_ylabel('Energy spectral density')
    plt.grid(ls = "--")
    plt.show()

    fig,ax = plt.subplots()
    ax.plot(T,F(T))
    ax.scatter(T,F(T))
    ax.set_title('Radiant flux (F) Vs Temperature (T)')
    ax.set_xlabel('Temperature (T)')
    ax.set_ylabel('Radiant flux (F)')
    plt.grid(ls = "--")
    plt.show()

    fig,ax = plt.subplots()
    ax.plot(np.log(T),np.log(F(T)))
    ax.scatter(np.log(T),np.log(F(T)))
    ax.set_title("Linear Regression")
    ax.set_xlabel('log(T)')
    ax.set_ylabel('log(F)')
    plt.grid(ls = "--")
    plt.show()
