import numpy as np
import matplotlib.pyplot as plt

# Maxwell Boltzmann Distribution
def boltz(x):
    return np.e**(-x)

# Bose-Einstein Distribution
def bose(x, alpha):
    X = []
    for i in x:
        if alpha> -i:
            X.append(i)
    X = np.array(X)
    return X, 1/(np.e**(X + alpha) - 1)

# Fermi-Dirac Distribution
def fermi(x, alpha):
    return 1/(np.e**(x + alpha) + 1)

x = np.linspace(-4, 4, 500) # E/kT

'''plot-1'''
fig, ax = plt.subplots(nrows = 1, ncols = 3)
ax[0].plot(x, boltz(x))
X, Y = bose(x, 0)
ax[1].plot(X, Y, linestyle = 'dashdot')
ax[2].plot(x, fermi(x, 0), linestyle = '--')
ax[0].set_xlabel("$\epsilon/kT$")
ax[1].set_xlabel("$\epsilon/kT$")
ax[2].set_xlabel("$\epsilon/kT$")
ax[0].set_ylabel("Occupancy")
ax[1].set_ylabel("Occupancy")
ax[2].set_ylabel("Occupancy")
ax[0].set_title("Maxwell-Boltzmann")
ax[1].set_title("Bose-Einstein")
ax[2].set_title("Fermi-Dirac")
ax[0].grid()
ax[1].grid()
ax[2].grid()
plt.title('Comparing Distribution Functions')

plt.figure("AD_Sub")
plt.plot(x,boltz(x), label = 'Maxwell-Boltzmann')
X, Y = bose(x, 0)
plt.plot(X, Y, linestyle = 'dashdot', label='Bose-Einstein')
plt.plot(x, fermi(x, 0), linestyle = '--', label='Fermi-Dirac')
plt.xlabel("ε/kT")
plt.ylabel("Occupancy")
plt.legend()
plt.grid()
plt.title('Comparing Distribution Functions')

fig, axn = plt.subplots(1, 3)
'''plot-2'''
T = np.array([10, 100, 1000, 5000])
for i in T:
    axn[0].plot(x*i*(8.62*10**(-5)), fermi(x, 0), label='T= '+str(i)+' K')
axn[0].legend()
axn[0].set_xlabel('Energy (in eV)')
axn[0].set_ylabel('Occupancy')
axn[0].set_title('Fermi-Dirac')
axn[0].grid()

'''plot-3'''
T = np.array([10, 100, 1000, 5000])
for i in T:
    axn[1].plot(bose(x, 0)[0]*i*(8.62*10**(-5)), bose(x, 0)[1], label='T= '+str(i)+' K')
axn[1].legend()
axn[1].set_xlabel('Energy')
axn[1].set_ylabel('Occupancy')
axn[1].set_title('Bose-Einstein')
axn[1].grid()

'''plot-4'''
T = np.array([10, 100, 1000, 5000])
for i in T:
    axn[2].plot(x*i*(8.62*10**(-5)), boltz(x), label='T= '+str(i)+' K')
axn[2].legend()
axn[2].set_xlabel('Energy')
axn[2].set_ylabel('Occupancy')
axn[2].set_title('Maxwell-Boltzmann')
axn[2].grid()

fig.suptitle("Varying Energy")

xT = np.linspace(0, 500000, 1000)
xE = np.array([0.1, 0.5, 1, 5]) # in eV

fig, ax = plt.subplots(1, 3)
'''plot-5'''
for i in xE:
    E = i/(8.62*10**(-5)*xT)
    ax[0].plot(xT, fermi(E, 0), label = 'E = '+str(i)+' eV')
ax[0].legend()
ax[0].set_xlabel('T (in K)')
ax[0].set_ylabel('Occupancy')
ax[0].set_title('Fermi-Dirac')
ax[0].grid()

'''plot-6'''
for i in xE:
    E = i/(8.62*10**(-5)*xT)
    ax[1].plot(xT, boltz(E), label = 'E = '+str(i)+' eV')
ax[1].legend()
ax[1].set_xlabel('T (in K)')
ax[1].set_ylabel('Occupancy')
ax[1].set_title('Maxwell-Boltzmann')
ax[1].grid()

'''plot-7'''
xT2 = np.linspace(0, 5000, 500)
for i in xE:
    E = i/(8.62*10**(-5)*xT)
    ax[2].plot(xT, bose(E, 0)[1], label='E= '+str(i)+' eV')
ax[2].legend()
ax[2].set_xlabel('T (in K)')
ax[2].set_ylabel('Occupancy')
ax[2].set_title('Bose-Einstein')
ax[2].grid()
fig.suptitle("Varying Temperature")

plt.show()









import numpy as np
import matplotlib.pyplot as plt

def max_bolt(x):
    return np.exp(-x)

def bose_einstein(x,alpha):
    return 1/(np.exp(x+alpha) - 1)

def fermi_dirac(x,alpha):
    return 1/(np.exp(x+alpha) + 1)


if __name__ == "__main__":
    x_range_max = np.linspace(-4,4,50)
    x_range_bose = np.linspace(0.1,4,50)
    x_range_fermi = np.linspace(-4,4,50)
    alpha = [0,1]
    k = 8.617333 * 10**(-5)
    T = np.array([10,100,1000,5000])
    fermi_total_fermi = [] ; fermi_total_bose = [] ; fermi_total_max = []
    
    f_max_bolt = max_bolt(x_range_max)
    f_bose_einstein = bose_einstein(x_range_bose,alpha[0])
    f_fermi_dirac = fermi_dirac(x_range_fermi,alpha[1])

    fig, ax1 = plt.subplots(1, 3, figsize=(10, 4))

    ax1[0].plot(x_range_max,f_max_bolt)
    ax1[0].scatter(x_range_max,f_max_bolt,marker = ".")
    ax1[1].plot(x_range_bose,f_bose_einstein,c = "g")
    ax1[1].scatter(x_range_bose,f_bose_einstein,c = "g",marker = ".")
    ax1[2].plot(x_range_fermi,f_fermi_dirac,c = "r")
    ax1[2].scatter(x_range_fermi,f_fermi_dirac,c = "r",marker = ".")
    for i in range(3):
        ax1[i].set(xlabel = "$\epsilon$/KT",ylabel = "f($\epsilon$)")
        ax1[i].grid(ls = "--")
    ax1[0].set(title = "Maxwell Boltzman Distribution")
    ax1[1].set(title = "Bose Einstein Distribution")
    ax1[2].set(title = "Fermi Dirac Distribution")
    plt.show()

    for i in range(len(T)):
        fermi_x = x_range_fermi * T[i]*k
        fermi_total_fermi.append(fermi_x)

    for i in range(len(T)):
        fermi_x = x_range_max * T[i]*k
        fermi_total_max.append(fermi_x)
        
    for i in range(len(T)):
        fermi_x = x_range_bose * T[i]*k
        fermi_total_bose.append(fermi_x)

    
    fig2,ax2 = plt.subplots()
    fig3,ax3 = plt.subplots()
    fig4,ax4 = plt.subplots()

    for i in range(len(fermi_total_fermi)):
        ax2.plot(fermi_total_fermi[i],f_fermi_dirac,label = "At T = "+str(T[i])+" K")
        ax2.scatter(fermi_total_fermi[i],f_fermi_dirac,marker = ".")
        ax3.plot(fermi_total_max[i],f_max_bolt,label = "At T = "+str(T[i])+" K")
        ax3.scatter(fermi_total_max[i],f_max_bolt,marker = ".")
        ax4.plot(fermi_total_bose[i],f_bose_einstein,label = "At T = "+str(T[i])+" K")
        ax4.scatter(fermi_total_bose[i],f_bose_einstein,marker = ".")
    ax2.set(xlabel = "$\epsilon$",ylabel = "f($\epsilon$)",title = "Fermi Dirac distribution at constant temperature")
    ax3.set(xlabel = "$\epsilon$",ylabel = "f($\epsilon$)",title = "Maxwell Boltzman distribution at constant temperature")
    ax4.set(xlabel = "$\epsilon$",ylabel = "f($\epsilon$)",title = "Bose Einstein distribution at constant temperature")
    ax2.grid(ls = "--")
    ax3.grid(ls = "--")
    ax4.grid(ls = "--")
    ax2.legend()
    ax3.legend()
    ax4.legend()
    plt.show()



