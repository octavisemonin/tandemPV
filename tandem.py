# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.interpolate import interp1d as interp

# <codecell>

sun = np.array(numpy.loadtxt(open("ASTMG173.csv","r"),delimiter=",",skiprows=2))

# <rawcell>

# Convert from W/m^2/nm to mW/cm^2/nm

# <codecell>

sun_nm = sun[:,0]
sun_space = sun[:,1]*1000/100**2
sun_GT = sun[:,2]*1000/100**2
sun_DC = sun[:,3]*1000/100**2

# <codecell>

plt.plot(sun_nm, sun_space, sun_nm, sun_GT, sun_nm, sun_DC)
plt.ylabel('Spectral Intensity (mW*cm-2*nm-1)')
plt.xlabel('Wavelength (nm)')
plt.legend(('AM0', 'AM1.5 Global Tilt', 'AM1.5 Direct and Circumsolar'))

# <codecell>

sun_eV = 1240/sun_nm

# <codecell>

def Jsc(Eg):
    integrand = sun_GT/sun_eV
    if Eg == 0: Eg = 0.31
    Eg_index = numpy.searchsorted(sun_nm, 1240/Eg)
    return numpy.trapz(integrand[0:Eg_index], sun_nm[0:Eg_index])

# <codecell>

kT = 0.0256   # eV
h = 4.136E-15 # eV-s
c = 2.998E8   # m/s
q = 1.602E-19 # C

# <codecell>

def planck(Eg):
    plancky = lambda Eph: 2*pi*Eph**2*q/(h**3*c**2)*1/(exp(Eph/kT) - 1)
    return quad(plancky, Eg, 4)[0]*1000/100**2 # in mA/cm^2

# <codecell>

def JV_Root(J, V, Jph, J0, Rs, Rsh, ERE):
    A = 1
    xbeta = 1
    return Jph - J0*(exp((V + J*A*Rs/1000)/kT/xbeta) -1) - (V + J*A*Rs/1000)/(A*Rsh)*1000 - J

# <codecell>

def JV(V, Jph, J0, Rs, Rsh, ERE): 
#    Jph = Jsc(Eg)
#    J0 = planck(Eg)
    if Rs == 0:
        return [Jph-J0*(exp(x/kT) - 1)/ERE for x in V]
    else:
        return [scipy.optimize.brentq(JV_Root, 70, -1E10, args=(x, Jph, J0, Rs, Rsh, ERE)) for x in V]

# <codecell>

def tandem_JV(V, Eg1, Eg2, series):

    if(Eg2 > Eg1):
        return 0
    
    Jph1 = Jsc(Eg1)
    Jph2 = Jsc(Eg2) - Jph1

    if (Jph1 > Jph2) and (series):
        Jph2 = (Jph1+Jph2)/2
        Jph1 = Jph2

    J0_1 = planck(Eg1)
    J0_2 = planck(Eg2)
    
    # sun
    # Eg1
    # Eg2
    # grnd
    
    Rs = 0
    Rsh = 0
    ERE = 1
    
    JV1 = JV(V, Jph1, J0_1, Rs, Rsh, ERE)
    JV2 = JV(V, Jph2, J0_2, Rs, Rsh, ERE)
    
#    plt.plot(V, JV1, V, JV2)
#    ylim(-10,60)

    if series:
        basis = JV1
        target = JV2

        f = interp(sorted(target), sorted(V, reverse=True), bounds_error = False, fill_value = 0)
        V_I = np.array([f(x) for x in basis])
        V_T = [x + y for x, y in zip(V, V_I)]
#        plt.plot(V_T,basis)

        return max([a*b for a,b in zip(V_T,basis)])
    else:
        return max(V*JV1)+max(V*JV2)

# <codecell>

V = np.linspace(-1, 3, 201)
p = tandem_JV(V, 0.75, 0.5, nonseries)
print p

# <codecell>

Eg1 = np.linspace(.3, 3, 136)
Eg2 = Eg1
series = 1
nonseries = 0

# <codecell>

PCE_S = np.array([[tandem_JV(V, x, y, series) for x in Eg1] for y in Eg2])

# <codecell>

PCE_NS = np.array([[tandem_JV(V, x, y, nonseries) for x in Eg1] for y in Eg2])

# <codecell>

PCE_S_1000 = np.array([[tandem_JV(V, x, y, series) for x in Eg1] for y in Eg2])

# <codecell>

PCE_NS_1000 = np.array([[tandem_JV(V, x, y, nonseries) for x in Eg1] for y in Eg2])

# <codecell>

Igor_PCE_NS = np.array(numpy.loadtxt(open("PCE_NS.txt","r"),delimiter="\t"))
Igor_PCE_S = np.array(numpy.loadtxt(open("PCE_ST.txt","r"),delimiter="\t"))

# <codecell>

matplotlib.pyplot.contourf(Eg1, Eg2, PCE_S, levels=arange(0,49,2))
xlabel('Top $E_g$ (eV)')
ylabel('Bottom $E_g$ (eV)')
cb = colorbar()

# <codecell>

matplotlib.pyplot.contourf(Eg1, Eg2, PCE_NS, levels=arange(0,49,2))
xlabel('Top $E_g$ (eV)')
ylabel('Bottom $E_g$ (eV)')
cb = colorbar()

# <codecell>

plt.pcolor(Eg1, Eg2, PCE_NS)
xlabel('Top $E_g$ (eV)')
ylabel('Bottom $E_g$ (eV)')
ylim([0.3, 3])
xlim([0.3, 3])
cb = colorbar()

# <codecell>

matplotlib.pyplot.contourf(Eg1, Eg2, PCE_NS-PCE_S+0.001, levels=linspace(0,35,15))
xlabel('Top $E_g$ (eV)')
ylabel('Bottom $E_g$ (eV)')
cb = colorbar()

# <codecell>

plt.imshow(PCE_NS-PCE_S, interpolation='nearest')
xlabel('Top $E_g$ (eV)')
ylabel('Bottom $E_g$ (eV)')
ylim([0,136])
cb = colorbar()

# <codecell>

print PCE_S.max()
print PCE_NS.max()

# <codecell>

thresholds = arange(34, 45.75, 0.25)
counts_S = np.array([float(sum(PCE_S > x)) for x in thresholds])
counts_NS = np.array([float(sum(PCE_NS > x)) for x in thresholds])
ratio = counts_NS/counts_S

# <codecell>

print ratio[4*(44-34)]

# <codecell>

plot(thresholds, ratio)
pyplot.ylim([0, 10])
ylabel('Phase Space Ratio')
xlabel('PCE Threshold (%)')

# <codecell>

subplot(2,2,1)
plt.pcolor(Eg1, Eg2, PCE_S)
ylim([0.3, 3])
xlim([0.3, 3])
subplot(2,2,2)
matplotlib.pyplot.contourf(PCE_S_1000, levels=arange(0,47,2))
cb = colorbar()
subplot(2,2,3)
matplotlib.pyplot.contourf(PCE_NS, levels=arange(0,47,2))
subplot(2,2,4)
matplotlib.pyplot.contourf(PCE_NS_1000, levels=arange(0,47,2))
cb = colorbar()

# <codecell>


