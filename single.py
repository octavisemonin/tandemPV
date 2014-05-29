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

%timeit planck(1.1)

# <codecell>

def JV_Root(J, V, J0, Jph, Eg, Rs, Rsh, ERE):
    A = 1
    xbeta = 1
    return Jph - J0*(exp((V + J*A*Rs/1000)/kT/xbeta) -1) - (V + J*A*Rs/1000)/(A*Rsh)*1000 - J

# <codecell>

def JV(V, Eg, Rs, Rsh, ERE): 
    Jph = Jsc(Eg)
    J0 = planck(Eg)
    if Rs == 0:
        return [Jph-J0*(exp(x/kT) - 1)/ERE for x in V]
    else:
        return [brentq(JV_Root, 70, -1E10, args=(x, J0, Jph, Eg, Rs, Rsh, ERE)) for x in V]

# <codecell>

V = np.arange(-1, 3, 0.001)
%prun J = JV(V, 1.1, 1, Inf, 1)

# <codecell>

f = interp(J, V, bounds_error = False)
f(-60)

# <codecell>

plt.plot(V, J)
xlabel('Voltage (V)')
ylabel('Current (mA/cm$^2$)')
ylim(-10,60)

# <codecell>

ERE = 1
Rs = 0
Rsh = Inf
bandgaps = linspace(0.3, 3, 50)
PCE = [max(V*JV(V, Eg, Rs, Rsh, ERE)) for Eg in bandgaps]

# <codecell>

plt.plot(bandgaps, PCE)
xlabel('Bandgap (eV)')
ylabel('PCE (%)')
matplotlib.pyplot.legend('ERE = 1')

# <codecell>

Rs = 1
PCE_1 = [max(V*JV(V, Eg, Rs, Rsh, ERE)) for Eg in bandgaps]
Rs = 2
PCE_2 = [max(V*JV(V, Eg, Rs, Rsh, ERE)) for Eg in bandgaps]

# <codecell>

Rs = 4
PCE_4 = [max(V*JV(V, Eg, Rs, Rsh, ERE)) for Eg in bandgaps]

# <codecell>

Rs = 8
PCE_8 = [max(V*JV(V, Eg, Rs, Rsh, ERE)) for Eg in bandgaps]

# <codecell>

plt.plot(bandgaps, PCE, bandgaps, PCE_1, bandgaps, PCE_2, bandgaps, PCE_4, bandgaps, PCE_8)
xlabel('Bandgap (eV)')
ylabel('PCE (%)')
matplotlib.pyplot.legend(('ERE = 1', 'ERE = 10$^{-1}$', 'ERE = 10$^{-2}$', 'ERE = 10$^{-4}$', 'ERE = 10$^{-8}$'))

# <codecell>

max(PCE)

# <codecell>


# <codecell>


