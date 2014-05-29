# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.interpolate import interp1d as interp

# <codecell>

bkgd = np.array(numpy.loadtxt(open("bkgd.txt","r"),skiprows=17))
bkgd_counts = bkgd[:,1]
nm = bkgd[:,0]
eV = 1240/nm
TotalIntensity = 0.025 # in mW/cm^2 room lights are = 25 uW/cm^2

# <codecell>

LED = np.array(numpy.loadtxt(open("LED.txt","r"),skiprows=17))
LED_counts = LED[:,1] - bkgd_counts
LED_I = LED_counts*eV
LED_I = TotalIntensity/numpy.trapz(LED_I, nm)*LED_I

# <codecell>

CFL = np.array(numpy.loadtxt(open("CFL.txt","r"),skiprows=17))
CFL_counts = CFL[:,1] - bkgd_counts
CFL_I = CFL_counts*eV
CFL_I = TotalIntensity/numpy.trapz(CFL_I, nm)*CFL_I

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

plt.plot(nm, LED_I, nm, CFL_I, sun_nm, sun_GT/1000)
plt.ylabel('Spectral Intensity (mW*cm-2*nm-1)')
plt.xlabel('Wavelength (nm)')
plt.legend(('LED', 'CFL', 'Sun'))
plt.xlim(400, 800)

# <codecell>

def Jsc(Eg, nm, spectrum):
    integrand = spectrum/(1240/nm)
    if Eg == 0: Eg = 0.31
    Eg_index = numpy.searchsorted(nm, 1240/Eg)
    return numpy.trapz(integrand[0:Eg_index], nm[0:Eg_index])

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

def JV(V, Eg, nm, spectrum, Rs, Rsh, ERE):  # ERE is in fractions, not %
    Jph = Jsc(Eg, nm, spectrum)
    J0 = planck(Eg)
    if Rs == 0:
        return [Jph-J0*(exp(x/kT) - 1)/ERE for x in V]
    else:
        return [brentq(JV_Root, 70, -1E10, args=(x, J0, Jph, Eg, Rs, Rsh, ERE)) for x in V]

# <codecell>

V = np.arange(-1, 3, 0.001)
J = JV(V, 1.1, nm, LED_I, 0, Inf, 1)

# <codecell>

f = interp(J, V, bounds_error = False)
f(-60)

# <codecell>

plt.plot(V, J)
xlabel('Voltage (V)')
ylabel('Current (mA/cm$^2$)')
ylim(-10,60)

# <codecell>

ERE = 1E-4 # this is in fractions, not %
Rs = 0
Rsh = Inf
bandgaps = linspace(0.3, 3, 50)
PCE_LED = [max(V*JV(V, Eg, nm, LED_I, Rs, Rsh, ERE))*100/TotalIntensity for Eg in bandgaps]
PCE_CFL = [max(V*JV(V, Eg, nm, CFL_I, Rs, Rsh, ERE))*100/TotalIntensity for Eg in bandgaps]
PCE_Sun = [max(V*JV(V, Eg, sun_nm, sun_GT, Rs, Rsh, ERE))*100/100 for Eg in bandgaps]

# <codecell>

plt.plot(bandgaps, PCE_LED, bandgaps, PCE_CFL, bandgaps, PCE_Sun)
xlabel('Bandgap (eV)')
ylabel('PCE (%)')
matplotlib.pyplot.legend(('LED','CFL', 'Sun'))

# <codecell>

max(PCE_Sun)

# <codecell>


# <codecell>


