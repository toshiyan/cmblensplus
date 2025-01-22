import numpy as np
from astropy import units as u
from scipy import constants as const


#Planck constant
h = const.h #m^2kg/s

#reduced Planck mass [GeV]
mpl = 2.435e18

# Speed of light
c = const.c #m/s
C = const.c*1e-3 #km/s

# Boltzman constant
k = const.k #m^2kg/s^2/K

# proton mass in kg
m_p = const.m_p

# Thomson cross-section for electron in m^2
sigmaT = const.value('Thomson cross section')  # [m^2]

# critical energy density in cosmology in h^2 kg/m^3
rho_c = 1.8788e-26

# Mpc/m
Mpc2m = 3.08567757e22

# CMB temperature
Tcmb  = 2.7255e6 #K

# eV -> Mpc
eV2Mpc = 806554.815354305*3.0856775814914e22*(2*np.pi)


# CMB 
#mtype = ['T','E','B']


# //// Unit Conversion //// #

# angle unit conversion
ac2rad  = np.pi/10800.
rad2ac  = 1./ac2rad
deg2rad = np.pi/180.
rad2deg = 1./deg2rad

# J/sr -> uk depending on frequency in GHz
def Jysr2uK(nu):
    return ( ( 1.*u.Jy/u.sr ).to( u.uK,equivalencies=u.thermodynamic_temperature(nu*u.GHz,Tcmb*1e-6*u.K)) ).value

def MJysr2uK(nu):
    return ( ( 1e6*u.Jy/u.sr ).to( u.uK,equivalencies=u.thermodynamic_temperature(nu*u.GHz,Tcmb*1e-6*u.K)) ).value

