
import matplotlib.pyplot as plt
import numpy as np
import os, sys

sys.path.append(r'C:\Users\User\Documents\Spectre\2020\fall')
sys.path.append(r'C:\Users\User\Documents\Spectre\2020\fall\CM1')
import utils as ut
import cm_utils as cm_ut

import spiceypy as spice


mu=398600.4415
a_e = 6378.1363 # km
w_e = 7.292115e-5 # rad/s
g = 9.81 # m/s^2
J2 = 1.082e-3

def get_LAN_dot(a, e, inc, J2=1.082e-3, mu=398600.4415, a_e=6378.1363, radians=True):
	# units in km, s
	n = np.sqrt(mu/a**3)

	if not radians:
		inc = np.radians(inc)

	LAN_dot = -3/2 * n * (a_e/a)**2 * J2 * (1/np.sqrt(1-e**2)) * np.cos(inc)
	return LAN_dot

def get_omega_dot(a, e, inc, J2=1.082e-3, mu=398600.4415, a_e=6378.1363, radians=True):
	# units in km, s
	n = np.sqrt(mu/a**3)

	if not radians:
		inc = np.radians(inc)

	omega_dot = -3/4 * n * (a_e/a)**2 * J2 * (1/(1-e**2)**2) * (1-5*np.cos(inc)**2)
	return omega_dot

def get_M0_dot(a, e, inc, J2=1.082e-3, mu=398600.4415, a_e=6378.1363, radians=True):
	# units in km, s
	n = np.sqrt(mu/a**3)

	if not radians:
		inc = np.radians(inc)	

	arg = 1 - 3/4 * (a_e/a)**2 * J2 * (1/(1-e**2)**(3/2)) * (1 - 3*np.cos(inc)**2)
	M0_dot = n*arg
	return M0_dot


def wrap(angle, radians=False):
	if not radians:
		angle = np.radians(angle)
	angle = (angle + np.pi) % (2 * np.pi) - np.pi
	if not radians:
		angle = np.degrees(angle)
	return angle

def sun_pos(JD0):
	# geocentric sun position
	#	units of km

	n = JD0 - 2451545.0
	L = 280.460 + 0.9856474*n # mean longitude of sun
	g = 357.528 + 0.9856003*n # mean anomaly of sun

	L, g = wrap(L), wrap(g)

	def sind(arg):
		return np.sin(np.radians(arg))
	def cosd(arg):
		return np.cos(np.radians(arg))
	def tand(arg):
		return sind(arg)/cosd(arg)

	lam = L + 1.915*sind(g) + 0.02*sind(2*g)
	# beta = 0.0
	eps = 23.439 - 0.0000004*n
	alpha = np.arctan( cosd(eps)*tand(lam) ) * 180/np.pi # deg

	# alpha2
	f = 180/np.pi
	t = tand(eps/2)**2
	RA = lam - f*t*sind(2*lam) + (f/2)*t**2 * sind(4*lam)
	DEC = np.arcsin(sind(eps)*sind(lam)) * 180/np.pi
	# print(RA,DEC)
	# sys.exit()

	R = 1.00014 - 0.01671*cosd(g) - 0.00014 * cosd(2*g)

	x = R*cosd(lam)
	y = R*cosd(eps)*sind(lam)
	z = R*sind(eps)*sind(lam)

	AU = 149597870.7 # km
	E = (L - alpha)*4 # eqn of time, in minutes

	r_sun = np.array([x,y,z])*AU
	# r_sun2 = cm_ut.sun_pos_approx(JD0)
	return r_sun


# angles = np.random.uniform(-321321, 321321, 1000)
# val1 = np.cos(np.radians(angles))
# for k in range(len(angles)):
# 	angles[k] = wrap(angles[k])
# val2 = np.cos(np.radians(angles))

# sys.exit()

# JD0 = cm_ut.date_to_jd(2006,4,2)
# r_sun = sun_pos(JD0)

JD0 = cm_ut.date_to_jd(2018,5,22)
JDf = cm_ut.date_to_jd(2018+12,5,22)
JD = np.arange(JD0,JDf,1)
t = (JD-JD0)/86400
t0 = t[0]

a = a_e + 500.0 # km
e = 0.0
inc = np.radians(89)
LAN_init = np.radians(0.0)
omega_init = np.radians(0.0)
M0_init = np.radians(0.0)

LAN_dot = get_LAN_dot(a, e, inc)
omega_dot = get_omega_dot(a, e, inc)
M0_dot = get_M0_dot(a, e, inc)

LAN = LAN_init + LAN_dot*(t-t0)
omega = omega_init + omega_dot*(t-t0)
M0 = M0_init + M0_dot*(t-t0)

N = len(t)
OE = np.transpose([[a]*N, [e]*N, [inc]*N, omega, LAN, M0])
r, v = np.zeros((N,3)), np.zeros((N,3))
r_sun = np.zeros((N,3))
for k in range(N):
	r0, v0 = cm_ut.OE2RV(OE[k], t[k], t0, mu)
	r0, v0 = r0.ravel(), v0.ravel()
	r[k] = r0
	v[k] = v0

	r_sun[k] = sun_pos(JD[k])


h = np.cross(r,v)
h_unit = (h.T / np.linalg.norm(h,axis=1)).T
s_unit = (r_sun.T / np.linalg.norm(r_sun,axis=1)).T

proj = h_unit.T[0]*s_unit.T[0] + h_unit.T[1]*s_unit.T[1] + h_unit.T[2]*s_unit.T[2]
beta = 90 - np.arccos(proj)*180/np.pi

# lon_eclip = cm_ut.lon_eclip(JD)
# beta2 = cm_ut.calc_beta_alt(lon_eclip, LAN, inc)



from scipy.signal import find_peaks

index, _ = find_peaks(beta)

fig, ax = ut.make_fig()
ax.plot(JD, beta)
ax.plot(JD[index], beta[index], 'r.')
fig.show()

