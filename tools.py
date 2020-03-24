# -*- coding: utf-8 -*-
########################################################################
# Copyright (C) 2013-2018 by Alexander Perminov
# perminov12@yandex.ru
#
########################################################################
# TRANSFORMATION FUNCTIONS
########################################################################
def dmy2jd(date):
	"""
	Converting day, month, year to JD.
	For dates after 23 november −4713 year.
	"""
	year, month, day, hour, min, sec = date
	a = (14 - month)/12
	y = year + 4800 - a
	m = month + 12*a - 3
	JD = day + (153*m + 2)/5 + 365*y + y/4 - y/100 + y/400 - 32045 + (hour - 12)/24. + min/1440. + sec/86400.
	return JD
########################################################################
def jd2dmy(JD):
	"""
	Converting JD to day, month, year.
	"""
	JDZ = int(JD)
	JDQ = JD - int(JD)
	# --------------------------------
	a = JDZ + 32044
	b = (4*a + 3)/146097
	c = a - 146097*b/4
	d = (4*c + 3)/1461
	e = c - 1461*d/4
	m = (5*e + 2)/153
	# --------------------------------
	day = e - (153*m + 2)/5 + 1
	month = m + 3 - 12*(m/10)
	year = 100*b + d - 4800 + m/10
	hourq = JDQ*24 + 12
	hour = int(hourq)
	minq = (hourq - hour)*60
	min = int(minq)
	secq = (minq - min)*60
	sec = int(secq)
	# --------------------------------
	m30 = [4, 6, 9, 11]
	m31 = [1, 3, 5, 7, 8, 10, 12]
	if hour >= 24:
		hour = hour - 24
		if month in m30:
			if day == 30: day, month = 1, month + 1
			else: day = day + 1
		elif month in m31:
			if day == 31:
				if month == 12: day, month, year = 1, 1, year + 1
				else: day, month = 1, month + 1
			else: day = day + 1
		elif month == 2:
			year = int(year)
			if year%400 == 0:	isLeapYear = True
			elif year%100 == 0:	isLeapYear = False
			elif year%4 == 0:	isLeapYear = True
			else:				isLeapYear = False
			if isLeapYear:
				if day == 29: day, month = 1, month + 1
				else: day = day + 1
			else: day = day + 1
	return year, month, day, hour, min, secq
########################################################################
def kepler2poincare(elem_k, kM):
	"""
	Transform Kepler orbital elements to Poincare elements.
	- 'kM' is sqrt(kappa2) * Mass,
	- 'aa' is semi-major axis,
	- 'ee'  is eccentricity,
	- 'ii'  is inclination,
	- 'om'  is argument of pericenter,
	- 'OM'  is longitude of ascending node,
	- 'MM'  is mean anomaly.
	"""
	from math import sin, cos, sqrt, fmod, pi
	aa, ee, ii, om, OM, MM = elem_k
	e2 = sqrt(1. - ee**2)
	LL = kM*sqrt(aa)
	PP = sqrt(2*LL * (1. - e2))
	QQ = sqrt(2*LL * e2 * (1. - cos(ii)))
	xx =  PP * cos(om + OM)
	yy = -PP * sin(om + OM)
	uu =  QQ * cos(OM)
	vv = -QQ * sin(OM)
	ll = om + OM + MM
	return [LL, xx, yy, uu, vv, ll]
########################################################################
def poincare2kepler(elem_p, kM):
	"""
	Transform Poincare elements to Kepler orbital elements.
	"""
	from math import acos, atan2, sqrt, fmod, pi
	LL, xx, yy, uu, vv, ll = elem_p
	aa = (LL/kM)**2
	ee = sqrt(1. - (1. - 0.5*(xx**2 + yy**2)/LL)**2)
	ii = acos(1. - 0.5*(uu**2 + vv**2)/LL/sqrt(1. - ee**2))
	OM = atan2(-vv, uu)
	UU = atan2(-yy, xx)
	om = UU - OM
	MM = ll - om - OM
	return [aa, ee, ii, om, OM, MM]
########################################################################
def cart2kepler(cart, k2):
	"""
	Calculation values of Kepler orbital elements.
	- 'k2' is gravitational parameter,
	- 'x, y, z' are coordinates,
	- 'vx, vy, vz' are velocities.
	"""
	from math import acos, atan, atan2, sin, cos, fabs, tan, sqrt, pi
	x, y, z, vx, vy, vz = cart
	r = sqrt(x**2 + y**2 + z**2)
	v = sqrt(vx**2 + vy**2 + vz**2)
	a = (2./r - v**2/k2)**-1
	if a > 0:
		s = x*vx + y*vy + z*vz
		sx = y*vz - z*vy
		sy = x*vz - z*vx
		sz = x*vy - y*vx
		p = (sx**2 + sy**2 + sz**2)/k2
		ra = 1. - r/a
		e = sqrt(s**2/(k2*a) + ra**2)
		t = sz/sqrt(k2*p)
		#print '%.20f'%t
		if fabs(fabs(t) - 1) <= 9e-16:
			i = 0 if t > 1 else pi*0
		else:
			i = acos(t)
		OM = atan2(sx, sy)
		th = atan2(sqrt(p)*s, sqrt(k2)*(p - r))
		uu = atan2((y*cos(OM)-x*sin(OM))/cos(i), x*cos(OM) + y*sin(OM))
		om = uu - th
		EE = 2.*atan(sqrt((1. - e)/(1. + e))*tan(th/2.))
		MM = EE - e*sin(EE)
		return [a, e, i, om, OM, MM]
	else:
		print a, 'WARNING: hyperbolic or parabolic orbit!'
		return [a, 0, 0, 0, 0, 0]
########################################################################
def kepler2cart(elem_k, k2):
	"""
	Calculation values of x, y, z.
	- 'k2' is gravitational parameter,
	- 'aa' is semi-major axis,
	- 'ee', 'ii', 'om', 'OM', 'MM' see above.
	"""
	from math import atan, tan, sin, cos, sqrt, fabs
	aa, ee, ii, om, OM, MM = elem_k
	EE, eps = MM, 1000
	while eps >= 1e-14:
		EM = MM + ee*sin(EE)
		eps = fabs(EM - EE)
		EE = EM
	th = 2.*atan(sqrt((1. + ee)/(1. - ee))*tan(EE/2.))
	uu = th + om
	pp = aa*(1. - ee**2)
	rr = pp/(1. + ee*cos(th))
	xx = rr * (cos(uu)*cos(OM) - sin(uu)*sin(OM)*cos(ii))
	yy = rr * (cos(uu)*sin(OM) + sin(uu)*cos(OM)*cos(ii))
	zz = rr * sin(uu)*sin(ii)
	vr = sqrt(k2/pp)*ee*sin(th)
	vn = sqrt(k2/pp)*(1. + ee*cos(th))
	vx = xx*vr/rr + (-sin(uu)*cos(OM) - cos(uu)*sin(OM)*cos(ii)) * vn
	vy = yy*vr/rr + (-sin(uu)*sin(OM) + cos(uu)*cos(OM)*cos(ii)) * vn
	vz = zz*vr/rr + cos(uu)*sin(ii) * vn
	return [xx, yy, zz, vx, vy, vz]
########################################################################
def barySolar(cart, mass):
	"""
	Calculation barycentric Cartesian coordinates of the Sun.
	- 'mass' is list of planets masses (by Solar mass);
	- 'cart' is list of Cartesian coordinates (or velocities).
	"""
	N = len(mass)
	solc = [sum([-mass[k]*cart[k][j] for k in range(0, N)]) for j in range(3)]
	return solc
########################################################################
def cart2jacobi(cart, mass):
	"""
	Transform barycentric Cartesian coordinates to Jacobi coordinates.
	- 'mass' is list of planets masses (by Solar mass);
	- 'cart' is list of Cartesian coordinates (or velocities).
	"""
	N = len(mass)
	sum_m = [sum([mass[i] for i in range(n+1)]) for n in range(N)]
	jc = []
	jc.append([sum([(mass[k]/sum_m[N-1])*cart[k][j] for k in range(1, N)], cart[0][j]/sum_m[N-1]) for j in range(3)])
	for i in range(1, N):
		jc.append([sum([-(mass[k]/sum_m[i-1])*cart[k][j] for k in range(1, i)], cart[i][j]-cart[0][j]/sum_m[i-1]) for j in range(3)])
	return jc
########################################################################
def jacobi2cart(jc, mass):
	"""
	Transform Jacobi coordinates to barycentric Cartesian coordinates.
	- 'mass' is list of planets masses (by Solar mass);
	- 'jc' is list of Jacobi coordinates (or velocities).
	"""
	N = len(mass)
	sum_m = [sum([mass[i] for i in range(n+1)]) for n in range(N)]
	cart = []
	cart.append([sum([-(mass[k]/sum_m[k])*jc[k][j] for k in range(1, N)], jc[0][j]) for j in range(3)])
	for i in range(1, N):
		cart.append([sum([-(mass[k]/sum_m[k])*jc[k][j] for k in range(i+1, N)], jc[0][j]+(sum_m[i-1]/sum_m[i])*jc[i][j]) for j in range(3)])
	return cart
########################################################################
# CALCULATION OF MASS PARAMETERS
########################################################################
def jacobiAMass(gmass, mu):
	"""
	Calculation mass-parameters in Jacobi coordinates.
	"""
	rlen = range(len(gmass))
	mp = [gmass[i]/gmass[0]/mu for i in rlen]
	sum_mp = [1.+mu*sum(mp[1:i+1]) for i in rlen]
	gmj = [gmass[0] if i == 0 else gmass[0]*sum_mp[i]/sum_mp[i-1] for i in rlen]
	mpj = [1. if i == 0 else gmass[0]*mp[i]/gmj[i] for i in rlen]
	km = [mpj[i]*gmj[i]**0.5 for i in rlen]
	return gmj[1:], mpj[1:], km[1:]
########################################################################
def jacobiGMass(gmass, mu):
	"""
	Calculation gravitational parameters in Jacobi coordinates.
	"""
	rlen = range(len(gmass))
	mp = [gmass[i]/gmass[0]/mu for i in rlen]
	sum_mp = [1.+mu*sum(mp[1:i+1]) for i in rlen]
	gmj = [gmass[0] if i == 0 else gmass[0]*sum_mp[i]/sum_mp[i-1] for i in rlen]
	return gmj[1:]
########################################################################
def jacobiNMass(gmass, mu):
	"""
	Calculation normalized masses in Jacobi coordinates.
	"""
	rlen = range(len(gmass))
	mp = [gmass[i]/gmass[0]/mu for i in rlen]
	sum_mp = [1.+mu*sum(mp[1:i+1]) for i in rlen]
	gmj = [gmass[0] if i == 0 else gmass[0]*sum_mp[i]/sum_mp[i-1] for i in rlen]
	mpj = [1. if i == 0 else gmass[0]*mp[i]/gmj[i] for i in rlen]
	return mpj[1:]
########################################################################
def jacobiKMass(gmass, mu):
	"""
	Calculation kappa*Mass in Jacobi coordinates.
	"""
	rlen = range(len(gmass))
	mp = [gmass[i]/gmass[0]/mu for i in rlen]
	sum_mp = [1.+mu*sum(mp[1:i+1]) for i in rlen]
	gmj = [gmass[0] if i == 0 else gmass[0]*sum_mp[i]/sum_mp[i-1] for i in rlen]
	mpj = [1. if i == 0 else gmass[0]*mp[i]/gmj[i] for i in rlen]
	km = [mpj[i]*gmj[i]**0.5 for i in rlen]
	return km[1:]
########################################################################
def massesRatio(gmass, mu):
	"""
	Calculation masses ratios.
	"""
	rlen = range(len(gmass)-1)
	mp = [gmass[i]/gmass[0]/mu for i in rlen]
	sum_mp = [1.+mu*sum(mp[1:i+1]) for i in rlen]
	mm = [mp[i]/sum_mp[i] for i in rlen]
	return mm
########################################################################
# TRANSFORMATION OF ORBITAL ELEMENTS
########################################################################
def cart2elem_jacobi(cart, gmass, mu, path_write = ''):
	"""
	The conversion of barycentric Cartesian coordinates to elements in Jacobi frame.
	"""
	from math import fmod, pi
	rlen1 = range(len(gmass))
	rlen2 = range(len(cart))
	km = jacobiKMass(gmass, mu)
	gm = jacobiGMass(gmass, mu)
	mass = [gmass[i]/gmass[0] for i in rlen1]
	cart_r = [cart[k][:3] for k in rlen2]
	cart_v = [cart[k][3:] for k in rlen2]
	barj_r = barySolar([[0. for i in range(3)]] + cart_r, mass)
	barj_v = barySolar([[0. for i in range(3)]] + cart_v, mass)
	jacb_r = cart2jacobi([barj_r] + cart_r, mass)[1:]
	jacb_v = cart2jacobi([barj_v] + cart_v, mass)[1:]
	kepj = [cart2kepler(jacb_r[k]+jacb_v[k], gm[k]) for k in rlen2]
	poij = [kepler2poincare(kepj[k], km[k]) for k in rlen2]
	for i in range(len(cart)):
		for j in [3, 4, 5]:
			if kepj[i][j] < 0.: kepj[i][j] += 2.*pi
			kepj[i][j] = fmod(kepj[i][j], 2*pi)
		if poij[i][5] < 0.: poij[i][5] += 2.*pi
		poij[i][5] = fmod(poij[i][5], 2*pi)
	if path_write != '':
		f = open(path_write + 'kepler_jacobi.txt', 'w')
		for planet in kepj:
			for elem in planet:
				f.write(str('%.16e'%elem).rjust(24, ' '))
		f.close()
		f = open(path_write + 'poincare_jacobi.txt', 'w')
		for planet in poij:
			for elem in planet:
				f.write(str('%.16e'%elem).rjust(24, ' '))
		f.close()
	return kepj, poij
########################################################################
def kepler_cart2jacobi(elements, gmass, mu, path_write = ''):
	"""
	The conversion of Keplerian elements in barycentric frame to elements in Jacobi frame.
	"""
	from math import fmod, pi
	rlen1 = range(len(gmass))
	rlen2 = range(len(elements))
	km = jacobiKMass(gmass, mu)
	gm = jacobiGMass(gmass, mu)
	mass = [gmass[i]/gmass[0] for i in rlen1]
	cart   = [kepler2cart(elements[k], gm[k]) for k in rlen2]
	cart_r = [cart[k][:3] for k in rlen2]
	cart_v = [cart[k][3:] for k in rlen2]
	barj_r = barySolar([[0. for i in range(3)]] + cart_r, mass)
	barj_v = barySolar([[0. for i in range(3)]] + cart_v, mass)
	jacb_r = cart2jacobi([barj_r] + cart_r, mass)[1:]
	jacb_v = cart2jacobi([barj_v] + cart_v, mass)[1:]
	kepj = [cart2kepler(jacb_r[k]+jacb_v[k], gm[k]) for k in rlen2]
	poij = [kepler2poincare(kepj[k], km[k]) for k in rlen2]
	for i in range(len(elements)):
		for j in [3, 4, 5]:
			if kepj[i][j] < 0.: kepj[i][j] += 2.*pi
			kepj[i][j] = fmod(kepj[i][j], 2*pi)
		if poij[i][5] < 0.: poij[i][5] += 2.*pi
		poij[i][5] = fmod(poij[i][5], 2*pi)
	if path_write != '':
		f = open(path_write + 'kepler_jacobi.txt', 'w')
		for planet in kepj:
			for elem in planet:
				f.write(str('%.16e'%elem).rjust(24, ' '))
		f.close()
		f = open(path_write + 'poincare_jacobi.txt', 'w')
		for planet in poij:
			for elem in planet:
				f.write(str('%.16e'%elem).rjust(24, ' '))
		f.close()
	return kepj, poij
########################################################################
def kepler_jacobi2cart(elements, gmass, mu, path_write = ''):
	"""
	The conversion of Keplerian elements from Jacobi frame to elements in barycentric frame.
	"""
	from math import fmod, pi
	rlen1 = range(len(gmass))
	rlen2 = range(len(elements))
	mp = [gmass[i]/gmass[0]/mu for i in rlen1]
	km = jacobiKMass(gmass, mu)
	gm = jacobiGMass(gmass, mu)
	mass = [gmass[i]/gmass[0] for i in rlen1]
	jacobi = [kepler2cart(elements[k], gm[k]) for k in rlen2]
	jacb_r = [jacobi[k][:3] for k in rlen2]
	jacb_v = [jacobi[k][3:] for k in rlen2]
	bary_r = jacobi2cart([[0. for i in range(3)]] + jacb_r, mass)[1:]
	bary_v = jacobi2cart([[0. for i in range(3)]] + jacb_v, mass)[1:]
	kepb = [cart2kepler(bary_r[k]+bary_v[k], gm[k]) for k in rlen2]
	poib = [kepler2poincare(kepb[k], km[k]) for k in rlen2]
	for i in range(len(elements)):
		for j in [3, 4, 5]:
			if kepb[i][j] < 0.: kepb[i][j] += 2.*pi
			kepb[i][j] = fmod(kepb[i][j], 2*pi)
		if poib[i][5] < 0.: poib[i][5] += 2.*pi
		poib[i][5] = fmod(poib[i][5], 2*pi)
	if path_write != '':
		f = open(path_write + 'kepler_bary.txt', 'w')
		for planet in kepb:
			for j in [3, 4, 5]:
				planet[j] = fmod(planet[j], 2*pi)
			for elem in planet:
				f.write(str('%.16e'%elem).rjust(24, ' '))
		f.close()
		f = open(path_write + 'poincare_bary.txt', 'w')
		for planet in poib:
			planet[5] = fmod(planet[5], 2*pi)
			for elem in planet:
				f.write(str('%.16e'%elem).rjust(24, ' '))
		f.close()
	return kepb, poib
########################################################################
# READING OF ORBITAL ELEMENTS AND COORDINATES
########################################################################
def read_coords(path_in, gmass, mu, num_pl = 4, isJacobi = True):
	"""
	Read data from a file with integration results.
	"""
	from math import pi
	from numpy import array
	f = open(path_in, 'r')
	s = f.read().split('\n')
	f.close()
	M = 2*(num_pl+1)+1
	N = len(s)/M
	data = [[] for i in range(5*num_pl+1)]
	for k in range(N):
		cart = [[0 for j in range(6)] for i in range(4)]
		time = float(s[k*M].strip().split()[0])/365.25
		for j in range(3, M, 2):
			cart[(j-1)/2-1][:3] = [float(item) for item in s[k*M+j].strip().split()]
			cart[(j-1)/2-1][3:] = [float(item) for item in s[k*M+j+1].strip().split()]
		kepler_j = cart2elem_jacobi(cart, gmass, mu)[0]
		kepler_b = kepler_jacobi2cart(kepler_j, gmass, mu)[0]
		for j in range(len(kepler_b)):
			for i in range(2, 6):
				kepler_b[j][i] *= 180./pi
				kepler_j[j][i] *= 180./pi
		for j in range(len(kepler_b)):
			for i in range(3, 6):
				if kepler_b[j][i] < 0.: kepler_b[j][i] += 360.
				if kepler_j[j][i] < 0.: kepler_j[j][i] += 360.
		data[0].append(time)
		ji = 1
		for j in range(num_pl):
			for i in range(5):
				if not isJacobi: data[ji].append(kepler_b[j][i])
				else: data[ji].append(kepler_j[j][i])
				ji += 1
	return (path_in, array(data), num_pl, True)
##################################################
def read_elems(path_in):
	"""
	Read data from a file with integration results.
	"""
	from numpy import array
	f = open(path_in, 'r')
	s = f.read().split('\n')
	f.close()
	N = len(s)-2
	is_kepler = True if '_ke' in path_in.split('/')[-1] else False
	t = 1 if '_osc' in path_in.split('/')[-1] else 0
	num_pl = (len(s[1].strip().split())-3+2*t)/5
	data = [[] for i in range(5*num_pl+3-2*t)]
	for k in range(N):
		line = s[k+1].strip().split()
		for j in range(5*num_pl+3-2*t):
			if is_kepler:
				if j in [4, 5, 9, 10, 14, 15, 19, 20] and float(line[j]) < 0.:
					data[j].append(float(line[j]) + 360.)
				else:
					data[j].append(float(line[j]))
			else:
				data[j].append(float(line[j]))
	return (path_in, array(data), num_pl, is_kepler)
########################################################################
def read_alpha(path_in):
	"""
	Read data from a file with longitudes.
	"""
	from numpy import array
	f = open(path_in, 'r')
	s = f.read().split('\n')
	f.close()
	N = len(s)-2
	num_pl = len(s[1].strip().split())-1
	data = [[] for i in range(num_pl+1)]
	for k in range(N):
		line = s[k+1].strip().split()
		for j in range(num_pl+1):
			data[j].append(float(line[j]))
	return (path_in, array(data), num_pl)
########################################################################
# OTHER UTILITY FUNCTIONS
########################################################################
def evl_f(elmass, series_input, series_type):
	"""
	Evaluation of the series.
	"""
	from mpmath import mpf
	from pyranha.math import evaluate
	# -------- coordinates & masses --------
	num_pl = len(elmass[0])
	gmj, mpj, km = jacobiAMass(elmass[1], elmass[2])
	mm = massesRatio(elmass[1], elmass[2])
	mp = [1./elmass[2]] + [elmass[1][i]/elmass[1][0]/elmass[2] for i in range(len(elmass[1]))][1:]
	sp = [1. + elmass[2]*sum(mp[1:i+1]) for i in range(len(elmass[1]))]
	# -------- dictionaries --------
	evaldict = {'K0': mpf(elmass[1][0])}
	evaldict.update({'m'+str(i+1): mpf(mp[i+1]) for i in range(num_pl)})
	evaldict.update({'s'+str(i+1): mpf(sp[i+1]) for i in range(num_pl)})
	evaldict.update({'L'+str(i+1): mpf(elmass[0][i][0]) for i in range(num_pl)})
	evaldict.update({'x'+str(i+1): mpf(elmass[0][i][1]) for i in range(num_pl)})
	evaldict.update({'y'+str(i+1): mpf(elmass[0][i][2]) for i in range(num_pl)})
	evaldict.update({'u'+str(i+1): mpf(elmass[0][i][3]) for i in range(num_pl)})
	evaldict.update({'v'+str(i+1): mpf(elmass[0][i][4]) for i in range(num_pl)})
	evaldict.update({'q'+str(i+1): mpf(elmass[0][i][5]) for i in range(num_pl)})
	evaldict.update({r'\nu_{q'+str(i+1)+'}': mpf(gmj[i]**2*mpj[i]**3*elmass[0][i][0]**-3) for i in range(num_pl)})
	if type(series_input) == str:
		series = series_type
		lf(series, series_input, df.boost_portable, cf.bzip2)
		series_eval = evaluate(series, evaldict)
	else:
		series_eval = evaluate(series_input, evaldict)
	return float(series_eval)
########################################################################
def rotateMatrix(xyz, axis, angle):
	"""
	The rotation matrix.
	"""
	from math import sin, cos
	if axis == 0:		# rotate around x-axis
		return [xyz[0], xyz[1]*cos(angle) - xyz[2]*sin(angle), xyz[1]*sin(angle) + xyz[2]*cos(angle)]
	elif axis == 1:		# rotate around y-axis
		return [xyz[0]*cos(angle) + xyz[2]*sin(angle), xyz[1], -xyz[0]*sin(angle) + xyz[2]*cos(angle)]
	elif axis == 2:		# rotate around z-axis
		return [xyz[0]*cos(angle) - xyz[1]*sin(angle), xyz[0]*sin(angle) + xyz[1]*cos(angle), xyz[2]]
	else:
		return [xyz[0], xyz[1], xyz[2]]
########################################################################
def zipping_dir(path_in, path_out):
	"""
	Zipping files in directory 'path_in' to 'path_out'.
	"""
	import os
	import tarfile
	if not os.path.exists(path_out): os.makedirs(path_out)
	for file in os.listdir(path_in):
		print path_in+file
		tar = tarfile.open(path_out + file[:-5] + '.tar.gz', mode='w:gz')
		tar.add(path_in+file, arcname=file, recursive=False)
		tar.close()
########################################################################
def latex_out(series, path_out):
	"""
	Print series in LaTex format.
	"""
	f = open(path_out + 'out.tex', 'w')
	latex = ''
	series_latex = series._repr_latex_()
	# replacing:
	series_latex = series_latex.replace('mu', '\mu')
	series_latex = series_latex.replace('K0', '\kappa_0^2')
	for item in ['K1' , 'K2','K3', 'K4', 'M1', 'M2', 'M3', 'M4', 'm1', 'm2', 'm3', 'm4']: series_latex = series_latex.replace(item, item[0]+'_{'+item[1]+'}')
	for item in ['s1', 's2', 's3', 's4']: series_latex = series_latex.replace(item, '\overline{m}_{'+item[1]+'}')
	for item in ['mm1', 'mm2', 'mm3', 'mm4']: series_latex = series_latex.replace(item, '\overline{'+item[0]+'}_{'+item[2]+'}')
	for item in ['L1', 'L2', 'L3', 'L4', 'q1', 'q2', 'q3', 'q4', 'x1', 'x2', 'x3', 'x4', 'y1', 'y2', 'y3', 'y4',\
				 'u1', 'u2', 'u3', 'u4', 'v1', 'v2', 'v3', 'v4']: series_latex = series_latex.replace(item, item[0]+'_{'+item[1]+'}')
	for item in [r'\nu_{q1}', r'\nu_{q2}', r'\nu_{q3}', r'\nu_{q4}']: series_latex = series_latex.replace(item, item[0:3]+'_{'+item[6]+'}')
	count1, count2, temp = 0, 0, ''
	for item in series_latex[9:-10]:
		temp1, temp2 = '', ''
		if item == '+' or item == '-':
			count1 += 1
			count2 += 1
		if count1 == 3:
			temp1 = '\\\\\\\\&' + item
			count1 = 0
		if count2 == 54:
			temp2 = '\n\end{align*}\n\pagebreak\n\\begin{align*}\n'
			count2 = 0
		temp = temp + item + temp2 + temp1
	series_latex = '&' + temp
	latex = latex + '\\begin{align*}\n'+series_latex+'\n\end{align*}\n'
	# writing:
	head = '\documentclass[8pt,russian]{report}\n\n\RequirePackage[utf8x]{inputenc}\n\RequirePackage[T2A]{fontenc}\n' + \
	'\RequirePackage[russian,english]{babel}\n\RequirePackage{extsizes}\n\RequirePackage{epsfig}\n' + \
	'\RequirePackage{amsmath}\n\RequirePackage{amssymb}\n\RequirePackage{amsfonts}\n\n' + \
	'\usepackage[left=1cm,right=1cm,top=1cm,bottom=1cm,bindingoffset=0cm]{geometry}\n\n' + \
	'\setlength{\\textwidth}{170mm}\n\setlength{\\textheight}{250mm}\n\n\selectlanguage{russian}\n\\footnotesize\n\n'
	f.write(head)
	f.write('\\begin{document}\n\n')
	f.write(latex)
	f.write('\n\end{document}\n\n')
	f.close()
########################################################################

