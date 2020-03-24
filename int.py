# -*- coding: utf-8 -*-
########################################################################
# Copyright (C) 2013-2018 by Alexander Perminov
# perminov12@yandex.ru
#
########################################################################
# THIS IS INTEGRATION TOOLSET
########################################################################
import os, time
from eps.config import *
########################################################################
# WRITING TEXT FILES FOR THE INTEGRATION
########################################################################
def sub_h(elmass = [[], [], 0], PROBLEM = '', prefix = '', postfix = '', part = 0, isDouble = False, isEval = False):
	'''
	The substitution mass parameters, denominators and L-elements (semi-major axes) into the Hamiltonian expansion.
	'''
	from mpmath import mpf
	from pyranha.math import degree, evaluate, subs
	from eps.tools import jacobiAMass
	# -------- elements & masses --------
	gmj, mpj, km = jacobiAMass(elmass[1], elmass[2])
	mp = [1./elmass[2]] + [elmass[1][i]/elmass[1][0]/elmass[2] for i in range(len(elmass[1]))][1:]
	sp = [1. + elmass[2]*sum(mp[1:i+1]) for i in range(len(elmass[1]))]
	list_elem = [jtem+item for item in '1234' for jtem in 'xyuv']
	num_pl = len(elmass[0])
	max_range = num_pl*4 if num_pl in [2, 3] else 16
	(epst2, type2) = (epsd, 'd') if isDouble else (epst, 't')
	typef = mpf if isEval else float
	# -------- dictionaries --------
	evaldict = {'K0': typef(elmass[1][0])}
	evaldict.update({'m'+str(i+1): typef(mp[i+1]) for i in range(num_pl)})
	evaldict.update({'s'+str(i+1): typef(sp[i+1]) for i in range(num_pl)})
	evaldict.update({'L'+str(i+1): typef(elmass[0][i][0]) for i in range(num_pl)})
	eval_dnu = {r'\nu_{q'+str(i+1)+'}': typef(gmj[i]**2*mpj[i]**3*elmass[0][i][0]**-3) for i in range(num_pl)}
	# -------- evaluation --------
	order = [0, 1, 2, 2, 3, 3, 3, 3, 3, 3, 3]
	#index = {0: '', 1: '', 2: '_t1h1_tor', 3: '_2', 4: '_t1h2', 5: '_t2h1_h2', 6: '_t2h1_t1a1', 7: '_t2h1_t1h1', 8: '_t1t1a1', 9: '_t1t1h1', 10: '_3'}
	const = [1., elmass[1][0]*elmass[2], (elmass[2]*elmass[1][0])**2, elmass[1][0]*elmass[2]**2]
	const += [elmass[1][0]**2*elmass[2]**3 for i in range(2)] + [(elmass[1][0]*elmass[2])**3 for i in range(4)] + [elmass[1][0]*elmass[2]**3]
	# -------- hamiltonian --------
	start = time.time()
	series_sum = 0
	PATH_DIR = PREFIX + 'SUM/HAM/HAM'+str(order[part])+prefix+'/'
	for file2 in os.listdir(PATH_DIR):
		if postfix == '' or postfix in file2[:-16]: 
			series = epst2(0)
			lf(series, PATH_DIR + file2, df.boost_portable, cf.bzip2)
			if num_pl < 4:
				for key in ['m'+str(i+1) for i in range(num_pl, 4)]: series = subs(series, key, 0)
				series = series.trim()
			if isEval:
				c, stemp = 0, 0
				for item in series.list:
					for jtem in item[0].list:
						c += 1
						stemp += jtem[0]*float(evaluate(jtem[1], eval_dnu))
				if type(stemp) == int:
					temp = stemp
				else:
					stemp = stemp.trim()
					print c, 'denominators'
					c, temp = 0, 0
					length = len(stemp)
					for item in stemp.list:
						c += 1
						poly = 1
						for j in range(max_range):
							poly *= pd(pq_list[2*num_pl+j])**degree(item[1], pq_list[2*num_pl+j:2*num_pl+1+j])
						temp += float(evaluate((item[0]*item[1]*poly**-1).trim(), evaldict))*poly
						if c%10000 == 0:
							print str(c)+'/'+str(length)
					temp = temp.trim()
			else:
				c, temp = 0, 0
				for item in series.list:
					for jtem in item[0].list:
						c += 1
						temp += jtem[0]*evaluate(jtem[1], eval_dnu)
				if type(temp) != int:
					for key in evaldict.keys():
						temp = subs(temp, key, evaldict[key])
						print key
					temp = temp.trim()
			series_sum += const[part]*temp
			print str(file2[:-16])+':', c, 'terms for', int(time.time() - start), 's'
	# -------- write to file --------
	PATH_SAVE = PREFIX + 'EVL/SUMS/SUMS_'+PROBLEM+'/'
	if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
	sf(series_sum.trim(), PATH_SAVE + '/h'+str(order[part])+prefix.lower()+postfix+'.epsd.boostp.bz2', df.boost_portable, cf.bzip2)
########################################################################
def txt_N(elmass = [[], [], 0], PROBLEM = '', order = 1, eps = 1):
	'''
	Writing text-files for the integration process (the slow variables).
	'''
	from math import fabs
	from pyranha.math import degree, evaluate, partial
	list_elem = [jtem+item for item in '1234' for jtem in 'xyuv']
	num_pl = len(elmass[0])
	max_range = num_pl*4 if num_pl in [2, 3] else 16
	series = {key: pd(0) for key in ['hm'] + pq_list[2*num_pl:]}
	postfix = ''
	# -------- loop --------
	for file2 in os.listdir(PREFIX + 'EVL/SUMS/SUMS_'+PROBLEM+'/'):
		if int(file2[1]) <= order:
			temp = pd(0)
			lf(temp, PREFIX + 'EVL/SUMS/SUMS_'+PROBLEM+'/'+file2, df.boost_portable, cf.bzip2)
			series['hm'] += temp
	if eps != 1:
		postfix = '_'+str(eps)
		evaldict = {'x1': elmass[0][0][1], 'y1': elmass[0][0][2], 'u1': elmass[0][0][3], 'v1': elmass[0][0][4],\
					'x2': elmass[0][1][1], 'y2': elmass[0][1][2], 'u2': elmass[0][1][3], 'v2': elmass[0][1][4],\
					'x3': elmass[0][2][1], 'y3': elmass[0][2][2], 'u3': elmass[0][2][3], 'v3': elmass[0][2][4],\
					'x4': elmass[0][3][1], 'y4': elmass[0][3][2], 'u4': elmass[0][3][3], 'v4': elmass[0][3][4]}
		temp_sum = 0
		for jtem in series['hm'].list:
			temp = jtem[0]*jtem[1]
			if degree(temp, pq_list[2*num_pl:]) <= 4:
				temp_sum += temp
				continue
			if fabs(evaluate(temp, evaldict)) > eps:
				temp_sum += temp
		series['hm'] = temp_sum
	# -------- motion equations --------
	for item in list_elem[:max_range]:
		dict_key = {'x': 'y', 'y': 'x', 'u': 'v', 'v': 'u'}
		for key in dict_key.keys():
			if item[0] == key:
				jtem = item.replace(item[0], dict_key[key])
		series[item] = -partial(series['hm'], jtem) if (jtem[0] in 'qyv') else partial(series['hm'], jtem)
	# -------- write to files hamiltonian and equations --------
	PATH_SAVE = PREFIX + 'EVL/TXTS/TXT_'+PROBLEM+postfix+'/'
	if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
	lines = {key: '' for key in ['hm'] + pq_list[2*num_pl:]}
	length = {key: 0 for key in ['hm'] + pq_list[2*num_pl:]}
	for what in ['hm'] + pq_list[2*num_pl:]:
		for item in series[what].list:
			if type(item[0]) != float:
				for jtem in item[0].list:
					for ktem in jtem[0].list:
						lines[what] += str('%.16e'%ktem[0]).rjust(24, ' ')
						length[what] += 1
						for elem in list_elem[:max_range]:
							deg_elem = degree(ktem[1], [elem])
							lines[what] += str(deg_elem).rjust(3, ' ')
						lines[what] += '\n'
			else:
				lines[what] += str('%.16e'%item[0]).rjust(24, ' ')
				length[what] += 1
				for elem in list_elem[:max_range]:
					deg_elem = degree(item[1], [elem])
					lines[what] += str(deg_elem).rjust(3, ' ')
				lines[what] += '\n'
		f2 = open(PATH_SAVE + what+'.txt', 'w')
		f2.write(str(length[what])+'\n'+lines[what])
		print what+':', length[what]
	# -------- write to files masses and elements --------
	f3 = open(PATH_SAVE + 'el.txt', 'w')
	for planet in elmass[0][:num_pl]:
		for elem in planet:
			f3.write(str('%.16e'%elem).rjust(24, ' '))
		f3.write('\n')
	f3.close()
	f4 = open(PATH_SAVE + 'gm.txt', 'w')
	for item in [elmass[2]]+elmass[1][:num_pl+1]:
		f4.write(str('%.16e'%item).rjust(24, ' ')+'\n')
	f4.close()
########################################################################
def txt_q(elmass = [[], [], 0], PROBLEM = '', paths = ['', '', ''], isDouble = False, isEval = False):
	'''
	Writing text-files for the integration process (the fast variables).
	'''
	from mpmath import mpf
	from pyranha.math import degree, evaluate, partial, subs
	from eps.tools import jacobiAMass
	# -------- elements & masses --------
	gmj, mpj, km = jacobiAMass(elmass[1], elmass[2])
	mp = [1./elmass[2]] + [elmass[1][i]/elmass[1][0]/elmass[2] for i in range(len(elmass[1]))][1:]
	sp = [1. + elmass[2]*sum(mp[1:i+1]) for i in range(len(elmass[1]))]
	list_elem = [jtem+item for item in '1234' for jtem in 'xyuv']
	num_pl = len(elmass[0])
	order = 2 if (paths[1] != '' and paths[2] != '') else 1
	max_range = num_pl*4 if num_pl in [2, 3] else 16
	epst2 = epsd if isDouble else epst
	typef = mpf if isEval else float
	series = [epst(0) for key in range(len(paths))]
	series2 = [0 for key in range(num_pl)]
	# -------- dictionaries & lists --------
	evaldict = {'K0': typef(elmass[1][0])}
	evaldict.update({'m'+str(i+1): typef(mp[i+1]) for i in range(num_pl)})
	evaldict.update({'s'+str(i+1): typef(sp[i+1]) for i in range(num_pl)})
	evaldict.update({'L'+str(i+1): typef(elmass[0][i][0]) for i in range(num_pl)})
	const = [elmass[1][0]*elmass[2], (elmass[1][0]*elmass[2])**2, elmass[1][0]*elmass[2]**2]
	# -------- hamiltonian --------
	for i in range(len(paths[:2*order-1])): lf(series[i], PREFIX + 'SUM/HAM/'+paths[i]+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
	if num_pl < 4:
		for key in paths[:2*order-1]:
			for key2 in ['m'+str(i+1) for i in range(num_pl, 4)]: series[key] = subs(series[key], key2, 0)
			series[key] = series[key].trim()
	# -------- evaluation of motion equations --------
	start = time.time()
	for i in range(num_pl):
		for key in range(len(paths[:2*order-1])):
			temp_q = partial(series[key], 'L'+str(i+1))
			if isEval:
				c, temp = 0, 0
				for item in temp_q.list:
					for jtem in item[0].list:
						for ktem in jtem[0].list:
							c += 1
							poly = 1
							for j in range(max_range):
								poly *= pt(pq_list[num_pl*2+j])**degree(ktem[1], pq_list[num_pl*2+j:num_pl*2+1+j])
							temp += float(evaluate((ktem[0]*ktem[1]*poly**-1).trim(), evaldict))*poly
							if c%100000 == 0: print c
				temp_q = const[key]*temp
				print 'q'+str(i+1), '->', c, 'terms for', int(time.time() - start), 's'
			else:
				for key2 in evaldict.keys():
					temp_q = subs(temp_q, key2, evaldict[key2])
				temp_q = const[key]*(temp_q.trim())
				print 'q'+str(i+1), 'for', int(time.time() - start), 's'
			series2[i] += temp_q
	# -------- write to files equations --------
	PATH_SAVE = PREFIX + 'EVL/TXTS/TXT_'+PROBLEM+'/'
	if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
	lines = ['' for key in range(num_pl)]
	length = [0 for key in range(num_pl)]
	for i in range(num_pl):
		for item in series2[i].list:
			if type(item[0]) != float:
				for jtem in item[0].list:
					for ktem in jtem[0].list:
						lines[i] += str('%.16e'%ktem[0]).rjust(24, ' ')
						length[i] += 1
						for elem in list_elem[:max_range]:
							deg_elem = degree(ktem[1], [elem])
							lines[i] += str(deg_elem).rjust(3, ' ')
						lines[i] += '\n'
			else:
				lines[i] += str('%.16e'%item[0]).rjust(24, ' ')
				length[i] += 1
				for elem in list_elem[:max_range]:
					deg_elem = degree(item[1], [elem])
					lines[i] += str(deg_elem).rjust(3, ' ')
				lines[i] += '\n'
		f2 = open(PATH_SAVE + 'q'+str(i+1)+'.txt', 'w')
		f2.write(str(length[i])+'\n'+lines[i])
		print 'q'+str(i+1)+':', length[i]
########################################################################
# THE INTEGRATION ALGORITHMS
########################################################################
def spline(path_to_files, num_pl = 4, korder = 5, prefix = '', is_day = True):
	'''
	The integration of the motion equations
	for the fast variables by B-splines interpolation.
	'''
	from numpy import fmod, log10, sqrt
	from scipy.interpolate import splint, splev, splrep
# set of step size in days or years
	SIZE_STEP = 365.25 if is_day else 1.
# length of the motion equations for the fast variables
	NTERMS = 25000
# names of files
	fnames = ['q'+str(i+1)+'.txt' for i in range(num_pl)]
# reading of planetary masses
	f = open(path_to_files + 'gm.txt', 'r')
	s = f.read().split()
	f.close()
	mu = float(s[0])
	gmass = [float(s[i]) for i in range(1, num_pl+2)]
# mass parameters
	mm, am, zo = [], [], []
	zm = [gmass[i]/gmass[0] for i in range(len(gmass))][1:]
	for i in range(num_pl): mm.append(1. + sum(zm[:i+1]))
	zz = [zm[0]/(mm[0]*mu)]			# reduced mass
	zk = [gmass[0]*mm[0]]			# gravitational parameter
	for i in range(1, num_pl): zz.append(zm[i]*mm[i-1]/(mm[i]*mu))
	for i in range(1, num_pl): zk.append(gmass[0]*mm[i]/mm[i-1])
# reading of initial values of averaged elements
	f = open(path_to_files + 'el.txt', 'r')
	s = f.read().split()
	f.close()
	lm, xm, ym, um, vm, qm = [], [], [], [], [], []
	for i in range(num_pl):
		lm.append(float(s[6*i+0]))
		xm.append(float(s[6*i+1]))
		ym.append(float(s[6*i+2]))
		um.append(float(s[6*i+3]))
		vm.append(float(s[6*i+4]))
		qm.append(float(s[6*i+5]))
	for i in range(num_pl):
		am.append(lm[i]**2*(zz[i]**-2)/zk[i])
		zo.append(sqrt(zk[i]*am[i]**-3))
# print of initial values
	print '*** initial longitudes and mean motions of planets ***'
	temp2 = 'l = '
	for i in range(num_pl): temp2 += str('%.16f'%qm[i]).rjust(20, ' ')
	print temp2
	temp2 = 'n = '
	for i in range(num_pl): temp2 += str('%.16f'%zo[i]).rjust(20, ' ')
	print temp2
# input data for the integration
	step = int(raw_input('> length of one step: '))
	step_ch = str(step)
	number = int(raw_input('> number of steps:    '))
	number_ch = str(int(log10(float(number*step))))
# reading series for right hands of the motion equations
	num = [0 for i in range(num_pl)]	
	rdot = [[0 for i in range(NTERMS)] for j in range(num_pl)]
	idot = [[[0 for k in range(4*num_pl)] for i in range(NTERMS)] for j in range(num_pl)]
	for j in range(num_pl):
		f = open(path_to_files + fnames[j], 'r')
		s = f.read().split('\n')
		f.close()
		for k in range(len(s)):
			if k == 0:
				num[j] = int(s[0])
				print fnames[j][:-4], num[j]
			else:
				r = s[k].split()
				for i in range(len(r)):
					if i == 0:
						rdot[j][k] = float(r[0])
					else:
						idot[j][k][i-1] = int(r[i])
# reading of averaged orbital elements (integration results)
	f = open(path_to_files + 'result_1e'+number_ch+'_'+step_ch+prefix+'_pe.txt', 'r')
	s = f.read().split('\n')
	f.close()
	t = []
	adot = [[] for j in range(num_pl)]
	for n in range(number+1):
		li, xi, yi, ui, vi = [], [], [], [], []
		r = s[n+1].split()
		t.append(float(r[0])*SIZE_STEP)
		for i in range(num_pl):
			li.append(float(r[5*i+1]))
			xi.append(float(r[5*i+2]))
			yi.append(float(r[5*i+3]))
			ui.append(float(r[5*i+4]))
			vi.append(float(r[5*i+5]))
		for j in range(num_pl):
			temp2 = 0
			for k in range(num[j]):
				temp = 1.
				for i in range(num_pl):
					temp *=  xi[i]**idot[j][k][0+4*i]*yi[i]**idot[j][k][1+4*i]\
							*ui[i]**idot[j][k][2+4*i]*vi[i]**idot[j][k][3+4*i]
				temp2 = temp2 + rdot[j][k]*temp
			adot[j].append(temp2)
	n = n + 1
# start of the integration process
	f = open(path_to_files + 'result_1e'+number_ch+'_'+step_ch+prefix+'_qe.txt', 'w')
	temp2 = 't'
	for i in range(num_pl): temp2 += ' alpha'+str(i+1)
	f.write(temp2 + '\n')
	temp2 = str(int(t[0]/SIZE_STEP)).rjust(12, ' ')
	for i in range(num_pl): temp2 += '\t'+str('%.16f'%(qm[i]*180./pi)).rjust(21, ' ')
	temp2 += '\n'
	f.write(temp2)
# find B-splines for the interpolation
	tck = [splrep(t, adot[j], k = korder) for j in range(num_pl)]
# the integration from t(1) to t(i)
	for i in range(1, n):
		temp = 2*i*pi
		qi = [splint(t[0], t[i], tck[j]) for j in range(num_pl)]
		for j in range(num_pl): qi[j] += (qm[j] + fmod(zo[j]*(t[i] - t[0]), 2*pi))
		for j in range(num_pl): qi[j] = fmod(qi[j] + temp, 2*pi)
		temp2 = str(int(t[i]/SIZE_STEP)).rjust(12, ' ')
		for j in range(num_pl): temp2 += '\t'+str('%.16f'%(qi[j]*180./pi)).rjust(21, ' ')
		temp2 += '\n'
		f.write(temp2)
	f.close()
########################################################################
def amb(paths, prefix, len_step, num_steps, num_pl = 4, order = 15, int_prec = 1e-11, num_iter_corr = 2, out_prec = 16, is_day = True):
	'''
	The integration of the motion equations
	for the slow variables by Adams-Moulton-Bashforth method.
	'''
	from eps.libs import amb
	path_gm, path_el, path_eqn, path_res = paths
	amb.init.ni = num_iter_corr									# number of iterations for corrector
	amb.init.etol = int_prec									# precision of integration
	amb.init.nor = order										# order of integrator
	amb.init.npl = num_pl										# number of planets
	amb.init.npldot = 4*num_pl									# number of motion equations
	amb.init.step = len_step									# length of one step
	amb.init.number_steps = num_steps							# number of steps
	amb.init.outprec = out_prec									# number of printed digits
	amb.init.isday = 1 if is_day else 0							# units: days or years
	amb.init.path_gm = path_gm+' '*(500 - len(path_gm))			# path to GM-file
	amb.init.path_el = path_el+' '*(500 - len(path_el))			# path to initial elements
	amb.init.path_eqn = path_eqn+' '*(500 - len(path_eqn))		# path to motion equations
	amb.init.path_res = path_res+' '*(500 - len(path_res))		# path to integration results
	amb.init.f_prefix = prefix+' '*(20 - len(prefix))			# file prefix
	amb.amb()
########################################################################
def gbs(paths, prefix, len_step, num_steps, num_pl = 4, order = 15, int_prec = 1e-11, initial_step = 1000, out_prec = 16, is_day = True):
	'''
	The integration of the motion equations
	for the slow variables by Gragg-Bulirsch-Stoer method.
	'''
	from eps.libs import gbs
	path_gm, path_el, path_eqn, path_res = paths
	gbs.init.hx = initial_step									# initial step size (if hx = 0: autoselection)
	gbs.init.etol = int_prec									# precision of integration (if etol <=0: constant step size)
	gbs.init.nor = order										# order of integrator
	gbs.init.npl = num_pl										# number of planets
	gbs.init.npldot = 4*num_pl									# number of motion equations
	gbs.init.step = len_step									# length of one step
	gbs.init.number_steps = num_steps							# number of steps
	gbs.init.outprec = out_prec									# number of printed digits
	gbs.init.isday = 1 if is_day else 0							# units: days or years
	gbs.init.path_gm = path_gm+' '*(500 - len(path_gm))			# path to GM-file
	gbs.init.path_el = path_el+' '*(500 - len(path_el))			# path to initial elements
	gbs.init.path_eqn = path_eqn+' '*(500 - len(path_eqn))		# path to motion equations
	gbs.init.path_res = path_res+' '*(500 - len(path_res))		# path to integration results
	gbs.init.f_prefix = prefix+' '*(20 - len(prefix))			# file prefix
	gbs.gbs()
########################################################################
def gauss_lg(paths, prefix, len_step, num_steps, num_pl = 4, order = 15, int_prec = 1e-11, num_iter_corr = 2, out_prec = 16, is_day = True):
	'''
	The integration of the motion equations
	for the slow variables by Gauss-Everhart method (Legendre Spacings).
	'''
	from eps.libs import gauss_lg
	path_gm, path_el, path_eqn, path_res = paths
	gauss_lg.init.ni = num_iter_corr								# number of iterations for corrector
	gauss_lg.init.etol = int_prec									# precision of integration
	gauss_lg.init.nor = order										# order of integrator
	gauss_lg.init.npl = num_pl										# number of planets
	gauss_lg.init.npldot = 4*num_pl									# number of motion equations
	gauss_lg.init.step = len_step									# length of one step
	gauss_lg.init.number_steps = num_steps							# number of steps
	gauss_lg.init.outprec = out_prec								# number of printed digits
	gauss_lg.init.isday = 1 if is_day else 0						# units: days or years
	gauss_lg.init.path_gm = path_gm+' '*(500 - len(path_gm))		# path to GM-file
	gauss_lg.init.path_el = path_el+' '*(500 - len(path_el))		# path to initial elements
	gauss_lg.init.path_eqn = path_eqn+' '*(500 - len(path_eqn))		# path to motion equations
	gauss_lg.init.path_res = path_res+' '*(500 - len(path_res))		# path to integration results
	gauss_lg.init.f_prefix = prefix+' '*(20 - len(prefix))			# file prefix
	gauss_lg.gauss_lg()
########################################################################
def gauss_rl(paths, prefix, len_step, num_steps, num_pl = 4, order = 15, int_prec = 1e-11, num_iter_corr = 2, out_prec = 16, is_day = True):
	'''
	The integration of the motion equations
	for the slow variables by Gauss-Everhart method (Radau & Lobatto Spacings).
	'''
	from eps.libs import gauss_rl
	path_gm, path_el, path_eqn, path_res = paths
	gauss_rl.init.ni = num_iter_corr								# number of iterations for corrector
	gauss_rl.init.etol = int_prec									# precision of integration
	gauss_rl.init.nor = order										# order of integrator
	gauss_rl.init.npl = num_pl										# number of planets
	gauss_rl.init.npldot = 4*num_pl									# number of motion equations
	gauss_rl.init.step = len_step									# length of one step
	gauss_rl.init.number_steps = num_steps							# number of steps
	gauss_rl.init.outprec = out_prec								# number of printed digits
	gauss_rl.init.isday = 1 if is_day else 0						# units: days or years
	gauss_rl.init.path_gm = path_gm+' '*(500 - len(path_gm))		# path to GM-file
	gauss_rl.init.path_el = path_el+' '*(500 - len(path_el))		# path to initial elements
	gauss_rl.init.path_eqn = path_eqn+' '*(500 - len(path_eqn))		# path to motion equations
	gauss_rl.init.path_res = path_res+' '*(500 - len(path_res))		# path to integration results
	gauss_rl.init.f_prefix = prefix+' '*(20 - len(prefix))			# file prefix
	gauss_rl.gauss_rl()
########################################################################
def everhart(paths, prefix, len_step, num_steps, num_pl = 4, order = 15, int_prec = 11, out_prec = 16, is_day = True):
	'''
	The integration of the motion equations
	for the slow variables by Everhart method.
	'''
	from eps.libs import ever
	path_gm, path_el, path_eqn, path_res = paths
	ever.init.ll = int_prec										# precision of integration
	ever.init.nor = order										# order of integrator
	ever.init.npl = num_pl										# number of planets
	ever.init.npldot = 4*num_pl									# number of motion equations
	ever.init.step = len_step									# length of one step
	ever.init.number_steps = num_steps							# number of steps
	ever.init.outprec = out_prec								# number of printed digits
	ever.init.isday = 1 if is_day else 0						# units: days or years
	ever.init.path_gm = path_gm+' '*(500 - len(path_gm))		# path to GM-file
	ever.init.path_el = path_el+' '*(500 - len(path_el))		# path to initial elements
	ever.init.path_eqn = path_eqn+' '*(500 - len(path_eqn))		# path to motion equations
	ever.init.path_res = path_res+' '*(500 - len(path_res))		# path to integration results
	ever.init.f_prefix = prefix+' '*(20 - len(prefix))			# file prefix
	ever.ever()
########################################################################
def rungekutta(paths, prefix, len_step, num_steps, num_pl = 4, order = 11, out_prec = 16, is_day = True):
	'''
	The integration of the motion equations
	for the slow variables by Runge-Kutta method.
	'''
	from eps.libs import rkutta
	path_gm, path_el, path_eqn, path_res = paths
	rkutta.init.nor = order											# order of integrator
	rkutta.init.npl = num_pl										# number of planets
	rkutta.init.npldot = 4*num_pl									# number of motion equations
	rkutta.init.step = len_step										# length of one step
	rkutta.init.number_steps = num_steps							# number of steps
	rkutta.init.outprec = out_prec									# number of printed digits
	rkutta.init.isday = 1 if is_day else 0							# units: days or years
	rkutta.init.path_gm = path_gm+' '*(500 - len(path_gm))			# path to GM-file
	rkutta.init.path_el = path_el+' '*(500 - len(path_el))			# path to initial elements
	rkutta.init.path_eqn = path_res+' '*(500 - len(path_eqn))		# path to motion equations
	rkutta.init.path_res = path_res+' '*(500 - len(path_res))		# path to integration results
	rkutta.init.f_prefix = prefix+' '*(20 - len(prefix))			# file prefix
	rkutta.rkutta()
########################################################################

