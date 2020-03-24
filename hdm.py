# -*- coding: utf-8 -*-
########################################################################
# Copyright (C) 2013-2018 by Alexander Perminov
# perminov12@yandex.ru
#
########################################################################
# THIS IS ALGORITHM FOR HORI-DEPRIT METHOD
########################################################################
import os, time
from eps.config import *
########################################################################
# CONSTRUCTION OF AVERAGED MOTION EQUATIONS
########################################################################
def ham_1(isAver = False, DIR = 'HAM', out_post = ''):
	"""
	"""
	from shutil import copyfile
	PATH_PHI = PREFIX + 'HAM/' + DIR+'/'
	NAME = ['PHI', 'HAM', 'TRG', 'INT']
	out_post = out_post if out_post == '' else '_'+out_post
	for file2 in os.listdir(PATH_PHI):
		order = int(file2[1])
		if order > 0:
			start = time.time()
			series = {NAME[i]: epst(0) for i in range(4)}
			lf(series['PHI'], PATH_PHI + file2, df.boost_portable, cf.bzip2)
			print 'read'
			if file2[0] == 'H':
				for item in series['PHI'].list:
					if item[1] == epst(1):
						series['HAM'] = item[0]*item[1]
				series['TRG'] = series['PHI'] - series['HAM']
				series['INT'] = series['TRG'].t_integrate()
			length = {}
			for key in series.keys():
				key2 = 'HAM' if (file2[0] == 'A' and key == 'PHI') else key
				PATH_SAVE = PREFIX + 'HDM/'+key2+str(order)+'/H'+str(order)+out_post+'/'
				if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
				if isAver:
					if file2[0] == 'A':
						is_write = True if key == 'PHI' else False
					else:
						is_write = True if key != 'HAM' else False
				else:
					is_write = True if file2[0] == 'H' else False
				if is_write:
					if key != 'PHI':
						sf(series[key], PATH_SAVE + key2[0].lower()+str(order)+'_'+file2[3:5]+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
					else:
						copyfile(PATH_PHI + file2, PATH_SAVE + key2[0].lower()+str(order)+'_'+file2[3:5]+'.epst.boostp.bz2')
				length.update({key: sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in series[key].list]])})
			print file2[:-16].ljust(16, ' '),
			for i in range(4): print NAME[i]+':', str(length[NAME[i]]).rjust(8, ' ')+',',
			print str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
########################################################################
def eqn_1(out_post = '', pair = ''):
	"""
	"""
	from pyranha.math import partial
	out_post = out_post if out_post == '' else '_'+out_post
	PATH = PREFIX + 'HDM/HAM1/H1'+out_post+'/'
	for file2 in os.listdir(PATH):
		if pair in file2 or pair == '':
			start = time.time()
			series = epst(0)
			lf(series, PATH + file2, df.boost_portable, cf.bzip2)
			for i in range(len(pq_list)):
				motion = -partial(series, pq_list[i]) if (pq_list[i][0] in 'qyv') else partial(series, pq_list[i])
				motion = motion.trim()
				if motion != epst(0):
					PATH_SAVE = PREFIX + 'HDM/EQN1/H1'+out_post+'/'+qp_list[i]+'/'
					if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
					sf(motion, PATH_SAVE + 'e'+file2[1:-16]+'_'+pq_list[i]+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
				print file2[:-16], qp_list[i], sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in motion.list]]), str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
########################################################################
def chn_1(out_post = '', pair = ''):
	"""
	"""
	from pyranha.math import partial
	out_post = out_post if out_post == '' else '_'+out_post
	PATH = PREFIX + 'HDM/INT1/H1'+out_post+'/'
	for file2 in os.listdir(PATH):
		if pair in file2 or pair == '':
			start = time.time()
			series = epst(0)
			lf(series, PATH + file2, df.boost_portable, cf.bzip2)
			for i in range(len(pq_list)):
				if i in range(4):
					s_prev = epst('s'+str(i)) if i != 0 else 1
					part_diff_nu = -3 * epst('K0')**2 * (epst('m'+str(i+1))**3 * s_prev * epst('s'+str(i+1))**-1) * epst('L'+str(i+1))**-4
					epst.register_custom_derivative(pq_list[i], lambda temp: temp.partial(pq_list[i]) + temp.partial(r'\nu_{'+qp_list[i]+r'}') * part_diff_nu)
					change = partial(series, pq_list[i])
					epst.unregister_all_custom_derivatives()
				else:
					change = -partial(series, pq_list[i]) if (pq_list[i][0] in 'qyv') else partial(series, pq_list[i])
				change = change.trim()
				if change != epst(0):
					PATH_SAVE = PREFIX + 'HDM/CHN1/H1'+out_post+'/'+qp_list[i]+'/'
					if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
					sf(change, PATH_SAVE + 'c'+file2[1:-16]+'_'+pq_list[i]+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
				print file2[:-16], qp_list[i], sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in change.list]]), str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
########################################################################
def d_phi1(out_post = '', pair = ''):
	"""
	"""
	from pyranha.math import partial
	out_post = out_post if out_post == '' else '_'+out_post
	PATH_PHI = PREFIX + 'HDM/PHI1/H1'+out_post+'/'
	for file2 in os.listdir(PATH_PHI):
		if pair in file2 or pair == '':
			start = time.time()
			series = epst(0)
			lf(series, PATH_PHI + file2, df.boost_portable, cf.bzip2)
			for i in range(len(pq_list)):
				diff = partial(series, pq_list[i])
				diff = diff.trim()
				if diff != epst(0):
					PATH_SAVE = PREFIX + 'HDM/D_PHI1/H1'+out_post+'/'+pq_list[i]+'/'
					if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
					sf(diff, PATH_SAVE + file2[:-16]+'_'+pq_list[i]+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
				print file2[:-16], pq_list[i], sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in diff.list]]), str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
########################################################################
def ham_s(max_deg, order = 2, num_item = 1, out_post = '', out_save = '', bracket = '', save_bracket = True):
	"""
	"""
	from pyranha.math import truncate_degree
	out_post = out_post if out_post == '' else '_'+out_post
	out_save = out_save if out_save == '' else '_'+out_save
	if order in range(1, 4) and num_item == 1:
		PATH = PREFIX + 'HDM/HAM'+str(order)+'/H'+str(order)+out_post+'/'
		index = '' if order == 1 else '_'+str(order)
	if order > 1 and num_item > 1:
		if order == 2: pbrackets = ['T1H1']
		if order == 3: pbrackets = ['T1H2', 'T2H1_H2', 'T2H1_T1A1', 'T2H1_T1H1', 'T1T1H1', 'T1T1A1']
		PATH = PREFIX + 'HDM/HAM'+str(order)+'/'+pbrackets[num_item-2]+out_post+'/'
		index = '_'+pbrackets[num_item-2].lower()+out_post.lower()
	PATH_SAVE = PREFIX + 'SUM/HAM/HAM'+str(order)+out_save+'/'
	if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
	print PATH[len(PREFIX):], '->', PATH_SAVE[len(PREFIX):]
	list_dirs = os.listdir(PATH) if num_item != 1 else ['']
	series = epst(0)
	for subdir in list_dirs:
		if (bracket == '') or (bracket == subdir):
			series_temp = epst(0)
			for file2 in os.listdir(PATH + subdir+'/'):
				temp = epst(0)
				lf(temp, PATH + subdir+'/'+file2, df.boost_portable, cf.bzip2)
				if max_deg > 0: temp = truncate_degree(temp, max_deg, pq_list[8:])
				if bracket == '':
					if save_bracket: series_temp += temp
					else:			 series 	 += temp
				if bracket == subdir:
					series_temp += temp
				print '>', file2[:-16]
			if (save_bracket) or (bracket != ''):
				sf(series_temp.trim(), PATH_SAVE + 'h'+str(order)+index+'_('+subdir+').epst.boostp.bz2', df.boost_portable, cf.bzip2)
	if (not save_bracket) and (bracket == ''):
		sf(series.trim(), PATH_SAVE + 'h'+str(order)+index+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
########################################################################
def eqn_s(max_deg, order = 2, out_save = ''):
	"""
	"""
	from pyranha.math import partial, truncate_degree
	out_save = out_save if out_save == '' else '_'+out_save
	PATH = PREFIX + 'SUM/HAM/HAM'+str(order)+out_save+'/'
	for file2 in os.listdir(PATH):
		if index in file2 and file2[0] == 'h':
			series = epst(0)
			lf(series, PATH + file2, df.boost_portable, cf.bzip2)
			for i in range(8, len(pq_list)):
				motion = -partial(series, pq_list[i]) if (pq_list[i][0] in 'qyv') else partial(series, pq_list[i])
				if max_deg > 0:
					motion = truncate_degree(motion, max_deg, pq_list[8:])
				PATH_SAVE = PREFIX + 'SUM/EQN/EQN'+str(order)+out_save+'/'
				if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
				sf(motion.trim(), PATH_SAVE + 'e'+file2[1:-16]+'_'+qp_list[i]+'.epst.boostp.bz2', df.boost_portable, cf.bzip2)
########################################################################
# CONSTRUCTION OF CHANGE-FUNCTIONS
########################################################################
def int_N(N, num_item = 1, t_order_max = 40, out_post = '', isDouble = False):
	"""
	"""
	if N == 2: BRACKETS = ['T1H1', 'T1A1']
	if N == 3: BRACKETS = ['T1H2', 'T1A2', 'T2H1_H2', 'T2H1_T1A1', 'T2H1_T1H1', 'T2A1_H2', 'T2A1_T1A1', 'T2A1_T1H1', 'T1T1H1', 'T1T1A1']
	PATH_PHI = PREFIX + 'HDM/PHI'+str(N)+'/'+BRACKETS[num_item-1]+out_post+'_TOR'+str(t_order_max).rjust(2, '0')+'/'
	NAME = ['PHI', 'HAM', 'TRG', 'INT']
	(epst2, type2) = (epsd, 'd') if isDouble else (epst, 't')
	for subdir in os.listdir(PATH_PHI):
		for file2 in os.listdir(PATH_PHI + subdir+'/'):
			start = time.time()
			order = int(file2[1])
			series = {NAME[i]: epst2(0) for i in range(4)}
			lf(series['PHI'], PATH_PHI + subdir+'/'+file2, df.boost_portable, cf.bzip2)
			for item in series['PHI'].list:
				if item[1] == epst2(1):
					series['HAM'] = item[0]*item[1]
			series['TRG'] = series['PHI'] - series['HAM']
			series['INT'] = series['TRG'].t_integrate()
			PATH_SAVE = PREFIX + 'HDM/INT'+str(N)+'/'+BRACKETS[num_item-1]+out_post+'_TOR'+str(t_order_max).rjust(2, '0')+'/'+subdir+'/'
			if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
			sf(series['INT'], PATH_SAVE + 'i'+str(order)+file2[2:-16]+'.eps'+type2+'.boostp.bz2', df.boost_portable, cf.bzip2)
			print subdir+'/'+file2[:-16].ljust(12, ' '), str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
########################################################################
def chn_N(N, num_item = 1, t_order_max = 40, out_post = '', isDouble = False):
	"""
	"""
	from pyranha.math import partial
	if N == 2: BRACKETS = ['T1H1', 'T1A1']
	if N == 3: BRACKETS = ['T1H2', 'T1A2', 'T2H1_H2', 'T2H1_T1A1', 'T2H1_T1H1', 'T2A1_H2', 'T2A1_T1A1', 'T2A1_T1H1', 'T1T1H1', 'T1T1A1']
	PATH_INT = PREFIX + 'HDM/INT'+str(N)+'/'+BRACKETS[num_item-1]+out_post+'_TOR'+str(t_order_max).rjust(2, '0')+'/'
	(type2, epst2) = ('d', epsd) if isDouble else ('t', epst)
	for subdir in os.listdir(PATH_INT):
		for file2 in os.listdir(PATH_INT + subdir+'/'):
				start = time.time()
				series = epst2(0)
				lf(series, PATH_INT + subdir+'/'+file2, df.boost_portable, cf.bzip2)
				for i in range(len(pq_list)):
					if i in range(4):
						s_prev = epst2('s'+str(i)) if i != 0 else 1
						part_diff_nu = -3 * epst2('K0')**2 * (epst2('m'+str(i+1))**3 * s_prev * epst2('s'+str(i+1))**-1) * epst2('L'+str(i+1))**-4
						epst2.register_custom_derivative(pq_list[i], lambda temp: temp.partial(pq_list[i]) + temp.partial(r'\nu_{'+qp_list[i]+r'}') * part_diff_nu)
						change = partial(series, pq_list[i])
						epst2.unregister_all_custom_derivatives()
					else:
						change = -partial(series, pq_list[i]) if (pq_list[i][0] in 'qyv') else partial(series, pq_list[i])
					if change != epst2(0):
						PATH_SAVE = PREFIX + 'HDM/CHN'+str(N)+'/'+BRACKETS[num_item-1]+out_post+'_TOR'+str(t_order_max).rjust(2, '0')+'/'+qp_list[i]+'/'+subdir+'/'
						if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
						sf(change.trim(), PATH_SAVE + 'c'+file2[1:-16]+'_'+pq_list[i]+'.eps'+type2+'.boostp.bz2', df.boost_portable, cf.bzip2)
				print PATH_INT+subdir+'/'+file2[:-16], str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
########################################################################
def int_s(max_deg, order = 1, num_item = 1, t_order_max = 15, prefix = '', save_bracket = True, isDouble = False):
	"""
	The summation of terms of the generating function.
	"""
	from pyranha.math import truncate_degree
	if order == 1:
		PATH, index = PREFIX + 'HDM/INT1/H1/', ''
	if order > 1:
		pbrackets = ['T1H1', 'T1A1'] if order == 2 else ['T1H2', 'T1A2', 'T2H1_H2', 'T2H1_T1A1', 'T2H1_T1H1', 'T2A1_H2', 'T2A1_T1A1', 'T2A1_T1H1', 'T1T1H1', 'T1T1A1']
		if num_item == 1: PATH, index = PREFIX + 'HDM/INT'+str(order)+'/H'+str(order)+'/', '_'+str(order)
		if num_item > 1:
			PATH = PREFIX + 'HDM/INT'+str(order)+'/'+pbrackets[num_item-2]+'_TOR'+str(t_order_max).rjust(2, '0')+'/'
			index = '_'+pbrackets[num_item-2].lower()+'_tor'+str(t_order_max).rjust(2, '0')
	(type2, epst2) = ('d', epsd) if isDouble else ('t', epst)
	PATH_SAVE = PREFIX + 'SUM/INT/INT'+str(order)+prefix+'/'
	if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
	print PATH[len(PREFIX):], '->', PATH_SAVE[len(PREFIX):]
	list_dirs = os.listdir(PATH) if num_item != 1 else ['']
	series = epst2(0)
	for subdir in list_dirs:
		series_temp = epst2(0)
		for file2 in os.listdir(PATH + subdir+'/'):
			temp = epst2(0)
			lf(temp, PATH + subdir+'/'+file2, df.boost_portable, cf.bzip2)
			if max_deg > 0: temp = truncate_degree(temp, max_deg, pq_list[8:])
			if save_bracket: series_temp += temp
			else:			 series 	 += temp
			print '>', file2[:-16]
		if save_bracket:
			sf(series_temp.trim(), PATH_SAVE + 'i'+str(order)+index+'_('+subdir+').eps'+type2+'.boostp.bz2', df.boost_portable, cf.bzip2)
	if not save_bracket:
		sf(series.trim(), PATH_SAVE + 'i'+str(order)+index+'.eps'+type2+'.boostp.bz2', df.boost_portable, cf.bzip2)
########################################################################
def chn_s(max_deg, order = 1, num_item = 1, t_order_max = 15, prefix = '', save_bracket = True, isDouble = False):
	"""
	The summation of terms of the change-functions.
	"""
	from pyranha.math import truncate_degree
	if order == 1:
		PATH, index = PREFIX + 'HDM/CHN1/H1'+prefix+'/', ''
	if order > 1:
		pbrackets = ['T1H1', 'T1A1'] if order == 2 else ['T1H2', 'T1A2', 'T2H1_H2', 'T2H1_T1A1', 'T2H1_T1H1', 'T2A1_H2', 'T2A1_T1A1', 'T2A1_T1H1', 'T1T1H1', 'T1T1A1']
		if num_item == 1: PATH, index = PREFIX + 'HDM/CHN'+str(order)+'/H'+str(order)+'/', '_'+str(order)
		if num_item > 1:
			PATH = PREFIX + 'HDM/CHN'+str(order)+'/'+pbrackets[num_item-2]+'_TOR'+str(t_order_max).rjust(2, '0')+'/'
			index = '_'+pbrackets[num_item-2].lower()+'_tor'+str(t_order_max).rjust(2, '0')
	(type2, epst2) = ('d', epsd) if isDouble else ('t', epst)
	PATH_SAVE = PREFIX + 'SUM/CHN/CHN'+str(order)+prefix+'/'
	if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
	print PATH[len(PREFIX):], '->', PATH_SAVE[len(PREFIX):]
	for i in range(len(pq_list)):
		list_dirs = os.listdir(PATH + qp_list[i]+'/') if num_item != 1 else ['']
		series = epst2(0)
		for subdir in list_dirs:
			series_temp = epst2(0)
			for file2 in os.listdir(PATH + qp_list[i]+'/'+subdir+'/'):
				temp = epst2(0)
				lf(temp, PATH + qp_list[i]+'/'+subdir+'/'+file2, df.boost_portable, cf.bzip2)
				if max_deg > 0: temp = truncate_degree(temp, max_deg, pq_list[8:])
				if save_bracket: series_temp += temp
				else:			 series 	 += temp
				print '>', file2[:-16]
			if save_bracket:
				sf(series_temp.trim(), PATH_SAVE + 'c'+str(order)+index+'_('+subdir+')_'+qp_list[i]+'.eps'+type2+'.boostp.bz2', df.boost_portable, cf.bzip2)
		if not save_bracket:
			sf(series.trim(), PATH_SAVE + 'c'+str(order)+index+'_'+qp_list[i]+'.eps'+type2+'.boostp.bz2', df.boost_portable, cf.bzip2)
########################################################################
# CONSTRUCTION OF POISSON BRACKETS
########################################################################
def pb_t1hN(N = 1, m = 6, pqi = 0, t_order_max = 40, in1_post = '', in2_post = '', out_post = '', isAver = True, isDouble = False, isPrint = True):
	"""
	Poisson brackets {T_1, H_1} and {T_1, H_2}.
	"""
	from pyranha.math import cos, sin, degree, ldegree, t_degree, t_ldegree, t_order, t_lorder, truncate_degree
	in1_post = in1_post if in1_post == '' else '_'+in1_post
	in2_post = in2_post if in2_post == '' else '_'+in2_post
	out_post = out_post if out_post == '' else '_'+out_post
	PATH_CHN = PREFIX + 'HDM/CHN1/H1'+in1_post+'/'
	PATH_PHI = PREFIX + 'HDM/D_PHI'+str(N)+'/H'+str(N)+in2_post+'/'
	(c, one, max_m, type2, pt2, epst2) = (0.5, 1., m, 'd', pd, epsd) if isDouble else (F(1,2), 1, F(m), 't', pt, epst)
	pairs = [['21_31', '31_21', '21_41', '41_21', '31_41', '41_31'], ['21_32', '32_21', '21_42', '42_21', '32_42', '42_32'],\
			 ['31_32', '32_31', '31_43', '43_31', '32_43', '43_32'], ['41_42', '42_41', '41_43', '43_41', '42_43', '43_42']]
	for i in range(len(pq_list))[pqi:pqi+1]:
		for file1 in os.listdir(PATH_CHN + qp_list[i]):
			diff1 = epst(0)
			lf(diff1, PATH_CHN + qp_list[i]+'/'+file1, df.boost_portable, cf.bzip2)
			diff1 = truncate_degree(diff1, F(m), pq_list[8:])
			new_diff1 = epst2(0)
			for item1 in diff1.list:
				if t_order(item1[1]) <= t_order_max:
					new_diff1 += one*(item1[0]*item1[1])
			if isPrint: print file1[:-16]#, '->', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in new_diff1.list]]), len(new_diff1.list)
			for file2 in os.listdir(PATH_PHI + qp_list[i]):
				TYPE = 'HAM' if isAver else 'PHI'
				PATH_SAVE = PREFIX + 'HDM/'+TYPE+str(N+1)+'/T1H'+str(N)+out_post+'_TOR'+str(t_order_max).rjust(2, '0')+'/'+pq_list[i]+'_'+qp_list[i]+'/'
				FILE_SAVE = '/'+TYPE[0].lower()+str(N+1)+'_'+file1[3:5]+'_'+file2[3:5]+'.eps'+type2+'.boostp.bz2'
				start = time.time()
				if not os.path.exists(PATH_SAVE + FILE_SAVE):
					diff2 = epst(0)
					lf(diff2, PATH_PHI + qp_list[i]+'/'+file2, df.boost_portable, cf.bzip2)
					diff2 = truncate_degree(diff2, F(m), pq_list[8:])
					new_diff2 = epst2(0)
					for item1 in diff2.list:
						if t_order(item1[1]) <= t_order_max:
							new_diff2 += one*(item1[0]*item1[1])
					if isPrint: print file2[:-16]#, '->', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in new_diff2.list]]), len(new_diff2.list)
					if isAver:
						count = 0
						new_result = epst2(0)
						pt2.set_auto_truncate_degree(max_m, pq_list[8:])
						for item in new_diff1.list:
							for jtem in new_diff2.list:
								trig = item[1]*jtem[1]
								for ktem in trig.list:
									if ktem[1] == 1:
										new_result = new_result + c*item[0]*jtem[0]*ktem[0]
							count += 1
							print count, item[1]
						pt2.unset_auto_truncate_degree()
					else:
						#count = 0
						#result = epst2(0)
						#pt2.set_auto_truncate_degree(max_m, pq_list[8:])
						#for item in new_diff1.list:
							#temp = c*item[0]*item[1]*new_diff2
							#new_temp = epst2(0)
							#for item1 in temp.list:
								#if t_order(item1[1]) <= t_order_max:
									#new_temp += item1[0]*item1[1]
							#result = result + new_temp
							#count += 1
							#print count, item[1]
						result = c*new_diff1*new_diff2
						pt2.unset_auto_truncate_degree()
					if not isAver:
						new_result = epst2(0)
						for item1 in result.list:
							if t_order(item1[1]) <= t_order_max:
								new_result += item1[0]*item1[1]
					if isPrint: print 'done'#, '->", sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in new_result.list]]), len(new_result.list)
					if new_result != epst2(0):						
						if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
						sf(new_result.trim(), PATH_SAVE + FILE_SAVE, df.boost_portable, cf.bzip2)
				print 'bracket:', pq_list[i]+'_'+qp_list[i], file1[:-16], file2[:-16], str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
########################################################################
def pb_t1a1(m = 6, pqi = 0, t_order_max = 40, in1_post = '', in2_post = '', out_post = '', isDouble = False, isPrint = True):
	"""
	Poisson bracket {T_1, A_1}.
	"""
	from pyranha.math import cos, sin, degree, ldegree, t_degree, t_ldegree, t_order, t_lorder, truncate_degree
	in1_post = in1_post if in1_post == '' else '_'+in1_post
	in2_post = in2_post if in2_post == '' else '_'+in2_post
	out_post = out_post if out_post == '' else '_'+out_post
	PATH_CHN = PREFIX + 'HDM/CHN1/H1'+in1_post+'/'
	PATH_EQN = PREFIX + 'HDM/EQN1/H1'+in2_post+'/'
	(c, one, max_m, type2, pt2, epst2) = (0.5, 1., m, 'd', pd, epsd) if isDouble else (F(1,2), 1, F(m), 't', pt, epst)
	for i in range(len(pq_list))[pqi:pqi+1]:
		if pq_list[i][0] != 'L':
			for file1 in os.listdir(PATH_CHN + qp_list[i]):
				diff1 = epst(0)
				lf(diff1, PATH_CHN + qp_list[i]+'/'+file1, df.boost_portable, cf.bzip2)
				diff1 = truncate_degree(diff1, F(m), pq_list[8:])
				new_diff1 = epst2(0)
				for item1 in diff1.list:
					if t_order(item1[1]) <= t_order_max:
						new_diff1 += one*(item1[0]*item1[1])
				if isPrint: print file1[:-16]+':', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in new_diff1.list]]), len(new_diff1.list)
				for file2 in os.listdir(PATH_EQN + pq_list[i]):
					start = time.time()
					diff2 = epst(0)
					lf(diff2, PATH_EQN + pq_list[i]+'/'+file2, df.boost_portable, cf.bzip2)
					diff2 = one*diff2
					if isPrint: print file2[:-16]+':', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in diff2.list]]), len(diff2.list)
					pt2.set_auto_truncate_degree(max_m, pq_list[8:])
					result = c*new_diff1*diff2 if (pq_list[i][0] in 'qyv') else -c*new_diff1*diff2
					pt2.unset_auto_truncate_degree()
					if isPrint: print 'result:', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in result.list]]), len(result.list)
					if result != epst2(0):
						PATH_SAVE = PREFIX + 'HDM/PHI2/T1A1'+out_post+'_TOR'+str(t_order_max).rjust(2, '0')+'/'+pq_list[i]+'_'+qp_list[i]+'/'
						if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
						sf(result.trim(), PATH_SAVE + '/p2'+'_'+file1[3:5]+'_'+file2[3:5]+'.eps'+type2+'.boostp.bz2', df.boost_portable, cf.bzip2)
					print 'bracket:', pq_list[i]+'_'+qp_list[i], file1[:-16], file2[:-16], str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
########################################################################
def pb_t1a2(): # empty
	"""
	Poisson bracket {T_1, A_2}.
	"""
	pass
##################################################
def pb_t2h1(m = 4, type_br = 0, num_br = 0, num_subr = 0, t_order_max = 20, in1_post = '', in2_post = '', out_post = '', isAver = True, isDouble = False, isPrint = True):
	"""
	Poisson bracket {T_2, H_1}.
	"""
	from pyranha.math import cos, sin, degree, ldegree, t_degree, t_ldegree, t_order, t_lorder, truncate_degree
	in1_post = in1_post if in1_post == '' else '_'+in1_post
	in2_post = in2_post if in2_post == '' else '_'+in2_post
	out_post = out_post if out_post == '' else '_'+out_post
	PBRACKETS = ['H2', 'T1H1_TOR'+str(t_order_max).rjust(2, '0'), 'T1A1_TOR'+str(t_order_max).rjust(2, '0')]
	PATHS_CHN = [PREFIX + 'HDM/CHN2/'+ITEM+in1_post+'/' for ITEM in PBRACKETS]
	PATH_PHI = PREFIX + 'HDM/D_PHI1/H1'+in2_post+'/'
	(c, one, max_m, type2, pt2, epst2) = (0.5, 1., m, 'd', pd, epsd) if isDouble else (F(1,2), 1, F(m), 't', pt, epst)
	for j in range(len(PATHS_CHN))[type_br:type_br+1]:
		for i in range(len(pq_list))[num_br:num_br+1]:
			listdir = os.listdir(PATHS_CHN[j] + qp_list[i]+'/') if not ('H2' in PATHS_CHN[j]) else ['']
			TYPE = 'HAM' if isAver else 'PHI'
			PATH_SAVE = PREFIX + 'HDM/'+TYPE+'3/T2H1_'+PBRACKETS[j]+out_post+'/'+pq_list[i]+'_'+qp_list[i]+'/'
			for subdir in listdir[num_subracket:num_subracket+1]:
				for file1 in os.listdir(PATHS_CHN[j] + qp_list[i]+'/'+subdir+'/'):
					#diff1 = epst2(0)
					#lf(diff1, PATHS_CHN[j] + qp_list[i]+'/'+subdir+'/'+file1, df.boost_portable, cf.bzip2)
					#diff1 = truncate_degree(diff1, F(m), pq_list[8:])
					#new_diff1 = epst2(0)
					#for item1 in diff1.list:
					#	if t_order(item1[1]) <= t_order_max:
					#		new_diff1 += one*(item1[0]*item1[1])
					new_diff1 = epst2(0)
					lf(new_diff1, PATHS_CHN[j] + qp_list[i]+'/'+subdir+'/'+file1, df.boost_portable, cf.bzip2)
					new_diff1 = truncate_degree(new_diff1, F(m), pq_list[8:])
					if isPrint: print file1[:-16]#, '->', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in new_diff1.list]]), len(new_diff1.list)
					for file2 in os.listdir(PATH_PHI + qp_list[i]+'/'):
						FILE_CHECK = PATH_SAVE + '/'+TYPE[0].lower()+'3_('+subdir+')_('+file1[3:8]+')_'+file2[3:5]+'.eps'+type2+'.boostp.bz2'
						if not os.path.isfile(FILE_CHECK):
							start = time.time()
							#diff2 = epst2(0)
							#lf(diff2, PATH_PHI + qp_list[i]+'/'+file2, df.boost_portable, cf.bzip2)
							#diff2 = truncate_degree(diff2, F(m), pq_list[8:])
							#new_diff2 = epst2(0)
							#for item1 in diff2.list:
								#if t_order(item1[1]) <= t_order_max:
									#new_diff2 += one*(item1[0]*item1[1])
							new_diff2 = epst2(0)
							lf(new_diff2, PATH_PHI + qp_list[i]+'/'+file2, df.boost_portable, cf.bzip2)
							new_diff2 = truncate_degree(new_diff2, F(m), pq_list[8:])
							if isPrint: print file2[:-16]#, '->', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in new_diff2.list]]), len(new_diff2.list)
							pt2.set_auto_truncate_degree(max_m, pq_list[8:])
							if isAver:
								count = 0
								new_result = epst2(0)
								for item in new_diff1.list:
									for jtem in new_diff2.list:
										trig = item[1]*jtem[1]
										for ktem in trig.list:
											if ktem[1] == 1:
												new_result = new_result + c*item[0]*jtem[0]*ktem[0]
									count += 1
									print count, item[1]
							else:
								result = c*new_diff1*new_diff2
							pt2.unset_auto_truncate_degree()
							if not isAver:
								new_result = epst2(0)
								for item1 in result.list:
									if t_order(item1[1]) <= t_order_max:
										new_result += item1[0]*item1[1]
							if isPrint: print 'done'#, '->', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in new_result.list]]), len(new_result.list)
							if new_result != epst2(0):
								if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
								sf(new_result.trim(), PATH_SAVE + '/'+TYPE[0].lower()+'3_('+subdir+')_('+file1[3:8]+')_'+file2[3:5]+'.eps'+type2+'.boostp.bz2', df.boost_portable, cf.bzip2)
							print 'bracket:', pq_list[i]+'_'+qp_list[i], file1[:-16], file2[:-16], str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
########################################################################
def pb_t2a1(): # empty
	"""
	Poisson bracket {T_2, A_1}.
	"""
	pass
##################################################
def pb_t1t1f1(f = 0, m = 4, num_br = 0, num_subr = 0, t_order_max = 20, in1_post = '', in2_post = '', out_post = '', isAver = True, isDouble = False, isPrint = True):
	"""
	Poisson brackets {T_1, {T_1, H_1}} and {T_1, {T_1, A_1}}.
	"""
	from pyranha.math import cos, sin, degree, ldegree, t_degree, t_ldegree, t_order, t_lorder, partial, truncate_degree
	in1_post = in1_post if in1_post == '' else '_'+in1_post
	in2_post = in2_post if in2_post == '' else '_'+in2_post
	out_post = out_post if out_post == '' else '_'+out_post
	if f == 0: c = -1./6. if isDouble else F(-1,6)
	if f == 1: c =  1./6. if isDouble else F( 1,6)
	TYPE2 = 'A' if f == 0 else 'H'
	PATH_CHN = PREFIX + 'HDM/CHN1/H1'+in1_post+'/'
	PATH_BRA = PREFIX + 'HDM/PHI2/T1'+TYPE2+'1'+in2_post+'_TOR'+str(t_order_max).rjust(2, '0')+'/'
	(one, max_m, type2, pt2, epst2) = (1., m, 'd', pd, epsd) if isDouble else (1, F(m), 't', pt, epst)
	for i in range(len(pq_list))[num_bracket:num_bracket+1]:
		TYPE = 'HAM' if isAver else 'PHI'
		PATH_SAVE = PREFIX + 'HDM/'+TYPE+'3/T1T1'+TYPE2+'1'+out_post+'_TOR'+str(t_order_max).rjust(2, '0')+'/'+pq_list[i]+'_'+qp_list[i]+'/'
		for file1 in os.listdir(PATH_CHN + qp_list[i]): # diff L,q (chn q,L)
			diff1 = epst(0)
			lf(diff1, PATH_CHN + qp_list[i]+'/'+file1, df.boost_portable, cf.bzip2)
			diff1 = truncate_degree(diff1, F(m), pq_list[8:])
			new_diff1 = epst2(0)
			for item1 in diff1.list:
				if t_order(item1[1]) <= t_order_max:
					new_diff1 += one*(item1[0]*item1[1])
			#new_diff1 = epst(0)
			#lf(new_diff1, PATH_CHN + qp_list[i]+'/'+file1, df.boost_portable, cf.bzip2)
			#new_diff1 = truncate_degree(new_diff1, F(m), pq_list[8:])
			if isPrint: print file1[:-16]#+':', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in new_diff1.list]]), len(new_diff1.list)
			for subdir in os.listdir(PATH_BRA)[num_subr:num_subr+1]:
				for file2 in os.listdir(PATH_BRA + subdir+'/'):
					CHECK_FILE = PATH_SAVE + '/'+TYPE[0].lower()+'3_'+file1[3:5]+'_('+subdir+'_'+file2[3:8]+').eps'+type2+'.boostp.bz2'
					if not os.path.isfile(CHECK_FILE):
						if '44' in file2[3:8]:
							num_list = range(1, 5)
						else:
							num_list = range(int(file2[4:5]), int(file2[3:4])+1)
							num_list += range(int(file2[7:8]), int(file2[6:7])+1)
						if int(pq_list[i][1]) in num_list:
							start = time.time()
							diff2 = epst(0)
							lf(diff2, PATH_BRA + subdir+'/'+file2, df.boost_portable, cf.bzip2)
							diff2 = truncate_degree(diff2, F(m), pq_list[8:])
							new_diff2 = epst2(0)
							for item1 in diff2.list:
								if t_order(item1[1]) <= t_order_max:
									new_diff2 += one*(item1[0]*item1[1])
							#new_diff2 = epst(0)
							#lf(new_diff2, PATH_BRA + subdir+'/'+file2, df.boost_portable, cf.bzip2)
							#new_diff2 = truncate_degree(new_diff2, F(m), pq_list[8:])
							if i in range(4, 8):
								s_prev = epst2('s'+str(i-4)) if i != 4 else 1
								part_diff_nu = -3 * epst2('K0')**2 * (epst2('m'+str(i-4+1))**3 * s_prev * epst2('s'+str(i-4+1))**-1) * epst2('L'+str(i-4+1))**-4
								epst2.register_custom_derivative(qp_list[i], lambda temp: temp.partial(qp_list[i]) + temp.partial(r'\nu_{'+pq_list[i]+r'}') * part_diff_nu)
								new_diff2 = partial(new_diff2, qp_list[i])
								epst2.unregister_all_custom_derivatives()
							else:
								new_diff2 = partial(new_diff2, qp_list[i])
							if isPrint: print subdir, file2[:-16]#+':', sum([len(item.list) for item in [jtem[0].list[0][0] for jtem in new_diff2.list]]), len(new_diff2.list)
							pt2.set_auto_truncate_degree(max_m, pq_list[8:])
							if isAver:
								count = 0
								new_result = epst2(0)
								for item in new_diff1.list:
									for jtem in new_diff2.list:
										trig = item[1]*jtem[1]
										for ktem in trig.list:
											if ktem[1] == 1:
												temp_result = c*(item[0]*jtem[0]*ktem[0])
												new_result = new_result + temp_result
									print c, item[1]
									count += 1
							else:
								result = c*new_diff1*new_diff2
							pt2.unset_auto_truncate_degree()
							if not isAver:
								new_result = epst2(0)
								for item1 in result.list:
									if t_order(item1[1]) <= t_order_max:
										new_result += item1[0]*item1[1]
							if new_result != epst2(0):
								if not os.path.exists(PATH_SAVE): os.makedirs(PATH_SAVE)
								sf(new_result.trim(), PATH_SAVE + '/'+TYPE[0].lower()+'3_'+file1[3:5]+'_('+subdir+'_'+file2[3:8]+').eps'+type2+'.boostp.bz2', df.boost_portable, cf.bzip2)
							print 'bracket:', pq_list[i]+'_'+qp_list[i], file1[:-16], file2[:-16], str('%.1f'%(time.time()-start)).rjust(4, ' '), 's'
########################################################################
# EVALUATION OF CONSTRUCTED SERIES
########################################################################
def evl_s(elmass = [[], [], 0], num_pl = 4, order = 2, max_n = 0, func = 'eqn', prefix = '', is_q = False):
	"""
	Evaluation of the constructed series.
	"""
	from mpmath import mpf
	from pyranha.math import evaluate, truncate_degree
	from eps.tools import jacobiAMass
	# -------- coordinates & masses --------
	gmj, mpj, km = jacobiAMass(elmass[1], elmass[2])
	mp = [1./elmass[2]] + [elmass[1][i]/elmass[1][0]/elmass[2] for i in range(len(elmass[1]))][1:]
	sp = [1. + elmass[2]*sum(mp[1:i+1]) for i in range(len(elmass[1]))]
	# -------- dictionaries --------
	evaldict = {'K0': mpf(elmass[1][0]),\
				'L1': mpf(elmass[0][0][0]), 'm1': mpf(mp[1]), 's1': mpf(sp[1]), r'\nu_{q1}': mpf(gmj[0]**2*mpj[0]**3*elmass[0][0][0]**-3),\
				'L2': mpf(elmass[0][1][0]), 'm2': mpf(mp[2]), 's2': mpf(sp[2]), r'\nu_{q2}': mpf(gmj[1]**2*mpj[1]**3*elmass[0][1][0]**-3),\
				'L3': mpf(elmass[0][2][0]), 'm3': mpf(mp[3]), 's3': mpf(sp[3]), r'\nu_{q3}': mpf(gmj[2]**2*mpj[2]**3*elmass[0][2][0]**-3),\
				'L4': mpf(elmass[0][3][0]), 'm4': mpf(mp[4]), 's4': mpf(sp[4]), r'\nu_{q4}': mpf(gmj[3]**2*mpj[3]**3*elmass[0][3][0]**-3),\
				'x1': mpf(elmass[0][0][1]), 'y1': mpf(elmass[0][0][2]), 'u1': mpf(elmass[0][0][3]), 'v1': mpf(elmass[0][0][4]), 'q1': mpf(elmass[0][0][5]),\
				'x2': mpf(elmass[0][1][1]), 'y2': mpf(elmass[0][1][2]), 'u2': mpf(elmass[0][1][3]), 'v2': mpf(elmass[0][1][4]), 'q2': mpf(elmass[0][1][5]),\
				'x3': mpf(elmass[0][2][1]), 'y3': mpf(elmass[0][2][2]), 'u3': mpf(elmass[0][2][3]), 'v3': mpf(elmass[0][2][4]), 'q3': mpf(elmass[0][2][5]),\
				'x4': mpf(elmass[0][3][1]), 'y4': mpf(elmass[0][3][2]), 'u4': mpf(elmass[0][3][3]), 'v4': mpf(elmass[0][3][4]), 'q4': mpf(elmass[0][3][5])}
	for i in range(num_pl, 4): evaldict.update({'m'+str(i+1): 0.})
	# -------- loop --------
	return_temp = {}
	for file2 in os.listdir(PREFIX + 'SUM/'+func.upper()+'/'+func.upper()+str(order)+prefix+'/'):
		if func[0] in file2[0]:
			if order == 1:
				c = (elmass[2]*elmass[1][0])
			if order == 2:
				c = elmass[1][0]*elmass[2]**order if file2[3] == str(order) else (elmass[1][0]*elmass[2])**order
			if order == 3:
				c = elmass[1][0]**2*elmass[2]**order if (('_t1h2' in file2[:-16]) or ('_t2h1_h2' in file2[:-16])) else (elmass[1][0]*elmass[2])**order
			series = epst(0)
			if is_q or (not is_q and file2.split('.')[0].split('_')[-1][0] != 'q'):
				lf(series, PREFIX + 'SUM/'+func.upper()+'/'+func.upper()+str(order)+prefix+'/'+file2, df.boost_portable, cf.bzip2)
				#if max_n > 0:
				#	series = truncate_degree(series, max_n, pq_list[8:])
			if num_pl in [2, 3]:
				for key3 in ['m1', 'm2', 'm3', 'm4']:
					series = series.subs(key3, F(evaldict[key3]))
				series = series.trim()
			eval_series = c*evaluate(series, evaldict)
			sum_items = 0
			for item in series.list:
				for jtem in item[0].list:
					sum_items += len(jtem[0].list)
			print file2[:-16]+':', '%.16e'%eval_series, sum_items
			return_temp.update({file2.split('.')[0]: eval_series})
	return return_temp
########################################################################

