# -*- coding: utf-8 -*-
########################################################################
# Copyright (C) 2013-2018 by Alexander Perminov
# perminov12@yandex.ru
#
########################################################################
# THESE ARE GLOBAL VARIABLES AND FUNCTIONS
########################################################################
from math import pi
from fractions import Fraction as F
from eps import pt, epst, pd, epsd, cf, df, sf, lf

# Settings:
NUM_PL = 4					# The number of planets
IS_DAY = True				# What is time unit? Day (True) or year (False)
PROBLEM = 'SIX_FULL'		# Working subdirectory

# NUM_PL = int(raw_input('Input number of planets (<= 4): '))
# IS_DAY = bool(raw_input('Time units (0 - year, 1 - day): '))
# PROBLEM = str(raw_input('Input name of problem:          '))

# Path to working directory:
PREFIX = '/media/alex/DATA/EPS/'+PROBLEM+'/'

# Used lists of orbital elements:
pq_list = [item+str(i+1) for item in 'Lqxyuv' for i in range(NUM_PL)]
qp_list = [item+str(i+1) for item in 'qLyxvu' for i in range(NUM_PL)]
names_ke = [item+str(i+1) for i in range(NUM_PL) for item in ['a', 'e', 'i', 'om', 'Om']]
names_pe = [item+str(i+1) for i in range(NUM_PL) for item in ['L', 'x', 'y', 'u', 'v']]

# Used units:
RHO, DAY, YEAR, AU = pi/180., 86400.0, 31557600.0, 149597870.700

# Time unit (day or year):
TIME_U = DAY if IS_DAY else YEAR

# Gravitational parameters of Sun and Jupiter:
GMS, GMJ = 132712440041.93938*TIME_U**2*AU**-3, 126686511.*TIME_U**2*AU**-3

