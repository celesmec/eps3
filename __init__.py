# -*- coding: utf-8 -*-
########################################################################
# Copyright (C) 2013-2018 by Alexander Perminov:)
# perminov12@yandex.ru
#
# This program is free software. It uses CAS Piranha
# (F. Biscani, https://github.com/bluescarni/piranha)
########################################################################
# THIS IS PROGRAM FOR MODELING OF PLANETARY SYSTEMS ORBITAL EVOLUTION
########################################################################
# piranha file compression and data format:
from pyranha import data_format as df, compression as cf

# piranha save/load functions:
from pyranha import save_file as sf, load_file as lf

# piranha series types:
from pyranha.types import divisor, divisor_series, monomial, polynomial, poisson_series, int16, rational, double

# Poisson series with rational coefficients and rational degrees:
pt = polynomial[rational,monomial[rational]]()

# Poisson series with double coefficients and rational degrees:
pd = polynomial[double,monomial[rational]]()

# echeloned Poisson series with rational coefficients and degrees, integer coefficients of divisor:
epst = poisson_series[divisor_series[polynomial[rational,monomial[rational]],divisor[int16]]]()

# echeloned Poisson series with double coefficients and rational degrees, integer coefficients of divisor:
epsd = poisson_series[divisor_series[polynomial[double,monomial[rational]],divisor[int16]]]()

del divisor, divisor_series, monomial, polynomial, poisson_series, int16, rational, double

# package modules:
__name__ = 'eps - evolution of planetary systems'
__all__ = ['ham', 'hdm', 'int', 'sol', 'keproc', 'tools', 'pt', 'pd', 'epst', 'epsd', 'cf', 'df', 'lf', 'sf']
__version__ = '3.1.4.16'

print '-'*64
print 'Package ' + __name__
print 'version ' + __version__
print 'modules: ' + str(__all__[:6])[1:-1]
print 'types:   ' + str(__all__[6:-4])[1:-1]
print '-'*64

