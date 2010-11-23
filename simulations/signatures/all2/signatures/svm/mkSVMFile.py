#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import sys

def transform(c):
	if c == ' ':
		return 0
	if c == '+':
		return 1
	if c == '-':
		return -1
	print >> sys.stderr, "*%s*" % c
	raise ValueError

for l in sys.stdin:
	print ' '.join(["%d:%d" % (i+1,transform(c)) for (i,c) in enumerate(l.replace('\n','')) if c != ' '])

