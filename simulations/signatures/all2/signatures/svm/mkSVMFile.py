#!/usr/bin/env python2

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

