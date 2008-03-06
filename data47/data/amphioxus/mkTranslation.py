#! /users/ldog/muffato/python -OO

import sys
import utils.myTools

dic = {}
for s in sys.argv[1:]:
	f = utils.myTools.myOpenFile(s, 'r')
	for l in utils.myTools.MySQLFileLoader(f):
		t = l.split('\t')
		dic[t[24]] = t[8]
	f.close()

def trans(s):
	if s in dic:
		return dic[s]
	else:
		return s.split('|')[3]

for l in sys.stdin:
	print " ".join([trans(x) for x in l.split()])

