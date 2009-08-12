#! /users/ldog/muffato/python

import sys

import utils.myFile

genes = set()

f = utils.myFile.openFile(sys.argv[1], 'r')
for l in f:
	t = l[:-1].split('\t')
	genes.add(t[4])

dicVert = {}
for s in sys.argv[2:]:
	f = utils.myFile.openFile(s, 'r')
	for l in f:
		t = l[:-1].split('\t')
		dicVert[t[1]] = t[0]
		assert len(t[0]) > 0
		assert len(t[1]) > 0
	f.close()

def trans(s):
	if s in dicVert:
		return dicVert[s]
	else:
		ll = s.split('|')
		assert len(ll) >= 4, s
		s = "JGI.Brafl1." + s.split('|')[3]
		if s in genes:
			return s
		else:
			return None

for l in sys.stdin:
	newNames = [trans(x) for x in l.split()]
	s = " ".join([x for x in newNames if x != None])
	print s
	if "JGI" not in s:
		print >> sys.stderr, l,

