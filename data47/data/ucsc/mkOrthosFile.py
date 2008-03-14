#! /users/ldog/muffato/python -OO

import sys
import collections

d = collections.defaultdict(list)

for l in sys.stdin:
	l = l.replace("\n", "")
	t = l.split("\t")
	for x in t[3:]:
		d[x].append(t[0])

f = open(sys.argv[1], "r")
for l in f:
	l = l.replace("\n", "")
	t = l.split("\t")
	res = set()
	for x in t[3:]:
		res.update(d[x])
	if len(res) > 0:
		res.add(t[0])
		print " ".join(res)
f.close()

