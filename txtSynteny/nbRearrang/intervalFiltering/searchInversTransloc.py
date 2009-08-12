#! /users/ldog/muffato/python

import sys
import collections

human = set()
dicHA = {}
ancestr = set()
dicAH = {}


for l in sys.stdin:
	t = l.split()
	if t[6] == "SAME_ORIENT":
		continue
	if t[6] == "TWO_ENDS":
		continue
	if t[6] == "DUPLICATES_SAME_ORIENT":
		continue
	if t[0] == "+":
		human.add((t[1],t[2]))
		human.add((t[2],t[1]))
		dicHA[t[1]] = t[3]
		dicHA[t[2]] = t[4]
	elif t[0] == "-":
		ancestr.add((t[1],t[2]))
		ancestr.add((t[2],t[1]))
		dicAH[t[1]] = t[3]
		dicAH[t[2]] = t[4]

print >> sys.stderr, len(human), len(ancestr)

for (ga1,ga2) in ancestr:
	gh1 = dicAH[ga1]
	gh2 = dicAH[ga2]
	for (ga3,ga4) in ancestr:
		gh3 = dicAH[ga3]
		gh4 = dicAH[ga4]
		if (gh1,gh3) in human and (gh2,gh4) in human:
			print "invert", " ".join(sorted((gh1,gh2,gh3,gh4)))

for (ga1,ga2) in ancestr:
	gh1 = dicAH[ga1]
	gh2 = dicAH[ga2]
	for (ga3,ga4) in ancestr:
		gh3 = dicAH[ga3]
		gh4 = dicAH[ga4]
		if (gh1,gh4) in human:
			for (ga5,ga6) in ancestr:
				gh5 = dicAH[ga5]
				gh6 = dicAH[ga6]
				if (gh5,gh2) in human and (gh3,gh6) in human:
					print "transloc", " ".join(sorted((gh1,gh2,gh3,gh4,gh5,gh6)))

