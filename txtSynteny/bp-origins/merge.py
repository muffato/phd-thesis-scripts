#! /users/ldog/muffato/python

import sys

f1 = open(sys.argv[1], "r")
f2 = open(sys.argv[2], "r")

l1 = None
l2 = None

while True:
	
	if l1 is None:
		l1 = f1.readline()
	if l2 is None:
		l2 = f2.readline()
	
	if (l1 == "") and (l2 == ""):
		break

	try:
		g1 = l1.split()[1]
		g2 = l2.split()[1]
	except IndexError:
		print l1
		l1 = None
		continue

	if g1 == g2:
		print l2
		l1 = None
		l2 = None
	else:
		print l1
		l1 = None

f1.close()
f2.close()

