#! /users/ldog/muffato/python -OO

import sys

nb = 0
for ligne in sys.stdin:
	champs = ligne.split()
	nb += 1
	print "sd%d" % nb, champs[0], champs[1], champs[1]
	print "sd%d" % nb, champs[2], champs[3], champs[3]

