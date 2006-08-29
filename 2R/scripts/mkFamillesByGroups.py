#! /usr/bin/python2.4

import sys

classes = [[] for i in range(4)]
noms = ['alpha', 'beta', 'gamma', 'epsilon']

for ligne in sys.stdin:
	champs = ligne.split()
	if champs[0] in noms:
		ind = noms.index(champs[0])
		classes[ind] = champs[1:]
	else:
		for g1 in classes[0]:
			print g1,
		for g2 in classes[1]:
			print g2,
		print
		
		for g1 in classes[2]:
			print g1,
		for g2 in classes[3]:
			print g2,
		print
		
		classes = [[] for i in range(4)]
	
