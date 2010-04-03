#! /users/ldog/muffato/python

# Parmi les sequences lues sur l'entree standard, conserve la derniere de chaque gene

import sys

res = {}
for l in sys.stdin:
	(name,_,seq) = l.replace('\n', '').rpartition("\t")
	res[name] = seq

for x in res.iteritems():
	print "\t".join(x)

