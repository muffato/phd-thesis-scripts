#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

# Parmi les sequences lues sur l'entree standard, conserve la derniere de chaque gene

import sys

res = {}
for l in sys.stdin:
	(name,_,seq) = l.replace('\n', '').rpartition("\t")
	res[name] = seq

for x in res.iteritems():
	print "\t".join(x)

