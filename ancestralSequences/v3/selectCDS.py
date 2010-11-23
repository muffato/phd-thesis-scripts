#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import sys

for l in sys.stdin:
	l = l.split()[1].upper()
	#print >> sys.stderr, l, len(l)%3
	if len(l) < 300:
		continue
	if len(l) > 3000:
		continue
	if (len(l) % 3) != 0:
		continue
	if "N" in l:
		continue
	cod = set(l[i:i+3] for i in xrange(0,len(l),3))
	#print >> sys.stderr, cod
	if "TAA" in cod:
		continue
	if "TAG" in cod:
		continue
	if "TGA" in cod:
		continue
	s3 = [l[i+2] for i in xrange(0,len(l),3)]
	print l, (s3.count("G")+s3.count("C"))/(len(s3)*0.01)

