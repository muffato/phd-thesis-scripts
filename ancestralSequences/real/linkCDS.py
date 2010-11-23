# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import sys

dic = {}
f = open(sys.argv[1], "r")
for l in f:
	dic[ l[:-4].split('/')[2] ] = l[:-1]
f.close()

notransc = {}
for l in sys.stdin:
	t = l.split("\t")
	if len(t) >= 7:
		notransc[t[4]] = t[6]

f = open(sys.argv[2], "r")
for (i,l) in enumerate(f):
	names = l.replace("(", " ").replace(")", " ").replace(",", " ").split()
	for s in names:
		s = s.split("/")[0]
		if s not in dic:
			s = notransc[s]
			#print >> sys.stderr, s
			#continue
		print "ln -s ../../../cds/%s fam/%d/seq/%s" % (dic[s], i+1, s)
f.close()


