#! /users/ldog/muffato/python -OO

import sys
import utils.myGenomes

genome = utils.myGenomes.Genome(sys.argv[1])

# On associe chaque gene aux elements alignes
d = utils.myTools.defaultdict(list)
for l in sys.stdin:
	l = l.replace("\n", "")
	t = l.split("\t")
	genes = list(genome.getGenesAt(t[0], beg = int(t[1]), end = int(t[2])))
	if len(genes) == 1:
		d[genes[0].names[0]].extend(t[3:])

# Affichage
for (name,lst) in d.iteritems():
	l = [name,"-","-"] + lst
	print utils.myTools.printLine([name,"-","-"] + lst)

