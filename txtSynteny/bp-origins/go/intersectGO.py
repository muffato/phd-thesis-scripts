#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Charge les liens GO et renvoie les niveaux d'intersection
"""

import sys
import math
import itertools
import collections

import utils.myFile
import utils.myTools
import utils.myMaths
import utils.myGenomes

# Arguments
arguments = utils.myTools.checkArgs( [("GO_data",file), ("gene_pairs",file)], [], __doc__)

GO_data = {}
f = utils.myFile.openFile(arguments["GO_data"], "r")
for l in f:
	t = l.split()
	GO_data[t[0]] = set(t[1:])
f.close()

f = utils.myFile.openFile(arguments["gene_pairs"], "r")
for l in f:
	(g1,g2) = l.split()
	if (g1 not in GO_data) or (g2 not in GO_data):
		continue
	inter = GO_data[g1].intersection(GO_data[g2])
	levels = [int(x.split(":")[2]) for x in inter]
	print g1, g2, max(levels)
f.close()


