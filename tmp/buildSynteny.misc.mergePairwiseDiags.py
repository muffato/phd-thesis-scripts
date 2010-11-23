#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import sys
import collections

import utils.myFile
import utils.myTools

arguments = utils.myTools.checkArgs([("diagList",utils.myTools.FileList(1))], [], "Rassemble les blocs pairwise de plusieurs fichiers en un seul")

data = collections.defaultdict(list)
for s in arguments["diagList"]:
	f = utils.myFile.openFile(s, "r")
	for l in f:
		t = l[:-1].split("\t")
		data["\t".join(t[1:8] + [t[9]])].append(t[0] + "=" + t[8])
	f.close()

for x in data:
	print x + "\t" + "|".join(data[x])

