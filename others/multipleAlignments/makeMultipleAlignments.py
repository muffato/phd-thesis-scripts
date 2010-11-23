#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Telecharge depuis le site d'Ensembl les fichiers des arbres de proteines
"""

import sys
import collections

import utils.myFile
import utils.myTools
import utils.myGenomes
import utils.myPhylTree
import utils.myProteinTree

arguments = utils.myTools.checkArgs( [("dataFile",file)], [("output",str,"output/align.%s.mfa")], __doc__)

def explodeSeq(seq, align):
	align = align.replace("M", "M ").replace("D", "D ").split()
	seq = seq
	res = []
	i = 0
	for x in align:
		l = int(x[:-1]) if len(x) > 1 else 1
		if x[-1] == "M":
			res.append(seq[i:(i+l)])
			i += l
		else:
			res.append("-" * l)
	#assert 0 <= len(seq)-3*i <= 5, (len(seq),i,3*i,align,seq)
	#if not (0 <= len(seq)-3*i <= 5):
	if len(seq) != i:
		print  (len(seq),i,align,seq,res)
	return "".join(res)

i = None
f = utils.myFile.openFile(arguments["dataFile"], "r")
for l in f:
	t = l.split()
	res = explodeSeq(t[5], t[4])
	if t[0] != i:
		if i != None:
			fo.close()
		i = t[0]
		fo = utils.myFile.openFile(arguments["output"] % i, "w")
	print >> fo, ">%s|%s|%s" % tuple(t[1:4])
	print >> fo, res
f.close()
if i != None:
	fo.close()
