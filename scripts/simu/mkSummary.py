#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2


import sys

import utils.myFile
import utils.myMaths
import utils.myTools

arguments = utils.myTools.checkArgs( [], [("withStddev",bool,True)], "Lit une liste de noms de fichiers sur l'entree standard et affiche le fichier des moyennes" )

alldata = []
for l in sys.stdin:
	data = []
	f = utils.myFile.openFile(l.replace('\n', ''), "r")
	for s in f:
		data.append(s.replace('\n', '').split(" "))
	f.close()
	alldata.append( data )

ref = data
for (i,l) in enumerate(ref):
	s = ""
	for (j,x) in enumerate(l):
		
		try:
			int(x)
			y = int
		except ValueError:
			try:
				float(x)
				y = float
			except ValueError:
				y = None
		
		if y is None:
			res = x
		else:
			lst = utils.myMaths.myStats.valSummary( [y(data[i][j]) for data in alldata] )
			if arguments["withStddev"]:
				res = "%.2f [%.2f]" % (lst[8], lst[9])
			else:
				res = "%.2f" % lst[8]

		s += " " + res
	print s[1:]

