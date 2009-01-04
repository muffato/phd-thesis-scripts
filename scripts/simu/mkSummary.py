#! /users/ldog/muffato/python


import sys

import utils.myFile
import utils.myMaths
import utils.myTools

arguments = utils.myTools.checkArgs( [], [], "Lit une liste de noms de fichiers sur l'entree standard et affiche le fichier des moyennes" )

alldata = []
for l in sys.stdin:
	data = []
	f = utils.myFile.openFile(l.replace('\n', ''), "r")
	for s in f:
		data.append(s.replace('\n', '').split("\t"))
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
			res = "%.2f [%.2f]" % (lst[8], lst[9])

		s += "\t" + res
	print s[1:]

