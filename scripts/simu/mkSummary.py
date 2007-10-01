#! /users/ldog/muffato/python -OO


# Librairies
import sys
import utils.myMaths
import utils.myTools

(noms_fichiers, options) = utils.myTools.checkArgs( [], [], "Lit une liste de noms de fichiers sur l'entree standard et affiche le fichier des moyennes" )

alldata = []
for l in sys.stdin:
	data = []
	f = open(l[:-1], "r")
	for s in f:
		data.append(s[:-1].split("\t"))
	f.close()
	alldata.append( data )

ref = data
for (i,l) in enumerate(ref):
	s = ""
	for (j,x) in enumerate(l):
		try:
			float(x)
			if "." not in x:
				0/0
			lst = utils.myMaths.myStats( [float(data[i][j]) for data in alldata] )
			res = "%.2f [%.2f]" % (lst.mean, lst.stddev)

		except Exception:
			res = x
		s += "\t" + res
	print s[1:]

