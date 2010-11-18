#!/usr/bin/env python2

__doc__ = """
	Lit des donnees dans un fichier et affiche les stats
"""

import utils.myFile
import utils.myMaths
import utils.myTools

arguments = utils.myTools.checkArgs([("files",utils.myTools.FileList(0))], [("txt",bool,True)], __doc__, showArgs=False)

# Les fichiers a regrouper (par defaut, on lit l'entree standard)
files = arguments["files"]
if len(files) == 0:
	files.append("-")

lst = []
for f in files:
	# Ouverture du fichier
	f = utils.myFile.openFile(f, 'r')
	# Lecture & regroupement
	for l in f:
		c = l.split()
		for x in c:
			try:
				x = int(x)
			except ValueError:
				x = float(x)
			lst.append(x)
	# Fermeture
	f.close()

# On affiche le resultat
if arguments["txt"]:
	print utils.myMaths.myStats.txtSummary(lst)
else:
	print " ".join(("%s" % x) for x in utils.myMaths.myStats.valSummary(lst))

