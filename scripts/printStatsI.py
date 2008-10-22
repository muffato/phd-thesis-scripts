#! /users/ldog/muffato/python

__doc__ = """
	Lit des donnees (int) dans un fichier et affiche les stats
"""

import sys
import utils.myTools
import utils.myMaths

# Les fichiers a regrouper (par defaut, on lit l'entree standard)
files = sys.argv[1:]
if len(files) == 0:
	files.append("-")

lst = []
for f in files:
	# Ouverture du fichier
	if f == "-":
		f = sys.stdin
	else:
		f = utils.myTools.myOpenFile(f, 'r')
	# Lecture & regroupement
	for l in f:
		c = l.split()
		lst.extend(int(x) for x in c)
	# Fermeture
	f.close()

# On affiche le resultat
print utils.myMaths.myStats.txtSummary(lst)

f.close()

