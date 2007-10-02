#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit des donnees dans un fichier et les regroupe
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myTools

########
# MAIN #
########

# Les fichiers a regrouper (par defaut, on lit l'entree standard)
files = sys.argv[1:]
if len(files) == 0:
	files.append("-")

comb = utils.myTools.myCombinator()
for f in files:
	# Ouverture du fichier
	if f == "-":
		f = sys.stdin
	else:
		f = utils.myTools.myOpenFile(f, 'r')
	# Lecture & regroupement
	for l in f:
		comb.addLink(l.split())
	# Fermeture
	f.close()

# On affiche le resultat
for x in comb:
	print " ".join(x)

