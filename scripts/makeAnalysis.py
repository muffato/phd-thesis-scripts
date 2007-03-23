#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit des donnees dans un fichier: (au choix)
		- les regroupe
		- affiche les stats
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myTools
import utils.myMaths

########
# MAIN #
########


# Arguments
modes = ["Stats", "Groups"]
(noms_fichiers, options) = utils.myTools.checkArgs(["fichier"], [("type",str,"float"), ("analysis",str,modes)], __doc__)

f = utils.myTools.myOpenFile(noms_fichiers["fichier"], 'r')

# Lit un flot de nombres et affiche les stats
if options.analysis == "Stats":

	lst = []
	t = eval(options.type)
	for l in f:
		c = l.split()
		lst.extend(t(x) for x in c)

	print utils.myMaths.myStats(lst)

elif options.analysis == "Groups":

	# On scanne toutes les lignes
	comb = utils.myTools.myCombinator([])
	for l in f:
		comb.addLink(l.split())

	# On affiche le resultat
	for x in comb:
		print " ".join(x)


f.close()
