#! /users/ldog/muffato/python -OO

__doc__ = """
Lit des donnees numeriques sur l'entree standard et affiche
Moyenne - Ecart type - Mediane - Min - Max
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
(noms_fichiers, options) = utils.myTools.checkArgs(["fichier"], [("type",str,"float")], __doc__)

lst = []
f = utils.myTools.myOpenFile(noms_fichiers["fichier"], 'r')
t = eval(options.type)
for l in f:
	c = l.split()
	for x in c:
		lst.append(t(x))
f.close()

if len(lst) > 0:
	res = (utils.myMaths.moyenne(lst), utils.myMaths.ecartType(lst), utils.myMaths.mediane(lst), min(lst), max(lst), len(lst))
else:
	res = (0, 0, 0, 0, 0, 0)

print "%.2f\t%.2f\t%d\t%d\t%d\t%d" % res
