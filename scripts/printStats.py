#! /users/ldog/muffato/python

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
(noms_fichiers, options) = utils.myTools.checkArgs([], [], __doc__)

lst = []
for l in sys.stdin:
	c = l.split()
	for x in c:
		lst.append(int(x))

if len(lst) > 0:
	res = (utils.myMaths.moyenne(lst), utils.myMaths.ecartType(lst), utils.myMaths.mediane(lst), min(lst), max(lst), len(lst))
else:
	res = (0, 0, 0, 0, 0, 0)

print "%.2f\t%.2f\t%d\t%d\t%d\t%d" % res
