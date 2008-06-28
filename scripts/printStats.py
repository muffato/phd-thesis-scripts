#! /users/ldog/muffato/python -OO

__doc__ = """
	Lit des donnees dans un fichier et affiche les stats
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
arguments = utils.myTools.checkArgs([("fichier",file)], [("type",str,"float")], __doc__)

f = utils.myTools.myOpenFile(arguments["fichier"], 'r')

# Lit un flot de nombres et affiche les stats
lst = []
t = eval(arguments["type"])
for l in f:
	c = l.split()
	lst.extend(t(x) for x in c)

print utils.myMaths.myStats(lst)

f.close()

