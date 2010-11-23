#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
Lit des familles de genes sur l'entree standard.
Rajoute a chaque famille les genes associes dans le genome ancestral passe en argument.
"""


# INITIALISATION #

# Librairies
import sys
import utils.myGenomes
import utils.myTools

arguments = utils.myTools.checkArgs((["genesAncestraux"],file), [], __doc__)

genesAnc = utils.myGenomes.Genome(arguments["genesAncestraux"])

lst = set()
nb = 0
for s in sys.stdin:
	t = set()
	for g in s.split():
		if g in genesAnc.dicGenes:
			t.add(genesAnc.dicGenes[g])
	
	for (c,i) in t:
		print " ".join(genesAnc.lstGenes[c][i].names),
	print s,

print >> sys.stderr, nb, "familles pour", len(lst), "couples conservees"
