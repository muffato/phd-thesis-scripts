#! /users/ldog/muffato/python -OO


__doc__ = """
Lit un genome sous forme de liste de numeros de genes ancestraux.
Convertit l'entree en un veritable fichier de genome ancestral.
"""

##################
# INITIALISATION #
##################

# Librairies
import utils.myGenomes
import utils.myTools


(noms_fichiers, options) = utils.myTools.checkArgs(["ancGenesList", "genome"], [], __doc__)

genesAnc = utils.myGenomes.loadGenome(noms_fichiers["ancGenesList"])

nb = 1
f = utils.myTools.myOpenFile(noms_fichiers["genome"], 'r')
for l in f:
	c = [int(x) for x in l.split()]
	for i in c:
		print nb, " ".join(genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr][i].names)
	nb += 1

f.close()

