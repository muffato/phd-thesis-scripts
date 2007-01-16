#! /users/ldog/muffato/python -OO

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""


##################
# INITIALISATION #
##################

# Librairies
import math
import sys
import utils.myGenomes
import utils.myTools
import utils.myMaths


#############
# FONCTIONS #
#############

# Charge le fichier de toutes les diagonales (supposees non chevauchantes)
def loadDiagsFile(nom, diagEntry):
	
	print >> sys.stderr, "Chargement du fichier de diagonales ...",
	f = utils.myTools.myOpenFile(nom, 'r')
	for l in f:

		ct = l.split('\t')
		anc = ct[0]
		l = int(ct[1])
		d = ct[2]
		esp = set([tuple(x.split('/')) for x in ct[3].split()])
		diagEntry[anc].append( (l, d, esp) )

	f.close()
	print >> sys.stderr, "OK"
	

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "diagsList"], \
	[("ancestr",str,"")], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
listEspeces = phylTree.species[phylTree.root] + ['Opossum']

# Les genes ancestraux
diagEntry = {}
for anc in phylTree.items:
	diagEntry[anc] = []
loadDiagsFile(noms_fichiers["diagsList"], diagEntry)


# Impression du resultat
nn = max([len(anc) for anc in phylTree.items])


anc = options["ancestr"]

for l in sys.stdin:
	for i in l.split():
		#i = int(l.split()[0])
		print diagEntry[anc][int(i)][1],
	print


