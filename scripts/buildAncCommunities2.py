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
		d = [int(x) for x in ct[2].split(' ')]
		esp = set([tuple(x.split('/')) for x in ct[3].split()])
		diagEntry[anc].append( (l, d, esp) )

	f.close()
	print >> sys.stderr, "OK"
	

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "genesList.conf", "diagsList"], \
	[("ancestr",str,""), ("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])

# Les genes ancestraux
#anc = options["ancestr"]
genesAnc = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % options["ancestr"])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr]

# On separe les especes en trois
(fils1, fils2) = phylTree.branchesSpecies[options["ancestr"]]
outgroup = phylTree.outgroupSpecies[options["ancestr"]]

print >> sys.stderr, fils1, fils2, outgroup
sys.exit(0)

# Les genomes modernes
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], phylTree.species[options["ancestr"]])

# Les diagonales
diagEntry = {}
for anc in phylTree.items:
	diagEntry[anc] = []
loadDiagsFile(noms_fichiers["diagsList"], diagEntry)
lstDiags = diagEntry[options["ancestr"]]
del diagEntry

# Test des diagonales
# Est-ce que des diagonales sont sur un seul chromosome sans qu'on le sache




#sys.exit(0)


print >> sys.stderr, len(lstDiags)
# Les genes seuls vont devenir des diagonales de 1
genesSeuls = set(xrange(len(lstGenesAnc)))
for (_,d,_) in lstDiags:
	genesSeuls.difference_update(d)
for i in genesSeuls:
	continue
	esp = [geneBank.dicGenes[s] for s in lstGenesAnc[i].names]
	if len(esp) >= 2:
		lstDiags.append( (1,[i],set([(e,c) for (e,c,_) in esp])) )
	#print "ajout de ", lstDiags[-1]

print >> sys.stderr, len(lstDiags)
#sys.exit(0)

combin = utils.myTools.myCombinator([])
dicAretes = dict([(i,{}) for i in xrange(len(lstDiags))])
for (i1,i2) in utils.myTools.myMatrixIterator(len(lstDiags), len(lstDiags), utils.myTools.myMatrixIterator.StrictUpperMatrix):
	(_,_,e1) = lstDiags[i1]
	(_,_,e2) = lstDiags[i2]
	commun = e1.intersection(e2)
	communEsp = set([e for (e,_) in commun])
	
	propF1 = float(len(communEsp.intersection(fils1))) / float(len(fils1))
	propF2 = float(len(communEsp.intersection(fils2))) / float(len(fils2))
	if len(outgroup) == 0:
		propOut = 0
	else:
		propOut = float(len(communEsp.intersection(outgroup))) / float(len(outgroup))
	
	score = propF1*propF2 + propF1*propOut + propF2*propOut
	
	if score > 0:
		combin.addLink([i1,i2])
		dicAretes[i1][i2] = score
		dicAretes[i2][i1] = score

#sys.exit(0)

# On traite chaque composante connexe
print >> sys.stderr, "impression"
indComp = 0
for g in combin:
	#print len(g)
	#continue
	indComp += 1
	nb = len(g)
	f = open('/users/ldog/muffato/work/temp/walktrap/%s/nodes.%d' % (options["ancestr"], indComp), 'w')
	for i in xrange(nb):
		print >> f, i, g[i]
	f.close()
	f = open('/users/ldog/muffato/work/temp/walktrap/%s/graph.%d' % (options["ancestr"], indComp), 'w')
	for (i1,i2) in utils.myTools.myMatrixIterator(nb, nb, utils.myTools.myMatrixIterator.StrictUpperMatrix):
		if g[i2] in dicAretes[g[i1]]:
			print >> f, i1, i2, dicAretes[g[i1]][g[i2]]
	f.close()
	
