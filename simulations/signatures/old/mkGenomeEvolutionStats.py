#!/usr/bin/env python2

__doc__ = """
	Compare les familles ancestrales pour renvoyer les stats des genes.
	Compare les genomes modernes pour les stats de rearrangements.
"""

import sys
import utils.myMaths
import utils.myTools
import utils.myGenomes
import utils.myPhylTree


# Renvoie la liste des chromosomes orthologues
def getOrthosChr(table, chr):
	
	res = {}
	for c1 in chr:
		# On ne garde que les chromosomes orthologues majoritaires
		lst = [x for (x,_) in utils.myMaths.flatten([x for (_,x) in table[c1]])]
		count = [(lst.count(x),x) for x in set(lst)]
		count.sort( reverse = True )
		nb = (len(lst)*arguments["minHomology"])/100
		tmp = []
		for (n,c2) in count:
			tmp.append( (c2,n) )
			nb -= n
			if nb <= 0:
				break
		res[c1] = tmp
	return res


def compare(esp1, esp2):

	genome1 = genes[esp1]
	chr1 = genome1.lstChr + genome1.lstScaff + genome1.lstRand
	genome2 = genes[esp2]
	chr2 = genome2.lstChr + genome2.lstScaff + genome1.lstRand
	genesAnc = genes[phylTree.dicParents[esp1][esp2]]

	# On a besoin des equivalences d'un genome a l'autre
	table12 = genome1.buildOrthosTable(chr1, genome2, chr2, False, genesAnc)
	res1 = getOrthosChr(table12, chr1)
	table21 = genome2.buildOrthosTable(chr2, genome1, chr1, False, genesAnc)
	res2 = getOrthosChr(table21, chr2)

	# Et d'en tirer les best-hits reciproques chromosomiques
	besthits1 = {}
	besthits2 = {}
	for c1 in chr1:
		for (c2,_) in res1[c1]:
			if len([(c,x) for (c,x) in res2[c2] if c == c1]) > 0:
				besthits1[c1] = besthits1.get(c1,[]) + [c2]
				besthits2[c2] = besthits2.get(c2,[]) + [c1]

	# Initialisation
	cassures = 0
	nbChr1 = 0
	pertes = 0
	nb1 = 0
	dupliques = {}
	tmp2 = {}

	# Parcours du genome 1
	for c1 in chr1:
		
		# Les calculs sur les nombres de genes
		nb1 += len(genome1.lstGenes[c1])
		pertes += len(genome1.lstGenes[c1]) - len(table12[c1])
		for (_,g2) in table12[c1]:
			# Les duplications du genome 2 vers le genome 1
			for x in g2:
				tmp2[x] = tmp2.get(x,0) + 1
			# Des pertes du genome 1 vers le genome 2
			if len(g2) == 0:
				pertes += 1
			# Les duplications du genome 1 vers le genome 2
			elif len(g2) > 1:
				dupliques[len(g2)] = dupliques.get(len(g2),0) + 1
		
		# Les calculs sur les rearrangements
		if len(besthits1.get(c1,[])) == 0:
			continue
		nbChr1 += 1
		cassures += len(besthits1[c1]) - 1
			
		
	# Parcours du genome 2
	nb2 = 0
	nbChr2 = 0
	fusions = 0
	for c2 in chr2:
		# On met a jour le nombre de genes du genome 1
		nb2 += len(genome2.lstGenes[c2])
		
		# Les calculs sur les rearrangements
		if len(besthits2.get(c2,[])) == 0:
			continue
			
		nbChr2 += 1
		fusions += len(besthits2[c2]) - 1

	nouveaux = nb2 - len(tmp2)

	# Les duplications dans l'autre sens
	dupliques2 = {}
	for x in tmp2.itervalues():
		dupliques2[x] = dupliques2.get(x,0) + 1
	dupliques2.pop(1, None)

	return ( (nbChr1,nb1, nbChr2,nb2), (nouveaux,dupliques,dupliques2,pertes), (cassures,fusions) )


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file)], \
	[("minHomology",int,90), ("genesFile",str,""), ("ancGenesFile",str,"")], \
	__doc__
)

# Chargement des tous les fichiers
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
genes = {}
for e in phylTree.listSpecies:
	genes[e] = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])
for a in phylTree.listAncestr:
	genes[a] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[a])


# 1. On charge chaque branche pour les stats de nb de genes
def do(node):
	for (e,l) in phylTree.items.get(node, []):
		print >> sys.stderr, "%s -> %s ..." % (node,e),
		(resume,statsGenes,statsChr) = compare(node, e)
		print utils.myFile.myTSV.printLine((node,e) + statsGenes)
		print >> sys.stderr, "OK"
		do(e)
do(phylTree.root)

