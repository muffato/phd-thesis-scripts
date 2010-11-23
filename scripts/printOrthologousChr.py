#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

__doc__ = """
	Compare deux genomes et renvoie les couples de chromosomes orthologues avec des stats sur l'evolution des genomes (nb genes, rearrangements)
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myMaths
import utils.myTools


#############
# FONCTIONS #
#############

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




########
# MAIN #
########

# Arguments
modeGenomeEvol = "GenomeEvolution"
arguments = utils.myTools.checkArgs( \
	[("studiedGenome",file), ("referenceGenome",file)], \
	[("orthologuesList",str,""), ("includeGaps",bool,False), ("includeScaffolds",bool,False), ("includeRandoms",bool,False), ("includeNones",bool,False), \
	("reverse",bool,False), ("minHomology",int,90)], \
	__doc__
)


# Chargement des fichiers
genome1 = utils.myGenomes.Genome(arguments["studiedGenome"])
genome2 = utils.myGenomes.Genome(arguments["referenceGenome"])
if arguments["reverse"]:
	x = genome1
	genome1 = genome2
	genome2 = x
if arguments["orthologuesList"] != "":
	genesAnc = utils.myGenomes.Genome(arguments["orthologuesList"])
else:
	genesAnc = genome2

# Les chromosomes a etudier
chr1 = []
chr2 = []
chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Chromosome]
chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Chromosome]
if arguments["includeScaffolds"]:
	chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Scaffold]
	chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Scaffold]
if arguments["includeRandoms"]:
	chr1.extend(genome1.chrList[utils.myGenomes.ContigType.Random]
	chr2.extend(genome2.chrList[utils.myGenomes.ContigType.Random]
if arguments["includeNones"]:
	chr1.extend(genome1.chrList[utils.myGenomes.ContigType.None]
	chr2.extend(genome2.chrList[utils.myGenomes.ContigType.None]


table12 = genome1.buildOrthosTable(chr1, genome2, chr2, arguments["includeGaps"], genesAnc)

# Liste de chromosomes orthologues

# On a besoin des equivalences d'un genome a l'autre
res1 = getOrthosChr(table12, chr1)
table21 = genome2.buildOrthosTable(chr2, genome1, chr1, arguments["includeGaps"], genesAnc)
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
#tt = 0
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

print >> sys.stderr, "Evolution de A (%d chr / %d genes) vers B (%d chr / %d genes)" % (nbChr1,nb1, nbChr2,nb2)
print >> sys.stderr, "Nb nouveaux genes", nouveaux
print >> sys.stderr, "Nb duplications", sum(dupliques.values()), dupliques
print >> sys.stderr, "Nb duplications (inverse)", sum(dupliques2.values()), dupliques2
print >> sys.stderr, "Nb genes perdus", pertes
print >> sys.stderr, "Nb points de cassures", cassures
print >> sys.stderr, "Nb points de fusions", fusions

for c1 in chr1:
	print "%s\t%s" % (c1, "\t".join(["%s (%d)" % (c2,n) for (c2,n) in res1[c1]]))
