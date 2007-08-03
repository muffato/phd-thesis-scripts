#! /users/ldog/muffato/python -OO

__doc__ = """
	Telecharge depuis le site de l'UCSC les fichiers utiles pour creer:
		- les listes de genes
		- les listes d'orthologues
		- les listes de paralogues
"""

# Librairies
import os
import sys
import utils.myPhylTree
import utils.myTools


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("species",str,""), ("root",str,""), ("OUT.directory",str,"")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

OUTgenesFile = options["OUT.directory"] + "/genes/genes.%s.list.bz2"
OUTorthosFile = options["OUT.directory"] + "/orthologs/orthos.%s.%s.list.bz2"
OUTparas2File = options["OUT.directory"] + "/orthologs/paras.%s.%s.list.bz2"
OUTparas1File = options["OUT.directory"] + "/paralogs/paras.%s.list.bz2"
OUTxrefFile = options["OUT.directory"] + "/genes/xref/xref.%s.list.bz2"
OUTancGenesFile = options["OUT.directory"] + "/fullAncGenes/ancGenes.%s.list.bz2"
OUTone2oneGenesFile = options["OUT.directory"] + "/ancGenes/one2one.%s.list.bz2"

# Les genes ancestraux
genesAnc = {}
lstAncestr = []
for anc in phylTree.dicLinks[options["species"]][options["root"]][1:]:
	try:
		genesAnc[anc] = utils.myGenomes.loadGenome(OUTancGenesFile % phylTree.fileName[anc])
		lstAncestr.insert(0, anc)
	except IOError:
		pass

# Le dictionnaire genes ancestraux -> xref
dicAncXRef = utils.myTools.defaultdict(set)
# Le dictionnaire xref -> genes ancestraux
dicXRefAnc = utils.myTools.defaultdict(set)
# Le dictionnaire genes UCSC -> xref
dicUCSCRef = utils.myTools.defaultdict(set)
for esp in phylTree.listSpecies:
	print >> sys.stderr, "Chargement des annotations xref de %s ..." % esp,
	try:
		for ligne in utils.myTools.myOpenFile(OUTxrefFile % phylTree.fileName[esp], 'r'):
			t = [intern(x) for x in ligne[:-1].split('\t')]
			if esp == options["species"]:
				dicUCSCRef[t[0]].update(t[3:])
			for (anc,genome) in genesAnc.iteritems():
				if t[0] not in genome.dicGenes:
					continue
				r = (anc,genome.dicGenes[t[0]][1])
				dicAncXRef[r].update(t[3:])
				for x in t[3:]:
					dicXRefAnc[x].add(r)
	except IOError:
		print >> sys.stderr, "-",
	print >> sys.stderr, "OK"

print >> sys.stderr, len(dicAncXRef), len(dicXRefAnc), len(dicUCSCRef)

nbPotentiels = 0
nbOK = 0
for (gene,xref) in dicUCSCRef.iteritems():
	
	# On fait la liste des genes ancestraux lies aux xref
	#   dans la limite de 1 / ancetre
	tab = {}
	#tab = utils.myTools.defaultdict(st)
	for x in xref:
		if x not in dicXRefAnc:
			continue
		for (anc,i) in dicXRefAnc[x]:
			#tab[anc].append(i)
			if (anc in tab) and (tab[anc] != i):
				tab[anc] = None
				continue
			tab[anc] = i
	
	# On fait le tri parmi les ancetres dont on est sur
	for anc in tab.keys():
		if tab[anc] == None:
			del tab[anc]

	#print gene, len(tab), tab
	#print gene, tab
	
	# Il faut qu'il en reste un pour placer la famille
	if len(tab) == 0:
		continue

	nbPotentiels += 1

	# On va verifier qu'on a defini une lignee sans ambiguite
	lastGene = None
	for anc in lstAncestr:
		if lastGene == None:
			if anc in tab:
				lastGene = genesAnc[anc].lstGenes[utils.myGenomes.Genome.defaultChr][tab[anc]]
			continue
		newGenes = [i for (_,i) in genesAnc[anc].getPosition(lastGene)]
		if anc in tab:
			# Si j'ai une definition, elle doit concorder
			if tab[anc] not in newGenes:
				break
		else:
			# Si pas de definition, il ne faut pas avoir de duplication
			if len(newGenes) > 1:
				break
	else:
		nbOK += 1
		print gene, len(tab), tab

print >> sys.stderr, nbPotentiels, nbOK


