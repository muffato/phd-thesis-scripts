#! /users/ldog/muffato/python -OO

__doc__ = """
	Insere les genes d'une espece dans les familles ancestrales grace aux annotations xref
"""

# Librairies
import os
import sys
import utils.myPhylTree
import utils.myTools


# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("ancestr",str,""), ("root",str,""), ("OUT.directory",str,"")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
anc = options["ancestr"]

OUTgenesFile = options["OUT.directory"] + "/genes/genes.%s.list.bz2"
OUTorthosFile = options["OUT.directory"] + "/orthologs/orthos.%s.%s.list.bz2"
OUTparas2File = options["OUT.directory"] + "/orthologs/paras.%s.%s.list.bz2"
OUTparas1File = options["OUT.directory"] + "/paralogs/paras.%s.list.bz2"
OUTxrefFile = options["OUT.directory"] + "/genes/xref/xref.%s.list.bz2"
OUTancGenesFile = options["OUT.directory"] + "/fullAncGenes/ancGenes.%s.list.bz2"
OUTone2oneGenesFile = options["OUT.directory"] + "/ancGenes/one2one.%s.list.bz2"


# Les genes ancestraux
genesAnc = {}
for anc in phylTree.dicLinks[anc][options["root"]][1:]:
	genesAnc[anc] = utils.myGenomes.Genome(OUTancGenesFile % phylTree.fileName[anc])
lstAncestr = genesAnc.keys()


especes = phylTree.species[anc]
ancSup = phylTree.parent[anc]

# 1. On se base sur les familles de l'ancetre juste au dessus
phylTree.loadAllSpeciesSince(ancSup, OUTgenesFile)
del phylTree.dicGenomes
comb = utils.myTools.myCombinator()
for g in genesAnc[ancSup]:
	t = [s for s in g.names if phylTree.dicGenes[s][0] in especes]
	comb.addLink(t)
del phylTree.dicGenes

# 2.a on charge les xref pour ancrer les genes sur les familles ancestrales de tous les ancetres disponibles

# Le dictionnaire genes ancestraux -> xref
dicAncXRef = utils.myTools.defaultdict(set)
# Le dictionnaire xref -> genes ancestraux
dicXRefAnc = utils.myTools.defaultdict(set)
# Le dictionnaire genes UCSC -> xref
dicUCSCRef = utils.myTools.defaultdict(set)
# Le dictionnaire gene moderne -> espece
dicGeneEsp = {}
for esp in phylTree.listSpecies:
	print >> sys.stderr, "Chargement des annotations xref de %s ..." % esp,
	for ligne in utils.myTools.myOpenFile(OUTxrefFile % phylTree.fileName[esp], 'r'):
		t = [intern(x) for x in ligne[:-1].split('\t')]
		if esp == options["species"]:
			dicUCSCRef[t[0]].update(t[3:])
		dicGenes[t[0]] = esp
		for (anc,genome) in genesAnc.iteritems():
			if t[0] not in genome.dicGenes:
				continue
			r = (anc,genome.dicGenes[t[0]][1])
			dicAncXRef[r].update(t[3:])
			for x in t[3:]:
				dicXRefAnc[x].add(r)
	print >> sys.stderr, "OK"
print >> sys.stderr, "%d genes ancestraux -> xref, %d xref -> genes ancestraux" % (len(dicAncXRef),len(dicXRefAnc))


# 2.b On selectionne celles coherentes avec les familles ancestrales
nbPotentiels = 0
nbOK = 0
lstOrthos = utils.myTools.defaultdict(set)
for (gene,xref) in dicUCSCRef.iteritems():
	
	# On fait la liste des genes ancestraux lies aux xref
	#   dans la limite de 1 / ancetre
	tab = {}
	for x in xref:
		if x not in dicXRefAnc:
			continue
		for (anc,i) in dicXRefAnc[x]:
			if (anc in tab) and (tab[anc] != i):
				tab[anc] = None
				continue
			tab[anc] = i
	
	# On fait le tri parmi les ancetres dont on est sur
	for anc in tab.keys():
		if tab[anc] == None:
			del tab[anc]

	# Il faut qu'il en reste un pour placer la famille
	if len(tab) == 0:
		continue

	nbPotentiels += 1

	# On va verifier qu'on a defini une lignee sans ambiguite
	lastGene = None
	for anc in lstAncestr:
		if lastGene == None:
			if anc in tab:
				lastGene = genesAnc[anc].lstGenes[None][tab[anc]]
			continue
		newGenes = [i for (_,i) in genesAnc[anc].getPosition(lastGene.names)]
		if anc in tab:
			# Si j'ai une definition, elle doit concorder
			if tab[anc] not in newGenes:
				break
		else:
			# Si pas de definition, il ne faut pas avoir de duplication
			if len(newGenes) > 1:
				break
			# La definition devient celle par defaut
			tab[anc] = newGenes[0]
	else:
		nbOK += 1
		# On remonte dans le temps voir si on n'a pas oublie un ancetre
		for anc in lstAncestr:
			if anc not in tab:
				oldGene = [i for (_,i) in genesAnc[anc].getPosition(lastGene.names)]
				if len(oldGene) > 1:
					# N'arrive jamais
					print >> sys.stderr, gene, xref, tab, anc, lastGene, oldGene
				elif len(oldGene) == 1:
					tab[anc] = oldGene[0]
		# On peut ecrire les nouveaux fichiers d'orthologues
		tmp = utils.myTools.defaultdict(set)
		for anc in lstAncestr:
			for g in genesAnc[anc].lstGenes[None][tab[anc]].names:
				tmp[dicGeneEsp[g]].append( (gene,g,phylTree.dicParents[dicGeneEsp[g]][options["species"]]) )
		for esp in tmp:
			lstOrthos[esp].add(gene)

#print >> sys.stderr, "%d genes a l'origine, %d potentiels, %d OK" % (len(dicUCSCRef),nbPotentiels,nbOK)

sys.exit(0)
# 3. On ecrit le d'orthologues adequat
for (esp,l) in lstOrthos.iteritems():
	print >> sys.stderr, "Ecriture des orthologues avec %d ..." % esp,
	fo1 = utils.myTools.myOpenFile(OUTorthosFile % (phylTree.fileName[esp],phylTree.fileName[options["species"]]), 'w')
	fo2 = utils.myTools.myOpenFile(OUTorthosFile % (phylTree.fileName[options["species"]],phylTree.fileName[esp]), 'w')
	anc = phylTree.dicParents[esp][options["species"]]
	for (gEsp,gOther) in l:
		obj1 = "\t".join( (gEsp,"-","-") )
		obj2 = "\t".join( (gOther,"-","-") )
		data = "\t".join( (anc,"0","0","0","0","ortholog_one2one") )
		print >> fo1, "\t".join( (obj2,obj1,data) )
		print >> fo2, "\t".join( (obj1,obj2,data) )
	fo1.close()
	fo2.close()
	print >> sys.stderr, "%d OK" % len(l)

