#! /usr/bin/python2.4

__doc__ = """
Extrait toutes les diagonales entre chaque paire d'especes.
Les diagonales apportent les genes qui etaient sur un meme chromosome
  depuis leur ancetre commun dans les deux lignees.
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myDiags


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf", "phylTree.conf"], \
	[("fusionThreshold",int,-1), ("minimalLength",int,2), ("sameStrand",bool,True), ("cutNodes",bool,True), ("combinSameChr",bool,False), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# 1. On lit tous les fichiers
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
listEspeces = phylTree.getSpecies(phylTree.root)
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], listEspeces)

# Pour sauver de la memoire
for esp in geneBank.dicEspeces:
	del geneBank.dicEspeces[esp].dicGenes

# La structure qui accueillera les diagonales
diagEntry = dict( [(anc, []) for anc in phylTree.items] )

# La fonction qui permet de traiter les diagonales
def combinDiag2(c1, c2, d1, d2):
	global diagEntry, toStudy, options
	global e1, e2, genomes

	if len(d1) < options["minimalLength"]:
		return
	
	#print >> sys.stderr, "reception de", c1, c2, d1, d2
	
	for (e, tmp) in toStudy:
		#print >> sys.stderr, "ecriture sur", e, tmp
		if e == e1:
			#print >> sys.stderr, tmp in genomes,
			#print >> sys.stderr, e in genomes[tmp],
			#print >> sys.stderr, c1 in genomes[tmp][e],
			#print >> sys.stderr, len(genomes[tmp][e][c1])
			dd = [genomes[tmp][e][c1][i][0] for i in d1]
		else:
			#print >> sys.stderr, tmp in genomes,
			#print >> sys.stderr, e in genomes[tmp],
			#print >> sys.stderr, c2 in genomes[tmp][e],
			#print >> sys.stderr, len(genomes[tmp][e][c2])
			dd = [genomes[tmp][e][c2][i][0] for i in d2]
		diagEntry[tmp].append( (dd,"%s.%s"%(e1,c1),"%s.%s"%(e2,c2)) )
			

# 2. On prepare tous les genomes ancestraux, les genomes traduits ...
genomes = {}
locations = {}
for anc in phylTree.items:
	
	genesAnc = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % anc, False, False)

	# Les listes des especes entre lesquelles on cherche des diagonales
	groupes = [phylTree.getSpecies(e) for (e,_) in phylTree.items[anc]]
	fils = utils.myMaths.flatten(groupes)
	
	# Traduction des genomes en liste des genes ancestraux
	print >> sys.stderr, "Traduction avec les genes de", anc, "",
	genomes[anc] = {}
	for e in fils:
		genomes[anc][e] = utils.myDiags.translateGenome(geneBank.dicEspeces[e], genesAnc)
		sys.stderr.write(".")
	print >> sys.stderr, " OK"
	
	# Liste des positions des genes ancestraux dans les genomes modernes
	print >> sys.stderr, "Extraction des positions des genes de %s ..." % anc,
	lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.AncestralGenome.defaultChr]
	tmp = dict( [(e,[[] for x in lstGenesAnc]) for e in fils] )
	for ianc in xrange(len(lstGenesAnc)):
		for g in lstGenesAnc[ianc].names:
			if g not in geneBank.dicGenes:
				continue
			(e,c,i) = geneBank.dicGenes[g]
			tmp[e][ianc].append( (c,i,geneBank.dicEspeces[e].lstGenes[c][i].strand) )
	locations[anc] = tmp
	print >> sys.stderr, "OK"
	
	del genesAnc

del geneBank

for anc in phylTree.items:

	groupes = [phylTree.getSpecies(e) for (e,_) in phylTree.items[anc]]
	print >> sys.stderr, "Extraction des diagonales de %s " % anc,

	for (i,j) in utils.myTools.myMatrixIterator(len(groupes), len(groupes), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		for e1 in groupes[i]:
			for e2 in groupes[j]:

				toStudy = [(e1,anc)]
				for tmp in phylTree.items:
					s = phylTree.getSpecies(tmp)
					if (e1 in s) and (e2 not in s):
						toStudy.append( (e1,tmp) )
					elif (e2 in s) and (e1 not in s):
						toStudy.append( (e2,tmp) )
			
				utils.myDiags.iterateDiags(genomes[anc][e1], locations[anc][e2], options["fusionThreshold"], options["sameStrand"], combinDiag2)
				sys.stderr.write(".")
	del locations[anc]
	print >> sys.stderr, " OK"

del genomes

for anc in diagEntry:

	print >> sys.stderr, "Traitement et impression de %s " % anc,
	if options["cutNodes"]:
		sys.stderr.write(".")
		voisins = {}
		for (d,_,_) in diagEntry[anc]:
			#print >> sys.stderr, d
			if len(d) == 1 and d[0] not in voisins:
				voisins[d[0]] = []
				continue
			for i in xrange(1,len(d)):
				x = d[i-1]
				y = d[i]
				#print >> sys.stderr, x, y
				voisins[x] = voisins.get(x,[]) + [y]
				voisins[y] = voisins.get(y,[]) + [x]
		sys.stderr.write(".")
		
		for x in voisins:
			#print >> sys.stderr, x
			voisins[x] = len(set(voisins[x]))
		sys.stderr.write(".")

		res = []
		for (d,o1,o2) in diagEntry[anc]:
			curr = []
			#print >> sys.stderr, d
			while len(d) > 0:
				x = d.pop(0)
				#print >> sys.stderr, x
				if voisins[x] < 3:
					curr.append(x)
				else:
					if len(curr) >= 2:
						res.append( (curr,o1,o2) )
					curr = []
			if len(curr) >= 2:
				res.append( (curr,o1,o2) )
		sys.stderr.write(".")
		del voisins
	else:
		res = diagEntry[anc]
	
	sys.stderr.write(".")

	combin = utils.myTools.myCombinator([])
	if options["combinSameChr"]:
		combin1 = utils.myTools.myCombinator([])
		for i in xrange(len(res)):
			(_,o1,o2) = res[i]
			combin1.addLink( [o1+o2,i] )
		for g in combin1:
			combin.addLink(utils.myMaths.flatten([res[i][0] for i in g if type(i) == int]))
	else:
		for (d,_,_) in res:
			combin.addLink(d)
	sys.stderr.write(".")
		
	for d in combin:
		print anc, " ".join([str(x) for x in d])
	print >> sys.stderr, " OK"

