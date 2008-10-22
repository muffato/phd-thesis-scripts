#! /users/ldog/muffato/python

__doc__ = """
Extrait toutes les diagonales entre chaque paire d'especes.
Les diagonales apportent les genes qui etaient sur un meme chromosome
  depuis leur ancetre commun dans les deux lignees.
"""


import sys
import collections

import utils.myPhylTree
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myDiags2



def getLongestDiags(oldDiags):

	# On construit la table "gene" -> "liste des diagonales qui le contiennent"
	# On s'en sert pour avoir les listes de diagonales chevauchantes
	dic = collections.defaultdict(list)
	diags = range(len(oldDiags))
	combin = utils.myTools.myCombinator()
	for (i,(_,_,da,ds)) in enumerate(oldDiags):
		diags[i] = [(x,ds[j]) for (j,x) in enumerate(da)]
		# Combinaison selon les genes ancestraux
		for j in da:
			dic[j].append(i)
		combin.addLink([i])
	
	for s in dic:
		combin.addLink(dic[s])
	
	del combin.dic

	for g in combin:
		# On casse les carrefours de diagonales les uns apres les autres
		gr = utils.myDiags2.WeightedDiagGraph([diags[i] for i in g])
		# Les diagonales resultat
		for (res,strand) in gr.getBestDiags():
			# Filtre de taille
			if len(res) < arguments["minimalLength"]:
				continue
			# On rajoute la liste des especes qui soutiennent la diagonale
			ok = set()
			# Test de chaque diagonale
			for i in g:
				d = [x for (x,_) in diags[i]]
				for j in xrange(len(d)-1):
					try:
						i1 = res.index(d[j])
						i2 = res.index(d[j+1])
					except ValueError:
						continue
					# Si la diagonale 'i' soutient la diagonale ancestrale, on rajoute les especes dont elle est tiree
					if abs(i1-i2) == 1:
						ok.add( (oldDiags[i][0][0],oldDiags[i][0][1]) )
						ok.add( (oldDiags[i][1][0],oldDiags[i][1][1]) )
						break
			yield (res,strand,tuple(ok))
	



# Cherche si la diagonale est localisee sur un chromosome d'une espece, non ordonnee
def findNewSpecies(d, esp, anc):
	
	# Les especes pour lesquelles la diagonale n'a pas ete detectee
	lstEsp = set(listSpecies).difference([e for (e,_) in esp])
	lstGenesAncAnc = genesAnc[anc].lstGenes[None]
	res = []
	for e in lstEsp:
		# L'ancetre commun
		a = phylTree.dicParents[anc][e]
		lstGenesAncA = genesAnc[a].lstGenes[None]
		genome = dicGenomes[e]

		# On va intersecter les chromosomes des orthologues de chaque gene
		poss = set()
		for i in d:
		
			# Les noms associes au gene ancestral
			names = lstGenesAncAnc[i].names
			if a != anc:
				names = genesAnc[a].getOtherNames(names[0])
			tmp = [c for (c,_) in genome.getPosition(names)]
			
			# On intersecte les ensembles de chromosomes entre eux
			if len(poss) == 0:
				# Il faut initialiser
				poss.update(tmp)
			elif len(tmp) > 0:
				poss.intersection_update(tmp)
				# Utile de continuer ?
				if len(poss) == 0:
					break
		
		# S'il en reste un, c'est gagne !
		if len(poss) != 0:
			res.append( (e,poss.pop()) )
	
	return res



# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("target",str)], \
	[("fusionThreshold",int,-1), ("minimalLength",int,2), ("sameStrand",bool,True), ("keepOnlyOrthos",bool,False),
	("useOutgroups",bool,False), \
	("showProjected",bool,False), ("showAncestral",bool,True), ("searchUndetectedSpecies",bool,True), ("cutLongestPath",bool,True), \
	("OUT.projDiags",str,"proj/diags.%s.list.bz2"), \
	("OUT.ancDiags",str,"anc/diags.%s.list.bz2"), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenomesFile",str,"~/work/ancestralGenomes/Genome.%s.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

# Les especes a utiliser
dicGenomes = {}
tmp = arguments["target"].split(',')
if len(arguments["target"]) == 0:
	print >> sys.stderr, "Aucune cible indiquee pour l'extraction des diagonales"
	sys.exit(1)
listSpecies = []
for x in tmp:
	if x[0] != '.':
		listSpecies.extend(phylTree.species[x])
		for e in phylTree.species[x]:
			dicGenomes[e] = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])
	else:
		listSpecies.append(x[1:])
		dicGenomes[x[1:]] = utils.myGenomes.Genome(arguments["ancGenomesFile"] % phylTree.fileName[x[1:]], withChr=True)

# Les outgroup du noeud le plus ancien
if arguments["useOutgroups"]:
	target = listSpecies[0]
	for e in listSpecies:
		target = phylTree.dicParents[target][e]
	listSpecies.extend(phylTree.outgroupSpecies[target])
	for e in phylTree.outgroupSpecies[target]:
		dicGenomes[e] = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])

# La liste des ancetres edites
dicLinks = [(e1,e2,set(phylTree.dicLinks[e1][e2][1:-1] + [phylTree.dicParents[e1][e2]])) for (e1,e2) in utils.myTools.myIterator.tupleOnStrictUpperList(listSpecies)]
tmp = set()
for (_,_,s) in dicLinks:
	tmp.update(s)
diagEntry = {}
genesAnc = {}
for anc in tmp:
	diagEntry[anc] = []
	genesAnc[anc] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[anc])


# On compare toutes les especes entre elles
for (e1,e2,toStudy) in dicLinks:
	statsDiags = []
	print >> sys.stderr, "Extraction des diagonales entre %s et %s ..." % (e1,e2),
	for ((c1,d1),(c2,d2),s) in utils.myDiags2.calcDiags(dicGenomes[e1], dicGenomes[e2], genesAnc[phylTree.dicParents[e1][e2]], arguments["minimalLength"], \
		arguments["fusionThreshold"], arguments["sameStrand"] and (e1 not in genesAnc) and (e2 not in genesAnc), arguments["keepOnlyOrthos"]):

		statsDiags.append(len(d1))
		pack1 = (e1,c1,tuple(d1))
		pack2 = (e2,c2,tuple(d2))
		dic1 = dicGenomes[e1].lstGenes[c1]
		dic2 = dicGenomes[e2].lstGenes[c2]
		for anc in toStudy:
			tmp = genesAnc[anc].dicGenes
			if dic1[d1[0]].names[0] in tmp:
				tmp = tuple(tmp[dic1[i1].names[0]][1] for i1 in d1)
			else:
				tmp = tuple(tmp[dic2[i2].names[0]][1] for i2 in d2)
			diagEntry[anc].append( (pack1,pack2,tmp,s) )
	
	print >> sys.stderr, utils.myMaths.myStats.txtSummary(statsDiags), "OK"


# Traitement final
for (anc,lst) in diagEntry.iteritems():

	if arguments["showProjected"]:
		f = utils.myTools.tsvWriter(arguments["OUT.projDiags"] % phylTree.fileName[anc])
		print >> sys.stderr, "Impression des diagonales projetees de %s ..." % anc,
		s = []
		for ((e1,c1,d1),(e2,c2,d2),da,ds) in lst:
			s.append( len(d1) )
			res = [anc,len(da), \
				e1,c1," ".join( [dicGenomes[e1].lstGenes[c1][i1].names[0] for i1 in d1] ), \
				e2,c2," ".join( [dicGenomes[e2].lstGenes[c2][i2].names[0] for i2 in d2] ), \
				utils.myTools.printLine(da, " "),utils.myTools.printLine(ds, " ")]
			f.csvobject.writerow(res)
		f.file.close()

		print >> sys.stderr, utils.myMaths.myStats.txtSummary(s), "OK"


	if arguments["showAncestral"]:
	
		if arguments["cutLongestPath"]:
			newlst = getLongestDiags(lst)
		else:
			newlst = ( (da, ds, ((e1,c1),(e2,c2))) for ((e1,c1,_),(e2,c2,_),da,ds) in lst )
		
		print >> sys.stderr, "Impression des diagonales ancestrales de %s ..." % anc,
		ft = utils.myTools.tsvWriter(arguments["OUT.ancDiags"] % phylTree.fileName[anc])
		s = []
		for (da,ds,esp) in newlst:
			s.append( len(da) )
			res = [anc,len(da),utils.myTools.printLine(da, " "),utils.myTools.printLine(ds, " ")]
			
			res.append( "|".join(["%s/%s" % x for x in esp]) )

			if arguments["searchUndetectedSpecies"]:
				res.append( "|".join(["%s/%s" % x for x in findNewSpecies(da, esp, anc)]) )
			
			f.csvobject.writerow(res)
		f.file.close()

		print >> sys.stderr, utils.myMaths.myStats.txtSummary(s), "OK"

