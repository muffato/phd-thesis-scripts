#! /users/ldog/muffato/python

__doc__ = """
Extrait toutes les diagonales entre chaque paire d'especes.
Les diagonales apportent les genes qui etaient sur un meme chromosome
  depuis leur ancetre commun dans les deux lignees.
"""


import sys
import itertools
import collections
import multiprocessing

import utils.myPhylTree
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myDiags

def projDiagsReader(filename):
	f = utils.myFile.myTSV.fileReader(filename)
	for (_,_,esp1,chr1,_,esp2,chr2,_,da,ds) in f.csvobject:
		yield ([int(x) for x in da.split()],[int(x) for x in ds.split()],((esp1,chr1),(esp2,chr2)))
	f.file.close()


def getLongestDiags(oldDiags):

	# On construit la table "gene" -> "liste des diagonales qui le contiennent"
	# On s'en sert pour avoir les listes de diagonales chevauchantes
	dic = collections.defaultdict(list)
	diags = list(oldDiags)
	combin = utils.myTools.myCombinator()
	for (i,(da,_,_)) in enumerate(diags):
		# Combinaison selon les genes ancestraux
		for j in da:
			dic[j].append(i)
		combin.addLink([i])
	
	for s in dic:
		combin.addLink(dic[s])
	
	del combin.dic

	for g in combin:
		# On casse les carrefours de diagonales les uns apres les autres
		gr = utils.myDiags.WeightedDiagGraph([zip(diags[i][0],diags[i][1]) for i in g])
		# Les diagonales resultat
		for (res,strand) in gr.getBestDiags():
			# Filtre de taille
			if len(res) < arguments["minimalLength"]:
				continue
			# On rajoute la liste des especes qui soutiennent la diagonale
			ok = set()
			# Test de chaque diagonale
			for i in g:
				for (g1,g2) in utils.myTools.myIterator.slidingTuple(diags[i][0]):
					try:
						i1 = res.index(g1)
						i2 = res.index(g2)
					except ValueError:
						continue
					# Si la diagonale 'i' soutient la diagonale ancestrale, on rajoute les especes dont elle est tiree
					if abs(i1-i2) == 1:
						ok.update(diags[i][2])
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
	[("searchUndetectedSpecies",bool,True), ("cutLongestPath",bool,True), ("minimalLength",int,2), \
	("MP.poolSize",int,None), ("MP.chunkSize",int,2), \
	("IN.projDiags",str,"proj/diags.%s.list.bz2"), \
	("OUT.ancDiags",str,"anc/diags.%s.list.bz2"), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
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

def processAncestr(anc):

	lst = projDiagsReader(arguments["IN.projDiags"] % phylTree.fileName[anc])

	if arguments["cutLongestPath"]:
		lst = getLongestDiags(lst)

	f = utils.myFile.myTSV.fileWriter(arguments["OUT.ancDiags"] % phylTree.fileName[anc])
	s = []
	for (da,ds,esp) in lst:
		s.append( len(da) )
		res = [anc,len(da),utils.myFile.myTSV.printLine(da, " "),utils.myTools.printLine(ds, " ")]
		
		res.append( "|".join(["%s/%s" % x for x in esp]) )

		if arguments["searchUndetectedSpecies"]:
			res.append( "|".join(["%s/%s" % x for x in findNewSpecies(da, esp, anc)]) )
		
		f.csvobject.writerow(res)
	f.file.close()

	return ("Statistiques de %s : %s" % (anc, utils.myMaths.myStats.txtSummary(s)))

goodAnc = [anc for anc in phylTree.listAncestr if utils.myFile.hasAccess(arguments["IN.projDiags"] % phylTree.fileName[anc])]
pool = multiprocessing.Pool(processes=arguments["MP.poolSize"])
for res in pool.imap_unordered(processAncestr, goodAnc, chunksize=arguments["MP.chunkSize"]):
	print >> sys.stderr, res

