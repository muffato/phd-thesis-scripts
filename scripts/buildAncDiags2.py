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
import utils.myFile
import utils.myTools
import utils.myMaths
import utils.myDiags


# Comparaison de deux especes
def compare(e1, e2, toStudy):
	statsDiags = []
	print >> sys.stderr, "Extraction des diagonales entre %s et %s ..." % (e1,e2),
	for ((c1,d1),(c2,d2),s) in utils.myDiags.calcDiags(dicGenomes[e1], dicGenomes[e2], genesAnc[phylTree.dicParents[e1][e2]], arguments["minimalLength"], \
		arguments["fusionThreshold"], arguments["sameStrand"] and (e1 not in genesAnc) and (e2 not in genesAnc), arguments["keepOnlyOrthos"]):

		l = len(d1)
		if l < arguments["minimalLength"]:
			continue

		statsDiags.append(l)
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


# Integration de toutes les diagonales d'un ancetre
def getLongestDiags(anc):

	# On construit la table "gene" -> "liste des diagonales qui le contiennent"
	# On s'en sert pour avoir les listes de diagonales chevauchantes
	dic = collections.defaultdict(list)
	oldDiags = diagEntry[anc]
	diags = range(len(oldDiags))
	combin = utils.myTools.myCombinator()
	lstGenesAncAnc = genesAnc[anc].lstGenes[None]

	# Etablit la liste des chromosomes de chaque espece, avec les comptes
	def getSpeciesCount(diag, synt):
		if not arguments["searchUndetectedSpecies"]:
			return []
		count = collections.defaultdict(int)
		countD = collections.defaultdict(int)
		for e in listSpecies:

			# L'ancetre commun
			a = phylTree.dicParents[anc][e]
			lstGenesAncA = genesAnc[a].lstGenes[None]
			genome = dicGenomes[e]

			for i in diag:
				# Les noms associes au gene ancestral
				names = lstGenesAncAnc[i].names
				if a != anc:
					names = set(genesAnc[a].getOtherNames(names[-1]) + names)
				tmp = set(c for (c,_) in genome.getPosition(names) if c not in genome.lstNoneS)
				tmp1 = set(utils.myGenomes.commonChrName(c.replace("_random","")) for c in tmp if c in genome.lstRandS)
				tmp.difference_update(genome.lstRandS)
				tmp.update(tmp1)

				tmpC = set(c for c in tmp if c in genome.lstChrS)
				tmpS = set(c for c in tmp if c in genome.lstScaffS)
				assert len(tmp) == (len(tmpC) + len(tmpS))

				if len(tmp) == 0:
					pass
				elif len(tmp) == 1:
					count[(e,tmp.pop())] += 1
				elif len(tmpC) == 1:
					count[(e,tmpC.pop())] += 1
				else:
					for x in tmp:
						(count if (e,x) in synt else countD)[(e,x)] += 1
		s = 1 if len(diag) == 1 else 2
		return [(e,c,i+countD[(e,c)]) for ((e,c),i) in count.iteritems() if i >= s]

	for (i,(_,_,da,ds)) in enumerate(oldDiags):
		diags[i] = zip(da, ds)
		# Combinaison selon les genes ancestraux
		for j in da:
			dic[j].append(i)
		combin.addLink([i])

	for j in xrange(len(lstGenesAncAnc)):
		if j in dic:
			combin.addLink(dic[j])
		else:
			# Envoi des genes singletons
			yield ([j], [0], [], getSpeciesCount([j], []))

	del combin.dic

	for g in combin:
		# Les paires de gene synteniques
		synt = collections.defaultdict(set)
		for i in g:
			for ((g1,_),(g2,_)) in utils.myTools.myIterator.slidingTuple(diags[i]):
				synt[(g1,g2)].add( oldDiags[i][0][:2] )
				synt[(g1,g2)].add( oldDiags[i][1][:2] )
				synt[(g2,g1)] = synt[(g1,g2)]

		ss = set()
		for i in g:
			for (d,_) in diags[i]:
				ss.add(d)
		ss = len(ss)
		tt = 0

		# On casse les carrefours de diagonales les uns apres les autres
		gr = utils.myDiags.WeightedDiagGraph([diags[i] for i in g])
		# Les diagonales resultat
		for (res,strand) in gr.getBestDiags():

			# On rajoute la liste des especes qui soutiennent la diagonale
			ok = collections.defaultdict(set)
			for x in utils.myTools.myIterator.slidingTuple(res):
				for y in synt[x]:
					ok[y].update(x)
			tt += len(res)
			yield (res, strand, [(e,c,len(s)) for ((e,c),s) in ok.iteritems()], getSpeciesCount(res, set(ok)))

		assert tt == ss


# Impression des diagonales projetees
def printProjected(anc):
	print >> sys.stderr, "Impression des diagonales projetees de %s ..." % anc,
	f = utils.myFile.myTSV.writer(arguments["OUT.projDiags"] % phylTree.fileName[anc])

	s = []
	for ((e1,c1,d1),(e2,c2,d2),da,ds) in diagEntry[anc]:
		s.append( len(d1) )
		res = [anc,len(da), \
			e1,c1," ".join("/".join(dicGenomes[e1].lstGenes[c1][i1].names) for i1 in d1), \
			e2,c2," ".join("/".join(dicGenomes[e2].lstGenes[c2][i2].names) for i2 in d2), \
			utils.myFile.myTSV.printLine(da, " "),utils.myFile.myTSV.printLine(ds, " ")]
		f.csvobject.writerow(res)
	f.file.close()

	print >> sys.stderr, utils.myMaths.myStats.txtSummary(s), "OK"


# Impression des diagonales (integrees) ancestrales
def printAncestral(anc):
	if arguments["cutLongestPath"]:
		newlst = getLongestDiags(anc)
	else:
		newlst = ( (da, ds, ((e1,c1),(e2,c2))) for ((e1,c1,_),(e2,c2,_),da,ds) in diagEntry[anc] )

	print >> sys.stderr, "Impression des diagonales ancestrales de %s ..." % anc,
	f = utils.myFile.myTSV.writer(arguments["OUT.ancDiags"] % phylTree.fileName[anc])
	s = []
	singles = 0
	for (da,ds,esp,esp2) in newlst:

		if len(da) > 1:
			s.append( len(da) )
		else:
			singles += 1

		res = [anc,len(da),utils.myFile.myTSV.printLine(da, " "),utils.myFile.myTSV.printLine(ds, " ")]

		res.append( "|".join("%s/%s/%d" % x for x in esp) )

		if arguments["searchUndetectedSpecies"]:
			res.append( "|".join("%s/%s/%d" % x for x in esp2) )

		f.csvobject.writerow(res)
	f.file.close()

	print >> sys.stderr, utils.myMaths.myStats.txtSummary(s), "+ %d singletons OK" % singles


# Arguments
arguments = utils.myTools.checkArgs( \
	[("phylTree.conf",file), ("target",str)], \
	[("fusionThreshold",int,-1), ("minimalLength",int,2), ("sameStrand",bool,True), ("keepOnlyOrthos",bool,False),
	("useOutgroups",bool,False), \
	("showProjected",bool,True), ("showAncestral",bool,True), ("searchUndetectedSpecies",bool,True), ("cutLongestPath",bool,True), \
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
assert len(arguments["target"]) > 0, "Aucune cible indiquee pour l'extraction des diagonales"

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
	compare(e1, e2, toStudy)

# Traitement final
for (anc,lst) in diagEntry.iteritems():

	if arguments["showProjected"]:
		printProjected(anc)

	if arguments["showAncestral"]:
		printAncestral(anc)

