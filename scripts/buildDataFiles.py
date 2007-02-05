#! /users/ldog/muffato/python -OO

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
import utils.myBioObjects
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myDiags


#############
# FONCTIONS #
#############

def calcDiags(e1, e2):

	# La fonction qui permet de stocker les diagonales sur les ancetres
	def combinDiag(c1, c2, d1, d2):
		global diagEntry

		if len(d1) < options["minimalLength"]:
			return
		
		dn1 = [g1.lstGenes[c1][i].names[0] for i in d1]
		dn2 = [g2.lstGenes[c2][i].names[0] for i in d2]

		for tmp in toStudy:
			diagEntry[tmp].append( ((e1,c1,dn1), (e2,c2,dn2)) )
	
	# Ecrire un genome en suite de genes ancestraux
	def translateGenome(genome):
		newGenome = {}
		for c in genome.lstChr:
			newGenome[c] = [(genesAnc.dicGenes.get(g.names[0], (0,-1))[1],g.strand) for g in genome.lstGenes[c]]
			if options["keepOnlyOrthos"]:
				newGenome[c] = [x for x in newGenome[c] if x[0] != -1]
		return newGenome

	# Les noeuds de l'arbre entre l'ancetre et les especes actuelles
	toStudy = [phylTree.getFirstParent(e1,e2)]
	for tmp in phylTree.items:
		s = phylTree.species[tmp]
		if (e1 in s) and (e2 not in s):
			toStudy.append( tmp )
		elif (e2 in s) and (e1 not in s):
			toStudy.append( tmp )

	# Chargement des orthologues
	#calcDiags(phylTree.commonNames[phylTree.listSpecies[i]][0], phylTree.commonNames[phylTree.listSpecies[j]][0])
	genesAnc = utils.myGenomes.GenomeFromOrthosList(options["orthosFile"] % (phylTree.commonNames[e1][0],phylTree.commonNames[e2][0]))
	newLoc = [[] for x in xrange(len(genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr]))]
	del genesAnc.lstGenes

	g1 = phylTree.dicGenomes[e1]
	g2 = phylTree.dicGenomes[e2]
	newGen = translateGenome(g1)
	tmp = translateGenome(g2)
	
	for c in g2.lstChr:
		for i in xrange(len(tmp[c])):
			(ianc,s) = tmp[c][i]
			if ianc != -1:
				newLoc[ianc].append( (c,i,s) )
	
	utils.myDiags.iterateDiags(newGen, newLoc, options["fusionThreshold"], options["sameStrand"], combinDiag)


# Cherche si la diagonale est localisee sur un chromosome d'une espece, non ordonnee
def findNewSpecies(d, esp, anc):
	
	# Les especes pour lesquelles la diagonale n'a pas ete detectee
	lstEsp = set(phylTree.listSpecies).difference([e for (e,_) in esp])
	res = []
	for e in lstEsp:
		# L'ancetre commun
		a = phylTree.getFirstParent(anc, e)
		# On va intersecter les chromosomes des orthologues de chaque gene
		poss = set([])
		for i in d:
			# Le gene chez l'ancetre commun
			g = genesAnc[a].getPosition(genesAnc[anc].lstGenes[utils.myGenomes.Genome.defaultChr][i])
			# Cas d'un gene specifique de la lignee
			if len(g) == 0:
				continue
			# Le gene dans l'autre espece
			tmp = [c for (c,_) in phylTree.dicGenomes[e].getPosition(genesAnc[a].lstGenes[utils.myGenomes.Genome.defaultChr][g[0][1]])]
			# Gene non trouve, on passe au suivant
			if len(tmp) == 0:
				continue
			# Sinon, on intersecte
			if len(poss) == 0:
				poss.update(tmp)
			else:
				poss.intersection_update(tmp)
				# Utile de continuer ?
				if len(poss) == 0:
					break
		# S'il en reste un, c'est gagne !
		if len(poss) != 0:
			res.append( (e,poss.pop()) )
	
	return res


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("OUT.genesFile",str,"~/work/data42/genes/genes.%s.list.bz2"), \
	("OUT.fullGenesFile",str,"~/work/data42/genes/full/genes.%s.list.bz2"), \
	("OUT.orthosFile",str,"~/work/data42/orthologs/orthos.%s.%s.list.bz2"), \
	("OUT.parasFile",str,"~/work/data42/paralogs/paras.%s.list.bz2"), \
	("IN.genesFile",str,""), \
	("IN.transcriptsFile",str,""), \
	("IN.orthosFile",str,""), \
	("IN.parasFile",str,"")], \
	__doc__ \
)


# L'arbre phylogenetique
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])

# Les noms utilises dans les fichiers "Homo Sapiens" -> "hsapiens"
nomReel = {}
for esp in phylTree.listSpecies:
	tmp = esp.lower().split()
	nomReel[esp] = tmp[0][0] + tmp[1]

# La table de traduction ENSGTXXX.N ENS
dicTranscrits = {}
for esp in phylTree.listSpecies:

	print >> sys.stderr, "Chargement de la table des transcrits de %s ..." % esp,
	f = utils.myTools.myOpenFile(options["IN.transcriptsFile"] % nomReel[esp], 'r')
	for ligne in f:
		c = ligne[:-1].split('\t')
		if len(c) <= 6:
			continue
		dicTranscrits[c[6]] = c[2]
	f.close()
	print >> sys.stderr, "OK"


# On genere les fichiers de listes de genes + Filtre scaffold/chromosome
forbiddenTokens = ["rand", "Un", "affold", "ont", "tig", "v6", "h", "MT", "Mt"]
forbiddenStartingTokens = ["NT", "U", "1099", "c", "E"]

for esp in phylTree.listSpecies:
	continue	
	print >> sys.stderr, "Mise en forme de la liste des genes de %s ..." % esp,
	fi = utils.myTools.myOpenFile(options["IN.genesFile"] % nomReel[esp], 'r')
	fo1 = utils.myTools.myOpenFile(options["OUT.fullGenesFile"] % phylTree.commonNames[esp][0], 'w')
	fo2 = utils.myTools.myOpenFile(options["OUT.genesFile"] % phylTree.commonNames[esp][0], 'w')
	
	nb1 = 0
	nb2 = 0
	for ligne in fi:
		c = ligne[:-1].split('\t')
		if c[1] == "\\N":
			continue
		gene = "\t".join( [c[11],c[7],c[8],c[9],c[1]] )
		print >> fo1, gene
		nb1 += 1

		for t in forbiddenTokens:
			if t in c[11]:
				break
		else:
			for t in forbiddenStartingTokens:
				if c[11].startswith(t):
					break
			else:
				print >> fo2, gene
				nb2 += 1

	fi.close()
	fo1.close()
	fo2.close()
	print >> sys.stderr, "%d/%d genes OK" % (nb2,nb1)

# On genere les fichiers d'ortologues
for esp1 in phylTree.listSpecies:
	for esp2 in phylTree.listSpecies:


		if esp1 == esp2:
			print >> sys.stderr, "Mise en forme de la liste des genes paralogues de %s ..." % esp1,
			fi = utils.myTools.myOpenFile(options["IN.parasFile"] % (nomReel[esp1],nomReel[esp1]), 'r')
			fo = utils.myTools.myOpenFile(options["OUT.parasFile"] % phylTree.commonNames[esp1][0], 'w')
			anc = 10000000
		else:
			print >> sys.stderr, "Mise en forme de la liste des genes orthologues entre %s et %s..." % (esp1,esp2),
			fi = utils.myTools.myOpenFile(options["IN.orthosFile"] % (nomReel[esp1],nomReel[esp2]), 'r')
			fo = utils.myTools.myOpenFile(options["OUT.orthosFile"] % (phylTree.commonNames[esp1][0],phylTree.commonNames[esp2][0]), 'w')
			anc = phylTree.ages[phylTree.getFirstParent(esp1, esp2)]

		nb = 0
		for ligne in fi:
			c = ligne[:-1].split('\t')
			if c[1] == "\\N":
				continue
			print phylTree.ages[c[16]], anc, (c[2],dicTranscrits[c[3]],c[20], c[4],dicTranscrits[c[5]],c[21], c[23],c[24],c[27],c[28], c[15],c[16])
			if phylTree.ages[c[16]] <= anc:
				nb += 1
				print >> fo, (c[2],dicTranscrits[c[3]],c[20], c[4],dicTranscrits[c[5]],c[21], c[23],c[24],c[27],c[28], c[15],c[16])

		fi.close()
		fo.close()
		print >> sys.stderr, "%d genes OK" % nb


