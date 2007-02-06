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


def fileIterator(nom):
	try:
		f = utils.myTools.myOpenFile(nom, 'r')
		tmp = ""
		for ligne in f:
			if ligne[-2] == '\\':
				tmp = ligne[:-2]
			else:
				yield tmp + ligne[:-1]
				tmp = ""
		f.close()
	except IOError:
		return



########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("OUT.genesFile",str,"~/work/data42/genes/genes.%s.list.bz2"), \
	("OUT.fullGenesFile",str,"~/work/data42/genes/full/genes.%s.list.bz2"), \
	("OUT.orthosFile",str,"~/work/data42/orthologs/orthos.%s.%s.list.bz2"), \
	("OUT.fullOrthosFile",str,"~/work/data42/orthologs/full/orthos.%s.%s.list.bz2"), \
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
dicFam = {}
for esp in phylTree.listSpecies:

	print >> sys.stderr, "Chargement de la table des transcrits de %s ..." % esp,
	for ligne in fileIterator(options["IN.transcriptsFile"] % nomReel[esp]):
		c = ligne.split('\t')
		dicTranscrits[c[6]] = c[2]
		dicFam[c[5]] = c[23]
	print >> sys.stderr, "OK"


# On genere les fichiers de listes de genes + Filtre scaffold/chromosome
forbiddenTokens = ["rand", "Un", "affold", "ont", "tig", "v6", "h", "MT", "Mt"]
forbiddenStartingTokens = ["NT", "U", "1099", "c", "E"]

for esp in phylTree.listSpecies:
	print >> sys.stderr, "Mise en forme de la liste des genes et des familles de %s ..." % esp,
	fo1 = utils.myTools.myOpenFile(options["OUT.fullGenesFile"] % phylTree.commonNames[esp][0], 'w')
	fo2 = utils.myTools.myOpenFile(options["OUT.genesFile"] % phylTree.commonNames[esp][0], 'w')
	
	nb1 = 0
	nb2 = 0
	for ligne in fileIterator(options["IN.genesFile"] % nomReel[esp]):
		c = ligne.split('\t')
		if c[1] == "\\N":
			continue
		gene = "\t".join( [c[11],c[7],c[8],c[9],dicFam[c[1]],c[1]] )
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

	fo1.close()
	fo2.close()
	print >> sys.stderr, "%d/%d genes OK" % (nb2,nb1)



# On genere les fichiers d'ortologues
for esp1 in phylTree.listSpecies:
	for esp2 in phylTree.listSpecies:


		if esp1 == esp2:
			print >> sys.stderr, "Mise en forme de la liste des genes paralogues de %s ..." % esp1,
			fi = options["IN.parasFile"] % (nomReel[esp1],nomReel[esp1])
			fo = utils.myTools.myOpenFile(options["OUT.parasFile"] % phylTree.commonNames[esp1][0], 'w')
			nb = 0
			for ligne in fileIterator(fi):
				c = ligne.split('\t')
				if c[2] == "\\N":
					continue
				nb += 1
				print >> fo, "\t".join([c[2],dicTranscrits[c[3]],c[20], c[4],dicTranscrits[c[5]],c[21], c[16], c[24],c[25],c[28],c[29]])

			fo.close()
			print >> sys.stderr, "%d genes OK" % nb

		else:
			print >> sys.stderr, "Extraction des genes orthologues entre %s et %s ..." % (esp1,esp2),
			fi = options["IN.orthosFile"] % (nomReel[esp1],nomReel[esp2])
			fo1 = utils.myTools.myOpenFile(options["OUT.orthosFile"] % (phylTree.commonNames[esp1][0],phylTree.commonNames[esp2][0]), 'w')
			fo2 = utils.myTools.myOpenFile(options["OUT.fullOrthosFile"] % (phylTree.commonNames[esp1][0],phylTree.commonNames[esp2][0]), 'w')
			anc = phylTree.ages[phylTree.getFirstParent(esp1, esp2)]
			
			nb1 = 0
			nb2 = 0
			for ligne in fileIterator(fi):
				c = ligne.split('\t')
				if c[2] == "\\N":
					continue
				s = "\t".join([c[2],dicTranscrits[c[3]],c[20], c[4],dicTranscrits[c[5]],c[21], c[15], c[24],c[25],c[28],c[29]])
				nb1 += 1
				print >> fo2, "%s\t%s" % (s, c[16])
				if phylTree.ages[c[16]] <= anc:
					nb2 += 1
					print >> fo1, s

			fo1.close()
			fo2.close()
			print >> sys.stderr, "%d/%d genes OK" %(nb2, nb1)


