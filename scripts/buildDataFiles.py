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
		f = utils.myTools.myOpenFile(options["IN.ensemblDirectory"] + "/" + nom, 'r')
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
	("OUT.orthosFile",str,"~/work/data42/orthologs/orthos.%s.%s.list.bz2"), \
	("OUT.parasFile",str,"~/work/data42/paralogs/paras.%s.list.bz2"), \
	("IN.ensemblDirectory",str,"~/workspace/ftp.ensembl.org/pub/release-41/mart_41/data/mysql/ensembl_mart_41/"), \
	("IN.genesFile",str,"%s_gene_ensembl__gene__main.txt.table.gz"), \
	("IN.transcriptsFile",str,"%s_gene_ensembl__transcript__main.txt.table.gz"), \
	("IN.orthosFile",str,"%s_gene_ensembl__homologs_%s__dm.txt.table.gz"), \
	("IN.parasFile",str,"%s_gene_ensembl__paralogs_%s__dm.txt.table.gz")], \
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


for esp in phylTree.listSpecies:
	print >> sys.stderr, "Mise en forme de la liste des genes et des familles de %s ..." % esp,
	fo = utils.myTools.myOpenFile(options["OUT.genesFile"] % phylTree.fileName[esp], 'w')
	
	nb = 0
	for ligne in fileIterator(options["IN.genesFile"] % nomReel[esp]):
		c = ligne.split('\t')
		if c[1] == "\\N":
			continue
		print >> fo, "\t".join( [c[11],c[7],c[8],c[9],dicFam[c[1]],c[1]] )
		
		nb += 1
	
	fo.close()
	print >> sys.stderr, "%d genes OK" % nb



# Ensembl n'utilise pas exactement le meme arbre que nous ...
def mkEnsemblPhylAdjustment(oldAnc, theoryAnc):

	if (theoryAnc == "Boreoeutheria") and (oldAnc == "Eutheria"):
		return "Boreoeutheria"

	if (theoryAnc == "FishInterm") and (oldAnc == "Percomorpha"):
		return "FishInterm"

	if oldAnc == "Smegmamorpha":
		return "Percomorpha"

	return oldAnc



# On genere les fichiers d'ortologues
for esp1 in phylTree.listSpecies:
	for esp2 in phylTree.listSpecies:

		theoryAnc = phylTree.getFirstParent(esp1, esp2)
		
		if esp1 == esp2:
			print >> sys.stderr, "Mise en forme de la liste des genes paralogues de %s ..." % esp1,
			fi = options["IN.parasFile"] % (nomReel[esp1],nomReel[esp1])
			fo = utils.myTools.myOpenFile(options["OUT.parasFile"] % phylTree.fileName[esp1], 'w')
			nb = 0
			for ligne in fileIterator(fi):
				c = ligne.split('\t')
				if c[2] == "\\N":
					continue
				nb += 1
				newAnc = mkEnsemblPhylAdjustment(c[16], theoryAnc)
				print >> fo, "\t".join([c[2],dicTranscrits[c[3]],c[20], c[4],dicTranscrits[c[5]],c[21], newAnc, c[24],c[25],c[28],c[29]])

			fo.close()
			print >> sys.stderr, "%d genes OK" % nb

		else:
			print >> sys.stderr, "Extraction des genes orthologues entre %s et %s ..." % (esp1,esp2),
			fi = options["IN.orthosFile"] % (nomReel[esp1],nomReel[esp2])
			fo = utils.myTools.myOpenFile(options["OUT.orthosFile"] % (phylTree.fileName[esp1],phylTree.fileName[esp2]), 'w')
			anc = phylTree.ages[theoryAnc]
			
			nb1 = 0
			nb2 = 0
			for ligne in fileIterator(fi):
				c = ligne.split('\t')
				if c[2] == "\\N":
					continue
				newAnc = mkEnsemblPhylAdjustment(c[16], theoryAnc)
				print >> fo, "\t".join([c[2],dicTranscrits[c[3]],c[20], c[4],dicTranscrits[c[5]],c[21], newAnc, c[24],c[25],c[28],c[29], c[15]])
				nb1 += 1
				if newAnc == theoryAnc:
					nb2 += 1

			fo.close()
			print >> sys.stderr, "%d/%d genes OK" % (nb2, nb1)


