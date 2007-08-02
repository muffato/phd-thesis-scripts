#! /users/ldog/muffato/python -OO

__doc__ = """
	Telecharge depuis le site d'Ensembl les fichiers utiles pour creer:
		- les listes de genes
		- les listes d'orthologues
		- les listes de paralogues
"""


##################
# INITIALISATION #
##################

# Librairies
import os
import sys
import utils.myPhylTree
import utils.myTools


#############
# FONCTIONS #
#############

#
# Permet de telecharger, decompresser et lire a la volee un fichier
#
def fileIterator(nom):
	stdout = utils.myTools.myOpenFile( ("%s/%s" % (options["IN.EnsemblURL"],nom)).replace("XXX", str(options["releaseID"])) , "r")
	#(stdin,stdout,stderr) = os.popen3( ("wget %s/%s -O - | gunzip" % (options["IN.EnsemblURL"],nom)).replace("XXX", str(options["releaseID"])) )
	#stdin.close()
	#stderr.close()
	tmp = ""
	for ligne in stdout:
		if ligne[-2] == '\\':
			tmp = ligne[:-2]
		else:
			yield tmp + ligne[:-1]
			tmp = ""
	stdout.close()


#
# Ensembl n'utilise pas exactement le meme arbre que nous ...
#
def mkEnsemblPhylAdjustment(oldAnc, theoryAnc):

	if (theoryAnc == "Boreoeutheria") and (oldAnc == "Eutheria"):
		return "Boreoeutheria"

	if (theoryAnc == "FishInterm") and (oldAnc == "Percomorpha"):
		return "FishInterm"

	if oldAnc == "Smegmamorpha":
		return "Percomorpha"

	return oldAnc


#
# Lit un fichier d'homologies d'Ensembl et le remet au format habituel
#
def proceedFile(fin, fout, foutR):
	fout = utils.myTools.myOpenFile(fout, 'w')
	if foutR != None:
		foutR = utils.myTools.myOpenFile(foutR, 'w')
	nb1 = 0
	nb2 = 0
	try:
		fin = fileIterator(fin)
		c = fin.next().split('\t')
		i1 = c.index(esp1B, 2) # Pour eviter le pb des genes paralogues pour lesquels anc=esp
		i2 = c.index(esp2B, i1+1)
		while True:
			nb1 += 1
			newAnc = mkEnsemblPhylAdjustment(c[1], theoryAnc)
			if newAnc == theoryAnc:
				nb2 += 1
			
			if options["releaseID"] in [44,45]:
				r = (c[7],c[4],c[24], c[i1+15],c[i1+12],c[i1+32], newAnc, c[i1+4],c[i1+5], c[i2+4],c[i2+5], c[0])
			else:
				r = (c[3],c[25],c[44], c[i2-8],c[i2+14],c[i2+33], newAnc, c[29],c[30], c[i2+18],c[i2+19], c[0])
			
			print >> fout, "\t".join(r)
			if foutR != None:
				print >> foutR, "\t".join([r[i] for i in (3,4,5, 0,1,2, 6, 9,10, 7,8, 11)])
			
			c = fin.next().split('\t')
	except StopIteration:
		pass
	fout.close()
	if foutR != None:
		foutR.close()
	return (nb1,nb2)



########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("ancestr",str,""), ("",str,""), \
	("OUT.genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("OUT.fullGenesFile",str,"~/work/data/genes/full/genes.%s.list.bz2"), \
	("OUT.orthosFile",str,"~/work/data/orthologs/orthos.%s.%s.list.bz2"), \
	("IN.UCSC-URL",str,"http://hgdownload.cse.ucsc.edu/goldenPath/XX/database/"), \
	("IN.genesFile",str,"genscan.txt.gz"), \
	("IN.parasFile",str,"xenoRefFlat.txt.gz"), \
	("IN.orthosFile",str,"compara_mart_homology_XXX/compara_%s_%s_orthologs__orthologs__main.txt.table.gz")], \
	__doc__ \
)



# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])


# Les fichiers de genes
for esp in phylTree.listSpecies:
	print >> sys.stderr, "Telechargement de la liste des genes de %s ..." % esp,
	fo1 = utils.myTools.myOpenFile(options["OUT.genesFile"] % phylTree.fileName[esp], 'w')
	fo2 = utils.myTools.myOpenFile(options["OUT.fullGenesFile"] % phylTree.fileName[esp], 'w')
	tmp = esp.lower().split()
	nb1 = 0
	nb2 = 0
	for ligne in fileIterator(options["IN.genesFile"] % (tmp[0][0] + tmp[1])):
		c = ligne.split('\t')
		if ("RNA" not in c[3]) and ("pseudogene" not in c[3]):
			print >> fo1, "\t".join( [c[11],c[7],c[8],c[9],c[1]] )
			nb1 += 1
		print >> fo2, "\t".join( [c[11],c[7],c[8],c[9],c[1]] )
		nb2 += 1
	
	fo1.close()
	fo2.close()
	print >> sys.stderr, "%d/%d genes OK" % (nb1,nb2)



# Les noms utilises dans les fichiers "Homo Sapiens" -> "hsap"
nomReel = []
dicNomsReels = {}
for esp in sorted(phylTree.listSpecies):
	tmp = esp.lower().split()
	tmp = tmp[0][0] + tmp[1][:3]
	nomReel.append(tmp)
	dicNomsReels[tmp] = esp


# On genere les fichiers d'homologues
for (esp1,esp2) in utils.myTools.myIterator.tupleOnUpperList(nomReel):

	esp1B = dicNomsReels[esp1]
	esp2B = dicNomsReels[esp2]

	theoryAnc = phylTree.dicParents[esp1B][esp2B]
		
	if esp1 == esp2:
		print >> sys.stderr, "Telechargement de la liste des genes paralogues de %s ..." % esp1B,
		nb,_ = proceedFile(options["IN.parasFile"] % (esp1,esp1), options["OUT.paras1File"] % phylTree.fileName[esp1B], None)
		print >> sys.stderr, "%d genes OK" % nb

	else:
		print >> sys.stderr, "Telechargement de la liste des genes orthologues entre %s et %s ..." % (esp1B,esp2B),
		nb1,nb2 = proceedFile(options["IN.orthosFile"] % (esp1,esp2), options["OUT.orthosFile"] % (phylTree.fileName[esp1B],phylTree.fileName[esp2B]), options["OUT.orthosFile"] % (phylTree.fileName[esp2B],phylTree.fileName[esp1B]))
		print >> sys.stderr, "%d/%d genes OK" % (nb2, nb1)

		print >> sys.stderr, "Telechargement de la liste des genes paralogues entre %s et %s ..." % (esp1B,esp2B),
		nb,_ = proceedFile(options["IN.parasFile"] % (esp1,esp2), options["OUT.paras2File"] % (phylTree.fileName[esp1B],phylTree.fileName[esp2B]), options["OUT.paras2File"] % (phylTree.fileName[esp2B],phylTree.fileName[esp1B]))
		print >> sys.stderr, "%d genes OK" % nb
	

