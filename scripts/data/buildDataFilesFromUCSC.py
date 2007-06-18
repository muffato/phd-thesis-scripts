#! /users/ldog/muffato/python -OO

__doc__ = """
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

def common(x):
	try:
		return int(x)
	except Exception:
		return x




# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesLizard","genesChicken", "chains"], \
	[], \
	__doc__ \
)

genomeChicken = utils.myGenomes.EnsemblGenome(noms_fichiers["genesChicken"])
genomeLizard = utils.myGenomes.EnsemblGenome(noms_fichiers["genesLizard"])

links = {}

lizChr = None
f = utils.myTools.myOpenFile(noms_fichiers["chains"], "r")
for l in f:
	c = l.split()
	if "net" in l:
		lizChr = c[1]
		continue
	if 'fill' not in l:
		continue
	#lizChr = c[3]
	liz1 = int(c[1])
	lizSize = int(c[2])

	lizGt = list(genomeLizard.getGenesNearB(lizChr, liz1+lizSize/2, lizSize/2))
	if len(lizGt) > 0:
		print lizGt
		print
f.close()

sys.exit(0)


f = utils.myTools.myOpenFile(noms_fichiers["chains"], "r")
for l in f:
	c = l.split()
	lizChr = c[1]
	liz1 = int(c[2])
	liz2 = int(c[3])
	chickChr = common(c[4][3:])
	chickStrand = c[7]
	chick1 = int(c[5])
	chick2 = int(c[6])
	lizGt = list(genomeLizard.getGenesNearB(lizChr, (liz1+liz2)/2, (liz2-liz1)/2))
	if len(lizGt) > 0:
		print lizGt
		print
	if chickStrand == "+":
		chickGt = list(genomeChicken.getGenesNearB(chickChr, (chick1+chick2)/2, (chick2-chick1)/2))
		if len(chickGt) > 0:
			print chickGt
			print
f.close()

sys.exit(0)


f = utils.myTools.myOpenFile(noms_fichiers["chains"], "r")
for l in f:
	c = l.split()
	lizChr = c[2]
	lizStrand = c[4]
	liz1 = int(c[5])
	liz2 = int(c[6])
	chickChr = common(c[7][3:])
	chickStrand = c[9]
	chick1 = int(c[10])
	chick2 = int(c[11])
	if chickStrand == "+":
		chickGt = list(genomeChicken.getGenesNearB(chickChr, (chick1+chick2)/2, (chick2-chick1)/2))
		if len(chickGt) > 0:
			print chickGt
			print
	if lizStrand == "+":
		lizGt = list(genomeLizard.getGenesNearB(lizChr, (liz1+liz2)/2, (liz2-liz1)/2))
		if len(lizGt) > 0:
			print lizGt
			print
	if lizStrand == "+" and chickStrand == "+":
		chickGt = list(genomeChicken.getGenesNearB(chickChr, (chick1+chick2)/2, (chick2-chick1)/2))
		lizGt = list(genomeLizard.getGenesNearB(lizChr, (liz1+liz2)/2, (liz2-liz1)/2))
		if len(chickGt) > 0 and len(lizGt) > 0:
			print chickGt
			print lizGt
			print
		
f.close()

sys.exit(0)

#
# Permet de telecharger, decompresser et lire a la volee un fichier
#
def fileIterator(nom):
	(stdin,stdout,stderr) = os.popen3( ("wget %s/%s -O - | gunzip" % (options["IN.EnsemblURL"],nom)).replace("XXX", str(options["releaseID"])) )
	stdin.close()
	stderr.close()
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
	[("releaseID",int,[42,43,44,45]),
	("OUT.genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("OUT.fullGenesFile",str,"~/work/data/genes/full/genes.%s.list.bz2"), \
	("OUT.orthosFile",str,"~/work/data/orthologs/orthos.%s.%s.list.bz2"), \
	("OUT.paras2File",str,"~/work/data/orthologs/paras.%s.%s.list.bz2"), \
	("OUT.paras1File",str,"~/work/data/paralogs/paras.%s.list.bz2"), \
	("IN.EnsemblURL",str,"ftp://ftp.ensembl.org/pub/release-XXX/mart_XXX/data/mysql/"), \
	("IN.genesFile",str,"ensembl_mart_XXX/%s_gene_ensembl__gene__main.txt.table.gz"), \
	("IN.parasFile",str,"compara_mart_homology_XXX/compara_%s_%s_paralogs__paralogs__main.txt.table.gz"), \
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
for (esp1,esp2) in utils.myTools.myMatrixIterator(nomReel, None, utils.myTools.myMatrixIterator.UpperMatrix):

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
	

