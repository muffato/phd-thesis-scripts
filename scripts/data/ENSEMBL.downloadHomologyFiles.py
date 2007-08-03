#! /users/ldog/muffato/python -OO

__doc__ = """
	Telecharge depuis le site d'Ensembl les fichiers d'homologie
"""


##################
# INITIALISATION #
##################

# Librairies
import os
import sys
import utils.myTools
import utils.myPhylTree


#############
# FONCTIONS #
#############

#
# Permet de telecharger, decompresser et lire a la volee un fichier
#
def fileIterator(nom):
	stdout = utils.myTools.myOpenFile( ("%s/%s" % (options["IN.EnsemblURL"],nom)).replace("XXX", str(options["releaseID"])) , "r")
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
	[("releaseID",int,[42,43,44,45]), ("OUT.directory",str,""), \
	("IN.EnsemblURL",str,"ftp://ftp.ensembl.org/pub/release-XXX/mart_XXX/data/mysql/"), \
	("IN.parasFile",str,"compara_mart_homology_XXX/compara_%s_%s_paralogs__paralogs__main.txt.table.gz"), \
	("IN.orthosFile",str,"compara_mart_homology_XXX/compara_%s_%s_orthologs__orthologs__main.txt.table.gz"), \
	("OUT.orthosFile",str,"orthologs/orthos.%s.%s.list.bz2"), \
	("OUT.paras2File",str,"orthologs/paras.%s.%s.list.bz2"), \
	("OUT.paras1File",str,"paralogs/paras.%s.list.bz2"), \
	], \
	__doc__ \
)



# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])

OUTorthosFile = os.path.join(options["OUT.directory"], options["OUT.orthosFile"])
OUTparas2File = os.path.join(options["OUT.directory"], options["OUT.paras2File"])
OUTparas1File = os.path.join(options["OUT.directory"], options["OUT.paras1File"])
for dir in [OUTorthosFile, OUTparas1File, OUTparas2File]:
	try:
		os.makedirs(os.path.dirname(dir))
	except OSError:
		pass


# Les noms utilises dans les fichiers "Homo Sapiens" -> "hsap"
nomReel = []
for esp in sorted(phylTree.listSpecies):
	tmp = esp.lower().split()
	tmp = tmp[0][0] + tmp[1][:3]
	nomReel.append( (tmp,esp) )

# On genere les fichiers d'homologues
for ((esp1,esp1B),(esp2,esp2B)) in utils.myTools.myIterator.tupleOnUpperList(nomReel):

	theoryAnc = phylTree.dicParents[esp1B][esp2B]
		
	if esp1 == esp2:
		print >> sys.stderr, "Telechargement de la liste des genes paralogues de %s ..." % esp1B,
		nb,_ = proceedFile(options["IN.parasFile"] % (esp1,esp1), OUTparas1File % phylTree.fileName[esp1B], None)
		print >> sys.stderr, "%d genes OK" % nb

	else:
		print >> sys.stderr, "Telechargement de la liste des genes orthologues entre %s et %s ..." % (esp1B,esp2B),
		nb1,nb2 = proceedFile(options["IN.orthosFile"] % (esp1,esp2), OUTorthosFile % (phylTree.fileName[esp1B],phylTree.fileName[esp2B]), OUTorthosFile % (phylTree.fileName[esp2B],phylTree.fileName[esp1B]))
		print >> sys.stderr, "%d/%d genes OK" % (nb2, nb1)

		print >> sys.stderr, "Telechargement de la liste des genes paralogues entre %s et %s ..." % (esp1B,esp2B),
		nb,_ = proceedFile(options["IN.parasFile"] % (esp1,esp2), OUTparas2File % (phylTree.fileName[esp1B],phylTree.fileName[esp2B]), OUTparas2File % (phylTree.fileName[esp2B],phylTree.fileName[esp1B]))
		print >> sys.stderr, "%d genes OK" % nb
	

