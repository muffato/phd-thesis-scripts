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
	

