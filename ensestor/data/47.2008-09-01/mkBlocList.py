#!/usr/bin/env python2


# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myPhylTree
import utils.walktrap


# Charge le fichier de toutes les diagonales (supposees non chevauchantes)
def loadDiagsFile(nom, ancName):
	
	print >> sys.stderr, "Chargement du fichier de diagonales ...",
	f = utils.myTools.myOpenFile(nom, 'r')
	lst = []
	for l in f:
		# Selection de l'ancetre
		if not l.startswith(ancName):
			continue
		# On enleve les "_random" et on extrait chaque colonne
		ct = l.replace('\n','').replace("_random", "").split('\t')
		# La diagonale
		d = [int(x) for x in ct[2].split(' ')]
		# On joint les especes qui ont vu la diagonale et celles qui n'apportent que le chromosome
		tmp = [y.split("/") for y in "|".join([x for x in ct[4:] if len(x) > 0]).split("|")]
		# Les chromosomes de ces especes
		espChr = frozenset( (phylTree.officialName[e],c) for (e,c) in tmp if ('Un' not in c) and (e in phylTree.officialName) )
		# On la garde en memoire
		lst.append( (d,espChr,ct[2],ct[3]) )

	f.close()
	print >> sys.stderr, "OK (%d diagonales)" % len(lst)
	return lst


#
# Rajoute les genes non presents dans des diagonales, en les considerant comme des diagonales de 1
#
def checkLonelyGenes():

	# Charge les genomes des ancetres outgroup
	genesAnc = {}
	for a in phylTree.listAncestr:
		if (phylTree.dicParents[arguments["ancestr"]][a] == a) and (a != arguments["ancestr"]):
			genesAnc[a] = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[a])
	
	print >> sys.stderr, "Ajout des genes solitaires ...",
	# Les genes seuls vont devenir des diagonales de 1
	genesSeuls = set(xrange(len(lstGenesAnc)))
	for (d,_,_,_) in lstDiagsIni:
		genesSeuls.difference_update(d)
	
	nb = 0
	new = []
	for i in genesSeuls:

		lst = []
		for e in phylTree.listSpecies:
			if e in lstEspOutgroup:
				a = phylTree.dicParents[arguments["ancestr"]][e]
				names = genesAnc[a].getOtherNames(lstGenesAnc[i].names)
			else:
				names = lstGenesAnc[i].names
			tmp = [dicGenes[x] for x in names if x in dicGenes]
			tmp = set( c for (x,c) in tmp if x == e )

			# Si les orthologues sont sur un unique chromosome
			if len(tmp) == 1:
				c = str(tmp.pop()).replace("_random","")
				if 'Un' not in c:
					lst.append( (e,c) )
		
		# Petit test: ne peuvent etre utilises que les genes avec au moins 2 especes
		# Pour etre exact, il faudrait avec 2 groupes parmi (fils1 + ... + filsN + outgroup)
		if len(lst) >= 2:
			new.append( ([i],frozenset(lst),str(i),"0") )
			nb += 1
	
	print >> sys.stderr, "%d genes OK" % len(new)
	return new


# Arguments
############
arguments = utils.myTools.checkArgs( [("phylTree.conf",file), ("diagsList",file)], [("ancestr",str,""), ("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), ("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], "Cree la liste des blocs de syntenie (y compris les genes singletons)" )


#  Chargement et initialisation
################################
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])
genesAnc = utils.myGenomes.Genome(arguments["ancGenesFile"] % phylTree.fileName[arguments.ancestr])
lstGenesAnc = genesAnc.lstGenes[None]

lstNoeudsFils = [a for (a,_) in phylTree.items[arguments["ancestr"]]]
lstEspParNoeudsFils = [phylTree.species[f] for f in lstNoeudsFils]
lstEspOutgroup = phylTree.outgroupSpecies[arguments["ancestr"]]

# Chargement des genomes si necessaire
#######################################
dicGenes = {}
for e in phylTree.listSpecies:
	genome = utils.myGenomes.Genome(arguments["genesFile"] % phylTree.fileName[e])
	for (g,(c,_)) in genome.dicGenes.iteritems():
		dicGenes[g] = (e,c)

# Les diagonales et les scores possibles
#########################################
lstDiagsIni = loadDiagsFile(arguments["diagsList"], arguments["ancestr"])

# On doit rajouter les genes non presents dans des diagonales
lstDiagsIni.extend(checkLonelyGenes())

for (i,(_,_,d,s)) in enumerate(lstDiagsIni):
	print "%d\t%s\t%s" % (i,d,s)


