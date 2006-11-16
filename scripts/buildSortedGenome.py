#! /users/ldog/muffato/python


# INITIALISATION #

# Librairies
import sys
import os
import utils.myGenomes
import utils.myTools
import utils.myMaths


# FONCTIONS #

def distInterGenes(tg1, tg2):

	def calcDist(anc):

		s = 0.
		n = 0.
		for (e,p) in phylTree.items[anc]:
			if e not in distAnc:
				if e in distEsp:
					val = distEsp[e]
				else:
					val = 0
			else:
				if distAnc[e] != 0:
					val = distAnc[e]
				else:
					val = calcDist(e)
			if val > 0  and (val <= options["seuilMaxDistInterGenes"] or options["seuilMaxDistInterGenes"] == 0):
				poids = 1./float(p)
				s += val*poids
				n += poids
		if n != 0:
			s /= n
		distAnc[anc] = s
		return s

	#print >> sys.stderr, 0, len(tg1), len(tg2)
	distEsp = {}
	for (e1,c1,i1) in tg1:
		for (e2,c2,i2) in tg2:
			if e1 == e2 and c1 == c2:
				x = abs(i1-i2)
				if e1 in distEsp:
					distEsp[e1] = min(distEsp[e1], x)
				else:
					distEsp[e1] = x
	
	#print >> sys.stderr, 1, distEsp
	distAnc = dict( [(a,0) for a in phylTree.items] )
	lst1 = [e for e in distEsp if distEsp[e] == 1]
	
	# On met les 1 dans l'arbre
	for i in xrange(len(lst1)):
		for j in xrange(i):
			e1 = lst1[i]
			e2 = lst1[j]
			anc = phylTree.getFirstParent(e1, e2)
			tmp = e1
			while tmp != anc:
				tmp = phylTree.parent[tmp]
				distAnc[tmp] = 1
			tmp = e2
			while tmp != anc:
				tmp = phylTree.parent[tmp]
				distAnc[tmp] = 1
			if distAnc[options["nomAncetre"]] == 1:
				return 1
	
	#print >> sys.stderr, 2, distAnc
			
	# On calcule par une moyenne les autres distances
	calcDist(options["nomAncetre"])
	#print >> sys.stderr, 3, distAnc
	return distAnc[options["nomAncetre"]]
	

# Initialisation & Chargement des fichiers
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf", "genomeAncestral", "phylTree.conf"], \
	[("nomAncetre",str,""), ("seuilMaxDistInterGenes",int,0), ("nbDecimales",int,2), ("penalite",int,1000000), \
	("nbConcorde",int,-1), ("withCondordeOutput",bool,False), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	"Trie les gens dans l'ordre indique par l'arbre phylogenetique" \
)

phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
if options["nomAncetre"] not in phylTree.items and options["nomAncetre"] not in phylTree.getSpecies(phylTree.root):
	print >> sys.stderr, "Can't retrieve the order of -%s- " % options["nomAncetre"]
	sys.exit(1)
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], phylTree.getSpecies(phylTree.root))
del geneBank.dicEspeces
genesAnc = utils.myGenomes.AncestralGenome(noms_fichiers["genomeAncestral"], True, False)

# On etend la liste des genes ancestraux pour utiliser les outgroup

anc = options["nomAncetre"]
dicOutgroupGenes = {}
while anc in phylTree.parent:
	anc = phylTree.parent[anc]
	tmpGenesAnc = utils.myGenomes.loadGenome(options["ancGenesFile"] % anc)
	del tmpGenesAnc.dicGenes
	for g in tmpGenesAnc:
		ianc = set([])
		newGenes = []
		for s in g.names:
			if s in genesAnc.dicGenes:
				ianc.add(genesAnc.dicGenes[s])
			elif s in geneBank.dicGenes:
				newGenes.append(geneBank.dicGenes[s])
		for i in ianc:
			if i in dicOutgroupGenes:
				dicOutgroupGenes[i].update(newGenes)
			else:
				dicOutgroupGenes[i] = set(newGenes)
	del tmpGenesAnc
	del ianc
	del newGenes

# On reecrit le genome en terme d'especes
genome = {}
for c in genesAnc.lstChr:
	genome[c] = []
	for i in xrange(len(genesAnc.lstGenes[c])):
		g = genesAnc.lstGenes[c][i]
		tmp = set([geneBank.dicGenes[s] for s in g.names if s in geneBank.dicGenes])
		if i in dicOutgroupGenes:
			tmp.extend(dicOutgroupGenes[i])
		genome[c].append(tmp)
del geneBank
del dicOutgroupGenes

nom = "mat"+str(os.getpid())
nbConcorde = max(1, options["nbConcorde"])
mult = pow(10, options["nbDecimales"])

# 2. On cree les blocs ancestraux tries et on extrait les diagonales
for c in genesAnc.lstChr:

	tab = genome[c]
	n = len(tab)
	
	print >> sys.stderr
	print >> sys.stderr, "- Chromosome %s (%d genes) -" % (c, n)
	
	print >> sys.stderr, "Ecriture de la matrice ... ",
	f = open(nom, 'w')
	
	print >> f, "NAME: CHRANC"
	print >> f, "TYPE: TSP"
	print >> f, "DIMENSION: %d" % (n+1)
	print >> f, "EDGE_WEIGHT_TYPE: EXPLICIT"
	print >> f, "EDGE_WEIGHT_FORMAT: UPPER_ROW"
	print >> f, "EDGE_WEIGHT_SECTION"
	print >> f, "0 " * n
	
	for i in xrange(n):
		for j in xrange(i+1,n):
			y = distInterGenes(tab[i], tab[j])
			if y == 0:
				print >> f, int(mult*options["penalite"]),
			elif y == 1:
				print >> f, 1,
			else:
				print >> f, int(mult*y),
		print >> f
	print >> f, "EOF"
	f.close()
	print >> sys.stderr, "OK"
	print >> sys.stderr, "Lancement de concorde ",
	lstTot = []
	for i in range(nbConcorde):
		comm = '~/work/scripts/concorde -m -x ' + nom
		if options["withCondordeOutput"]:
			os.system(comm + ' >&2')
		else:
			os.system(comm + ' > /dev/null')
		lstTot.append(utils.myBioObjects.ConcordeFile(nom + ".sol"))
		os.system('rm -f 0%s* %s.*' % (nom,nom) )
		sys.stderr.write(".")

	# On remet chaque liste dans le meme sens que la premiere
	for i in range(1, nbConcorde):
		if not lstTot[i].isMemeSens(lstTot[0]):
			lstTot[i].reverse()

	for i in xrange(n):
		q = set([s.res[i] for s in lstTot])
		if options["nbConcorde"] < 1:
			print c,
		else:
			print c, len(q),
		print " ".join(genesAnc.lstGenes[c][lstTot[0].res[i]-1].names)
	
	print >> sys.stderr, len(utils.myMaths.unique(lstTot)), "solutions"


os.system('rm -f *%s*' % nom )
