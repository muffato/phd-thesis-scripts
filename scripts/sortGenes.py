#! /usr/bin/python2.4


# INITIALISATION #

# Librairies
import sys
import os
import utils.myGenomes
import utils.myTools
import utils.myMaths


# FONCTIONS #

def buildDistTree(arbre, distEsp):
		
	def calcDist(anc):

		s = 0.
		n = 0.
		for (e,p) in arbre.items[anc]:
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


	distAnc = dict( [(a,0) for a in arbre.items] )
	
	# On met les 1 dans l'arbre
	lst1 = set([e for e in distEsp if distEsp[e] == 1])
	while len(lst1) >= 2:
		lst1b = set([arbre.parent[e] for e in lst1 if e in arbre.parent])
		for e in lst1b:
			distAnc[e] = 1
		lst1 = lst1b
		
	# On calcule par une moyenne les autres distances
	calcDist(arbre.root)
	return distAnc
	


def distInterGenes(tg1, tg2):
	
	dic = {}
	for (e1,c1,i1) in tg1:
		for (e2,c2,i2) in tg2:
			if e1 == e2 and c1 == c2:
				x = abs(i1-i2)
				if e1 in dic:
					dic[e1] = min(dic[e1], x)
				else:
					dic[e1] = x
	
	return buildDistTree(phylTree, dic)[options["nomAncetre"]]

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genesList.conf", "genomeAncestral", "phylTree.conf"], \
	[("nomAncetre", str, ""), ("seuilMaxDistInterGenes", int, 0), ("nbConcorde", int, -1), ("nbDecimales", int, 2), ("penalite", int, 100000000)], \
	"Trie les gens dans l'ordre indique par l'arbre phylogenetique" \
)

phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], phylTree.getSpecies(phylTree.root))
genesAnc = utils.myGenomes.AncestralGenome(noms_fichiers["genomeAncestral"], True)

if options["nomAncetre"] not in phylTree.items and options["nomAncetre"] not in geneBank.dicEspeces:
	print >> sys.stderr, "Can't retrieve the order of -%s- " % options["nomAncetre"]
	sys.exit(1)


nom = "mat"+str(os.getpid())
nbConcorde = max(1, options["nbConcorde"])

# 2. On cree les blocs ancestraux tries et on extrait les diagonales
for c in genesAnc.lstChr:

	tab = [[geneBank.dicGenes[s] for s in g.names if s in geneBank.dicGenes] for g in genesAnc.lstGenes[c]]
	n = len(tab)
	
	print >> sys.stderr
	print >> sys.stderr, "- Chromosome %s (%d genes) -" % (c, n)
	(blocsAnc, dico) = genesAnc.splitChr(geneBank, c)
	
	print >> sys.stderr, "Ecriture de la matrice ... ",
	f = open(nom, 'w')
	
	print >> f, "NAME: CHRANC"
	print >> f, "TYPE: TSP"
	print >> f, "DIMENSION: %d" % (n+1)
	print >> f, "EDGE_WEIGHT_TYPE: EXPLICIT"
	print >> f, "EDGE_WEIGHT_FORMAT: UPPER_ROW"
	print >> f, "EDGE_WEIGHT_SECTION"
	print >> f, "0 " * n
	
	for i in xrange(n-1):
		for j in xrange(i+1,n):
			y = distInterGenes(tab[i], tab[j])
			if y == 0:
				print >> f, options["penalite"],
			elif y == 1:
				print >> f, 1,
			else:
				print >> f, int(pow(10, options["nbDecimales"])*y),
		print >> f
	print >> f, "EOF"
	f.close()
	print >> sys.stderr, "OK"
	print >> sys.stderr, "Lancement de concorde ",
	lstTot = []
	for i in range(nbConcorde):
		os.system('/users/ldog/muffato/work/scripts/concorde -m -x ' + nom + ' >&2')
		lstTot.append(utils.myBioObjects.ConcordeFile(nom + ".sol"))
		os.system('rm -f 0%s* %s.*' % (nom,nom) )
		sys.stderr.write(".")

	print >> sys.stderr

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


os.system('rm -f *%s*' % nom )
