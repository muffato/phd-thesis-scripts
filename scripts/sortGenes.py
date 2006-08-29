#! /usr/bin/python2.4


# INITIALISATION #

# Librairies
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myOrthos
import myTools
import myMaths


# FONCTIONS #

def calcDist(dic, lst1, anc):

	if anc in dic:
		return dic[anc]
	if anc not in phylTree.items:
		return 0
	
	# Est-ce que les genes sont cote a cote dans deux especes de deux
	#  sous-arbres differents
	desc = phylTree.getDescendants(anc)
	cote = [len(set(x) & lst1) != 0 for x in desc]
	if cote.count(True) >= 2:
		return 1
	
	s = 0.
	n = 0.
	for (e,p) in phylTree.items[anc]:
		val = calcDist(dic, lst1, e)
		if val > 0  and (val <= options["seuilMaxDistInterGenes"] or options["seuilMaxDistInterGenes"] == 0):
			poids = 1./float(p)
			s += val*poids
			n += poids
	if n == 0:
		return 0.
	return s/n


def distInterGenes(tg1, tg2):
	
	lst1 = set([])
	dic = {}
	for g1 in tg1:
		(e1,c1,i1) = geneBank.dicGenes[g1]
		for g2 in tg2:
			(e2,c2,i2) = geneBank.dicGenes[g2]
			if e1 == e2 and c1 == c2:
				dic[e1] = float(abs(i1-i2))
				if abs(i1-i2) == 1:
					lst1.add(e1)
		
	return calcDist(dic, lst1, options["nomAncetre"])

# Arguments
(noms_fichiers, options) = myTools.checkArgs( \
	["genesList.conf", "genomeAncestral", "phylTree.conf"], \
	[("nomAncetre", str, ""), ("seuilMaxDistInterGenes", int, 0), ("nbConcorde", int, 1), ("nbDecimales", int, 2), ("penalite", int, 1000000)], \
	"Trie les gens dans l'ordre indique par l'arbre phylogenetique" \
)

geneBank = myOrthos.MyGeneBank(noms_fichiers[0])
genesAnc = myOrthos.AncestralGenome(noms_fichiers[1], True)
phylTree = myOrthos.PhylogeneticTree(noms_fichiers[2])

if options["nomAncetre"] not in phylTree.items and options["nomAncetre"] not in geneBank.dicEspeces:
	print >> sys.stderr, "Can't retrieve the order of -%s- " % options["nomAncetre"]
	sys.exit(1)


nom = "mat"+str(os.getpid())

# 2. On cree les blocs ancestraux tries et on extrait les diagonales
for c in genesAnc.lstGenes:

	tab = [[g for g in t if g in geneBank.dicGenes] for t in genesAnc.lstGenes[c]]
	n = len(tab)
	
	print >> sys.stderr
	print >> sys.stderr, "- Chromosome %s (%d genes) -" % (c, n)
	(blocsAnc, dico) = genesAnc.splitChr(geneBank, c)
	
	print >> sys.stderr, "Ecriture de la matrice ... ",
	f = open(nom, 'w')
	
	print >> f, "NAME: CHRANC_%s" % c
	print >> f, "TYPE: TSP"
	print >> f, "DIMENSION: %d" % (n+1)
	print >> f, "EDGE_WEIGHT_TYPE: EXPLICIT"
	print >> f, "EDGE_WEIGHT_FORMAT: UPPER_ROW"
	print >> f, "EDGE_WEIGHT_SECTION"
	
	for i in range(n):
		print >> f, 0,
	print >> f
	
	for i in range(n-1):
		for j in range(i+1,n):
			y = distInterGenes(tab[i], tab[j])
			if y == 0:
				print >> f, options["penalite"],
			elif y == 1:
				print >> f, 1,
			else:
				print >> f, int(pow(10, options["nbDecimales"])*y),
		print >> f
	print >> f, "EOF"
	print >> sys.stderr, "OK"
	print >> sys.stderr, "Lancement de concorde ",
	f.close()
	
	lstTot = []
	for i in range(options["nbConcorde"]):
		#os.system('/users/ldog/muffato/work/scripts/concorde -x ' + nom + ' > /dev/null')
		os.system('/users/ldog/muffato/work/scripts/concorde -x ' + nom + ' >&2')
		lstTot.append(myOrthos.ConcordeFile(nom + ".sol"))
		os.remove(nom + ".sol")
		if os.access(nom + ".res", os.W_OK):
			os.remove(nom + ".res")
		if os.access("O" + nom + ".res", os.W_OK):
			os.remove("O" + nom + ".res")
		sys.stderr.write(".")

	print >> sys.stderr

	# On remet chaque liste dans le meme sens que la premiere
	for i in range(1, options["nbConcorde"]):
		if not lstTot[i].isMemeSens(lstTot[0]):
			lstTot[i].reverse()

	for i in range(n):
		q = set([s.res[i] for s in lstTot])
		print "%c %d " % (c, len(q)),
		for x in tab[lstTot[0].res[i]-1]:
			print x,
		print

os.remove(nom)


