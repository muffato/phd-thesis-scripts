#! /users/ldog/muffato/python -OO


# INITIALISATION #

# Librairies
import os
import sys
import random
import utils.myGenomes
import utils.myTools
import utils.myMaths


# FONCTIONS #

def distInterGenes(tg1, tg2):

	def calcDist(anc):

		d = 0.
		s = 0.
		nb = 0
		# Moyenne sur tous les fils du noeud
		for (e,p) in phylTree.items[anc]:
			# Le fils est un nouveau noeud
			if e in phylTree.items:
				# Distance de 1 ou appel recursif
				if e in lst1Anc:
					val = 1
				else:
					(val,x) = calcDist(e)
					p += x
			# C'est une espece representee par les deux genes ancestraux
			elif e in distEsp:
				val = distEsp[e]
			# Une espece non representee
			else:
				continue
			
			# Si on a une vraie valeur de distance, on continue la moyenne
			if val != 0:
				poids = 1./float(p)
				d += val*poids
				s += poids
				nb += 1

		# Test final
		if nb != 0:
			d /= s
			if (options["seuilMaxDistInterGenes"] == 0) or (d <= options["seuilMaxDistInterGenes"]):
				if nb == 1:
					return (d,1./s)
				else:
					return (d,0)
		return (0,0)

	# Les distances chez chaque espece
	distEsp = {}
	for (e1,c1,i1) in tg1:
		for (e2,c2,i2) in tg2:
			if e1 == e2 and c1 == c2:
				x = abs(i1-i2)
				if e1 in distEsp:
					distEsp[e1] = min(distEsp[e1], x)
				else:
					distEsp[e1] = x
	
	# On fait la liste des especes qui presentent une distance de 1
	lst1Esp = [e for e in distEsp if distEsp[e] == 1]
	
	# On met les 1 dans les noeuds de l'arbre entre les especes
	lst1Anc = set()
	for i in xrange(len(lst1Esp)):
		for j in xrange(i):
			e1 = lst1Esp[i]
			e2 = lst1Esp[j]
			lst1Anc.update(dicNodesBetween[(e1,e2)])
			if options["ancestr"] in lst1Anc:
				return 1
	
			
	# On calcule par une moyenne les autres distances
	return calcDist(options["ancestr"])[0]
	

# Initialisation & Chargement des fichiers
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["genomeAncestral", "phylTree.conf"], \
	[("ancestr",str,""), ("seuilMaxDistInterGenes",float,0), ("nbDecimales",int,2), ("penalite",int,1000000), \
	("nbConcorde",int,-1), ("withConcordeOutput",bool,False), ("withConcordeStats",bool,False),\
	("concordeExec",str,"~/work/scripts/concorde"),\
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	"Trie les gens dans l'ordre indique par l'arbre phylogenetique" \
)


# L'arbre phylogentique
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
if options["ancestr"] not in (phylTree.listAncestr + phylTree.listSpecies):
	print >> sys.stderr, "Can't retrieve the order of -%s- " % options["ancestr"]
	sys.exit(1)

# Les liens entre les paires d'especes
dicNodesBetween = {}
for e1 in phylTree.listSpecies:
	for e2 in phylTree.listSpecies:
		dicNodesBetween[(e1,e2)] = phylTree.getNodesBetween(e1,e2)

# On charge les genomes
phylTree.loadAllSpeciesSince(None, options["genesFile"])
del phylTree.dicGenomes
genesAnc = utils.myGenomes.loadGenome(noms_fichiers["genomeAncestral"])

# On etend la liste des genes ancestraux pour utiliser les outgroup
anc = options["ancestr"]
dicOutgroupGenes = {}
while anc in phylTree.parent:
	anc = phylTree.parent[anc]
	tmpGenesAnc = utils.myGenomes.loadGenome(options["ancGenesFile"] % anc)
	del tmpGenesAnc.dicGenes
	for g in tmpGenesAnc:
		ianc = set()
		newGenes = []
		for s in g.names:
			if s in genesAnc.dicGenes:
				ianc.add(genesAnc.dicGenes[s])
			elif s in phylTree.dicGenes:
				newGenes.append(phylTree.dicGenes[s])
		for i in ianc:
			if i in dicOutgroupGenes:
				dicOutgroupGenes[i].update(newGenes)
			else:
				dicOutgroupGenes[i] = set(newGenes)
	del tmpGenesAnc

# On reecrit le genome en terme d'especes
genome = {}
for c in genesAnc.lstChr:
	genome[c] = []
	for i in xrange(len(genesAnc.lstGenes[c])):
		g = genesAnc.lstGenes[c][i]
		tmp = set([phylTree.dicGenes[s] for s in g.names if s in phylTree.dicGenes])
		if i in dicOutgroupGenes:
			tmp.extend(dicOutgroupGenes[i])
		genome[c].append(tmp)
del phylTree.dicGenes
del dicOutgroupGenes

nom = "mat%08d" % ((os.getpid() ^ os.getppid() ^ random.randint(1,16777215)) & 16777215)
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
			#print >> sys.stderr, tab[i], tab[j], y
			if y == 0:
				print >> f, int(mult*options["penalite"]),
			elif y == 1:
				print >> f, 0,
			else:
				print >> f, int(mult*y),
		print >> f
	print >> f, "EOF"
	f.close()
	print >> sys.stderr, "OK"
	print >> sys.stderr, "Lancement de concorde ",
	lstTot = []
	for i in range(nbConcorde):
		comm = options["concordeExec"] + ' -x ' + nom
		if options["withConcordeOutput"]:
			os.system(comm + ' >&2')
		else:
			os.system(comm + ' > /dev/null')
		if os.access(nom + ".sol", os.R_OK):
			lstTot.append(utils.myBioObjects.ConcordeFile(nom + ".sol"))
		os.system('rm -f 0%s* %s*' % (nom,nom) )
		sys.stderr.write(".")

	# On remet chaque liste dans le meme sens que la premiere
	for i in range(1, len(lstTot)):
		if not lstTot[i].isMemeSens(lstTot[0]):
			lstTot[i].reverse()

	if len(lstTot) == 0:
		lstTot.append(range(1, n+1))
	
	for i in xrange(n):
		print c,
		if (options["nbConcorde"] > 1) and options["withConcordeStats"]:
			print len(set([s.res[i] for s in lstTot])),
		#if len(lstTot) == 0:
		#	print " ".join(genesAnc.lstGenes[c][i].names)
		#else:
		print " ".join(genesAnc.lstGenes[c][lstTot[0].res[i]-1].names)
	
	solUniq = utils.myMaths.unique([l.res for l in lstTot])
	print >> sys.stderr, len(solUniq), "solutions"
	if options["withConcordeStats"]:
		for sol in solUniq:
			print ".%s" % c, " ".join([str(i) for i in sol])


os.system('rm -f *%s*' % nom )
