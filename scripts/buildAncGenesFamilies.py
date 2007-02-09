#! /users/ldog/muffato/python -OO

__doc__ = """
Lit l'arbre phylogenetique et construit les fichiers de genes ancestraux
 pour tous les noeuds en utilisant les fichiers d'orthologues de Biomart.
Affine la liste en exhibant les genes specifiques d'une lignee et ceux qui
 n'ont plus qu'un seul gene moderne correspondant
"""

##################
# INITIALISATION #
##################

# Librairies
import sys
import operator
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myCommunities
import utils.myPsOutput


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf"], \
	[("homologyLevels",str,"ortholog_one2many,ortholog_many2many,apparent_ortholog_one2one,ortholog_one2one"), \
	("ancestr",str,""),\
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2"), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("one2oneFile",str,"~/work/data/one2one/one2one.%s.list.bz2"), \
	("orthoFile",str,"~/work/data/orthologs/orthos.%s.%s.list.bz2")], \
	__doc__ \
)

phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])
phylTree.loadAllSpeciesSince(options["ancestr"], options["genesFile"])
homologies = options["homologyLevels"].split(",")
utils.myPsOutput.initColor()



def buildAncFile(anc, lastComb):

	def doLoad(s):

		f = utils.myTools.myOpenFile(s, 'r')
		
		# On lit chaque ligne
		for ligne in f:
			champs = ligne.split()
			gA = champs[0]
			gB = champs[3]
			
			# Si les genes ont ete separes au noeud superieur, on les laisse separes
			if lastComb.dic.get(gA,0) != lastComb.dic.get(gB,0):
				continue
			
			comb.addLink([gA, gB])
			if gA not in aretes:
				aretes[gA] = dict([])
			
			# Les apparent_one2one comptent pour 1/2 dans les poids des aretes
			if champs[6].startswith("apparent"):
				aretes[gA][gB] = 0.5
			else:
				aretes[gA][gB] = 1

		f.close()
		


	# 1. on combine tous les fichiers d'orthologues
	comb = utils.myTools.myCombinator([])
	esp = [phylTree.commonNames[x][0] for x in phylTree.species[anc]]
	aretes = {}
	
	print >> sys.stderr, "Construction des familles d'orthologues de", anc, ":", "-".join(esp), "",
	for (i,j) in utils.myTools.myMatrixIterator(len(esp), len(esp), utils.myTools.myMatrixIterator.StrictUpperMatrix):
		doLoad(options["orthoFile"] % (esp[i],esp[j]))
		utils.myTools.stderr.write('.')
	print >> sys.stderr, " OK"
	
	# 2. On affiche les groupes d'orthologues
	print >> sys.stderr, "Construction des fichiers de", anc, "..."
	
	s = anc.replace('/', '_').replace(' ', '_')
	f = utils.myTools.myOpenFile(options["ancGenesFile"] % s, 'w')
	ff = utils.myTools.myOpenFile(options["one2oneFile"] % s, 'w')
	nbA = 0
	nbO = 0
	res = utils.myTools.myCombinator([])
	for x in comb:
	
		score = phylTree.findFamilyComposition(x)
		poidsBranches = [max([len(score[e]) for e in espGrp]) for espGrp in phylTree.branchesSpecies[anc]]
	
		# A. Famille specifique d'une unique sous-branche
		if len([e for e in poidsBranches if e > 0]) == 1:
			
			# A.1/ Si elle est heritee, on doit la garder
			if x[0] in lastComb.dic:
				print >> sys.stderr, "Famille heritee"
				clusters = [x]
			# A.2/ Sinon, elle correspond a une famille qui sera cree plus tard
			else:
				continue
		
		
		# B. Dans une des sous-branches, tous les genes sont en 1 copie
		elif min(poidsBranches) == 1:
			print >> sys.stderr, "Branche sans duplication"
			clusters = [x]

		# C. On calcule les aretes du graphe
		else:

			# C.1/ Construction des aretes
			nbAretesTot = 0
			for (i1,i2) in utils.myTools.myMatrixIterator(len(x), len(x), utils.myTools.myMatrixIterator.StrictUpperMatrix):
				g1 = x[i1]
				g2 = x[i2]
				if (g1 in aretes.get(g2,[])) or (g2 in aretes.get(g1,[])):
					nbAretesTot += 1
				if (phylTree.dicGenes[g1][0] == phylTree.dicGenes[g2][0]):
					nbAretesTot += 1

			
			# C.2/ Graphe complet -> 1 cluster
			if nbAretesTot == (len(x)*(len(x)-1))/2:
				print >> sys.stderr, "Graphe complet"
				clusters = [x]
				
			elif nbAretesTot >= (len(x)*(len(x)-3))/2:
				print >> sys.stderr, "Graphe quasi-complet"
				clusters = [x]
				
			# C.3/ On est oblige de clusteriser
			else:
				lstCommunitiesOrig = utils.myCommunities.launchCommunitiesBuild1(x, aretes)
				
				lstCommunities = []
				for comm in lstCommunitiesOrig:
				
					print >> sys.stderr, "Test de [alpha=%f relevance=%f parts=%d N/A=%d]" % (comm[0],comm[1],len(comm[2]),len(comm[3]))
					
					# Ne nous interessent que les clusterisations completes
					if len(comm[-1]) != 0:
						print >> sys.stderr, "genes oublies"
						continue
					
					# Ne nous interessent que les clusterisations avec une relevance suffisante
					if comm[1] < 0.35:
						print >> sys.stderr, "relevance insuffisante"
						continue
						
					# Une clusterisation est correcte si tous les clusters sont repartis sur les deux sous-branches
					ok = True
					for c in comm[2]:
						score = phylTree.findFamilyComposition(c)
						poidsBranches = [max([len(score[e]) for e in espGrp]) for espGrp in phylTree.branchesSpecies[anc]]
						if len([e for e in poidsBranches if e > 0]) == 1:
							ok = False
							print >> sys.stderr, "Clusterisation incoherente"
							break
					if ok:
						print >> sys.stderr, "Communaute recevable"
						lstCommunities.append(comm)
			
				# On trie suivant la meilleure relevance
				lstCommunities.sort(key = operator.itemgetter(1), reverse = True)

				# C.3.a/ Aucune clusterisation satisfaisante
				if len(lstCommunities) == 0:
					print >> sys.stderr, "Aucune communaute satisfaisante"
					clusters = [x]
					
				# C.3.b/ Unique solution
				elif len(lstCommunities) == 1:

					# TODO: Verifier le score de relevance !!
					clusters = lstCommunities[0][2]
					print >> sys.stderr, "Resultat unique [alpha=%f relevance=%f parts=%d]" % (lstCommunities[0][0],lstCommunities[0][1],len(clusters))
				
				# C.3.c/ Le choix s'impose
				else:
					print >> sys.stderr, "Hesitations !"
					clusters = lstCommunities[0][2]
				
			
		nbA += len(clusters)

		# E. Ecriture
		for c in clusters:
			res.addLink(c)
			print >> f, " ".join(c)
			
			score = dict( [(e,[]) for e in phylTree.dicGenomes] )
			for i in c:
				score[phylTree.dicGenes[i][0]].append(i)
			l = [score[e][0] for e in score if len(score[e]) == 1]
			if len(l) >= 1:
				print >> ff, " ".join(l)
				nbO += 1
	
	comb = None
	
	# 3. On rajoute les genes qui n'ont plus qu'une seule copie
	print >> sys.stderr, "Rajout des genes solitaires ...",
	for e in phylTree.species[anc]:
		for g in phylTree.dicGenomes[e].dicGenes:
			if (g not in lastComb.dic) or (g in res.dic):
				continue
			print >> f, g
			nbA += 1
			print >> ff, g
			nbO += 1
			res.addLink([g])

	f.close()
	ff.close()
	
	print >> sys.stderr, "OK (%d/%d)" % (nbA, nbO)

	del res.grp
	for (esp,_) in phylTree.items[anc]:
		buildAncFile(esp, res)
		#if esp not in phylTree.dicGenomes:
		#	buildAncFile(esp, res)

buildAncFile(options["ancestr"], utils.myTools.myCombinator([]))
