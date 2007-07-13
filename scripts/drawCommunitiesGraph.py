#! /users/ldog/muffato/python -OO

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""


##################
# INITIALISATION #
##################

# Librairies
import os
import sys
import math
import operator
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myPhylTree
import utils.walktrap
import utils.myPsOutput

#############
# FONCTIONS #
#############

# Charge le fichier de toutes les diagonales (supposees non chevauchantes)
def loadDiagsFile(nom, ancName):
	
	print >> sys.stderr, "Chargement du fichier de diagonales ...",
	f = utils.myTools.myOpenFile(nom, 'r')
	lst = []
	for l in f:
		
		if not l.startswith(ancName):
			continue
		ct = l[:-1].replace("_random", "").split('\t')
		d = [int(x) for x in ct[2].split(' ')]
		esp = set()
		if len(ct[3]) > 0:
			esp.update( set([tuple(x.split('/')) for x in ct[3].split('|')]) )
		if len(ct) == 5 and len(ct[4]) > 0:
			esp.update( set([tuple(x.split('/')) for x in ct[4].split('|')]) )
		esp = set([(phylTree.officialName[e],c) for (e,c) in esp if (e in phylTree.officialName) and ('Un' not in c)])
		lst.append( (d, esp) )

	f.close()
	print >> sys.stderr, "OK (%d diagonales)" % len(lst)
	return lst


#
# Rajoute les genes non presents dans des diagonales, en les considerant comme des diagonales de 1
#
def checkLonelyGenes():

	# Charge les genomes des especes modernes si ca n'a pas encore ete fait
	if len(phylTree.dicGenomes) == 0:
		phylTree.loadAllSpeciesSince(None, options["genesFile"])
	
	# Charge les genomes des ancetres outgroup
	genesAnc = {}
	for a in phylTree.listAncestr:
		if phylTree.dicParents[options["ancestr"]][a] == a:
			genesAnc[a] = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % phylTree.fileName[a])
	
	print >> sys.stderr, "Ajout des genes solitaires ...",
	# Les genes seuls vont devenir des diagonales de 1
	genesSeuls = set(xrange(len(lstGenesAnc)))
	for (d,_) in lstDiags:
		genesSeuls.difference_update(d)
	nb = 0
	new = {}
	for i in genesSeuls:

		lst = []
		# Les especes filles
		esp = [phylTree.dicGenes[s] for s in lstGenesAnc[i].names if s in phylTree.dicGenes]

		for e in phylTree.listSpecies:
			if e in outgroup:
				# Le gene chez l'ancetre commun
				a = phylTree.dicParents[options["ancestr"]][e]
				g = genesAnc[a].getPosition(lstGenesAnc[i])
				# Cas d'un gene specifique de la lignee
				if len(g) == 0:
					continue
				# Le gene dans l'autre espece
				tmp = [c for (c,_) in phylTree.dicGenomes[e].getPosition(genesAnc[a].lstGenes[utils.myGenomes.Genome.defaultChr][g[0][1]])]
			else:
				# Le gene dans l'autre espece
				tmp = [c for (x,c,_) in esp if x == e]
			
			# Au final
			if len(tmp) == 1:
				c = str(tmp[0]).replace("_random","")
				if 'Un' not in c:
					lst.append( (e,c) )
		
		if len(lst) >= 2:
			lst = tuple(sorted(set(lst)))
			new[lst] = new.get(lst,[]) + [i]
			nb += 1
	
	lstDiags.extend([(d,set(e)) for (e,d) in new.iteritems()])
	print >> sys.stderr, "%d (%d) OK" % (nb,len(new))


# Charge les ancetres deja construits et rajoute les infos dans les diagonales
def checkAlreadyBuildAnc():
	genAlready = {}
	for f in filsAnc:
		s = options["alreadyBuiltAnc"] % phylTree.fileName[f]
		print >> sys.stderr, "Checking %s ..." % f,
		if os.access(s, os.R_OK):
			genAlready[f] = utils.myGenomes.AncestralGenome(s, chromPresents=True)
			s = len(genAlready[f].lstChr)
			if options["weightNbChr+"]:
				espCertitude[f] = 1. - 1. / float(s)
			if options["weightNbChr-"]:
				espIncertitude[f] = s * alphaIncertitude
		else:
			print >> sys.stderr, "Not found"

	print >> sys.stderr, "Mise a jour des chromosomes des diagonales ...",
	for (d,e) in lstDiags:
		g = utils.myMaths.flatten([lstGenesAnc[i].names for i in d])
		for f in genAlready:
			e.update([(f,genAlready[f].dicGenes[s][0]) for s in g if s in genAlready[f].dicGenes])
	del genAlready
	print >> sys.stderr, "OK"




def calcPoids(node):

	# Les fils a egalite avec un poids de 1
	for f in phylTree.branches[node]:
		calcPoidsFils(f, 1.)
	if not options["useOutgroups"]:
		return
	outgroup = []
	anc = node
	while anc in phylTree.parent:
		par = phylTree.parent[anc]
		outgroup.extend([(e,1./float(2*phylTree.ages[par]-phylTree.ages[node])) for e in phylTree.branches[par] if e != anc])
		anc = par
	s = sum([a for (_,a) in outgroup])
	for (e,a) in outgroup:
		calcPoidsFils(e, a/s)


########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "diagsList"], \
	[("ancestr",str,""), ("alreadyBuiltAnc",str,""), \
	("useOutgroups",bool,True), ("useLonelyGenes",bool,False), ("weightNbChr+",bool,False), ("weightNbChr-",bool,False), ("newScoring",bool,False), \
	("walktrapLength",int,5), ("qualityFunction",int,[2,1,3]), ("cutoff",float,1), \
	("genesFile",str,"~/work/data/genes/genes.%s.list.bz2"), \
	("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

#  Chargement et initialisation
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])
genesAnc = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % phylTree.fileName[options.ancestr])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr]
lstDiags = loadDiagsFile(noms_fichiers["diagsList"], options["ancestr"])

lstToutesEspeces = set(phylTree.listSpecies)
lstEspOutgroup = set(phylTree.outgroupSpecies[options["ancestr"]])
lstEspFilles = set(phylTree.species[options["ancestr"]])
lstNoeudsFils = phylTree.branches[options["ancestr"]]
lstEspParNoeudsFils = [set(s) for s in phylTree.branchesSpecies[options["ancestr"]]]

filsAnc = phylTree.branches[options["ancestr"]]
filsEsp = [set(s) for s in phylTree.branchesSpecies[options["ancestr"]]]
outgroup = set(phylTree.outgroupSpecies[options["ancestr"]])


# L'apport en terme de couverture de l'arbre phylogenetique
def calcPoidsFils(node, calc):
	dicPoidsEspeces[node] = calc
	if node in phylTree.tmpItems:
		poids = calc / float(len(phylTree.tmpItems[node]))
		for (f,_) in phylTree.tmpItems[node]:
			calcPoidsFils(f, poids)
dicPoidsEspeces = dict.fromkeys(phylTree.listSpecies, 0.)
phylTree.initCalcDist(options["ancestr"], options["useOutgroups"])
calcPoidsFils(options["ancestr"], float(len(phylTree.tmpItems[0])))



# Le taux de precision de chaque espece
espCertitude = dict.fromkeys(phylTree.commonNames, 1.)
espIncertitude = dict.fromkeys(phylTree.commonNames, 0.)
if options["weightNbChr+"] or options["weightNbChr-"]:
	phylTree.loadAllSpeciesSince(None, options["genesFile"])
	tmp = []
	for e in phylTree.listSpecies:
		s = len(phylTree.dicGenomes[e].lstChr)
		if s > 0:
			tmp.append(s)
			if options["weightNbChr+"]:
				espCertitude[e] = 1. - 1. / float(s)
	if options["weightNbChr-"]:
		alphaIncertitude = 1. / float(min(tmp) * max(tmp))
		for e in phylTree.listSpecies:
			espIncertitude[e] = len(phylTree.dicGenomes[e].lstChr) * alphaIncertitude


# On doit rajouter les genes non presents dans des diagonales
if options["useLonelyGenes"]:
	checkLonelyGenes()

# On doit noter les chromosomes des diagonales sur ces ancetres deja construits
if options["alreadyBuiltAnc"] != "":
	checkAlreadyBuildAnc()

phylTree.dicGenomes.clear()

print >> sys.stderr, "%-25s\t%s\t%s\t%s" % ("Ancetre", "Poids", "Cert", "Incert")
for e in dicPoidsEspeces:
	print >> sys.stderr, "%-25s\t%.3f\t%.3f\t%.3f" % (e, dicPoidsEspeces[e], espCertitude[e], espIncertitude[e])


def calcScore(i1, i2):

	(_,e1) = lstDiags[i1]
	(_,e2) = lstDiags[i2]
	comparedEsp = set([e for (e,_) in e1]).intersection([e for (e,_) in e2])
	communEsp = set([e for (e,_) in e1.intersection(e2)])

	propF = []
		
	for (i,f) in enumerate(filsAnc):
		# Chez l'ancetre du dessous - meme chr
		if f in communEsp:
			propF.append( dicPoidsEspeces[f]*espCertitude[f] )
		# Chez l'ancetre du dessous - diff chr
		elif f in comparedEsp:
			propF.append( dicPoidsEspeces[f]*espIncertitude[f] )
		# Sinon, on revient aux genomes modernes
		else:
			f = filsEsp[i]
			s = sum([dicPoidsEspeces[e]*espCertitude[e] for e in f.intersection(communEsp)])
			# On peut rajouter l'incertitude des autres especes
			#   - si au moins une espece a valide la fusion
			#   - si la branche est constituee d'une unique espece
			if (s > 0) or (len(f) == 1):
				s += sum([dicPoidsEspeces[e]*espIncertitude[e] for e in f.difference(communEsp)])
			propF.append(s)
	
	#for i in xrange(len(filsEsp)):
	#	# Chez l'ancetre du dessous - meme chr
	#	if filsAnc[i] in communEsp:
	#		propF.append( dicPoidsEspeces[filsAnc[i]]*espCertitude[filsAnc[i]] )
	#	# Chez l'ancetre du dessous - diff chr
	#	elif filsAnc[i] in comparedEsp:
	#		propF.append( dicPoidsEspeces[filsAnc[i]]*espIncertitude[filsAnc[i]] )
	#	# Sinon, on revient aux genomes modernes
	#	else:
	#		s = sum([dicPoidsEspeces[e]*espCertitude[e] for e in filsEsp[i].intersection(communEsp)])
	#		# On peut rajouter l'incertitude des autres especes
	#		#   - si au moins une espece a valide la fusion
	#		#   - si la branche est constituee d'une unique espece
	#		if (s > 0) or (len(filsEsp[i]) == 1):
	#			s += sum([dicPoidsEspeces[e]*espIncertitude[e] for e in filsEsp[i].difference(communEsp)])
	#		propF.append(s)

		
	propOut = sum([dicPoidsEspeces[e]*espCertitude[e] for e in communEsp.intersection(outgroup)])
	if (propOut > 0):
		propOut += sum([dicPoidsEspeces[e]*espIncertitude[e] for e in communEsp.difference(outgroup)])
	
	s = sum(propF) * propOut
	for (f1,f2) in utils.myTools.myMatrixIterator(propF, None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
		s += f1 * f2

	return s
	
def calcScore2(i1, i2):

	(_,e1) = lstDiags[i1]
	(_,e2) = lstDiags[i2]
	comparedEsp = set([e for (e,_) in e1]).intersection([e for (e,_) in e2])
	communEsp = set([e for (e,_) in e1.intersection(e2)])

	values = {}
	for e in comparedEsp:
		values[e] = espIncertitude[e]
	for e in communEsp:
		values[e] = espCertitude[e]

	propF = [phylTree.calcDist(values, f) for f in filsEsp]
	if len(outgroup) == 0:
		propOut = 0
	else:
		propOut = phylTree.calcDist(values, phylTree.parent[options["ancestr"]])
	
	s = sum(propF) * propOut
	for (f1,f2) in utils.myTools.myMatrixIterator(propF, None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
		s += f1 * f2

	return s


walktrapInstance = utils.walktrap.WalktrapLauncher()
print >> sys.stderr, "Calcul de la matrice ...",
if options["newScoring"]:
	walktrapInstance.updateFromFunc(range(len(lstDiags)), calcScore2)
else:
	walktrapInstance.updateFromFunc(range(len(lstDiags)), calcScore)
walktrapInstance.doWalktrap(randomWalksLength=options["walktrapLength"], qualityFunction=options["qualityFunction"])


clusters = []
# Chaque composante connexe
for (nodes,cuts,_,dend) in walktrapInstance.res:
	print >> sys.stderr, cuts
	# Un score de relevance > 0.2
	interessant = [(alpha,score,dend.cut(alpha)) for (alpha,score) in cuts if score > 0.2]
	# Les noeuds seuls representent < de l'ensemble des noeuds
	interessant = [(alpha,score,clust) for (alpha,score,(clust,lonely)) in interessant if len(lonely) < len(nodes)/2]
	if len(interessant) == 0:
		print >> sys.stderr, "-",
		clusters.append(nodes)
	else:
		# Au choix, on prend la version la moins fusionnee
		print >> sys.stderr, "+%d/%f/%f" % (len(interessant[-1][-1]), interessant[-1][0], interessant[-1][1]),
		clusters.extend(interessant[-1][-1])
print >> sys.stderr



print >> sys.stderr, "Impression des %d chromosomes ancestraux ..." % len(clusters),
dic = {}
for (i,c) in enumerate(clusters):
	for d in c:
		dic[d] = i

print >> sys.stderr, "Affichage du graphe ...",

refSpecies = "Homo sapiens"

utils.myPsOutput.initColor()

print "graph {"

for i in xrange(len(lstDiags)):
	(d,esp) = lstDiags[i]
	l = 0.25 * math.log(len(d))/math.log(2)
	lst = [c for (e,c) in esp if e == refSpecies]
	print '%d  [height=%.2f,width=%.2f,fixedsize=true,label="",style=filled' % (i,l,l),
	if len(lst) > 0:
		(r,g,b) = utils.myPsOutput.colorTableUNIX2RGB[utils.myPsOutput.colorTransl.get(str(lst[0]),"black")]
		#print ',color="#%02x%02x%02x"' % (r,g,b),
	print "]"


for i1 in xrange(len(lstDiags)):
	#if i1 not in dic:
	#	continue
	for i2 in xrange(i1+1, len(lstDiags)):
		#if i2 not in dic:
		#	continue
		#if dic[i1] != dic[i2]:
		#	continue
		if options["newScoring"]:
			s = calcScore2(i1, i2)
		else:
			s = calcScore(i1, i2)

		if s >= options["cutoff"]:
			(_,e1) = lstDiags[i1]
			(_,e2) = lstDiags[i2]
			communEsp = set([e for (e,_) in e1.intersection(e2)])
			
			#s = 3 - (3.-s)/100.
			#if dic[i1] != dic[i2]:
			#	s = 0

			if dic.get(i1,-1) == dic.get(i2,-2):
				style = "dotted"
			else:
				style = "invis"
			
			print '%d -- %d [style=%s,w=%.3f,weight=%.3f' % (i1,i2,style,s,s),

			#if dic[i1] != dic[i2]:
			#	print ',constraint=false',

			if dic.get(i1,-1) == dic.get(i2,-2):
				(r,g,b) = utils.myPsOutput.colorTableUNIX2RGB[utils.myPsOutput.colorTransl[str(dic[i1]+1)]]
			#if refSpecies in communEsp:
			#	lst = [c for (e,c) in e1 if e == refSpecies]
			#	(r,g,b) = utils.myPsOutput.colorTableUNIX2RGB[utils.myPsOutput.colorTransl.get(str(lst[0]),"black")]
				print ',color="#%02x%02x%02x"' % (r,g,b),
			print "]"

print "}"

print >> sys.stderr, "OK"


