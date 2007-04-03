#! /users/ldog/muffato/python -OO

__doc__ = """
A partir de toutes les diagonales extraites entre les especes,
  reconstruit les chromosomes (ou scaffold) de chaque ancetre.
"""


##################
# INITIALISATION #
##################

# Librairies
import sys
import math
import random
import operator
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myDiags
import utils.myCommunities
#import utils.myCommunities2


genesAnc = utils.myGenomes.AncestralGenome(sys.argv[1])
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr]

random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)
random.shuffle(lstGenesAnc)

for i in xrange(1000):
	print "A", " ".join(lstGenesAnc[i].names)
sys.exit(0)

w = utils.myCommunities2.WalktrapWrapper()

#f = open('/workspace/muffato/proteines/reducedGraphs/graph.0', 'r')
f = open('/users/ldog/muffato/work/graph.3109', 'r')

w.updateFromFile(f)
w.doWalktrap()

#print w.res
#sys.exit(0)

for (comp,best,d) in w.res:

	print >> sys.stderr, "comp connexe"
	for (scale,rel) in best:
		print >> sys.stderr, "scale", scale
		print >> sys.stderr, "relevance", rel
		print >> sys.stderr, d.cut(scale)

sys.exit(0)

comb = utils.myTools.myCombinator([])
random.seed(sys.argv[1])

ll = []
for i in xrange(1000):
	
	l = [random.randint(1,10000000) for i in xrange(random.randint(2,1000))]

	ll.extend(l)
	#comb.addLink(l)

#print comb.getNbGrp()
#print len(utils.myMaths.flatten(ll))
#print len(ll)
print utils.myMaths.myStats(ll)

sys.exit(0)


for i in xrange(100000):
	dic[random.random()] = random.random()
	dic[random.random()] = str(random.random())
	dic[random.random()] = 1000000*random.random()
	dic[str(random.random())] = random.random()
	dic[str(random.random())] = str(random.random())
	dic[str(random.random())] = 1000000*random.random()
	dic[1000000*random.random()] = random.random()
	dic[1000000*random.random()] = str(random.random())
	dic[1000000*random.random()] = 1000000*random.random()

for (k,v) in dic.iteritems():
	print k, v

sys.exit(0)


def countAltern(lst):
	
	#lst = [x[eD] for (_,_,x) in lstDCS]

	def countChr(c):
		nb = 0
		for x in lst:
			if c not in x:
				break
			nb += 1
		#print >> sys.stderr, "count", c, nb, lst
		return nb
	
	count = {}
	last = {}
	while len(lst) > 0:
		#print >> sys.stderr, "iteration", lst
		curr = lst.pop(0)
		#print >> sys.stderr, "count", count
		#print >> sys.stderr, "last", last
		for x in curr:
			for y in last:
				if y == x:
					continue
				s = (countChr(x)+1) * last[y] + count.get((x,y), 0)
				count[(x,y)] = count[(y,x)] = s
			for y in curr:
				if y >= x:
					continue
				s = 1 + count.get((x,y), 0)
				count[(x,y)] = count[(y,x)] = s
				
				
		for y in last:
			if y not in curr:
				last[y] = 0
		for x in curr:
			last[x] = last.get(x,0) + 1
		#print >> sys.stderr

	return count

#print countAltern([[5],[13],[5],[5],[5,13],[5,13],[5],[13]])
print countAltern([[5],[6],[5],[5],[5,13],[7,13],[13],[5],[1],[5],[5],[13],[13],[1,5]])


sys.exit(0)



# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf" ], [], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])


print phylTree.parent['Homo sapiens']
print phylTree.parent['Human']
print phylTree.species['Boreoeutheria']
print phylTree.species['Glires']
print phylTree.branchesSpecies['Vertebreta']

print phylTree.fileName




#############
# FONCTIONS #
#############

# Charge le fichier de toutes les diagonales (supposees non chevauchantes)
def loadDiagsFile(nom):
	
	print >> sys.stderr, "Chargement du fichier de diagonales ...",
	f = utils.myTools.myOpenFile(nom, 'r')
	res = []
	for l in f:

		ct = l[:-1].split('\t')
		anc = ct[0]
		#l = int(ct[1])
		d = [int(x) for x in ct[2].split(' ')]
		if len(d) <= 2:
			continue
		#esp = set()
		#if len(ct[3]) > 0:
		#	esp.update( set([tuple(x.split('/')) for x in ct[3].split('|')]) )
		#if len(ct) == 5 and len(ct[4]) > 0:
		#	#print >> sys.stderr, ct[4], ct[4].split('|')
		#	esp.update( set([tuple(x.split('/')) for x in ct[4].split('|')]) )
		#print >> sys.stderr, esp
		#esp = set([(phylTree.officialName[e],c) for (e,c) in esp])
		#diagEntry[anc].append( (l, d, esp) )
		res.append( d )

	f.close()
	print >> sys.stderr, "OK"
	return res
	

########
# MAIN #
########

# Arguments
(noms_fichiers, options) = utils.myTools.checkArgs( \
	["phylTree.conf", "diagsList"], \
	[("ancestr",str,""), ("ancGenesFile",str,"~/work/data/ancGenes/ancGenes.%s.list.bz2")], \
	__doc__ \
)

# L'arbre phylogenetique
phylTree = utils.myBioObjects.PhylogeneticTree(noms_fichiers["phylTree.conf"])

lstDiags = loadDiagsFile(noms_fichiers["diagsList"])

combin = utils.myTools.myCombinator([])
for d in lstDiags:
	combin.addLink(d)

	
gr = utils.myDiags.DiagGraph(lstDiags)
print >> sys.stderr, "[%d]" % len(lstDiags),


ind = 0
for g in combin:
	ind += 1
	sys.stdout = open('/users/ldog/muffato/work/temp/tmp/graphA.%d' % ind, 'w')
	gr.printGraph(g)
	sys.stdout.close()

print >> sys.stderr, "!",
last = None
while last != len(gr.sommets):
	print >> sys.stderr, "{%d}" % len(gr.sommets),
	last = len(gr.sommets)
	gr.reduceGraph()
print >> sys.stderr

ind = 0

for g in combin:
	ind += 1
	sys.stdout = open('/users/ldog/muffato/work/temp/tmp/graphB.%d' % ind, 'w')
	gr.printGraph(g)
	print >> utils.myTools.stdout, len(set(g).intersection(gr.sommets))
	sys.stdout.close()



sys.exit(0)
# Les genes ancestraux
genesAnc = utils.myGenomes.AncestralGenome(options["ancGenesFile"] % options["ancestr"].replace('/', '_').replace(' ', '_'))
lstGenesAnc = genesAnc.lstGenes[utils.myGenomes.Genome.defaultChr]

# On separe les especes en trois
(fils1, fils2) = phylTree.branchesSpecies[options["ancestr"]]
outgroup = phylTree.outgroupSpecies[options["ancestr"]]

print >> sys.stderr, fils1, fils2, outgroup

# Les genomes modernes
#geneBank = utils.myGenomes.GeneBank(noms_fichiers["genesList.conf"], phylTree.species[options["ancestr"]])

# Les diagonales
diagEntry = {}
for anc in phylTree.items:
	diagEntry[anc] = []
loadDiagsFile(noms_fichiers["diagsList"], diagEntry)
lstDiags = diagEntry[options["ancestr"]]
del diagEntry


print >> sys.stderr, len(lstDiags)
# Les genes seuls vont devenir des diagonales de 1
#genesSeuls = set(xrange(len(lstGenesAnc)))
#for (_,d,_) in lstDiags:
#	genesSeuls.difference_update(d)
#for i in genesSeuls:
#	esp = [geneBank.dicGenes[s] for s in lstGenesAnc[i].names]
#	if len(esp) >= 2:
#		lstDiags.append( (1,[i],set([(e,c) for (e,c,_) in esp])) )
#	#print "ajout de ", lstDiags[-1]
#
#print >> sys.stderr, len(lstDiags)
#sys.exit(0)


# Il faut calculer l'apport de chaque espece
dicPoidsEspeces = {}
def calcPoidsFils(node, calc):
	if node in phylTree.listSpecies:
		dicPoidsEspeces[node] = calc
	else:
		poids = calc / float(len(phylTree.items[node]))
		for f in phylTree.branches[node]:
			calcPoidsFils(f, poids)

def calcPoids(node):
	# Les fils a egalite avec un poids de 1
	for f in phylTree.branches[node]:
		calcPoidsFils(f, 1.)
	outgroup = []
	anc = node
	while anc in phylTree.parent:
		par = phylTree.parent[anc]
		outgroup.extend([(e,2*phylTree.ages[par]-phylTree.ages[node]) for (e,_) in phylTree.items[par] if e != anc])
		anc = par
	s = sum([1./math.log(a) for (_,a) in outgroup])
	for (e,a) in outgroup:
		calcPoidsFils(e, 1. / (math.log(a)*s))

calcPoids(options["ancestr"])

print >> sys.stderr, dicPoidsEspeces

def calcScore(i1, i2):

	global lstDiags, dicPoidsEspeces, fils1, fils2, outgroup
	
	(_,_,e1) = lstDiags[i1]
	(_,_,e2) = lstDiags[i2]
	communEsp = set([e for (e,_) in e1.intersection(e2)])

	propF1 = sum([dicPoidsEspeces[e] for e in communEsp.intersection(fils1)])
	propF2 = sum([dicPoidsEspeces[e] for e in communEsp.intersection(fils2)])
	propOut = sum([dicPoidsEspeces[e] for e in communEsp.intersection(outgroup)])
	
	return propF1*propF2 + propF1*propOut + propF2*propOut
	

lstLstComm = utils.myCommunities.launchCommunitiesBuildB(len(lstDiags), calcScore)
clusters = []

# Chaque composante connexe
for lst in lstLstComm:
	# relevance >= 0.3 && noeudsOublies = 0
	interessant = [comm for comm in lst if (len(comm[3]) == 0) and (comm[1] >= 0.3)]
	interessant.sort(key = operator.itemgetter(1), reverse = True)
	if len(interessant) == 0:
		clusters.append(utils.myMaths.flatten(lst[0][2])+lst[0][-1])
	else:
		clusters.extend(interessant[0][2])


print >> sys.stderr, "Impression des chromosomes ancestraux ...",
lstChr = []
for c in clusters:
	lst = set()
	for i in c:
		lst.update(lstDiags[i][1])
	lstChr.append(lst)

for (i1,i2) in utils.myTools.myMatrixIterator(len(lstChr), len(lstChr), utils.myTools.myMatrixIterator.StrictUpperMatrix):
	inter = lstChr[i1].intersection(lstChr[i2])
	lstChr[i1].difference_update(inter)
	lstChr[i2].difference_update(inter)

chrIndex = 0
for c in lstChr:
	chrIndex += 1
	for i in c:
		print chrIndex, " ".join(lstGenesAnc[i].names)


print >> sys.stderr, "OK"
