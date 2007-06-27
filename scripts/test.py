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
import time
import random
import operator
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myDiags
import utils.myPsOutput
import utils.myPhylTree
#import utils.walktrap.myCommunitiesProxy
import utils.walktrap





#genomes = []
#for i in sys.argv[1:]:
#	genomes.append(utils.myGenomes.loadGenome(i))
#	print len(genomes[-1].lstChr)

s = utils.walktrap.myCommunitiesProxy.WalktrapProxy()
f = open("/users/ldog/muffato/heimdall/proteines/reduced/graph.4741", "r")

#for i in xrange(100):
#	for j in xrange(100):
#		s.addEdge(i, j, random.random())
#		#pass

s.updateFromFile(f)
f.close()

s.doWalktrap(internal = False)
print s.res
sys.exit(0)

def linear(a, b, ratio):
	return int(a + ratio*(b-a))

def linear3((ra,ga,ba), (rb,gb,bb), ratio):
	return (linear(ra,rb,ratio), linear(ga,gb,ratio), linear(ba,bb,ratio))

def printLine(start, end, nbsteps, y):
	for i in xrange(nbsteps):
		col = linear3(start, end, i/(nbsteps-1.))
		utils.myPsOutput.drawBox(i*21./nbsteps, y, 21/nbsteps, 1, col, col)

chrC = {}
chrC["1"] = (153,102,0)
chrC["2"] = (102,102,0)
chrC["3"] = (153,153,30)
chrC["4"] = (204,0,0)
chrC["5"] = (255,0,0)
chrC["6"] = (255,0,204)
chrC["7"] = (255,204,204)
chrC["8"] = (255,153,0)
chrC["9"] = (255,204,0)
chrC["10"] = (255,255,0)
chrC["11"] = (204,255,0)
chrC["12"] = (0,255,0)
chrC["13"] = (53,128,0)
chrC["14"] = (0,0,204)
chrC["15"] = (102,153,255)
chrC["16"] = (153,204,255)
chrC["17"] = (0,255,255)
chrC["18"] = (204,255,255)
chrC["19"] = (153,0,204)
chrC["20"] = (204,51,255)
chrC["21"] = (204,153,255)
chrC["22"] = (102,102,102)
chrC["23"] = (153,153,153)
chrC["24"] = (204,204,204)


utils.myPsOutput.printPsHeader()


#printLine( (10,10,100), (230,230,255), 3, 1)
#printLine( (10,10,100), (10,10,255), 3, 2)

#printLine( (10,100,10), (230,255,230), 3, 4)
#printLine( (100,10,10), (255,230,230), 3, 3)

#printLine( (10,10,255), (230,230,255), 3, 5)
#printLine( (10,255,10), (230,255,230), 3, 6)
#printLine( (255,10,10), (255,230,230), 3, 7)
#printLine( (10,255,255), (230,255,255), 3, 8)
##printLine( (255,10,255), (255,230,255), 3, 9)
#printLine( (255,255,10), (255,255,230), 3, 10)

#for i in xrange(1,25+1):
#	utils.myPsOutput.drawBox(i*.5, 15, .5, 1, i, i)

#for i in chrC:
#	utils.myPsOutput.drawBox(int(i)*.5, 16, .5, 1, chrC[i], chrC[i])

c5 = []
c5 += ["red4", "coral4", "firebrick", "red2", "coral2", "darkorange", "gold"]
c5 += ["yellow", "khaki1","wheat1","peachpuff","lightsalmon", "hotpink2", "magenta2", "darkorchid2", "purple2", "darkorchid4"]
c5 += ["blue2", "royalblue2", "blue4"]
c5 += ["turquoise4", "darkseagreen4", "chartreuse4","mediumaquamarine",(121,204,61), "chartreuse2", "olivedrab2", "darkolivegreen1"]
c5 += ["darkseagreen1", "paleturquoise2", "lightblue", "skyblue1", "turquoise2", "lavender","thistle2"]
c5 += [(204,204,153),"lightgoldenrod3","ivory2","honeydew3","slategray"]
for i in xrange(len(c5)):
	utils.myPsOutput.drawBox(i*.5, 21, .5, 1, c5[i], c5[i])




utils.myPsOutput.printPsFooter()

sys.exit(0)

def func(x):
	return x*x

for i in xrange(100):

	test = [random.random() for j in xrange(100000)]
	#todo = [x*x for x in test]
	#todo = map(lambda x: x*x, test)
	#todo = map(func, test)
	todo = range(len(test))
	for j in xrange(len(test)):
		todo[j] = test[j]*test[j]


sys.exit(0)



phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1])

for anc in phylTree.listAncestr:
	
	nbO = len(phylTree.outgroupSpecies[anc])
	f = [len(x) for x in phylTree.branchesSpecies[anc]]
	s = nbO * sum(f)
	for (n1,n2) in utils.myTools.myMatrixIterator(f, None, utils.myTools.myMatrixIterator.StrictUpperMatrix):
		s += n1*n2

	print anc, s


sys.exit(0)
nbb = 10

x = "1"
x = x.zfill(int(sys.argv[1]))
xx = [int(x[i:i+nbb],2) for i in xrange(0,len(x)+1-nbb,nbb)]

y = "1"
y = y.zfill(int(sys.argv[2]))
yy = [int(y[i:i+nbb],2) for i in xrange(0,len(y)+1-nbb,nbb)]

dic = {}
for i in xrange(2**nbb):
	dic[i] = {}
	for j in xrange(2**nbb):
		n = i ^ j
		tmp = n - ((n >> 1) & 033333333333) - ((n >> 2) & 011111111111)
		dic[i][j] = ((tmp + (tmp >> 3)) & 030707070707) % 63
tyy = [dic[i] for i in yy]

for i in xrange(0, len(xx)-len(yy)):
	s = sum([tyy[j][xx[i+j]] for j in xrange(len(yy))])
	print i, 2*s-len(y)

sys.exit(0)

n = len(y)
while len(x) >= n:

	s = len([None for j in xrange(-1,-n-1,-1) if x[j] == y[j]])
	print len(x)-n, 2*s-n
	x.pop()

sys.exit(0)

#genome = utils.myGenomes.EnsemblGenome("~/work/data43/genes/genes.Gallus.gallus.list.bz2")
#genome = utils.myGenomes.EnsemblGenome(sys.argv[1])
phylTree = utils.myPhylTree.PhylogeneticTree(sys.argv[1])
phylTree.initCalcDist("Boreoeutheria", True)
#phylTree.buildPhylLinks()
#print phylTree.dicLinks["Homo sapiens"]
#print phylTree.dicLinks["Boreoeutheria"]
#print phylTree.dicParents["Homo sapiens"]
#print phylTree.dicParents["Boreoeutheria"]

print phylTree.tmpItems

sys.exit(0)

for e in phylTree.listSpecies:

	genome = utils.myGenomes.EnsemblGenome("~/work/data43/genes/genes.%s.list.bz2" % phylTree.fileName[e])
	
	s = 0
	lst = {}
	for c in genome.lstChr + genome.lstScaff:
		s += len(genome.lstGenes[c])
		lst[len(genome.lstGenes[c])] = lst.get(len(genome.lstGenes[c]), 0) + 1
		

	lst2 = lst.items()
	lst2.sort(reverse=True)

	p = s
	nbc = 0
	for (l,n) in lst2:
		nbc += n
		p -= l
		if p*100 < s*1:
			break

	print nbc, "\t", e

sys.exit(0)

genome = utils.myGenomes.EnsemblGenome("~/work/data43/genes/genes.Gallus.gallus.list.bz2")
#for g in genome.lstGenes[4]:
#	print g.beginning, g.strand, g.names

genesAnc = utils.myGenomes.EnsemblOrthosListGenome("~/work/data43/orthologs/orthos.Gallus.gallus.Homo.sapiens.list.bz2", ancFilter = ["Amniota"])
(c,i) = genesAnc.dicGenes["ENSGALG00000008006"]
print genesAnc.lstGenes[c][i].names
#print genesAnc.dicGenes["ENSGALG00000007993"].names
#print genesAnc.dicGenes["ENSGALG00000008006"].names
#print genesAnc.dicGenes["ENSG00000077279"].names
#print genesAnc.dicGenes["ENSG00000077274"].names

sys.exit(0)

utils.myPsOutput.printPsHeader()
utils.myPsOutput.drawLine(10, 10, 1, 1, utils.myPsOutput.getColor(1, "black"))
utils.myPsOutput.drawLine(10, 10, -1, 1, utils.myPsOutput.getColor(2, "black"))
utils.myPsOutput.drawLine(10, 10, -1, -1, utils.myPsOutput.getColor(3, "black"))
utils.myPsOutput.drawLine(10, 10, 1, -1, utils.myPsOutput.getColor(4, "black"))
utils.myPsOutput.printPsFooter()


sys.exit(0)
for i in xrange(8):
	for j in xrange(3):
		s = -1
		if (i & (1 << j)) > 1:
			s = 1
		print i, j, j, s, "gene.%d.%d" % (i,j)



sys.exit(0)
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
phylTree = utils.myPhylTree.PhylogeneticTree(noms_fichiers["phylTree.conf"])


print phylTree.parent['Homo sapiens']
print phylTree.parent['Human']
print phylTree.species['Boreoeutheria']
print phylTree.species['Glires']
print phylTree.branchesSpecies['Vertebreta']

print phylTree.fileName


