#!/usr/bin/env python2

__doc__ = """
	Dessine les densites de coorientation
"""

import sys
import math
import collections

import utils.myFile
import utils.myTools
import utils.myMaths
import utils.myGenomes
import utils.myPsOutput

# Arguments
arguments = utils.myTools.checkArgs( [("dataFile",file), ("cpglist",file)], [], __doc__)

nx = 15.
ny = 8.
stepx = 1./nx
stepy = 10./ny
dx = (29.7/9.)/nx
dy = (21./11.)/ny

valcpg = {}
f = utils.myFile.openFile(arguments["cpglist"], "r")
for l in f:
	t = l[:-1].split()
	valcpg[t[0]] = 1. if float(t[5]) > 0.375 else 0.
f.close()

def mkTripleDict(t):
	dConv = collections.defaultdict(t)
	dDiv = collections.defaultdict(t)
	dCoor = collections.defaultdict(t)
	return  {"1/1": dCoor, "-1/-1": dCoor, "1/-1": dConv, "-1/1": dDiv}

toy = {"1/1": 49, "-1/-1": 49, "1/-1": 50, "-1/1": 51}
toy = {"1/1": 25, "-1/-1": 25, "1/-1": 50, "-1/1": 75}
toyS = {"1/1": 1, "-1/-1": 1, "1/-1": 21, "-1/1": 41, None: 61}
toyA = dict([(str(x),i+10) for (i,x) in enumerate([0, 5, 10, 15, 30, 55, 60, 90, 95, 105, 130, 175, 320, 340, 420, 510, 550, 570, 580])])
count = mkTripleDict(int)
tags = mkTripleDict(list)

f = utils.myFile.openFile(arguments["dataFile"], "r")
for l in f:
	
	if "None" in l:
		continue
	

	t = l[:-1].split('\t')

	if t[21] == "1":
		continue

	# x = log10(size)
	x = int(math.log10(int(t[4])-int(t[3])+1)/stepx) + nx
	# x = GC%
	#x = int((max(float(t[9])-20,0.)/10)/stepx) + nx
	# x = Age
	#x = int((float(t[6])/200+1)/stepx) + nx
	
	# y = GC%
	y = int(float(t[9])/stepy) + ny
	# y = Age
	#y = int((float(t[6])/6)/stepy) + ny+2
	#y = toyA[t[6]] + toyS[t[5]]
	# y = orientation
	#y = toy[t[5]]

	#if int(t[6]) < 100:
	#	continue
	
	count[t[5]][(x,y)] += 1
	#count[t[5]][(x,toyS[None]+toyA[t[6]])] += 1

	# tag = Age
	tags[t[5]][(x,y)].append(int(t[6]))
	# tag = log10(size)
	#tags[t[5]][(x,y)].append(int(math.log10(int(t[4])-int(t[3])+1)/stepx))
	#tags[t[5]][(x,y)].append(math.log10(int(t[4])-int(t[3])+1))
	# tag = GC%
	#tags[t[5]][(x,y)].append( float(t[9]) )
	#tags[t[5]][(x,toyS[None]+toyA[t[6]])].append( float(t[9]) )
	# tag = CpG
	#(g1,g2) = t[1].split('/')
	#(s1,s2) = t[5].split('/')
	#if (s1 == "-1") and (g1 in valcpg):
	#	tags[t[5]][(x,y)].append(valcpg[g1])
	#if (s2 == "1") and (g2 in valcpg):
	#	tags[t[5]][(x,y)].append(valcpg[g2])
	# tag = CNE
	#tags[t[5]][(x,y)].append(float(t[11]))
	# tag = %repeat
	#tags[t[5]][(x,y)].append(float(t[13]))
	
	#tags[t[5]][(x,y)].append(float(t[22]))

f.close()

#print >> sys.stderr, utils.myMaths.myStats.txtSummary(cpg.values())

keys = set()
keys.update(count["1/1"].keys())
keys.update(count["1/-1"].keys())
keys.update(count["-1/1"].keys())

# Age
(minTagValue,maxTagValue) = (15., 420.)
# GC
#(minTagValue,maxTagValue) = (0.3, 1.)
#(minTagValue,maxTagValue) = (1., 6.)
# GC
#(minTagValue,maxTagValue) = (30., 90.)
# Repeat (30/40/50)
#(minTagValue,maxTagValue) = (0., 30.)
# CNE
#(minTagValue,maxTagValue) = (2., 20.)
# CNE chien
#(minTagValue,maxTagValue) = (0., 20.)
#(minTagValue,maxTagValue) = (0., 1.)

U = -0.3
V = 0.3

utils.myPsOutput.printPsHeader(landscape=True)

for x in xrange(1, 10):
	x = int(x/stepx)
	utils.myPsOutput.drawLine(x*dx, 0, 0, 21, (0,0,0))
utils.myPsOutput.drawBox(0,ny*dy, 31,dy, (0,0,0), (0,0,0))

for y in xrange(10, 101, 10):
	y = int(y/stepy)
	utils.myPsOutput.drawLine(0, y*dy, 29.7, 0, (0,0,0))
utils.myPsOutput.drawBox(nx*dx,0, dx,31, (0,0,0), (0,0,0))

totx = mkTripleDict(float)
toty = mkTripleDict(float)
nbx = mkTripleDict(float)
nby = mkTripleDict(float)
nbx = collections.defaultdict(float)
nby = collections.defaultdict(float)
ax = collections.defaultdict(list)
ay = collections.defaultdict(list)



def color(nbConv, nbDiv, nbCoor):
	
	nbTot = float(nbConv+nbDiv+nbCoor)
	propConv = nbConv / nbTot
	propDiv = nbDiv / nbTot
	propCoor = nbCoor / nbTot

	R = int(propConv*255)
	G = int(propDiv*255)
	B = int(propCoor*255)
	
	return (R,G,B)



for (x,y) in keys:
	nbCoor = count["1/1"][(x,y)]
	nbConv = count["1/-1"][(x,y)]
	nbDiv = count["-1/1"][(x,y)]

	ltags = []
	ltags.extend(tags["1/1"][(x,y)])
	#ltags.extend(tags["1/-1"][(x,y)])
	#ltags.extend(tags["-1/1"][(x,y)])
	
	#nbTot = float(nbConv+nbDiv+nbCoor)
	nbTot = len(ltags)
	
	totx["1/1"][x] += count["1/1"][(x,y)]
	totx["1/-1"][x] += count["1/-1"][(x,y)]
	totx["-1/1"][x] += count["-1/1"][(x,y)]
	
	toty["1/1"][y] += count["1/1"][(x,y)]
	toty["1/-1"][y] += count["1/-1"][(x,y)]
	toty["-1/1"][y] += count["-1/1"][(x,y)]
	
	if nbTot == 0:
	#if (2*nbDiv) == 0:
		continue
	
	#Y = propConv
	Y = ((sum(ltags)/len(ltags))-minTagValue)/(maxTagValue-minTagValue)
	Y = max(Y, 0.)
	Y = min(Y, 1.)
	#U = propCoor * .872 - .436
	#V = propConv * 1.23 - .615
	#(R,G,B) = utils.myPsOutput.YUV2RGB((Y,.1,.5))
	(R,G,B) = utils.myPsOutput.YUV2RGB((Y,U,V))
	
	#(R,G,B) = color(nbConv, nbDiv, nbCoor)

	nbx[x] += len(ltags)
	nby[y] += len(ltags)
	ax[x].extend(ltags)
	ay[y].extend(ltags)
	
	if max(R,G,B) == 0:
		continue

	#utils.myPsOutput.drawBox(x*dx,y*dy, dx,dy, (R,G,B), (R,G,B))
	utils.myPsOutput.drawBox(x*dx,y*dy, dx,dy, (0,0,0), (R,G,B))



mnbx = max(nbx.values())
mnby = max(nby.values())
print >> sys.stderr, "maxX, maxY:", mnbx, mnby
for x in nbx:
	(R,G,B) = color(totx["1/-1"][x], totx["-1/1"][x], totx["1/1"][x])
	Y = ((sum(ax[x])/len(ax[x]))-minTagValue)/(maxTagValue-minTagValue)
	Y = max(Y, 0.)
	Y = min(Y, 1.)
	(R,G,B) = utils.myPsOutput.YUV2RGB((Y,U,V))
	utils.myPsOutput.drawBox(x*dx,ny*dy, dx,-(ny*nbx[x]*dy)/mnbx, (R,G,B), (R,G,B))
	utils.myPsOutput.drawBox(x*dx,ny*dy, dx,dy, (0,0,0), (R,G,B))

for y in nby:
	(R,G,B) = color(toty["1/-1"][y], toty["-1/1"][y], toty["1/1"][y])
	Y = ((sum(ay[y])/len(ay[y]))-minTagValue)/(maxTagValue-minTagValue)
	Y = max(Y, 0.)
	Y = min(Y, 1.)
	(R,G,B) = utils.myPsOutput.YUV2RGB((Y,U,V))
	utils.myPsOutput.drawBox(nx*dx,y*dy, -(nx*nby[y]*dx)/mnby,dy, (R,G,B), (R,G,B))
	utils.myPsOutput.drawBox(nx*dx,y*dy, dx,dy, (0,0,0), (R,G,B))


utils.myPsOutput.printPsFooter()

