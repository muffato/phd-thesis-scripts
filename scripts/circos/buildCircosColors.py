#! /usr/bin/python2.4

##################
# INITIALISATION #
##################

# Librairies
import sys
import os

sys.path.append(os.environ['HOME'] + "/work/scripts/utils")
import myTools
import myPsOutput


########
# MAIN #
########

myTools.checkArgs( [], [], "Affiche sur stderr le fichier color.conf de circos")


myPsOutput.initColor()

print >> sys.stderr, "white = 255,255,255"
for c in myPsOutput.color:
	if c.isalpha() or c.isdigit():
		(r,g,b) = myPsOutput.colorTable[myPsOutput.color[c]]
		print >> sys.stderr, "%s = %d,%d,%d" % (c, 255*r, 255*g, 255*b)



