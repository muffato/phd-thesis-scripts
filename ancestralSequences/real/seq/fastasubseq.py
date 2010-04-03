#! /users/ldog/muffato/python

import sys

for l in sys.stdin:
	t = l.split()
	#print "fastasubseq cdna/%s --start %s --length %s > cds/%s" % (t[0], t[7], int(t[8])-int(t[7]), t[0])
	print "fastasubseq.sh %s %s %s" % (t[0], t[7], int(t[8])-int(t[7]))

#Rattus_norvegicus/ENSRNOT00000030329.fa vulgar  ENSRNOT00000030329 0 583 . ENSRNOT00000030329 0 1749 + 3203 M 583 1749

