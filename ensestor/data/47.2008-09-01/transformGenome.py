#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2


# Librairies
import sys
import utils.myGenomes
import utils.myTools
import utils.myMaths
import utils.myPhylTree

# Arguments
############
arguments = utils.myTools.checkArgs( [("blocsList",file), ("genome",file)], [], "transforme un genome avec les indices des blocs")

lstBlocs = {}
f = utils.myTools.myOpenFile(arguments["blocsList"], "r")
for l in f:
	(i,_,bloc) = l.replace('\n','').partition('\t')
	lstBlocs[bloc] = i
f.close()

f = utils.myTools.myOpenFile(arguments["genome"], "r")
for l in f:
	t = l.replace('\n','').split('\t')
	bloc = t[2] + "\t" + t[3]
	if len(t[2].split()) == 1:
		t[1] = "0"
	print t[0] + "\t" + t[1] + "\t" + lstBlocs[bloc]
f.close()

