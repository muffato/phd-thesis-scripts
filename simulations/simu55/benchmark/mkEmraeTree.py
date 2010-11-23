#!/usr/bin/env python2
# Copyright (c) 2010 IBENS/Dyogen : Matthieu MUFFATO, Alexandra LOUIS, Hugues ROEST CROLLIUS
# Mail : hrc@ens.fr or alouis@biologie.ens.fr
# Licences GLP v3 and CeCILL v2

import os
import sys
import math
import random
import itertools
import collections

import utils.myMaths
import utils.myTools
import utils.myPhylTree


# Arguments
arguments = utils.myTools.checkArgs([("phylTree.conf",file)], [], "")

# L'arbre phylogenetique
phylTree = utils.myPhylTree.PhylogeneticTree(arguments["phylTree.conf"])

