#!/usr/bin/env python2

__doc__ = """
	Convertit un genome (suite de diagonales) en genome (suite de genes)
"""

import sys

import utils.myTools
import utils.myGenomes

arguments = utils.myTools.checkArgs( [("diagsFile",file), ("ancGenesFile",file)], [], __doc__)

ancGenes = utils.myGenomes.Genome(arguments["ancGenesFile"])
genome = utils.myGenomes.Genome(arguments["diagsFile"], ancGenes=ancGenes)
genome.printEnsembl(sys.stdout)

