#! /users/ldog/muffato/python -OO

import utils.myPsOutput
import utils.myPhylTree

(noms,_) = utils.myTools.checkArgs(["phylTree.conf"],[],"")

phylTree = utils.myPhylTree.PhylogeneticTree(noms["phylTree.conf"])
print phylTree.convertToFlatFile(phylTree.root), ";"

