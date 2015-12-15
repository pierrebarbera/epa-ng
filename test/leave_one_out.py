#!/usr/bin/python

import sys
import os.path
from tempfile import mkdtemp
from dendropy import *

def help():
    print "USAGE:\traxml_path epa_path tree_file MSA_file output_dir"
def wrng(msg):
    print msg
    help()
    exit()
def err(msg):
    print msg
    print "Aborting"
    exit()

utils.Logging.log_to_stdout()

if __name__ != "__main__":
    print "Not running as main program."
    exit()

if len(sys.argv) != 6:
    print "incorrect number of arguments"
    help()
    exit()

raxml = sys.argv[1]
epa = sys.argv[2]
tree_file = sys.argv[3]
MSA_file = sys.argv[4]
output_dir = sys.argv[5]

if not os.path.isfile(raxml)
    wrng("raxml doesn't exist or isn't a file")
if not os.path.isfile(epa)
    wrng("epa doesn't exist or isn't a file")
if not os.path.isfile(tree_file)
    wrng("tree_file doesn't exist or isn't a file")
if not os.path.isfile(MSA_file)
    wrng("MSA_file doesn't exist or isn't a file")
if not os.path.isdir(output_dir)
    wrng("output_dir doesn't exist or isn't a directory")

# read in tree and MSA files
tree = Tree.get(path=tree_file, schema="newick")
msa = DnaCharacterMatrix.get(path=MSA_file, schema="fasta")

# for every tip:

for node in tree.leaf_iter():
    tmp = tempfile.mkdtemp()
    # trim from the tree and MSA
    lou_tree = tree.clone()
    lou_msa = msa.clone()
    lou_tree.prune_taxa_with_labels([node.taxon])
    if not lou_msa.contains(node.taxon)
        err("MSA does not contain a pruned taxon")
    lou_msa.delitem(node.taxon)

    # write both to tmp tmp folder

    # call raxml with trimmed files

    # call epa with trimmed files

    # call validation script on both jplace files, log to log file

    # clear tmp folder
