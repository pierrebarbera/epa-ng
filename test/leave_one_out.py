#!/usr/bin/python

import sys
import os
from tempfile import mkdtemp
from dendropy import *
from subprocess import call
import glob

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

if len(sys.argv) != 6:
    print "incorrect number of arguments"
    help()
    exit()

raxml = sys.argv[1]
epa = sys.argv[2]
tree_file = os.path.abspath(sys.argv[3])
MSA_file = os.path.abspath(sys.argv[4])
output_dir = os.path.abspath(sys.argv[5])

if not os.path.isfile(raxml):
    wrng("raxml doesn't exist or isn't a file")
if not os.path.isfile(epa):
    wrng("epa doesn't exist or isn't a file")
if not os.path.isfile(tree_file):
    wrng("tree_file doesn't exist or isn't a file")
if not os.path.isfile(MSA_file):
    wrng("MSA_file doesn't exist or isn't a file")
if not os.path.isdir(output_dir):
    wrng("output_dir doesn't exist or isn't a directory")

# read in tree and MSA files
tree = Tree.get(path=tree_file, schema="newick", rooting="force-unrooted")
msa = DnaCharacterMatrix.get(path=MSA_file, schema="fasta")

num_failed = 0
num_run = 0

# for every tip:
for node in tree.leaf_node_iter():
    # trim taxon from the cloned tree
    lou_tree = tree.clone()
    to_prune = TaxonNamespace()
    to_prune.add_taxon(node.taxon)
    lou_tree.prune_taxa(to_prune)
    lou_tree.deroot()

    # write tree to tmp folder
    cur_outdir = os.path.join(output_dir, str(node.taxon))
    if not os.path.exists(cur_outdir):
        os.makedirs(cur_outdir)

    cur_treefile = os.path.join(cur_outdir, "tree.newick")
    lou_tree.write(path=cur_treefile, schema="newick", suppress_rooting=True)

    # call raxml with trimmed files
    print "calling raxml:"
    params = [raxml, "-f", "v", "-s", MSA_file, "-t", cur_treefile, "-n", "leave_one_out", "-m","GTRGAMMA",
    "-w", cur_outdir]
    print params
    ret = call(params, stdout=open(os.devnull, 'wb'))

    # call epa with trimmed files
    print "calling epa: "
    params = [epa, cur_treefile, MSA_file, "-oO", "-w", cur_outdir]
    print params
    ret = call(params, stdout=open(os.devnull, 'wb'))

    # call validation script on both jplace files, log to log file
    jplace_files = glob.glob(os.path.join(cur_outdir, "*.jplace"))
    assert len(jplace_files) == 2

    print "calling compare script:"
    params = ["./jplace_compare.py", "-v"] + jplace_files
    print params
    ret = call(params)

    num_failed += ret
    num_run += 1

print "Failed " + str(num_failed) + " out of " + str(num_run) + " tests"
