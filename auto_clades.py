#!/usr/bin/env python3

import os, sys
import dendropy as dp
import matplotlib.pyplot as plt
from sklearn import linear_model
from optparse import OptionParser
import numpy as np

##########################################
# Function to Get command line arguments #
##########################################

def parsecommandline():

	usage = "usage: %prog [options]"
	parser = OptionParser(usage=usage)
	parser.add_option("-g", "--goldtree", action="store", dest="gold", help="Gold standard tree in newick format", default="")
	parser.add_option("-o", "--outfile", action="store", dest="outfile", help="Output file for clades in csv format", default="")
	parser.add_option("-f", "--force", action="store_true", dest="force", help="Force overwrite of output file", default=False)
	parser.add_option("-n", "--maxclades", action="store", type="int", dest="maxclades", help="Maximum umber of clades to produce (if this is > median branch length then only clades up to the median will be output)", default=100)
	parser.add_option("-m", "--mintaxa", action="store", type="int", dest="mintaxa", help="Minimum number of taxa permitted in a clade", default=10)
	return parser.parse_args()

def checkoptions(options):
	if options.gold=="":
		print("Gold tree file is required (-g)")
		sys.exit()
	elif not os.path.isfile(options.gold):
		print("Cannot find gold tree file", options.gold)
		sys.exit()
	if options.outfile=="":
		print("Output file is required (-m)")
		sys.exit()
	if not options.force and os.path.isfile(options.outfile):
		print("Output file", options.outfile, "already exists. Either choose a new filename or use the force (-f) option to overwrite")
		sys.exit()


def clades_from_tree(goldTreeFile, outfile, maxclades=100, mintaxa=3):
	
	#Open the tree file
	tns = dp.TaxonNamespace()
	try:
		goldTree = dp.Tree.get(file=open(goldTreeFile, 'r'), schema="newick", preserve_underscores=True, taxon_namespace=tns)
	except:
		print("Failed to open tree file: "+goldTreeFile)
		return 1
		
	goldTree.encode_bipartitions()
	
	numtaxa=float(len(goldTree.leaf_nodes()))
	maxtaxa=int(numtaxa/2)
	
	nodebls=[]
	lengths=[]
	for node in goldTree.postorder_internal_node_iter():
		labels=[]
		for leaf in node.leaf_iter():
			labels.append(leaf.taxon.label)
		if len(labels)>mintaxa and len(labels)<maxtaxa:
			nodebls.append([node.edge.length, labels])
			lengths.append(node.edge.length)
	
	#print(np.mean(lengths), np.median(lengths), np.std(lengths))
	median=np.median(lengths)
	nodebls.sort()
	nodebls.reverse()
	output=open(outfile, "w")
	print("ID,Clade", file=output)
	for count, i in enumerate(nodebls):
		if i[0]>median and count<maxclades:
			for label in i[1]:
				print(label+",clade"+str(count+1), file=output)
		else:
			break
			
	output.close()
	print(count, "clades added to output file")
	
	

def main():
	(options, args) = parsecommandline()
	checkoptions(options)
	clades_from_tree(options.gold, options.outfile, options.maxclades, options.mintaxa)
	



if __name__ == "__main__":
	main()


