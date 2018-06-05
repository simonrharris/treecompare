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

	usage = "usage: %prog [options] <input file>"
	parser = OptionParser(usage=usage)
	parser.add_option("-t", "--testtree", action="store", dest="test", help="Test tree in newick format", default="")
	parser.add_option("-g", "--goldtree", action="store", dest="gold", help="Gold standard tree in newick format", default="")
	parser.add_option("-m", "--metadata", action="store", dest="metadata", help="Metadata file in csv format", default="")
	parser.add_option("-c", "--cladecolumn", action="store", dest="cladecolumn", help="Name of column in metadata containing clades of interest[default = %default]", default="clade")
	parser.add_option("-i", "--idcolumn", action="store", dest="idcolumn", help="Name of column in metadata containing id of isolates. This must match both trees. [default = %default]", default="id")
	parser.add_option("-p", "--png", action="store_true", dest="printpng", help="Print pngs of rtt plots. [Default = false]", default=False)
	parser.add_option("-r", "--remove_mismatches", action="store_true", dest="remove", help="Remove mismatching taxa between gold and test trees. [Default = false]", default=False)
	parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="Be more verbose. [Default = false]", default=False)
	return parser.parse_args()

def checkoptions(options):
	if options.test=="":
		print("Test tree file is required (-t)")
		sys.exit()
	elif not os.path.isfile(options.test):
		print("Cannot find test tree file", options.test)
		sys.exit()
	if options.gold=="":
		print("Gold tree file is required (-g)")
		sys.exit()
	elif not os.path.isfile(options.gold):
		print("Cannot find gold tree file", options.gold)
		sys.exit()
	if options.metadata=="":
		print("Metadata file is required (-m)")
		sys.exit()
	elif not os.path.isfile(options.metadata):
		print("Cannot find metadata file", options.metadata)
		sys.exit()
	if options.cladecolumn==options.idcolumn:
		print("ID and clade must be in different columns")
		sys.exit()

def get_rtt_distances(tree, taxa, scaling=1.0):
	#Calculate root to tip distances and scale by multiplying by the scaling factor if necessary (mainly so the plots show properly)
	rtt=[]
	tree.seed_node.edge_length=0
	for taxon in taxa:
		leafNode=tree.find_node_with_taxon_label(taxon)
		rtt.append([leafNode.distance_from_root()*scaling])
	
	return np.array(rtt)


def find_tree_mismatches(tree1, tree2):
	tree1taxa=set([])
	tree2taxa=set([])
	for node in tree1.leaf_node_iter():
		tree1taxa.add(node.taxon.label)
	for node in tree2.leaf_node_iter():
		tree2taxa.add(node.taxon.label)

	missing=tree1taxa.difference(tree2taxa)
	missing2=tree2taxa.difference(tree1taxa)
	return missing.union(missing2)


def transferIndex(tree,taxa):
	
	l=len(tree.leaf_nodes())
	if len(taxa)<(float(l)/2):
		p=len(taxa)
		lightside=True
	else:
		lightside=False
		p=totaltaxa-len(taxa)
	
	minti=[l+1,None]

	for node in tree.postorder_node_iter():
		if node.is_leaf():
			if node.annotations.get_value('l0')!=None:
				print(l, p, node.annotations.get_value('l0'), node.annotations.get_value('l1'))
			
			if node.taxon.label in taxa:
				if lightside:
					node.annotations.add_new(name='l0', value=1)
					node.annotations.add_new(name='l1', value=0)
				else:
					node.annotations.add_new(name='l0', value=0)
					node.annotations.add_new(name='l1', value=1)
			else:
				if lightside:
					node.annotations.add_new(name='l0', value=0)
					node.annotations.add_new(name='l1', value=1)
				else:
					node.annotations.add_new(name='l0', value=1)
					node.annotations.add_new(name='l1', value=0)
			
		else:
			l0=0
			l1=0
			for child in node.child_node_iter():
				l0+=child.annotations.get_value('l0')
				l1+=child.annotations.get_value('l1')
			node.annotations.add_new(name='l0', value=l0)
			node.annotations.add_new(name='l1', value=l1)
		
		ti=min(p-node.annotations.get_value('l0')+node.annotations.get_value('l1'), l-p-node.annotations.get_value('l0')+node.annotations.get_value('l1'))
		#print(l, p, node.annotations.get_value('l0'), node.annotations.get_value('l1'), p-node.annotations.get_value('l0')+node.annotations.get_value('l1'), l-p-node.annotations.get_value('l0')+node.annotations.get_value('l1'), ti)
		
		if ti<minti[0]:
			minti=[ti, round(1-(ti/(p-1)),4), node]
	
	#print(minti)
	interlopers=[]
	if minti[0]>0:
		missing=taxa[:]
		for leaf in minti[2].leaf_iter():
			if (lightside and leaf.annotations.get_value('l0')==1) or (not lightside and leaf.annotations.get_value('l1')==1):
				missing.remove(leaf.taxon.label)
			else:
				interlopers.append(leaf.taxon.label)
	else:
		missing=[]
	return minti[:2], missing, interlopers
		
	


def validate_tree(testTree, goldTree, clusterTaxa=[], cladeName="my_clade", printpng=False, verbose=False):
	
	if len(clusterTaxa)==1:
		print(prefix, 1, "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", sep=',')
		return
		
	outputline=[cladeName,len(clusterTaxa)]
	
	#calculate tranfser index for trees
	testTreeClone=testTree.clone(1)
	goldTreeClone=goldTree.clone(1)
	
	ti, testMissing, testInterlopers=transferIndex(testTreeClone,clusterTaxa)
	outputline+=ti
	ti, goldMissing, goldInterlopers=transferIndex(goldTreeClone,clusterTaxa)
	outputline+=ti
	
	#Calculate a scaling factor to scale branch lengths by the total tree length for the rtt calculations. This mostly helps pyplot which doesn't like the tiny distances from RAxmL for some reason
	#if testTreeLength>goldTreeLength:
	#	goldScaling=testTreeLength/goldTreeLength
	#	testScaling=1.0
	#	testLabel=testTreeFile+" root to tip distances"
	#	goldLabel="rescaled "+goldTreeFile+" root to tip distances"
	#else:
	#	goldScaling=1.0
	#	testScaling=goldTreeLength/testTreeLength
	#	testLabel="rescaled "+testTreeFile+" root to tip distances"
	#	goldLabel=goldTreeFile+" root to tip distances"
	
	to_prune=[]
	for l in testTreeClone.leaf_nodes():
		if not l.taxon.label in clusterTaxa:
			to_prune.append(l.taxon.label)
	
	testTreeClone.prune_taxa_with_labels(to_prune)
	goldTreeClone.prune_taxa_with_labels(to_prune)
	
	#Calculate the (scaled) root to tip distances for the clade, but only for the taxa in the cluster list
	testScaling=1
	goldScaling=1
	testRtt=get_rtt_distances(testTreeClone, clusterTaxa, testScaling)
	goldRtt=get_rtt_distances(goldTreeClone, clusterTaxa, goldScaling)
	
	#outputline.append(round(np.mean(testRtt),0))
	#outputline.append(round(np.mean(goldRtt),0))
	#outputline.append(round(np.median(testRtt),0))
	#outputline.append(round(np.median(goldRtt),0))
	
	#do the linear regression
	regr = linear_model.LinearRegression()
	regr.fit(goldRtt, testRtt)
	#print("Regression:")
	#print("R-square =", regr.score(goldRtt, testRtt))
	outputline.append(round(regr.score(goldRtt, testRtt),4))
	
	min_samples=int(3*(len(clusterTaxa)/4))
	min_samples=len(clusterTaxa)-1
	#print(min_samples, len(clusterTaxa))
	if min_samples<3:
		min_samples=3
	
	# Robustly fit linear model with RANSAC algorithm
	#ransac = linear_model.RANSACRegressor(min_samples=min_samples)
	#ransac = linear_model.RANSACRegressor(min_samples=3, max_trials=1000)
	#ransac.fit(goldRtt, testRtt)
	#inlier_mask = ransac.inlier_mask_
	#outlier_mask = np.logical_not(inlier_mask)
	#print("R-square =", ransac.score(goldRtt[inlier_mask], testRtt[inlier_mask]), sum(outlier_mask), sum(inlier_mask), min_samples)
	#print(inlier_mask, ransac.n_trials_)#, ransac.n_skips_no_inliers_, ransac.n_skips_invalid_data_, ransac.n_skips_invalid_model_)
	
	
	if printpng:
		#Plot the regression plot
		plt.figure()
		plt.scatter(goldRtt, testRtt, color='black')
		#plt.scatter(goldRtt[inlier_mask], testRtt[inlier_mask], color='black')
		#plt.scatter(goldRtt[outlier_mask], testRtt[outlier_mask], color='grey')
		plt.plot(goldRtt, regr.predict(goldRtt), color='red',linewidth=1)
		#plt.plot(goldRtt, ransac.predict(goldRtt), color='cornflowerblue',linewidth=1)
		plt.xlabel(goldLabel)
		plt.ylabel(testLabel)
		plt.savefig(prefix+"_rtt.png")
		plt.close()
	
	
	outputline.append(dp.calculate.treecompare.symmetric_difference(testTreeClone,goldTreeClone))
	missingbiparts=dp.calculate.treecompare.find_missing_bipartitions(testTreeClone,goldTreeClone)
	allbiparts=testTreeClone.encode_bipartitions(suppress_unifurcations=True)
	outputline.append(round(float(len(missingbiparts))/len(allbiparts),4))
	
	outputline.append('; '.join(testMissing))
	outputline.append('; '.join(testInterlopers))
	outputline.append('; '.join(goldMissing))
	outputline.append('; '.join(goldInterlopers))
	
	print(','.join(map(str,outputline)))


def main():
	(options, args) = parsecommandline()
	checkoptions(options)
	cladecolumnnumber=-1
	idcolumnnumber=-1
	clades={}
	for i, line in enumerate(open(options.metadata, "rU")):
		words=line.strip().split(",")
		if i==0:
			for x, word in enumerate(words):
				if word.lower()==options.cladecolumn.lower():
					cladecolumnnumber=x
				elif word.lower()==options.idcolumn.lower():
					idcolumnnumber=x
		else:
			if cladecolumnnumber==-1:
				print("Failed to find column called", options.cladecolumn, "in headler line of metadata file")
				sys.exit()
			elif idcolumnnumber==-1:
				print("Failed to find column called", options.idcolumn, "in headler line of metadata file")
				sys.exit()
			try:
				clade=words[cladecolumnnumber]
			except:
				print("Misformed metadata file on line", i+1)
				sys.exit()
			try:
				name=words[idcolumnnumber]
			except:
				print("Misformed metadata file on line", i+1)
				sys.exit()
			try:
				clades[clade].append(name)
			except:
				
				clades[clade]=[name]

	#Open the two tree files
	tns = dp.TaxonNamespace()
	try:
		testTree = dp.Tree.get(file=open(options.test, 'r'), schema="newick", preserve_underscores=True, taxon_namespace=tns)
	except:
		print("Failed to open tree file: "+testTreeFile)
		return 1
	
	
	try:
		goldTree = dp.Tree.get(file=open(options.gold, 'r'), schema="newick", preserve_underscores=True, taxon_namespace=tns)
	except:
		print("Failed to open tree file: "+goldTreeFile)
		return 1
		
	#testTree.encode_bipartitions()
	#goldTree.encode_bipartitions()

	print("Clade", "Size", "Test TI", "Test Scaled TI", "Gold TI", "Gold Scaled TI", "R^2", "RF distance", "Mismatching bipartitions", "Test Moved In", "Test Moved Out", "Gold Moved In", "Gold Moved Out", sep=',')
	for clade in clades:
		#if clade in ["clade14", "clade15", "2-4"]:
		validate_tree(testTree, goldTree, clades[clade], clade, options.printpng, options.verbose)
	#validate_tree(options.test, options.gold, [], "all")



if __name__ == "__main__":
	main()
		
