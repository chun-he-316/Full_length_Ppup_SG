#!/usr/bin/python

from goatools.obo_parser import GODag
from collections import defaultdict
import argparse
import sys
from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS


class read_bg :
	def __init__(self,bg):
		self.bg= bg

	def pocessAssosiation(self) :
		def trans(a):
			if a == 'biological_process' :
				return('BP')
			elif a == 'molecular_function' :
				return('MF')
			else :
				return('CC')
		adict = defaultdict(dict)
		obodag = GODag("go-basic.obo")
		bgflie = self.bg
		for line in open(bgflie,'r') :
			tdict = defaultdict(set)
			line = line.strip()
			gene = line.split('\t')[0]
			if ',' in line :
				gos = line.split('\t')[1].split(',')
			else :
				gos = [line.split('\t')[1]]
			for go in gos :
				r = obodag.query_term(go)
				if r is None :
					pass
				else :
					namespace = trans(obodag.query_term(go).namespace)
					tdict[namespace].update({go})
			for key,value in tdict.items() :
				adict[key][gene] = value
		return (adict)

	def pocessBggenes(self) :
		genes = []
		for line in open(self.bg,'r') :
			line = line.strip()
			genes.append(line.split('\t')[0])
		return (genes)

class enrich :
	def __init__(self,associations,bggenes,studygenes,out):
		self.associations = associations
		self.bggenes = bggenes
		self.studygenes = studygenes
		self.out = out

	def go(self) :
		obodag = GODag("go-basic.obo")
		goeaobj = GOEnrichmentStudyNS(
		self.bggenes, # List of mouse protein-coding genes
		self.associations, # geneid/GO associations
		obodag, # Ontologies
		alpha = 0.05,
		methods = ['fdr_bh'],
		propagate_counts = False)
		goea_results_all = goeaobj.run_study(self.studygenes)
		for t in ['BP','CC','MF'] :
			goea_results = [r for r in goea_results_all if r.NS == t]
			goeaobj.wr_xlsx(self.out + 'go_' + t + '.xlsx', goea_results)
#	def kegg(self) :

def main() :
	if len(sys.argv) == 1:
		print("Please input arguments. Use -h or --help to see manual.")
	else :
		parser = argparse.ArgumentParser()
		parser.add_argument("-bg", "--background", help="Input BackgroundFile in .txt format",default='')
		parser.add_argument("-gl", "--studygenes", help="Gene list to be enrichment", default='')
		parser.add_argument("-o", "--out",help="Output name",default="out")
		args = parser.parse_args()
	
	# read bg for enrichgo
		r = read_bg(args.background)
		assosiation = r.pocessAssosiation()
		bggenes = r.pocessBggenes()
	
	# get study genes
		studygenes=[]
		for line in open(args.studygenes,'r') :
			line = line.strip()
			studygenes.append(line)
	# run go
		e = enrich(assosiation,bggenes,studygenes,args.out)
		e.go()
	
	# run kegg
if __name__ == '__main__':
	main()
