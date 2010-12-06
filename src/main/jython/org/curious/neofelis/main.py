import os
import sys
from getopt import getopt
from org.curious.neofelis import genemark
from org.curious.neofelis import extend
from org.curious.neofelis import intergenic
from org.curious.neofelis import promoters
from org.curious.neofelis import terminators

if __name__ == "__main__":
	opts, args = getopt(sys.argv[1:], "ht:c:i:m:p:d:b:g:e:n:lx:a:", ["help", 
                                                                         "transterm=", 
                                                                         "cut-offs=",
                                                                         "interval=",
                                                                         "matrix=", 
                                                                         "prefix=", 
                                                                         "database=", 
                                                                         "blastall=", 
                                                                         "genemark=", 
                                                                         "evalue=", 
                                                                         "min-length=", 
                                                                         "local", 
                                                                         "excel-suffix=", 
                                                                         "max-e="])
	documentation = """
-m --matrix        :Matrix with which to run genemark
-u --heuristically :Run genemark heuristically
-g --genemark      :Location of genemark
-b --blast         :Location of blast
-e --eValue        :Minimal evalue for any genes detected
-h --help          :Print help documentation
"""
	blastall = "~/blast"
	database = "nr"
	genemark = "~/genemark"
	evalue = 0.1
	matrix = None
	
	queries = []
	for arg in sys.argv[1:]:
		if arg.find('-') != 0:
			queries.append(arg)

        for opt, arg in opts:
                if arg in queries:
                        quieries.remove(arg)
                if opt in ("-h", "--heuristically"):
                        heuristically = True
                if opt in ("-g", "--genemark"):
                        genemark = arg
                if opt in ("-b", "--blast"):
                        blast = arg
                if opt in ("-d", "--database"):
                        database = arg
                if opt in ("-m", "--matrix"):
                        matrix = arg
                if opt in ("-e", "--eValue"):
                        e_value = arg
                        
	for query in queries:
		if not os.path.split(query)[1]:
			queries.remove(query)
			queries.extend([os.path.split(query)[0] + '/' + file for file in os.listdir(query) if file[0] != '.'])
			
	for query in queries:
                name = os.path.splitext(query)[0]
                name = os.path.split(name)[1]

                if matrix:
                        initialGenes = genemark.findGenes(query, name, genemark, matrix, blast, e_value)
                elif:
                        initialGenes = genemark.findGenes(query, name, genemark, blast, e_value)
                extendedGenes = extend.extendGenes(initial_genes, name, blast, e_value)
                intergenicGenes = intergenics.findIntergenics(query, extendedGenes, name, minLength, eValue)
                genes = {}
                for k, v in extendedGenes.items() + intergenicGenes.items():
                        genes[k] = v
                
                initialTerminators = terminators.findTerminators(query, transterm)
                initialPromoters = promoters.findPromoters(query)
                
                writeArtemisFile(genes, promoters, terminators)
                
		
		#remove all signals inside orfs
		for f in os.listdir("artemis_complete/"):
			for c in cut_offs:
				print name + "_" + c
				print f
				if name + "_" + c in f:
					signal_command_now = signal_command_now + "artemis_complete/" + f + " "
					
		print signal_command_now
					
		os.system(signal_command_now)
		
		#build fasta file of all orfs
		fasta_command_now = "python make_fasta.py " + q
		print fasta_command_now
		os.system(fasta_command_now)
		
