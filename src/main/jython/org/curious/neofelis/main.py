"""
This Module is the main file for the Neofelis Genome Annotator
"""

"""
Copyright 2010 Jarl Haggerty
This file is part of Neofelis.

Neofelis is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Neofelis is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with Neofelis.  If not, see <http://www.gnu.org/licenses/>.
"""

import os
import sys
from getopt import getopt
from org.curious.neofelis import genemark
from org.curious.neofelis import extend
from org.curious.neofelis import intergenic
from org.curious.neofelis import promoters
from org.curious.neofelis import terminators

from org.python.util import PythonInterpreter

interpreter = PythonInterpreter()
interpreter.exec("print 3+6");
sys.exit(0)

if __name__ == "__main__":
  print sys.argv
  opts, args = getopt(sys.argv, "m:d:g:b:e:h", ["matrix=", "database=", "genemark=", "blast=", "eValue=", "help"])
  documentation = """
-m --matrix        :Matrix with which to run genemark
-d --database      :Database to use when running blast
-g --genemark      :Location of genemark
-b --blast         :Location of blast
-e --eValue        :Minimal evalue for any genes detected
-h --help          :Print help documentation
"""
  
  blastLocation = "~/blast"
  database = "nr"
  genemarkLocation = "~/genemark"
  eValue = 0.1
  matrix = None
	
  queries = []
  for arg in sys.argv:
    if arg.find('-') != 0:
      queries.append(arg)
    for opt, arg in opts:
      if arg in queries:
        quieries.remove(arg)
      if opt in ("-g", "--genemark"):
        genemarkLocation = arg
      if opt in ("-b", "--blast"):
        blastLocation = arg
      if opt in ("-d", "--database"):
        database = arg
      if opt in ("-m", "--matrix"):
        matrix = arg
      if opt in ("-e", "--eValue"):
        eValue = arg
                        
  for query in queries:
		if not os.path.split(query)[1]:
			queries.remove(query)
			queries.extend([os.path.split(query)[0] + '/' + file for file in os.listdir(query) if file[0] != '.'])

  for query in queries:
    name = os.path.splitext(query)[0]
    name = os.path.split(name)[1]
    
    initialGenes = genemark.findGenes(query, name, blastLocation, database, eValue, genemarkLocation, matrix)

    for k, v in initialGenes.items():
      print k, v
    sys.exit(0)
    
    extendedGenes = extend.extendGenes(initial_genes, name, blastLocation, e_value)
    intergenicGenes = intergenics.findIntergenics(query, extendedGenes, name, minLength, eValue)
    genes = {}
    for k, v in extendedGenes.items() + intergenicGenes.items():
      genes[k] = v
        
    initialTerminators = terminators.findTerminators(query, transterm)
    initialPromoters = promoters.findPromoters(query)
                
    writeArtemisFile(genes, promoters, terminators)
                
"""
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
"""	
