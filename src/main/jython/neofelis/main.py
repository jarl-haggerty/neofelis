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
from neofelis import genemark
from neofelis import extend
from neofelis import intergenic
from neofelis import promoters
from neofelis import terminators
from neofelis import artemis
from neofelis import utils

if __name__ == "__main__":
  print sys.argv
  opts, args = getopt(sys.argv, "m:d:g:b:e:rq:h", ["matrix=", "database=", "genemark=", "blast=", "eValue=", "remote", "query=", "help"])
  documentation = """
-m --matrix        :Matrix with which to run genemark
-d --database      :Database to use when running blast
-g --genemark      :Location of genemark
-b --blast         :Location of blast
-l --min-length    :Minimum length of any genes discovered
-e --eValue        :Minimal evalue for any genes detected
-r --remote        :Run blast with remote NCBI servers.
-h --help          :Print help documentation
-q --query         :Genome or directory of genomes to run pipeline on
"""
  
  blastLocation = "~/blast"
  database = "nr"
  genemarkLocation = "~/genemark"
  eValue = 0.1
  matrix = None
  minLength = 100
  remote = False
	
  for opt, arg in opts:
    if opt in ("-q", "--query"):
      sources = [arg]
    elif opt in ("-g", "--genemark"):
      genemarkLocation = arg
    elif opt in ("-b", "--blast"):
      blastLocation = arg
    elif opt in ("-d", "--database"):
      database = arg
    elif opt in ("-m", "--matrix"):
      matrix = arg
    elif opt in ("-e", "--eValue"):
      eValue = float(arg)
    elif opt in ("-l", "--min-length"):
      minLength = int(arg)
    elif opt in ("-r", "--remote"):
      remote = True
    elif opt in ("-h", "--help"):
      print documentation
      sys.exit(0)

  queries = []
  while sources:
    source = sources.pop()
    if os.path.isdir(source):
      newSources = filter(lambda x: x[0] != ".", os.listdir(query))
      newSources = map(lambda x: os.path.join(source, x), newSources)
      sources.extend(newSources)
    else:
      queries.append(source)
  print queries

  for query in queries:
    name = os.path.splitext(query)[0]
    name = os.path.split(name)[1]
    genome = utils.loadGenome(query)
    
    initialGenes = genemark.findGenes(query, name, blastLocation, database, eValue, genemarkLocation, matrix, remote)
    artemis.writeArtemisFile(name + ".art", genome, initialGenes.values())
    extendedGenes = extend.extendGenes(query, initialGenes, name, blastLocation, database, eValue, remote)
    artemis.writeArtemisFile(name + "extended.art", genome, extendedGenes.values())
    intergenicGenes = intergenic.findIntergenics(query, extendedGenes, name, minLength, blastLocation, database, eValue, remote)
    """
    for k, v in intergenicGenes.items():
      print k
      print v
    """
    genes = {}
    for k, v in extendedGenes.items() + intergenicGenes.items():
      genes[k] = v
    artemis.writeArtemisFile(name + "intergenic.art", genome, genes.values())
    scaffolds = scaffolds.extractScaffolds(genes)
    artemis.writeArtemisFile(name + "scaffolds.art", genome, scaffolds.values())

    initialPromoters = promoters.findPromoters(query)
    initialTerminators = terminators.findTerminators(query, transterm)
    
                
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
