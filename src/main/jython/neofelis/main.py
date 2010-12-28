"""
This Module is the main file for the Neofelis Genome Annotator
"""

"""
Copyright 2010 Jarl Haggerty

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0
       
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import os
import sys
from getopt import getopt
from neofelis import pipeline

if __name__ == "__main__":
  opts, args = getopt(sys.argv, "m:d:g:b:e:rq:h", ["matrix=", "database=", "genemark=", "blast=", "eValue=", "remote", "query=", "help"])
  documentation = """
-m --matrix               :Matrix with which to run genemark
-d --database             :Database to use when running blast
-g --genemark             :Location of genemark
-b --blast                :Location of blast
-l --min-length           :Minimum length of any genes discovered
-e --eValue               :Minimal evalue for any genes detected
-r --remote               :Run blast with remote NCBI servers.
-s --scaffolding-distance :Distance to allow between genes when determining scaffolds
-h --help                 :Print help documentation
-q --query                :Genome or directory of genomes to run pipeline on
"""
  
  blastLocation = "~/blast"
  database = "nr"
  genemarkLocation = "~/genemark"
  eValue = 0.1
  matrix = None
  minLength = 100
  scaffoldingDistance = 100
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
    elif opt in ("-s", "--scaffolding-distance"):
      scaffoldingDistance = int(arg)
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

  pipeline.run(blastLocation, genemarkLocation, database, eValue, matrix, minLength, scaffoldingDistance, remote, queries)

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
