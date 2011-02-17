"""
This module use used to predict genes with Genemark and then annotate them using
the function findGenes.
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

from neofelis import utils
import sys
import re
import os
import subprocess


#Have to make sure that a directory to store the blasts this module creates exists.
if not os.path.isdir("initialBlasts"):
  os.mkdir("initialBlasts")

def modifyFastaHeader(fileName, name):
  """
  Modify the headers of the genes in the fasta file named by fileName so that
  they contain name.
  """
  input = open(fileName, "r")
  swap = ""
  for line in input:
    matches = re.search(">orf_(\d+).*, (\d+ - \d+)", line)
    if matches:
      swap += ">" + name + "~" + ":".join(matches.groups()).replace(" ", "") + "\n"
    else:
      swap += line
  input.close()

  output = open(fileName, "w")
  output.write(swap)
  output.close()

def removeInvalidGenes(fileName, genomeLength):
  """
  This function exists because there was a case where genemark predicted a gene that terminated
  2 base pairs after the end of the genome.  This function will remove any genes that
  start or stop outside the genome.
  """
  input = open(fileName, "r")
  swap, reading = "", True
  for line in input:
    location = re.search(r">orf_\d+.*, (\d+) - (\d+)", line)
    if location:
      location = [int(q) for q in location.groups()]
      reading = location[0] >= 1 and location[1] >= 1 and location[0] <= genomeLength and location[1] <= genomeLength
    if reading:
      swap += line
  input.close()

  output = open(fileName, "w")
  output.write(swap)
  output.close()

def findGenes(query, name, blastLocation, database, eValue, genemark, matrix, remote):
  """
  Uses genemark to predict genes in query and then uses blast with the given eValue
  to find annotations for those genes.  If a matrix is not specified the GC program in
  genemark will be used to select a heuristic matrix.  Name is a name used
  to refer to the genome(typically query without the file extension), and for retrieving any
  cached blast results.
  """
  genome = utils.loadGenome(query)
  if not matrix:
    gc = int(utils.getGCContent(genome))
    matrix = genemark + "/" + "heuristic_mat/heu_11_" + str(gc) + ".mat"
  subprocess.Popen([genemark + "/gm", "-opq", "-m", matrix, query]).wait()
  print
  removeInvalidGenes(query + ".orf", len(genome))
  modifyFastaHeader(query + ".orf", name)
  
  result = utils.cachedBlast("initialBlasts/" + name + ".blastp.xml", blastLocation, database, eValue, query + ".orf", remote)
  print len(utils.loadGenome(query))
  #os.remove(query + ".orf")
  #os.remove(query + ".lst")
  return result
