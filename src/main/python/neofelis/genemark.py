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
  fileName: The name of the file to modify.
  name:     The name of the genome that will be included in the header of each sequence
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
  fileName:     Name of the file to modify.
  genomeLength: Length of the genome that the sequences in fileName came from
  
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

def findGenes(query, name, blastLocation, database, eValue, genemark, matrix, pipeline):
  """
  query:         File name of the query.
  name:          Name of the genome in the query.
  blastLocation: Location of blast installation.
  database:      Name of the database to search.
  eValue:        E value to use when searching.
  genemark:      Location of the genemark installation.
  matrix:        Name of the matrix to use, or None
  
  
  Uses genemark to predict genes in query and then uses blast with the given eValue
  to find annotations for those genes.  If a matrix is not specified the GC program in
  genemark will be used to select a heuristic matrix.
  """
  genome = utils.loadGenome(query)
  if not matrix:
    gc = int(utils.getGCContent(genome))
    matrix = genemark + "/" + "heuristic_mat/heu_11_" + str(min(max(30, gc), 70)) + ".mat"
  genemarkProcess = subprocess.Popen([genemark + "/gm", "-opq", "-m", matrix, query], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
  while genemarkProcess.poll() == None:
    genemarkProcess.stdout.read()
    genemarkProcess.stderr.read()
  removeInvalidGenes(query + ".orf", len(genome))
  modifyFastaHeader(query + ".orf", name)
  
  result = utils.cachedBlast("initialBlasts/" + name + ".blastp.xml", blastLocation, database, eValue, query + ".orf", pipeline)
  os.remove(query + ".orf")
  os.remove(query + ".lst")
  return result
