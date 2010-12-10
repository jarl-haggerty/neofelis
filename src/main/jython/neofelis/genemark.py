"""
This Module contains findGenes, the function used to predict genes inside a genome
using genemark and then use blast to seach for annotations for those ganes.
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

from java.lang import Double
from java.lang import Runtime
from org.python.core import PyFile
from org.python.core.util import FileUtil
from neofelis import utils
import sys
import re
import os


"""Have to make sure that a directory to store the blasts this module creates exists."""
if not os.path.isdir("initialBlasts"):
  os.mkdir("initialBlasts")

def getGC(query, genemark):
  """
  Calculate the GC content of the query(a fasta file) with the gc program in the
  specified genemark directory.
  """
  process = Runtime.getRuntime().exec(genemark + "/gc " + query)
  process.waitFor()
  result = FileUtil.wrap(process.getInputStream()).read()
  return float(re.search("GC% = (\d+)", result).group(1))

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

def findGenes(query, name, blast, database, eValue, genemark, matrix = None):
  """
  Uses genemark to predict genes in query and then uses blast with the given eValue
  to find annotations for those genes.  If a matrix is not specified the GC program in
  genemark will be used to select a heuristic matrix.  Name is simply a name used
  to refer to the genome(typically query without the file extension).
  """
  if not matrix:
    gc = int(getGC(query, genemark))
    matrix = genemark + "/" + "heuristic_mat/heu_11_" + str(gc) + ".mat"
  Runtime.getRuntime().exec(genemark + "/gm -opq -m " + matrix + " " + query).waitFor()
  modifyFastaHeader(query + ".orf", name)
  utils.cachedBlast("initialBlasts/" + name + ".blastp.xml", blast, database, eValue, query)
  #Runtime.getRuntime().exec("rm " + query + ".orf")
  return utils.parseBlast("initialBlasts/" + name + ".blastp.xml")
