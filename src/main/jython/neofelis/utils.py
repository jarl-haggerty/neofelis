"""
Utilities for the Neofelis Annotation Pipeline.
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

from javax.xml.parsers import DocumentBuilderFactory
from org.python.core.util import FileUtil
from java.lang import Double
from org.xml.sax.helpers import XMLReaderFactory
from org.xml.sax import XMLReader
from org.xml.sax.helpers import DefaultHandler
import os
import subprocess
import re

"""Start and stop codons."""
startCodons = ("ATG", "GTG", "TTG")
stopCodons = ("TGA", "TAA", "TAG")

"""Dictionary containing the complements of Nucleotides"""
complementDictionary = {}
complementDictionary["A"] = "T"
complementDictionary["T"] = "A"
complementDictionary["G"] = "C"
complementDictionary["C"] = "G"

"""Dictionary containing the mappings from codons to proteins"""
translationDictionary = {}
translationDictionary["TTT"] = "F"
translationDictionary["TTC"] = "F"
translationDictionary["TTA"] = "L"
translationDictionary["TTG"] = "L"
translationDictionary["CTT"] = "L"
translationDictionary["CTC"] = "L"
translationDictionary["CTA"] = "L"
translationDictionary["CTG"] = "L"
translationDictionary["TCT"] = "S"
translationDictionary["TCC"] = "S"
translationDictionary["TCA"] = "S"
translationDictionary["TCG"] = "S"
translationDictionary["TAT"] = "Y"
translationDictionary["TAC"] = "Y"
translationDictionary["TGA"] = "*"
translationDictionary["TAA"] = "*"
translationDictionary["TAG"] = "*"
translationDictionary["TGT"] = "C"
translationDictionary["TGC"] = "C"
translationDictionary["TGG"] = "W"
translationDictionary["CCA"] = "P"
translationDictionary["CCC"] = "P"
translationDictionary["CCG"] = "P"
translationDictionary["CCT"] = "P"
translationDictionary["CAC"] = "H"
translationDictionary["CAT"] = "H"
translationDictionary["CAA"] = "Q"
translationDictionary["CAG"] = "Q"
translationDictionary["CGA"] = "R"
translationDictionary["CGC"] = "R"
translationDictionary["CGG"] = "R"
translationDictionary["CGT"] = "R"
translationDictionary["AGA"] = "R"
translationDictionary["AGG"] = "R"
translationDictionary["ATT"] = "I"
translationDictionary["ATC"] = "I"
translationDictionary["ATA"] = "I"
translationDictionary["ATG"] = "M"
translationDictionary["ACA"] = "T"
translationDictionary["ACC"] = "T"
translationDictionary["ACG"] = "T"
translationDictionary["ACT"] = "T"
translationDictionary["AAT"] = "N"
translationDictionary["AAC"] = "N"
translationDictionary["AAA"] = "K"
translationDictionary["AAG"] = "K"
translationDictionary["AGT"] = "S"
translationDictionary["AGC"] = "S"
translationDictionary["GTA"] = "V"
translationDictionary["GTT"] = "V"
translationDictionary["GTG"] = "V"
translationDictionary["GTC"] = "V"
translationDictionary["GCA"] = "A"
translationDictionary["GCT"] = "A"
translationDictionary["GCG"] = "A"
translationDictionary["GCC"] = "A"
translationDictionary["GAT"] = "D"
translationDictionary["GAC"] = "D"
translationDictionary["GAA"] = "E"
translationDictionary["GAG"] = "E"
translationDictionary["GGA"] = "G"
translationDictionary["GGT"] = "G"
translationDictionary["GGG"] = "G"
translationDictionary["GGC"] = "G"

def translate(input):
  """
  Returns a Neucleotide sequence translated into Proteins.
  """
  result = ""
  for i in xrange(len(input)-2):
    result += translationDictionary[input[i:i+3]]
  return result

def reverseComplement(input):
  """
  Returns the reverse complement of a Neucleotide sequence.
  """
  result = map(lambda x: complementDictionary[x], input)
  result.reverse()
  return "".join(result)

def loadGenome(fileName):
  """
  Loads the genome from a fasta file containing a single genome.
  """
  input = open(fileName, "r")
  result = ""
  for line in input:
    match = re.match("([ACGT]+)\n", line.upper())
    if match:
      result += match.group(1)
  return result

def getGeneLocations(genes):
  """
  Takes a a map with GeneStructs as values and returns a two Dictionaries.
  These Dictionies will contain as tuples the left and right ends of genes
  such that calling genomes.[left:right] will return the gene.  the first
  dictionary is for genes on the forward strand and the second for the
  reverse
  """
  forward = {}
  reverse = {}
  for k, v in genes.items():
    if v.location[0] < v.location[1]:
      forward[k] = [v.location[0]-1, v.location[1]]
    else:
      reverse[k] = [v.location[1]-1, v.location[0]]
  return forward, reverse

def cachedBlast(fileName, blastLocation, database, eValue, query):
  """
  Performs a blast search using the blastp executable and database in blastLocation on
  the query with the eValue.  The result is an XML file saved to fileName.  If fileName
  already exists the search is skipped.
  """
  if not os.path.isfile(fileName):
    print "Caching"
    output = open(fileName, "w")
    temp = [blastLocation + "/bin/blastp",
                      "-db", blastLocation + "/db/" + database,
                      "-num_threads", str(Runtime.getRuntime().availableProcessors()),
                      "-evalue", str(eValue),
                      "-outfmt", "5",
                      "-query", query + ".orf"]
    print temp
    subprocess.Popen(temp,
                     stdout = output).wait()
    output.close()
    print "Cached"

class Iteration:
  """
  A structure for holding information about a gene's blast result.
  """
  query =           None
  location =        []
  numHits =         0
  bitScore =        0
  eValue =          Double.POSITIVE_INFINITY
  identity =        0
  alignmentLength = 0
  id =              None
  title =           None
  organism =        None

  def __str__(self):
    result = "<"
    result += "Query = " + str(self.query) + ", "
    result += "Location = " + str(self.location) + ", "
    result += "NumHits = " + str(self.numHits) + ", "
    result += "BitScore = " + str(self.bitScore) + ", "
    result += "EValue = " + str(self.eValue) + ", "
    result += "Identity = " + str(self.identity) + ", "
    result += "AlignmentLength = " + str(self.alignmentLength) + ", "
    result += "ID = " + str(self.id) + ", "
    result += "Title = " + str(self.title) + ", "
    result += "Organism = " + str(self.organism)
    result += ">"
    return result

class Hit:
  """
  A structure for holding information about a hit.
  """
  eValue = Double.POSITIVE_INFINITY
  bitScore = 0
  identity = 0
  alignmentLength = 0
  id = None
  title = None
  organism = None

  def __str__(self):
    result = "<"
    result += "BitScore = " + str(self.bitScore) + ", "
    result += "EValue = " + str(self.eValue) + ", "
    result += "Identity = " + str(self.identity) + ", "
    result += "AlignmentLength = " + str(self.alignmentLength) + ", "
    result += "ID = " + str(self.id) + ", "
    result += "Title = " + str(self.title) + ", "
    result += "Organism = " + str(self.organism)
    result += ">"
    return result

class Hsp:
  """
  A structure for holding information about a Hsp.
  """
  eValue = Double.POSITIVE_INFINITY
  bitScore = 0
  identity = 0
  alignmentLength = 0

class BlastHandler(DefaultHandler):
  """
  A SAX handler for parsing Blast XML output.
  """
  iterations = []
  hits = []
  hsps = []
  tag = None
  
  def startElement(self, uri, tag, name, attributes):
    """
    Records the tag of the current node and generates a new
    object to store the information in the iteration,
    hit, and hsp nodes.
    """
    if name == "Iteration":
      self.iterations += [Iteration()]
    elif name == "Hit":
      self.hits += [Hit()]
    elif name == "Hsp":
      self.hsps += [Hsp()]
    self.tag = tag

  def endElement(self, uri, tag, name):
    """
    Calculates the contents of a Hit structure once the end of a Hit node has been reached,
    and calculates the contents of a Iteration structure once the end of an Iteration node
    has been reached.
    """
    if tag == "Iteration" and self.hits:
      bestHit = min(self.hits, key = lambda hit:hit.eValue)
      self.iterations[-1].eValue = bestHit.eValue
      self.iterations[-1].bitScore = bestHit.bitScore
      self.iterations[-1].identity = bestHit.identity
      self.iterations[-1].alignmentLength = bestHit.alignmentLength
      self.iterations[-1].id = bestHit.id
      self.iterations[-1].title = bestHit.title
      self.iterations[-1].organism = bestHit.organism
      self.iterations[-1].numHits = len(self.hits)
      self.hits = []
    elif tag == "Hit":
      bestHsp = min(self.hsps, key = lambda hsp:hsp.eValue)
      self.hits[-1].eValue = bestHsp.eValue
      self.hits[-1].bitScore = bestHsp.bitScore
      self.hits[-1].identity = bestHsp.identity
      self.hits[-1].alignmentLength = bestHsp.alignmentLength
      self.hsps = []
    self.tag = ""

  def characters(self, raw, start, length):
    """
    Pulls the character information from the current node depending on the
    tag of the parent.
    """
    text = raw[start:start+length].tostring()
    if self.tag == "Iteration_query-def":
      self.iterations[-1].query, location = text.split(":")
      self.iterations[-1].location = [int(l) for l in location.split("-")]
    elif self.tag == "Hit_id":
      self.hits[-1].id = text
    elif self.tag == "Hit_def" and not self.hits[-1].title:
      match = re.search("(.+)\\[(.+)\\]", text)
      if match:
        self.hits[-1].title, self.hits[-1].organism = match.group(1), match.group(2)
      else:
        self.hits[-1].title, self.hits[-1].organism = text, ""
    elif self.tag == "Hsp_bit-score":
      self.hsps[-1].bitScore = float(text)
    elif self.tag == "Hsp_evalue":
      self.hsps[-1].eValue = float(text)
    elif self.tag == "Hsp_identity":
      self.hsps[-1].identity = float(text)
    elif self.tag == "Hsp_align-len":
      self.hsps[-1].alignmentLength = int(text)

def parseBlast(fileName):
  """
  A function for parsing XML blast output.
  """
  reader = XMLReaderFactory.createXMLReader()
  reader.setContentHandler(BlastHandler())
  reader.parse(fileName)

  return dict(map(lambda iteration: (iteration.query, iteration), reader.getContentHandler().iterations))

def getGCContent(genome):
  """
  A function for calculating the GC content of a genome.
  """
  return reduce(lambda x, y: x+int(y in ("G", "C")), genome, 0)/float(len(genome))*100
