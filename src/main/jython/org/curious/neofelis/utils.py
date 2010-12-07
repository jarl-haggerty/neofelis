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

from javax.xml.xpath import XPathFactory
from javax.xml.xpath import XPath
from javax.xml.xpath import XPathConstants
from javax.xml.parsers import DocumentBuilderFactory
from org.python.core.util import FileUtil
from java.lang import Double
from java.lang import Runtime
import os

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

"""Helper objects for parsing blast results."""
xPath = XPathFactory.newInstance().newXPath()
documentBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder()
iterationString = "/BlastOutput/BlastOutput_iterations/Iteration[%d]"
hitString = "/BlastOutput/BlastOutput_iterations/Iteration[%d]/Iteration_hits/Hit[%d]"
hspString = "/BlastOutput/BlastOutput_iterations/Iteration[%d]/Iteration_hits/Hit[%d]/Hit_hsps/Hsp[%d]"

def getDocument(fileName):
  """
  Returns the Document of a blast XML file.
  """
  return documentBuilder.parse(fileName)

def countIterations(document):
  """
  Counts the number of iterations in a blast Document.
  """
  return xPath.evaluate(iterationString[:-4], document, XPathConstants.NODESET).getLength()

def countIterationHits(document, n):
  """
  Counts the number of hits in the nth iteration of a blast Document.
  """
  return xPath.evaluate(hitString[:-4] % (n), document, XPathConstants.NODESET).getLength()

def countHsps(document, n, i):
  """
  Count the number of Hsps for the ith hit of the nth iteration of a blast Document.
  """
  return xPath.evaluate(hspString[:-4] % (n, i), document, XPathConstants.NODESET).getLength()

def getIterationValue(document, n, value):
  """
  Retreives the text contained in the XML node with the tag specified by value under the nth
  iteration in a blast Document.
  """
  return xPath.evaluate(iterationString % (n) + "/" + value, document)

def getHitValue(document, n, i, value):
  """
  Retreives the text contained in the XML node with the tag specified by value under the nth
  iteration and ith hit in a blast Document.
  """
  return xPath.evaluate(hitString % (n, i) + "/" + value, document)

def getHspValue(document, n, i, j, value):
  """
  Retreives the text contained in the XML node with the tag specified by value under the nth
  iteration, ith hit, jth hsp in a blast Document.
  """
  return xPath.evaluate(hspString % (n, i, j) + "/" + value, document)

def getMinHsp(document, n, i, measure):
  """
  Retreives the index of the smallest hsp in the nth iteration and ith hit
  as defined by the text in the node with the tag specified by measure under the hsps.
  """
  minValue, minHsp = Double.POSITIVE_INFINITY, 0
  for hsp in xrange(countHsps(document, iteration, hit)):
    value = float(xPath.evaluate(hspString % (iteration, hit, hsp) + "/" + measure, document))
    if value < minValue:
      minValue, minHsp = value, hsp
  return minHsp

def getMaxHsp(document, iteration, hit, measure):
  """
  Retreives the index of the largest hsp in the nth iteration and ith hit
  as defined by the text in the node with the tag specified by measure under the hsps.
  """
  maxValue, maxHsp = Double.NEGATIVE_INFINITY, 0
  for hsp in xrange(countHsps(document, iteration, hit)):
    value = float(xPath.evaluate(hspString % (iteration, hit, hsp) + "/" + measure, document))
    if value > maxValue:
      maxValue, maxHsp = value, hsp
  return maxHsp

class GeneStruct:
  """
  A structure for holding information about a gene's blast result.
  """
  location =        []
  numHits =         0
  bitScore =        0
  eValue =          Double.POSITIVE_INFINITY
  identity =        0
  alignmentLength = 0
  hitId =           -1
  title =           "None"
  organism =        "None"

def loadGenome(fileName):
  """
  Loads the genome from a fasta file containing a single genome.
  """
  input = open(fileName, "r")
  result = ""
  for line in input:
    match = re.match("([ACGT]+)\n", line)
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
    process = Runtime.getRuntime().exec([blastLocation + "/bin/blastp",
                                        "-db", blastLocation + "/db/" + database,
                                        "-num_threads", str(Runtime.getRuntime().availableProcessors()),
                                        "-evalue", str(eValue),
                                        "-outfmt", "5",
                                        "-query", query + ".orf"])
    process.waitFor()
    input = FileUtil.wrap(process.getInputStream())
    output = open(fileName, "w")
    output.write(input.read())
    output.close()
    input.close()
