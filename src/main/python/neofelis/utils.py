"""
Utilities for the Neofelis Annotation Pipeline.
"""

from org.python.core.util import FileUtil
from java.lang import Double
from java.lang import Runtime
from java.lang import ClassLoader
from org.xml.sax.helpers import XMLReaderFactory
from org.xml.sax import XMLReader
from org.xml.sax import InputSource
from org.xml.sax.helpers import DefaultHandler
from org.xml.sax import SAXParseException
import os
import os.path
import subprocess
import re

"""Start and stop codons."""
startCodons = ("ATG", "GTG", "TTG")
stopCodons = ("TGA", "TAA", "TAG")

"""Dictionary containing the complements of Nucleotides"""
complementDictionary = {"A":"T", "T":"A", "G":"C", "C":"G"}

"""Dictionary containing the mappings from codons to proteins"""
translationDictionary = {"TTT":"F","TTC":"F",
                         "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
                         "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "AGT":"S", "AGC":"S",
                         "TAT":"Y", "TAC":"Y",
                         "TGA":"*", "TAA":"*", "TAG":"*",
                         "TGT":"C", "TGC":"C",
                         "TGG":"W",
                         "CCA":"P", "CCC":"P", "CCG":"P", "CCT":"P",
                         "CAC":"H", "CAT":"H",
                         "CAA":"Q", "CAG":"Q",
                         "CGA":"R", "CGC":"R", "CGG":"R", "CGT":"R", "AGA":"R", "AGG":"R",
                         "ATT":"I", "ATC":"I", "ATA":"I",
                         "ATG":"M",
                         "ACA":"T", "ACC":"T", "ACG":"T", "ACT":"T",
                         "AAT":"N", "AAC":"N",
                         "AAA":"K", "AAG":"K",
                         "GTA":"V", "GTT":"V", "GTG":"V", "GTC":"V",
                         "GCA":"A", "GCT":"A", "GCG":"A", "GCC":"A",
                         "GAT":"D", "GAC":"D",
                         "GAA":"E", "GAG":"E",
                         "GGA":"G", "GGT":"G", "GGG":"G", "GGC":"G"}

def translate(input):
  """
  Returns a Neucleotide sequence translated into Proteins.
  """
  result = ""
  for i in xrange(0, len(input)-2, 3):
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
    match = re.match("([ACGT]+)", line.upper())
    if match:
      result += match.group(1)
  input.close()
  return result

def getHeader(fileName):
  input = open(fileName, "r")
  for line in input:
    match = re.match(">(.+)", line.upper())
    if match:
      input.close()
      return match.group(1)
  input.close()
  return None

def isGenome(fileName):
  """
  Returns true if the file represented by fileName is a fasta file containing one genome.
  """
  if os.path.isdir(fileName):
    return False
  input = open(fileName, "r")
  if not re.match(">.+", input.next()):
    return False
  for line in input:
    if line.strip() and not re.match("[ACGT]+", line.upper()):
      return False
  return True

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

class Iteration:
  """
  A structure for holding information about a gene's blast result.
  """
  def __init__(self):
    self.query =           None
    self.location =        []
    self.numHits =         0
    self.bitScore =        0
    self.eValue =          Double.POSITIVE_INFINITY
    self.identity =        0
    self.alignmentLength = 0
    self.id =              "None"
    self.title =           "None"
    self.organism =        "None"
    self.note =            ""
    self.color =           "0 255 255"
    self.intergenic =      False
    
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
    result += "Organism = " + str(self.organism) + ", "
    result += "Intergenic = " + str(self.intergenic)
    result += ">"
    return result

class Hit:
  """
  A structure for holding information about a hit.
  """
  def __init__(self):
    self.eValue = Double.POSITIVE_INFINITY
    self.bitScore = 0
    self.identity = 0
    self.alignmentLength = 0
    self.id = None
    self.title = None
    self.organism = None

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
  def __init__(self):
    self.eValue = Double.POSITIVE_INFINITY
    self.bitScore = 0
    self.identity = 0
    self.alignmentLength = 0

class BlastHandler(DefaultHandler):
  """
  A SAX handler for parsing Blast XML output.
  """
  def __init__(self):
    self.iterations = []
    self.hits = []
    self.hsps = []
    self.tag = None
    self.text = ""
  
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
    self.tag, self.text = tag, ""

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
    elif self.tag == "Iteration_query-def":
      self.iterations[-1].query, location = self.text.split(":")
      self.iterations[-1].location = [int(l) for l in location.split("-")]
    elif self.tag == "Hit_id":
      self.hits[-1].id = self.text
    elif self.tag == "Hit_def":
      match = re.search(r"([^\[]+)\[([^\]]+)", self.text)
      if match:
        self.hits[-1].title, self.hits[-1].organism = match.group(1).strip(), match.group(2).strip()
      else:
        self.hits[-1].title, self.hits[-1].organism = self.text.strip(), ""
    elif self.tag == "Hsp_bit-score":
      self.hsps[-1].bitScore = float(self.text)
    elif self.tag == "Hsp_evalue":
      self.hsps[-1].eValue = float(self.text)
    elif self.tag == "Hsp_identity":
      self.hsps[-1].identity = float(self.text)
    elif self.tag == "Hsp_align-len":
      self.hsps[-1].alignmentLength = int(self.text)

  def characters(self, raw, start, length):
    """
    Pulls the character information from the current node depending on the
    tag of the parent.
    """
    self.text += raw[start:start+length].tostring()

  def resolveEntity(self, publicId, systemId):
    return InputSource(ClassLoader.getSystemResourceAsStream("dtds/" + os.path.split(systemId)[1]))

def parseBlast(fileName):
  """
  A function for parsing XML blast output.
  """
  reader = XMLReaderFactory.createXMLReader()
  reader.entityResolver = reader.contentHandler = BlastHandler()
  reader.parse(fileName)

  return dict(map(lambda iteration: (iteration.query, iteration), reader.getContentHandler().iterations))

def cachedBlast(fileName, blastLocation, database, eValue, query, pipeline, force = False):
  """
  Performs a blast search using the blastp executable and database in blastLocation on
  the query with the eValue.  The result is an XML file saved to fileName.  If fileName
  already exists the search is skipped.
  """
  if not os.path.isfile(fileName) or force:
    output = open(fileName, "w")
    command = [blastLocation + "/bin/blastp",
               "-evalue", str(eValue),
               "-outfmt", "5",
               "-query", os.getcwd() + "/" + query,
               "-num_threads", str(Runtime.getRuntime().availableProcessors()),
               "-db", os.path.split(database)[1]]

    blastProcess = subprocess.Popen(command,
                                    stdout = subprocess.PIPE,
                                    cwd = database)
    while blastProcess.poll() == None:
      output.write(blastProcess.stdout.read())
      if pipeline.exception:
        psProcess = subprocess.Popen(["ps", "aux"], stdout = subprocess.PIPE)
        awkProcess = subprocess.Popen(["awk", "/" + " ".join(command).replace("/", "\\/") + "/"], stdin = psProcess.stdout, stdout = subprocess.PIPE)
        for line in awkProcess.stdout:
          subprocess.Popen(["kill", "-9", re.split(r"\s+", line)[1]])
        output.close()
        raise pipeline.exception
    remaining = blastProcess.stdout.read()
    while remaining:
      output.write(remaining)
      remaining = blastProcess.stdout.read()

    output.close()

  try:
    return parseBlast(fileName)
  except SAXParseException:
    return cachedBlast(fileName, blastLocation, database, eValue, query, pipeline, True)

def getGCContent(genome):
  """
  A function for calculating the GC content of a genome.
  """
  return reduce(lambda x, y: x+int(y in ("G", "C")), genome, 0)/float(len(genome))*100

def isNaN(number):
  """
  Returns true if number actually is a number.
  """
  return number != number
