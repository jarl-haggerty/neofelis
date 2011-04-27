"""
This module is for writing files that summarize the results of the pipeline.
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
import re
import copy
from xml.dom import minidom
from org.xml.sax.helpers import XMLReaderFactory
from org.xml.sax import XMLReader
from org.xml.sax import InputSource
from org.xml.sax import EntityResolver
from org.xml.sax.helpers import DefaultHandler
from javax.xml.parsers import DocumentBuilderFactory
from java.lang import ClassLoader
from javax.xml.transform import OutputKeys
from javax.xml.transform import TransformerFactory
from javax.xml.transform.dom import DOMSource
from javax.xml.transform.stream import StreamResult
from java.io import File

xmlDictionary = {"&" : "&amp;",
                 "\"" : "&quot;"}

def handleXMLCharacters(input):
    result = ""
    for q in input:
        result += xmlDictionary[q] if q in xmlDictionary else q
    return result

class BlastMerger(DefaultHandler):
    """
    A SAX handler for Deleting Genes.
    """
    def __init__(self, sources, genes, output, printing = False, parent = None):
      print sources
      self.sources = sources
      self.output = output
      self.genes = genes
      self.text = []
      self.iterationQueryDef = True
      self.iterationQueryDefString = ""
      self.holding = False
      self.parent = parent
      self.printing = printing
      self.whitespace = []

    def startDocument(self):
      if isinstance(self.output, str):
        self.output = open(self.output, "w")
        self.output.write("<?xml version=\"1.0\"?>\n")
        self.output.write("<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"NCBI_BlastOutput.dtd\">")

    def endDocument(self):
      self.output.write("\n")
      self.output.close()

    def startElement(self, uri, tag, name, attributes):
        if tag == "Iteration":
            self.printing = self.holding = True
        self.iterationQueryDef = tag == "Iteration_query-def"
        self.iterationQueryDefString = ""
            
        if self.printing:
            if self.holding:
                self.text += re.sub("\n\s+\n", "\n", "".join(self.whitespace)) + "<" + tag + ">"
            else:
                self.output.write(re.sub("\n\s+\n", "\n", "".join(self.whitespace)) + "<" + tag + ">")
            self.whitespace = []
    
    def endElement(self, uri, tag, name):
      if tag == "BlastOutput_iterations":
        if self.sources:
          self.output.write("  ")
          reader = XMLReaderFactory.createXMLReader()
          reader.entityResolver = reader.contentHandler = BlastMerger(self.sources[1:], self.genes, self.output, False, self)
          try:
            reader.parse(self.sources[0])
          except BreakParsingException:
            if self.parent:
              raise BreakParsingException()
            else:
              pass
        else:
          raise BreakParsingException()
        
      if self.iterationQueryDef:
        self.iterationQueryDef = False
        if self.iterationQueryDefString[:self.iterationQueryDefString.rfind(":")] in self.genes:
          self.printing = True
          self.holding = False
          self.output.write("".join(self.text))
          self.text = []
        else:
          self.printing = False
          self.holding = False
          self.text = []
              
      if self.printing:
        if self.holding:
          self.text += re.sub("\n\s+\n", "\n", "".join(self.whitespace)) + "</" + tag + ">"
        else:
          self.output.write(re.sub("\n\s+\n", "\n", "".join(self.whitespace)) + "</" + tag + ">")
        self.whitespace = []

      if tag == "Iteration":
        self.holding = False
        self.printing = True

    def characters(self, raw, start, length):
        if self.iterationQueryDef:
            self.iterationQueryDefString += raw[start:start+length].tostring()
        if self.printing:
            if self.holding:
                self.text += handleXMLCharacters(raw[start:start+length].tostring())
            else:
                self.output.write(handleXMLCharacters(raw[start:start+length].tostring()))

    def ignorableWhitespace(self, raw, start, length):
      self.whitespace += raw[start:start+length].tostring()

    def resolveEntity(self, publicId, systemId):
        return InputSource(ClassLoader.getSystemResourceAsStream("dtds/" + os.path.split(systemId)[1]))

class BreakParsingException(Exception):
    pass

def writeSpreadsheet(genes, output):
  """
  genes:  A list of Iteration objects.
  output: File object to write to.

  Writes a summary of genes into output.
  """
  genes = sorted(filter(lambda x: x.numHits != 0, genes), key = lambda x: min(x.location)) + sorted(filter(lambda x: x.numHits == 0, genes), key = lambda x: min(x.location))
  output.write("Id\tLocation\tHits\tBest bit score\tBest evalue\tBest identity\tAlignment Length\tBest hit gi\tDefinition\tOrganism\n")
  for gene in genes:
    output.write(gene.query + "\t")
    output.write("-".join(map(str, gene.location)) + "\t")
    output.write(str(gene.numHits) + "\t")
    output.write(str(gene.bitScore) + "\t")
    output.write(str(gene.eValue) + "\t")
    output.write(str(gene.identity) + "\t")
    output.write(str(gene.alignmentLength) + "\t")
    if gene.id == "None":
      output.write(gene.id + "\t")
    else:
      output.write(re.search(r"gi\|(\d+)\|", gene.id).groups()[0] + "\t")
    output.write(gene.title + "\t")
    output.write(gene.organism + "\n")

def report(name, genes, output):
  """
  name:   Name of the genome.
  genes:  A dictionary that maps query names to Iteration objects.
  output: Output file name without an extension.

  Writes a report of the contents of the blast searchs for the queries in
  genes into "<name>.dat" and "<name>.xls".
  """
  reader = XMLReaderFactory.createXMLReader()
  reader.entityResolver = reader.contentHandler = BlastMerger(["extendedBlasts/" + name + ".blastp.xml", "intergenicBlasts/" + name + ".blastp.xml"], genes.keys(), output + ".blastp.xml", True)
  reader.parse("initialBlasts/" + name + ".blastp.xml")

  #writeSpreadsheet(genes.values(), output)
