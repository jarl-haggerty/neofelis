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
from org.xml.sax.helpers import XMLReaderFactory
from org.xml.sax import XMLReader
from org.xml.sax import InputSource
from org.xml.sax.helpers import DefaultHandler
from java.lang import ClassLoader

class Node():
  def __init__(self, tag, attributes, parent):
    self.tag = tag
    self.parent = parent
    self.children = {}
    self.childrenList = []
    if parent:
      self.parent.childrenList.append(self)
      if tag in parent.children:
        parent.children[tag].append(self)
      else:
        parent.children[tag] = [self]
    self.text = None

  def __str__(self):
    content  = self.content+">" if isinstance(self.content, str) else str(map(str, self.content))+">"
    return "<tag = " + self.tag + ", content = " + content

  def __getitem__(self, index):
    return self.children[index]

  def clearChildren(self):
    self.children = {}
    self.childrenList = []

def writeNode(root, output, depth = 0):
  maxLength = reduce(lambda x, y: max(x, len(y)), root.children.keys(), 0)
  for node in root.childrenList:
    if node.text:
      output.write("    "*depth + node.tag + " "*(maxLength-len(node.tag)) + " = " + node.text + "\n")
    else:
      output.write("    "*depth + node.tag + "\n")
      writeNode(node, output, depth+1)

def writeSpreadsheet(genes, output):
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

class Collector(DefaultHandler):
  """
  A SAX handler for parsing Blast XML output.
  """
  def __init__(self, genes, output):
    self.genes = genes
    self.output = output
    self.root = None
    self.text = ""
  
  def startElement(self, uri, tag, name, attributes):
    """
    Records the tag of the current node and generates a new
    object to store the information in the iteration,
    hit, and hsp nodes.
    """
    self.root = Node(tag, attributes, self.root)
    self.text = ""

  def endElement(self, uri, tag, name):
    """
    Calculates the contents of a Hit structure once the end of a Hit node has been reached,
    and calculates the contents of a Iteration structure once the end of an Iteration node
    has been reached.
    """
    if not self.root.children:
      self.root.text = self.text

    if tag == "Iteration":
      query = self.root.children["Iteration_query-def"][0].text
      if query[:query.find(":")] in self.genes:
        writeNode(self.root.parent, self.output)
      self.root.parent.clearChildren()
      
    if self.root.parent:
      self.root = self.root.parent

  def characters(self, raw, start, length):
    """
    Pulls the character information from the current node depending on the
    tag of the parent.
    """
    self.text += raw[start:start+length].tostring()

  def resolveEntity(self, publicId, systemId):
    return InputSource(ClassLoader.getSystemResourceAsStream("dtds/" + os.path.split(systemId)[1]))

def report(name, genes, output):
  dataOutput = open(output + ".dat", "w")
  spreadsheetOutput = open(output + ".xls", "w")
  
  reader = XMLReaderFactory.createXMLReader()
  reader.setContentHandler(Collector(genes.keys(), dataOutput))
  reader.setEntityResolver(reader.getContentHandler())

  reader.parse("initialBlasts/" + name + ".blastp.xml")
  reader.parse("extendedBlasts/" + name + ".blastp.xml")
  reader.parse("intergenicBlasts/" + name + ".blastp.xml")

  writeSpreadsheet(genes.values(), spreadsheetOutput)
  
  dataOutput.close()
  spreadsheetOutput.close()
