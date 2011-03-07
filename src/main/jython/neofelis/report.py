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
from org.xml.sax.helpers import XMLReaderFactory
from org.xml.sax import XMLReader
from org.xml.sax import InputSource
from org.xml.sax.helpers import DefaultHandler
from java.lang import ClassLoader

class Node():
  """
  Class for storing the information of an XML node.
  """
  def __init__(self, tag, attributes, parent):
    """
    tag:        A string.
    attributes: A dictionary mapping strings to strings.
    parent:     A Node or None.
    
    Initializes the node with the tag as it's tag, attributes as attribute, and parent as the parent.
    If parent is not None then this Node is added to the parent's list of children
    """
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
    """
    Clears the list and dictionary of children.
    """
    self.children = {}
    self.childrenList = []

def writeNode(root, output, depth = 0):
  """
  root:   A Node object to write.
  output: A file object to write into.
  depth:  Amount of indentation.
  
  Writes root into output.
  """
  maxLength = reduce(lambda x, y: max(x, len(y)), root.children.keys(), 0)
  for node in root.childrenList:
    if node.text:
      output.write("    "*depth + node.tag + " "*(maxLength-len(node.tag)) + " = " + node.text + "\n")
    else:
      output.write("    "*depth + node.tag + "\n")
      writeNode(node, output, depth+1)

def writeSpreadsheet(genes, output):
  """
  genes:  A list of Iteration objects.
  output: File object to write to.

  Writes a summary of genes into output.
  """
  genes = sorted(filter(lambda x: x.numHits != 0, genes), key = lambda x: min(x.location)) + sorted(filter(lambda x: x.numHits == 0, genes), key = lambda x: min(x.location))
  output.write("Id\tStart\tStop\tHits\tBest bit score\tBest evalue\tBest identity\tAlignment Length\tBest hit gi\tDefinition\tOrganism\n")
  for gene in genes:
    identification = int(re.search(r"\w+~(\d+)", gene.query).group(1))
    identification = identification+(len(genes)) if gene.intergenic else identification
    output.write(str(identification) + "\t")
    output.write(str(gene.location[0]) + "\t")
    output.write(str(gene.location[1]) + "\t")
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
    """
    genes:  List of gene names.
    output: File object for writing to.
    """
    self.genes = genes
    self.output = output
    self.root = None
    self.text = ""
  
  def startElement(self, uri, tag, name, attributes):
    """
    Move down the tree and add the new root to the children of the old root,
    and clear the text.
    """
    self.root = Node(tag, attributes, self.root)
    self.text = ""

  def endElement(self, uri, tag, name):
    """
    If the root is a leaf node then the text is add to the root and then the
    root moves to the current root's parent.  Also, if the current root is an
    Iteration node then it is written to a file if it's name(without the coordinates)
    is in genes and removed from the tree.
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
    Pulls the character information from the current node.
    """
    self.text += raw[start:start+length].tostring()

  def resolveEntity(self, publicId, systemId):
    return InputSource(ClassLoader.getSystemResourceAsStream("dtds/" + os.path.split(systemId)[1]))

def report(name, genes, output):
  """
  name:   Name of the genome.
  genes:  A dictionary that maps query names to Iteration objects.
  output: Output file name without an extension.

  Writes a report of the contents of the blast searchs for the queries in
  genes into "<name>.dat" and "<name>.xls".
  """
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
