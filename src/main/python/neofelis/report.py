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
import math
from org.xml.sax.helpers import XMLReaderFactory
from org.xml.sax import XMLReader
from org.xml.sax import InputSource
from org.xml.sax import EntityResolver
from org.xml.sax.helpers import DefaultHandler
from javax.xml.parsers import DocumentBuilderFactory
from java.lang import ClassLoader
from java.io import File

xmlDictionary = {"&" : "&amp;",
                 "\"" : "&quot;"}

def handleXMLCharacters(input):
    """
    Turns special characters into their xml form.
    """
    result = ""
    for q in input:
        result += xmlDictionary[q] if q in xmlDictionary else q
    return result

class HTMLWriter(DefaultHandler):
    """
    An xml parser that transforms an xml file into html.  The data is presented
    as iterations containing a list of hits, which then contain a list of hsps that
    show the alignment, evalue, bit score, and length of the hsp.
    """
    def __init__(self, output):
        self.output = output
        self.htmlDepth = 0
        self.querySequence = self.hitSequence = self.midline = self.tag = ""

    def writeHTMLStartTag(self, tag):
        self.output.write("  "*self.htmlDepth + "<" + tag + ">\n")
        self.htmlDepth += 1

    def writeHTMLEndTag(self, tag):
        self.htmlDepth -= 1
        self.output.write("  "*self.htmlDepth + "</" + tag + ">\n")

    def writeHTMLParagraph(self, text):
        self.output.write("  "*self.htmlDepth + "<p>" + text + "</p>\n")

    def writeHTMLLine(self, text):
        self.output.write("  "*self.htmlDepth + text + "\n")

    def writeHTMLDefinitionItem(self, text):
        self.output.write("  "*self.htmlDepth + "<dt>" + text + "</dt>\n")

    def startDocument(self):
        self.output = open(self.output, "w")
        self.writeHTMLStartTag("html")
        self.writeHTMLStartTag("body  style=\"font-family:monospace;white-space:nowrap;\"")

    def endDocument(self):
        self.writeHTMLEndTag("html")
        self.writeHTMLEndTag("body")
        self.output.write("\n")
        self.output.close()

    def startElement(self, uri, tag, name, attributes):
        self.text, self.tag = "", tag
        if tag == "BlastOutput_iterations":
            self.writeHTMLStartTag("dl")
            self.writeHTMLDefinitionItem("Iterations:")
            self.writeHTMLStartTag("dd")
        elif tag == "Iteration_hits":
            self.writeHTMLStartTag("dl")
            self.writeHTMLDefinitionItem("Hits:")
            self.writeHTMLStartTag("dd")
        elif tag == "Hit_hsps":
            self.writeHTMLStartTag("dl")
            self.writeHTMLDefinitionItem("Hsps:")
            self.writeHTMLStartTag("dd")

    def endElement(self, uri, tag, name):
        if tag == "BlastOutput_iterations":
            self.writeHTMLEndTag("dd")
            self.writeHTMLEndTag("dl")
        elif tag == "Iteration_hits":
            self.writeHTMLEndTag("dd")
            self.writeHTMLEndTag("dl")
        elif tag == "Hit_hsps":
            self.writeHTMLEndTag("dd")
            self.writeHTMLEndTag("dl")
        elif tag == "Iteration_query-def":
            self.writeHTMLStartTag("dl")
            self.writeHTMLDefinitionItem("Iteration: " + self.text)
            self.writeHTMLStartTag("dd")
        elif tag == "Iteration":
            self.writeHTMLEndTag("dd")
            self.writeHTMLEndTag("dl")
        elif tag == "Hit_def":
            self.writeHTMLStartTag("dl")
            self.writeHTMLDefinitionItem("Hit: " + self.text)
            self.writeHTMLStartTag("dd")
        elif tag == "Hit":
            self.writeHTMLEndTag("dd")
            self.writeHTMLEndTag("dl")
        elif tag == "Hsp_qseq":
            self.querySequence = self.text
        elif tag == "Hsp_hseq":
            self.hitSequence = self.text
        elif tag == "Hsp_evalue":
            self.eValue = self.text
        elif tag == "Hsp_midline":
            self.midline = self.text.replace(" ", "&nbsp;")
        elif tag == "Hsp_align-len":
            self.alignmentLength = self.text
        elif tag == "Hsp_identity":
            self.identity = self.text
        elif tag == "Hsp":
            self.writeHTMLStartTag("p")  
            self.writeHTMLLine("E Value = " + self.eValue + "<br/>")
            self.writeHTMLLine("Alignment Length = " + self.alignmentLength + "<br/>")
            self.writeHTMLLine("Identity = " + self.identity + "<br/>")
            self.writeHTMLLine(self.querySequence + "<br/>")
            self.writeHTMLLine(self.midline + "<br/>")
            self.writeHTMLLine(self.hitSequence + "<br/>")
            self.writeHTMLEndTag("p")  

    def characters(self, raw, start, length):
        self.text += handleXMLCharacters(raw[start:start+length].tostring())

    def resolveEntity(self, publicId, systemId):
        return InputSource(ClassLoader.getSystemResourceAsStream("dtds/" + os.path.split(systemId)[1]))

class BlastMerger(DefaultHandler):
    """
    A SAX handler for Deleting Genes.

    xml contents are echoed to the output depending on the printing and holding flags.  If printing is True then
    text may be written depending on holding.  If holding is true then chunks of text is saved in a list, which may
    simply be flushed, otherwise the text is written to the output.

    This parser starts with a root parse and a list of xml files.  Once the ending Iteration tag has been reached a
    new parser will be spawned to append the iterations in the next xml file.  The parent will start with printing
    set to True but the children won't so the contents of the xml outside the iterations are only written once.

    When an Iteration element is encountered printing and holding is set to True, so text is buffered.  Once the query definition
    is reached the buffered text and the rest of the iteration element will be written to the file if the definition is in the list
    of genes, otherwise it will be dumped and the rest of the iteration ignored.  Also, so iterations have unique numbers the iteration
    numbers are replaced with an incrementing counter.  Also, whitespace is buffered until a new element is reached so that empty
    lines can be removed before being written to a file.
    """
    def __init__(self, sources, genes, output, printing = False, parent = None, iteration = 1):
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
      self.iteration = iteration
      self.iterationNumberTag = False

    def startDocument(self):
        """
        If output is a string then nothing has been written, so it needs to be made into a file object
        and a header needs to be written.
        """
        if isinstance(self.output, basestring):
            self.output = open(self.output, "w")
            self.output.write("<?xml version=\"1.0\"?>\n")
            self.output.write("<!DOCTYPE BlastOutput PUBLIC \"-//NCBI//NCBI BlastOutput/EN\" \"NCBI_BlastOutput.dtd\">\n")
        
    def endDocument(self):
        self.output.write("\n")
        self.output.close()

    def startElement(self, uri, tag, name, attributes):

        if tag == "Iteration":
            self.printing = self.holding = True
        self.iterationNumberTag = tag == "Iteration_iter-num"
        self.iterationQueryDef = tag == "Iteration_query-def"
        self.iterationQueryDefString = ""

        if self.printing:
            if self.holding:
                self.text += re.sub("\n\s*\n", "\n", "".join(self.whitespace)) + "<" + tag + ">"
            else:
                self.output.write(re.sub("\n\s*\n", "\n", "".join(self.whitespace)) + "<" + tag + ">")
            self.whitespace = []

    def endElement(self, uri, tag, name):
        if tag == "BlastOutput_iterations":
            if self.sources:
                self.output.write("  ")
                reader = XMLReaderFactory.createXMLReader()
                reader.entityResolver = reader.contentHandler = BlastMerger(self.sources[1:], self.genes, self.output, False, self, self.iteration)
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
                self.text += re.sub("\n\s*\n", "\n", "".join(self.whitespace)) + "</" + tag + ">"
            else:
                self.output.write(re.sub("\n\s*\n", "\n", "".join(self.whitespace)) + "</" + tag + ">")
            self.whitespace = []

        if tag == "Iteration":
            self.holding = False
            self.printing = True

    def characters(self, raw, start, length):
        if self.iterationQueryDef:
            self.iterationQueryDefString += raw[start:start+length].tostring()
        if self.printing:
            if self.holding:
                if self.iterationNumberTag:
                    self.text += str(self.iteration)
                    self.iteration += 1
                else:
                    self.text += handleXMLCharacters(raw[start:start+length].tostring())
            else:
                self.output.write(handleXMLCharacters(raw[start:start+length].tostring()))

    def ignorableWhitespace(self, raw, start, length):
        self.whitespace += raw[start:start+length].tostring()

    def resolveEntity(self, publicId, systemId):
        return InputSource(ClassLoader.getSystemResourceAsStream("dtds/" + os.path.split(systemId)[1]))

class BreakParsingException(Exception):
    """
    Exception for halting xml parsing.
    """
    pass

def writeSpreadsheet(genes, output):
  """
  genes:  A list of Iteration objects.
  output: File name to write iterations to.

  Writes a spreadsheet of genes into output.
  """
  output = open(output + ".xls", "w")
  i = 0
  genes = sorted(filter(lambda x: x.numHits != 0, genes), key = lambda x: min(x.location)) + sorted(filter(lambda x: x.numHits == 0, genes), key = lambda x: min(x.location))
  output.write("Id\tStart\tStop\tHits\tBest bit score\tBest evalue\tBest identity\tAlignment Length\tBest hit gi\tDefinition\tOrganism\n")
  for gene in genes:
    i += 1
    output.write(str(i) + "\t")
    output.write("\t".join(map(str, gene.location)) + "\t")
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
  output.close()

def report(name, genes, output):
  """
  name:   Name of the genome.
  genes:  A dictionary that maps query names to Iteration objects.
  output: Output file name without an extension.

  Writes a report of the contents of the blast searchs for the queries in
  genes into "name.html", "name.blastp.xml", and "name.xls".
  """
  reader = XMLReaderFactory.createXMLReader()
  reader.entityResolver = reader.contentHandler = BlastMerger(["extendedBlasts/" + name + ".blastp.xml", "intergenicBlasts/" + name + ".blastp.xml"], genes.keys(), output + ".blastp.xml", True)
  reader.parse("initialBlasts/" + name + ".blastp.xml")

  reader = XMLReaderFactory.createXMLReader()
  reader.entityResolver = reader.contentHandler = HTMLWriter(output + ".blastp.html")
  reader.parse(output + ".blastp.xml")

  writeSpreadsheet(genes.values(), output)
