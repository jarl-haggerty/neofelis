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
from org.xml.sax.helpers import XMLReaderFactory
from org.xml.sax import XMLReader
from org.xml.sax import InputSource
from org.xml.sax.helpers import DefaultHandler
from java.lang import ClassLoader

class BlastTranslater(DefaultHandler):
  def __init__(self, output):
    self.output = open(output, "w")
    self.depth = 0
    self.lastDepth = -1
    self.buffer = ""
    self.texts = []
    self.previousStart = ""
    self.previousStop = ""
    self.nodes = []
    self.endTags = 0

  def endDocument(self):
    self.output.close()
  
  def startElement(self, uri, tag, name, attributes):
    if self.previousStart:
      maxLength = max(map(lambda x: len(x[0]), self.nodes)) if self.nodes else None
      for key, value in self.nodes:
        self.output.write((self.depth-1)*"    " + key + (maxLength-len(key))*" " + " = " + value + "\n")
      self.nodes = []
      self.output.write((self.depth-1)*"    " + self.previousStart + "\n")
      
    self.previousStart = tag
    self.previousStop = ""
    self.buffer = ""
    self.depth += 1

  def endElement(self, uri, tag, name):
    if self.previousStop:
      maxLength = max(map(lambda x: len(x[0]), self.nodes)) if self.nodes else None
      for key, value in self.nodes:
        self.output.write(self.depth*"    " + key + (maxLength-len(key))*" " + " = " + value + "\n")
      self.nodes = []
    else:
      self.nodes.append([tag, self.buffer])

    self.previousStart = ""
    self.previousStop = tag
    self.buffer = ""
    self.depth -= 1

  def characters(self, raw, start, length):
    self.buffer += raw[start:start+length].tostring()

  def resolveEntity(self, publicId, systemId):
    return InputSource(ClassLoader.getSystemResourceAsStream("dtds/" + os.path.split(systemId)[1]))

def xmlToText(input, output):
  reader = XMLReaderFactory.createXMLReader()
  translater = BlastTranslater(output)
  reader.setContentHandler(translater)
  reader.setEntityResolver(translater)
  reader.parse(input)
