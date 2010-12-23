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

from neofelis import utils
from HTMLParser import HTMLParser
import re

class Promoter():
  def __init__(self):
    self.LDF = 0
    self.position = None
    self.signal10Location = None
    self.signal35Location = None
    

class MyHTMLParser(HTMLParser):
  def __init__(self):
    self.recording = False
    self.text = ""
    self.promoters = []
  def handle_starttag(self, tag, attrs):
    self.recording = tag == "pre"
  def handle_endtag(self, tag):
    if tag == "pre":
      self.recording = False
      promoters = re.findall(r"\s+Promoter\s+Pos:\s+(\d+)\s+LDF-\s+(\d+\.?\d*)", self.text)
      tenBoxes = re.findall(r"\s+-10\s+box\s+as\s+pos.\s+(\d+)\s+([ACGT]+)\s+Score\s+\d+", self.text)
      thiryFiveBoxes = re.findall(r"\s+-35\s+box\s+as\s+pos.\s+(\d+)\s+([ACGT]+)\s+Score\s+\d+", self.text)
      def parsePromoter(promoter, tenBox, thiryFiveBox):
        result = Promoter()
        result.LDF = float(promoter[1])
        result.position = int(promoter[0])
        result.signal10Location = [int(tenBox[0]), int(tenBox[0])+len(tenBox[1])]
        result.signal35Location = [int(thirtyFiveBox[0]), int(thirtyFiveBox[0])+len(thirtyFiveBox[1])]
      self.promoters = map(parsePromoter, promoters, tenBoxes, thirtyFiveBoxes)
  def handle_data(self, text):
    if self.recording:
      self.text += text
    

def findPromoters(query, LDFCutoff):
  results = urllib.urlopen("http://linux1.softberry.com/cgi-bin/programs/gfindb/bprom.pl",
                           urllib.urlencode({"DATA" : utils.loadGenome(query)}))
  parser = BPROMParser()
  parser.feed(results.read())
  results.close()
  return filter(lambda x: x.LDF > LDFCutoff, parser.promoters)
