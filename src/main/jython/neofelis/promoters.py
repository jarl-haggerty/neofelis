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
