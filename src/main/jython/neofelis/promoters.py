"""
This module contains classes and functions for predicting promoters in a genome.
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

from neofelis import utils
import re
import urllib
import functools
import os.path

#Make sure the directory for BPROM results exists.
if not os.path.isdir("promoterPredictions"):
  os.mkdir("promoterPredictions")

class Promoter():
  """
  Class for storing information about a promoter.
  """
  def __init__(self, start, stop, score):
    self.location = [start, stop]
    self.score = score

  def __str__(self):
    result = "<"
    result += "location = " + str(map(str, self.location)) + ", "
    result += "score = " + str(self.score)
    result += ">"
    return result

def parseFruitfly(fruitflyData):
  """
  bpromData: Contents of the body of a Fruitfly prediction.

  return:    A list of Promoter objects representing the results of the BPROM prediction
  """
  fruitflyData = fruitflyData.replace("<font size=\"+2\">", "").replace("</font>", "")
  forwardPromoters = re.findall(r"\s*(\d+)\s+\d+\s+(\d?\.?\d+)\s+([ACGT]+)<br>", fruitflyData[:fruitflyData.find("Predictions for the reverse strand")])
  print "Hi Joe", len(forwardPromoters)
  reversePromoters = re.findall(r"\s*\d+\s+(\d+)\s+(\d?\.?\d+)\s+([ACGT]+)<br>", fruitflyData[fruitflyData.find("Predictions for the reverse strand"):])
  result = map(lambda x: Promoter(int(x[0]), int(x[0])+len(x[2])-1, float(x[1])), forwardPromoters)
  result += map(lambda x: Promoter(int(x[0])+len(x[2]), int(x[0])+1, float(x[1])), reversePromoters)
  return result

def cachedPrediction(genome, fileName):
  """
  genome:   Genome as a string.
  fileName: File to save the prediction results in.

  return:   Results of the promoter prediction stored in a list of Promoter objects.

  If the file specified by fileName already exists then this function simply parses the file
  already there.
  """
  if not os.path.isfile(fileName):
    results = urllib.urlopen("http://www.fruitfly.org/cgi-bin/seq_tools/promoter.pl",
                             urllib.urlencode({"organism" : "prokaryote", "reverse" : "yes", "threshold" : "0.9", "text" : genome}))
    resultString = results.read()
    results.close()
    output = open(fileName, "w")
    output.write(resultString)
    output.close()
  input = open(fileName, "r")
  results = parseFruitfly(input.read())
  input.close()
  return results

def findPromoters(query, name):
  """
  query:     Name of the query file.
  name:      Name of the genome.

  return:    A list of Promoter objects for the forward and reverse strands.

  This function uses the Berkeley Drosophila Genome Project website to predict promoters and parses the results into the list of Promoter objects
  that are returned.
  """
  genome = utils.loadGenome(query)
  return cachedPrediction(genome, "promoterPredictions/" + name + ".html")
