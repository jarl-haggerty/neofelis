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
from javax.swing import JFrame, JLabel, JButton, WindowConstants
from java.awt import GridBagLayout, GridBagConstraints

#Make sure the directory for BPROM results exists.
if not os.path.isdir("promoterPredictions"):
  os.mkdir("promoterPredictions")

class Promoter():
  """
  Class for storing information about a promoter.
  """
  def __init__(self):
    self.score = 0
    self.position = None
    self.location = None
    self.signal10Location = None
    self.signal35Location = None

  def __str__(self):
    result = "<"
    result += "Score = " + str(self.score) + ", "
    result += "Position = " + str(self.position) + ", "
    result += "Signal10Location = " + str(self.signal10Location) + ", "
    result += "Signal35Location = " + str(self.signal35Location)
    result += ">"
    return result

def parseBPROM(bpromData):
  """
  bpromData: Contents of the body of a BPROM prediction.
  
  return: A list of Promoter objects representing the results of the BPROM prediction
  """
  promoters = re.findall(r"\s+Promoter\s+Pos:\s+(\d+)\s+LDF-\s+(\d+\.?\d*)", bpromData)
  tenBoxes = re.findall(r"\s+-10\s+box\s+at\s+pos.\s+(\d+)\s+([ACGT]+)\s+Score\s+-?\d+", bpromData)
  thirtyFiveBoxes = re.findall(r"\s+-35\s+box\s+at\s+pos.\s+(\d+)\s+([ACGT]+)\s+Score\s+-?\d+", bpromData)
  def parsePromoter(promoter, tenBox, thirtyFiveBox):
    result = Promoter()
    result.score = float(promoter[1])
    result.position = int(promoter[0])
    result.signal10Location = [int(tenBox[0]), int(tenBox[0])+len(tenBox[1])]
    result.signal35Location = [int(thirtyFiveBox[0]), int(thirtyFiveBox[0])+len(thirtyFiveBox[1])]
    result.location = [min(result.signal10Location + result.signal35Location), max(result.signal10Location + result.signal35Location)]
    return result
  return map(parsePromoter, promoters, tenBoxes, thirtyFiveBoxes)

def cachedBPROM(genome, fileName, frame):
  """
  genome: Genome as a string.
  fileName: File to save the BPROM results in.
  swing: A JFrame or None.  If this is None then messages will be printed, if it isn't then
         they will also be put in a dialog box.

  return: Results of the BPROM prediction stored in a list of Promoter objects.

  If the file Specified by fileName already exists then this function simply parses the file
  already there.  Also, if a request is made to BPROM and nothing is returned, no file is
  created, the user is warned, and an empty list is returned.
  """
  offset = 25 if ".forward.bprom" in fileName else 50
  direction = "forward" if offset == 50 else "reverse"
    
  def getPromoters():
    input = open(fileName, "r")
    results = parseBPROM(input.read())
    input.close()
    return results

  if not os.path.isfile(fileName):
    results = urllib.urlopen("http://linux1.softberry.com/cgi-bin/programs/gfindb/bprom.pl",
                             urllib.urlencode({"DATA" : genome}))
    resultString = results.read()
    results.close()
    resultString = resultString[resultString.find("<pre>"):resultString.find("</pre>")]
    resultString = re.sub("<+.+>+", "", resultString).strip()
    if resultString:
      output = open(fileName, "w")
      output.write(resultString)
      output.close()
      return getPromoters()
    else:
      if frame:
        messageFrame = JFrame("BPROM Error",
                              defaultCloseOperation = WindowConstants.DISPOSE_ON_CLOSE)
        messageFrame.setLocation(frame.location().x + offset, frame.location().y + offset)
        messageFrame.contentPane.layout = GridBagLayout()
        constraints = GridBagConstraints()
        constraints.gridx, constraints.gridy = 0, 0
        constraints.gridwidth, constraints.gridheight = 1, 1
        constraints.fill = GridBagConstraints.BOTH
        constraints.weightx, constraints.weighty = 1, 1
        messageFrame.contentPane.add(JLabel("<html>The pipeline will continue to run but BPROM<br/>did not process the request for promoters on the<br/>" + direction + " strand.  Try again tomorrow.</html>"), constraints)
        constraints.gridx, constraints.gridy = 0, 1
        constraints.fill = GridBagConstraints.NONE
        constraints.weightx, constraints.weighty = 1, 1
        constraints.anchor = GridBagConstraints.LINE_END
        messageFrame.contentPane.add(JButton("Ok", actionPerformed = lambda e: messageFrame.dispose()), constraints)
        messageFrame.pack()
        messageFrame.visible = True
      print "BPROM Error:", "The pipeline will continue to run but BPROM did not process the request for promoters on the " + direction + " strand.  Try again tomorrow"
      return []
  else:
    return getPromoters()

def reverseCoordinates(genomeLength, promoter):
  """
  genomeLength: Length of the genome.
  promoter: Promoter object.

  return: Original promoter but with it's coordinates altered so it appears to be on the reverse strand.
  """
  newPromoter = Promoter()
  newPromoter.score = promoter.score
  newPromoter.position = genomeLength+1 - promoter.position
  newPromoter.location = map(lambda x: genomeLength+1 - x, promoter.location)
  newPromoter.signal10Location = map(lambda x: genomeLength+1 - x, promoter.signal10Location)
  newPromoter.signal35Location = map(lambda x: genomeLength+1 - x, promoter.signal35Location)
  return newPromoter

def findPromoters(query, name, scoreCutoff, frame):
  """
  query: Name of the query file.
  name: Name of the genome.
  scoreCutoff: Minimum promoter score value for any promoters.
  frame: A JFrame that may be used as the parent for a JDialog to display messages.  If it is none then messages
         are just printed.

  return: A list of Promoter objects for the forward and reverse strands.
  
  This function uses BPROM to predict promoters and parses the results into the list of Promoter objects
  that are returned. Promoters with a score lower than scoreCutoff are filtered out.
  """
  genome = utils.loadGenome(query)
  forwardResults = cachedBPROM(genome, "promoterPredictions/" + name + ".forward.bprom", frame)
  reverseResults = cachedBPROM(utils.reverseComplement(genome), "promoterPredictions/" + name + ".reverse.bprom", frame)
  reverseResults = map(functools.partial(reverseCoordinates, len(genome)), reverseResults)
  return filter(lambda x: x.score > scoreCutoff, forwardResults + reverseResults)
