"""
This module contains classes and functions for predicting promoters in a genome.
"""

import subprocess
from neofelis import utils
import re
import os

class Terminator():
  """
  Class used to represent a terminator.
  """
  def __init__(self):
    self.location = None
    self.confidence = None
    self.hpScore = None
    self.tailScore = None

def writeCoords(fileName, name, genes):
  """
  name:  Name of the genome.
  ganes: List of Iteration objects.

  Writes the genes into a .crd file.
  """
  output = open(fileName + ".crd", "w")
  for gene in genes:
    output.write("gene\t%d\t%d\t%s\n" %  (gene.location[0], gene.location[1], name))
  output.close()

def parseTransterm(input):
  """
  input:  Contents of a transterm file.

  return: Contents of the file parsed into a list of Terminator objects.
  """
  matches = re.findall(r"\s+TERM\s+\d+\s+(\d+)\s+-\s+(\d+)\s+[+-]\s+\S+\s+(\d+)\s+(-?\d+.?\d+)\s+(-?\d+.?\d+)", input)
  def buildTerminator(pieces):
    """
    Turns a list of strings into a Terminator object.
    """
    result = Terminator()
    result.location = [int(pieces[0]), int(pieces[1])]
    result.confidence = int(pieces[2])
    result.hpScore = float(pieces[3])
    result.tailScore = float(pieces[4])
    return result

  return map(buildTerminator, matches)

def findTerminators(query, name, genes, transterm):
  """
  query:     File name of the query.
  name:      Name of the genome.
  genes:     List of Iteration objects.
  transterm: Location of the transterm installation.

  return:    A list of Terminator objects.

  This function runs transterm with the query and the genes and parses the results into the return value.
  """
  fileName = os.path.splitext(os.path.split(query)[1])[0]
  writeCoords(fileName, name, genes)
  output = ""
  transtermProcess = subprocess.Popen([transterm + "/transterm", "-p", transterm + "/expterm.dat", query, fileName + ".crd"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  
  while transtermProcess.poll() == None:
    output += transtermProcess.stdout.read()
    transtermProcess.stderr.read()
    
  remaining = transtermProcess.stdout.read()
  
  while remaining:
    output += remaining
    remaining = transtermProcess.stdout.read()

  result = parseTransterm(output)
  os.remove(fileName + ".crd")
  return result
