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

import subprocess
from neofelis import utils
import re

class Terminator():
  def __init__(self):
    self.location = None
    self.confidence = None
    self.hpScore = None
    self.tailScore = None

def writeCoords(name, genes):
  output = open(name + ".crd", "w")
  for gene in genes:
    output.write("gene\t%d\t%d\t%s\n" %  (gene.location[0], gene.location[1], name))
  output.close()

def parseTransterm(input):
  matches = re.findall(r"\s+TERM\s+\d+\s+(\d+)\s+-\s+(\d+)\s+[+-]\s+\S+\s+(\d+)\s+(-?\d+.?\d+)\s+(-?\d+.?\d+)", input)
  def buildTerminator(pieces):
    result = Terminator()
    result.location = [int(pieces[0]), int(pieces[1])]
    result.confidence = int(pieces[2])
    result.hpScore = float(pieces[3])
    result.tailScore = float(pieces[4])
    return result

  return map(buildTerminator, matches)

def findTerminators(query, name, genes, transterm):
  writeCoords(name, genes)
  output = open(name + ".tt", "w")
  subprocess.Popen([transterm + "/transterm", "-p", transterm + "/expterm.dat", query, name + ".crd"], stdout=output).wait()
  output.close()
  input = open(name + ".tt", "r")
  result = parseTransterm(input.read())
  input.close()
  return result
