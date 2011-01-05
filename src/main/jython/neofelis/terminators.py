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

def writeCoords(name, genes):
  output = open(name + ".crd", "w")
  for gene in genes:
    output.write("gene\t%d\t%d\t%s\n" %  (gene.location[0], gene.location[1], name))
  output.close()

def findTerminators(query, name, genes, transterm):
  writeCoords(name, genes)
  output = open("terms", "w")
  subprocess.Popen([transterm + "/transterm", "-p", transterm + "/expterm.dat", query, name + ".crd"], stdout=output).wait()
  output.close()
