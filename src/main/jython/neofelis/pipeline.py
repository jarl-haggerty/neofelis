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
import sys
from neofelis import genemark
from neofelis import extend
from neofelis import intergenic
from neofelis import promoters
from neofelis import terminators
from neofelis import artemis
from neofelis import utils
from neofelis import scaffolds

def run(blastLocation, genemarkLocation, database, eValue, matrix, minLength, scaffoldingDistance, remote, queries):
  for query in queries:
    name = os.path.splitext(query)[0]
    name = os.path.split(name)[1]
    genome = utils.loadGenome(query)
    
    initialGenes = genemark.findGenes(query, name, blastLocation, database, eValue, genemarkLocation, matrix, remote)
    artemis.writeArtemisFile(name + ".art", genome, initialGenes.values())
    extendedGenes = extend.extendGenes(query, initialGenes, name, blastLocation, database, eValue, remote)
    artemis.writeArtemisFile(name + "extended.art", genome, extendedGenes.values())
    intergenicGenes = intergenic.findIntergenics(query, extendedGenes, name, minLength, blastLocation, database, eValue, remote)
    genes = {}
    for k, v in extendedGenes.items() + intergenicGenes.items():
      genes[k] = v
    artemis.writeArtemisFile(name + "intergenic.art", genome, genes.values())
    scaffolded = scaffolds.refineScaffolds(genes, scaffoldingDistance)
    artemis.writeArtemisFile(name + "scaffolds.art", genome, scaffolded.values())

    initialPromoters = promoters.findPromoters(query)
    initialTerminators = terminators.findTerminators(query, transterm)
    
                
    writeArtemisFile(genes, promoters, terminators)
