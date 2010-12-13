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

import os
import sys
from neofelis import utils

"""Have to make sure that a directory to store the blasts this module creates exists."""
if not os.path.isdir("extendedBlasts"):
  os.mkdir("extendedBlasts")
  
def getStops(genes):
  return map(lambda x: x.location[1], filter(lambda x: x.location[0] < x.location[1], genes)), map(lambda x: x.location[0], filter(lambda x: x.location[1] < x.location[0], genes))

def getExtensions(genome, genes):
  forwardStops, reverseStops = getStops(genes)
  forwardStops.append(1)
  reverseStops.append(len(genome))
  results = {}
  for gene in genes:
    results[gene] = []
    if gene.location[0] < gene.location[1]:
      bound = max(filter(lambda x: x < gene.location[0], forwardStops))
      for i in xrange(gene.location[0]-1, bound-1, -3):
        if genome[i-3:i] in utils.startCodons:
          results[gene].append(i-3)
        elif genome[i-3:i] in utils.stopCodons:
          break
    else:
      bound = max(filter(lambda x: x > gene.location[1], reverseStops))
      for i in xrange(gene.location[0], bound-1, 3):
        if utils.reverseComplement(genome[i:i+3]) in utils.startCodons:
          results[gene].append(i+3)
        elif utils.reverseComplement(genome[i:i+3]) in utils.stopCodons:
          break
  return results

def writeExtensions(genome, extensions):
  output = open("extensions.fas", "w")
  for gene, extensionList in extensions.items():
    for extension in extensionList:
      if gene.location[0] < gene.location[1]:
        output.write(">" + gene.query + ":" +
                     "-".join(map(str, [extension, gene.location[1]])) + "\n")
        proteins = utils.translate(genome[extension:gene.location[1]])
        for i in xrange(0, len(proteins), 50):
          output.write(proteins[i:min(i+50, len(proteins))] + "\n")
      else:
        output.write(">" + gene.query + ":" +
                     "-".join(map(str, [gene.location[0], extension])) + "\n")
        proteins = utils.translate(utils.reverseComplement(genome[gene.location[1]-1:extension]))
        for i in xrange(0, len(proteins), 50):
          output.write(proteins[i:min(i+50, len(proteins))] + "\n")
  output.close()

def applyExtensions(genes, extendedGenes):
  forwardStops, reverseStops = getStops()
  def forwardReduce(gene, x, y):
    if y.find(gene) == 0:
      if gene.location[0] < gene.location[1]:
        gapSize = y.location[0] - max(filter(lambda z: z < gene.location[0], forwardStops))
      else:
        gapSize = min(filter(lambda z: z > gene.location[0], reverseStops)) - y.location[0]
      if gapSize < 100:
        return max(x, y, key = lambda z: abs(z.location[1] - z.location[0]))
      else:
        return min(x, y, key = lambda z: z.eValue)
    else:
      return x
        
  result = {}
  for gene, geneData in genes.items():
    result[gene] = reduce(reduceFunction, extendedGenes.items(), geneData)
  return result

def extendGenes(query, genes, name, blast, database, eValue):
        genome = utils.loadGenome(query)
        extensions = getExtensions(genome, genes.values())
        writeExtensions(genome, extensions)
        print "extending"
        extendedGenes = utils.cachedBlast("extendedBlasts/" + name + ".blastp.xml", blast, database, eValue, "extensions.fas")
        sys.exit(0)
        return applyExtensions(genes, extendedGenes)
