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
from java.lang import Runtime
from neofelis import utils

def getStops(genes):
  return map(lambda x: x.location[1]-2, filter(lambda x: x.location[0] < x.location[1], genes)), map(lambda x: x.location[0]-2, filter(lambda x: x.location[1] < x.location[0], genes))

def getExtensions(genome, genes):
  forwardStops, reverseStops = getStops(genes)
  result = {}
  for gene, geneData in locations.items():
    results[gene] = []
    if geneData.location[0] < geneData.location[1]:
      bound = max(filter(lambda x: x < geneData.location[0], forwardStops))
      for i in xrange(geneData.location[0] - 3, bound, -3):
        if genome[i:i+3] in startCodons:
          result[geneData].append(i)
        elif genome[i:i+3] in stopCodons:
          break
    else:
      bound = max(filter(lambda x: x > geneData.location[1], reverseStops))
      for i in xrange(geneData.location[1], bound, 3):
        if reverseComplement(genome[i:i+3]) in startCodons:
          result[geneData].append(i)
        elif reverseComplement(genome[i:i+3]) in stopCodons:
          break
  return result

def writeExtensions(genome, extensions):
  output = open("extensions.fas", "w")
  for gene, extension in extension.items():
    if gene.location[0] < gene.location[1]:
      output.write(gene.query + "~" +
                   "-".join([gene.location[0], gene.location[1]]) + "\n")
      proteins = translate(genome[gene.location[0]-1:gene.location[1]])
      for i in xrange(0, len(proteins), 50):
        output.write(proteins[i:min(i+50, len(proteins))] + "\n")
    else:
      output.write(gene.query + "~" +
                   "-".join([gene.location[0], gene.location[1]]) + "\n")
      proteins = translate(reverseComplement(genome[gene.location[1]-1:gene.location[0]]))
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
        return max(x, y, key = lambda z: z.location[1] - z.location[0])
      else:
        return min(x, y, key = lambda z: z.eValue)
    else:
      return x
        
  result = {}
  for gene, geneData in genes.items():
    result[gene] = reduce(reduceFunction, extendedGenes.items(), geneData)
  return result

def extendGenes(query, genes, name, blast, database, eValue):
        genome = loadGenome(query)
        extensions = getExtensions(genome, locations)
        writeExtensions(genome, extensions)
        extendedGenes = cachedBlast("extendedBlasts/" + name + ".blastp.xml")
        return applyExtensions(genes, extendedGenes)
