"""
This module use used to extend genes predicted with Genemark and then annotate them using
the function extendGenes.
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

import os
import re
import sys
import functools
from neofelis import utils

#Have to make sure that a directory to store the blasts this module creates exists.
if not os.path.isdir("extendedBlasts"):
  os.mkdir("extendedBlasts")
  
def getStops(genes):
  """
  genes:  A list of Iteration objects.

  return: A 2-tuple, first object is a list of where all the forward coding genes stop,
          second is a list of where all the reverse coding genes stop.
  """
  forwardStops = map(lambda x: x.location[1], filter(lambda x: x.location[0] < x.location[1], genes))
  reverseStops = map(lambda x: x.location[1], filter(lambda x: x.location[1] < x.location[0], genes))
  return forwardStops, reverseStops

def getExtensions(genome, genes):
  """
  genome: The genome as a string.
  genes:  A list of Iteration objects.

  return: A dictionary mapping genes(Iteration objects) to alternative locations where that gene could start.
  
  The alternate starts are calculated by starting at the original start of the gene and iterating backwards.
  When a start codon is found the start of that start codon is added to the list of alternate starts.  If this start
  codon comes before the start of the previous gene then is it still added to the list but the search terminates.
  """
  forwardStops, reverseStops = getStops(genes)
  forwardStops.append(1)
  reverseStops.append(len(genome))
  results = {}
  for gene in genes:
    results[gene] = []
    if gene.location[0] < gene.location[1]:
      bound = max(filter(lambda x: x < gene.location[1], forwardStops))
      for i in xrange(gene.location[0]-1, 0, -3):
        if genome[i-3:i] in utils.startCodons:
          results[gene].append(i-3)
          if i <= bound-1:
            break
        elif genome[i-3:i] in utils.stopCodons:
          break
    else:
      bound = min(filter(lambda x: x > gene.location[1], reverseStops))
      for i in xrange(gene.location[0], len(genome), 3):
        if utils.reverseComplement(genome[i:i+3]) in utils.startCodons:
          results[gene].append(i+3)
          if i >= bound-1:
            break
        elif utils.reverseComplement(genome[i:i+3]) in utils.stopCodons:
          break
  return results

def writeExtensions(genome, extensions):
  """
  genome: The genome as a string.
  extensions: A dictionary mapping genes(Iteration objects) to alternative locations where that gene could start.
  
  This function will write the translation of each possible extension to the file, "extensions.fas".
  """
  output = open("extensions.fas", "w")
  q = 0
  for gene, extensionList in extensions.items():
    for extension in extensionList:
      q += 1
      if gene.location[0] < gene.location[1]:
        ext = extension+1
        proteins = utils.translate(genome[extension:gene.location[1]])
      else:
        ext = extension
        proteins = utils.translate(utils.reverseComplement(genome[gene.location[1]-1:extension]))
      output.write(">" + gene.query + "~" + str(q) + ":" +
                   "-".join(map(str, [ext, gene.location[1]])) + "\n")
      for i in xrange(0, len(proteins), 50):
        output.write(proteins[i:min(i+50, len(proteins))] + "\n")
      
  output.close()

def applyExtensions(genome, genes, extendedGenes):
  """
  genome:        The genome as a string.
  genes:         A dictionary that maps query names to Iteration objects
  extendedGenes: A dictionary that maps query names to Iteration objects, extended versions of genes

  return:        A merging of genes with extendedGenes consisting of the, "better" gene in the event of a conflict
  
  The merging is done by iterating over the dictionary genes, for each entry in genes extendedGenes
  is iterated over.  If an entry in extendedGenes has a query name that starts with the query name
  of the original gene then that entry is an extension of the original gene.  This extension will replace
  the gene in the new dictionary if it either has an eValue that is lower than the original gene or the extension places
  it within 100 bps of the preceeding gene and is closer to the stop of the preceding gene.
  """
  forwardStops, reverseStops = getStops(genes.values())
  forwardStops.append(1)
  reverseStops.append(len(genome))
  
  def reduceFunction(gene, x, y):
    if re.sub(r"(~\d+)~\d+", r"\1", y.query) == gene.query:
      if gene.location[0] < gene.location[1]:
        stop = max(filter(lambda z: z < gene.location[1], forwardStops))
        gapSize = y.location[0] - stop
      else:
        stop = min(filter(lambda z: z > gene.location[1], reverseStops))
        gapSize = stop - y.location[0]
      if gapSize < 0:
        return min(x, y, key = lambda z: abs(z.location[0] - stop))
      elif gapSize < 100 or abs(x.eValue - y.eValue) < 10e-5 or utils.isNaN(x.eValue-y.eValue):
        return max(x, y, key = lambda z: abs(z.location[1] - z.location[0]))
      else:
        return min(x, y, key = lambda z: z.eValue)
    else:
      return x
        
  result = {}
  for gene, geneData in genes.items():
    result[gene] = reduce(functools.partial(reduceFunction, geneData), extendedGenes.values(), geneData)
    if result[gene] != geneData:
      result[gene].color = "0 255 0"
      result[gene].note = "Extended"
  return result

def extendGenes(query, genes, name, blast, database, eValue, pipeline):
  """
  query:    File name of the query.
  ganes:    A dictionary that maps query names to Iteration objects
  name:     Name of the genome
  blast:    Location of the installation of blast.
  database: The database to use with blast.
  eValue:   The E Value to use with blast.

  return:   A new dictionary mapping query names to Iteration objects with any better extensions replacing the originals.
  
  This function will search for any possible extensions of the genes in the query.  An extension will replace the original gene in the resulting
  dictionary if it either brings the start of the gene sufficiently close to the end of a previous gene or it has
  a lower eValue.
  """
  genome = utils.loadGenome(query)
  extensions = getExtensions(genome, genes.values())
  
  writeExtensions(genome, extensions)
  extendedGenes = utils.cachedBlast("extendedBlasts/" + name + ".blastp.xml", blast, database, eValue, "extensions.fas", pipeline)
  os.remove("extensions.fas")
  return applyExtensions(genome, genes, extendedGenes)
