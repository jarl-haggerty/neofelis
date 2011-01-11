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

"""Have to make sure that a directory to store the blasts this module creates exists."""
if not os.path.isdir("extendedBlasts"):
  os.mkdir("extendedBlasts")
  
def getStops(genes):
  """
  Given a list of genes this function returns two lists, one for the forward strand and one for the
  reverse strand, consisting of number marking where all the genes stop encoding.
  """
  forwardStops = map(lambda x: x.location[1], filter(lambda x: x.location[0] < x.location[1], genes))
  reverseStops = map(lambda x: x.location[1], filter(lambda x: x.location[1] < x.location[0], genes))
  return forwardStops, reverseStops

def getExtensions(genome, genes):
  """
  Given a genome and a list of genes in that genome this function returns a map.  The keys of the map
  are the genes, and the values are lists of all the possible alternative starts in the genome for that
  gene.
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
  Given a genome and a dictionary mapping genes to lists of their alternative starts(extensions)
  this function will write the translation of each possible extension to the file extensions.fas.
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
  This function takes a genome, and two dictionaries of genes.  The result of this function
  is essentially a merging of the two dictionaries of genomes.  The merging is done by iterating
  over the dictionary genes, for each entry in genes extendedGenes is iterated over.  If an entry
  in extendedGenes has a query name that starts with the query name of the original gene then that
  entry is an extension of the original gene.  This extension will replace the gene in the new
  dictionary if it either has an eValue that is lower than the original gene or the extension places
  it within 100 bps of the preceeding gene.
  """
  forwardStops, reverseStops = getStops(genes.values())
  forwardStops.append(1)
  reverseStops.append(len(genome))
  def reduceFunction(gene, x, y):
    if re.sub(r"(~\d+)~\d+", r"\1", y.query) == gene.query:
      if gene.location[0] < gene.location[1]:
        gapSize = abs(y.location[0] - max(filter(lambda z: z < gene.location[1], forwardStops)))
      else:
        gapSize = abs(min(filter(lambda z: z > gene.location[1], reverseStops)) - y.location[0])
      if gapSize < 100 or abs(x.eValue - y.eValue) < 10e-5 or utils.isNaN(x.eValue-y.eValue):
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

def extendGenes(query, genes, name, blast, database, eValue, remote):
  """
  This function will search for any possible extensions of the genes in the fasta file query.  These extensions will
  be blasted using the blast installation located at blast with the database database and the eValue.  If remote is
  true then the blast search is run remotely.  An extension will replace the original gene in the resulting
  dictionary if it either brings the start of the gene sufficiently close to the end of a previous gene or it has
  a lower eValue.
  """
  genome = utils.loadGenome(query)
  extensions = getExtensions(genome, genes.values())
  
  writeExtensions(genome, extensions)
  extendedGenes = utils.cachedBlast("extendedBlasts/" + name + ".blastp.xml", blast, database, eValue, "extensions.fas", remote)
  os.remove("extensions.fas")
  return applyExtensions(genome, genes, extendedGenes)
