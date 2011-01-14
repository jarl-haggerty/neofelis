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

import copy
import sys
import os
from neofelis import utils

"""Have to make sure that a directory to store the blasts this module creates exists."""
if not os.path.isdir("intergenicBlasts"):
  os.mkdir("intergenicBlasts")

class Region():
  def __init__(self, startGene, start, stop, stopGene):
    self.startGene = startGene
    self.start = start
    self.stop = stop
    self.stopGene = stopGene

  def __str__(self):
    result = ">"
    result += "StartGene = " + str(self.startGene) + ", "
    result += "Start = " + str(self.start) + ", "
    result += "StopGene = " + str(self.stopGene) + ", "
    result += "Stop = " + str(self.stop) + ">"
    return result
  
def deleteGenes(genomeLength, genes, minLength = 3):
    def filterFunction(region):
        inset = abs(region.stopGene.location[1] - region.stopGene.location[0])/2 if region.stopGene else 0
        return region.stop + inset - region.start > minLength

    forwardResult, reverseResult = [Region(None, 0, genomeLength, None)], [Region(None, 0, genomeLength, None)]
    for gene in genes:
        intermediate = []
        if gene.location[0] < gene.location[1]:
            result = forwardResult
            bottom, top = gene.location[0], gene.location[1]
        else:
            result = reverseResult
            bottom, top = genomeLength-gene.location[0]+1, genomeLength-gene.location[1]+1
        for region in result:
            if bottom < region.stop and region.stop <= top:
                intermediate += [Region(region.startGene, region.start, bottom, gene)]
            elif bottom <= region.start and region.start < top:
                intermediate += [Region(gene, top, region.stop, region.stopGene)]
            elif region.start < bottom and top < region.stop:
                intermediate += [Region(region.startGene, region.start, bottom, gene), Region(gene, top, region.stop, region.stopGene)]
            else:
                intermediate += [region]
        if gene.location[0] < gene.location[1]:
            forwardResult = filter(filterFunction, intermediate)
        else:
            reverseResult = filter(filterFunction, intermediate)
    
    forwardResult = filter(filterFunction, forwardResult)
    reverseResult = filter(filterFunction, reverseResult)
    return forwardResult, reverseResult

def findPotentialGenes(genome, openLocations, minLength = 3):
    result = []
    stop = None
    for region in openLocations:
        debug = False
        if (region.start, region.stop) == (79616, 79693):
            debug = True
        inset = abs(region.stopGene.location[1] - region.stopGene.location[0])/2 if region.stopGene else 0
        if debug:
            print "Debug"
            print region.stop, inset
            print region.stopGene
        for frame in xrange(3):
            if debug: print ""
            for start in xrange(region.stop+frame+inset, 0, -3):
                if debug:
                  print start - 3, genome[start-3:start]
                if genome[start-3:start] in utils.stopCodons:
                    stop = start
                    if debug: print "stop", stop
                    if start <= region.start:
                        break
                if genome[start-3:start] in utils.startCodons and stop and start < region.stop:
                    if debug: print "start", start-3
                    result.append((start-3, stop))
                    if start <= region.start:
                        break
    return filter(lambda x: x[1]-x[0] > minLength, result)
                    
def writePotentialGenes(genome, locations):
    output = open("intergenics.fas", "w")
    q = 0
    for location in locations:
        q += 1
        if location[0] < location[1]:
            output.write(">intergenic~" + str(q) + ":" + str(location[0]+1) + "-" + str(location[1]) + "\n")
            proteins = utils.translate(genome[location[0]:location[1]])
        else:
            output.write(">intergenic~" + str(q) + ":" + str(location[0]) + "-" + str(location[1]+1) + "\n")
            proteins = utils.translate(utils.reverseComplement(genome[location[1]:location[0]]))
        for i in xrange(0, len(proteins), 50):
            output.write(proteins[i:min(i+50, len(proteins))] + "\n")
    output.close()

def removeCommonStops(genes):
    stopDictionary = {}
    for gene in genes.values():
        if gene.location[1] in stopDictionary:
            stopDictionary[gene.location[1]].append(gene)
        else:
            stopDictionary[gene.location[1]] = [gene]

    def reduceFunction(x, y):
        if utils.isNaN(x.eValue - y.eValue):
            return max(x, y, key = lambda z: abs(z.location[1] - z.location[0]))
        else:
            return min(x, y, key = lambda z: z.eValue)            
    result = {}
    for stoppedGenes in stopDictionary.values():
        temp = reduce(reduceFunction, stoppedGenes)
        result[temp.query] = temp
    return result

def findIntergenics(query, genes, name, minLength, blast, database, eValue, remote):
    genome = utils.loadGenome(query)
    reverseComplementGenome = utils.reverseComplement(genome)
    openForwardLocations, openReverseLocations = deleteGenes(len(genome), genes.values(), minLength)

    potentialGenes = findPotentialGenes(genome, openForwardLocations, minLength)
    reversePotentialGenes = findPotentialGenes(reverseComplementGenome, openReverseLocations, minLength)
    potentialGenes += map(lambda x: (len(genome)-x[0], len(genome)-x[1]), reversePotentialGenes)

    writePotentialGenes(genome, potentialGenes)

    result = utils.cachedBlast("intergenicBlasts/" + name + ".blastp.xml", blast, database, eValue, "intergenics.fas", remote)
    #os.remove("intergenics.fas")
    result = removeCommonStops(result)
    for r in result.values():
        r.intergenic = True
        r.note = "Intergenic"
        r.color = "160 32 240"
    return result
