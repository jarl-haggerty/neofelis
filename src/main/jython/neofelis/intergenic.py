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

import copy
from java.lang import Runtime
from neofelis import utils

def deleteGenes(genomeLength, genes, minLength = 3):
    forwardResult, reverseResult = [0, genomeLength], [0, genomeLength]
    for gene in genes:
        intermediate = []
        if gene.location[0] < gene.location[1]:
            result = forwardResult
            bottom, top = gene.location[0]-1, gene.location[1]
        else:
            result = reverseResult
            bottom, top = gene.location[1]-1, gene.location[0]
        for region in result:
            if bottom < region[1] and region[1] < top:
                intermediate += [(region[0], bottom)]
            elif bottom < region[0] and region[0] < top:
                intermediate += [(top, region[1])]
            elif region[0] < bottom and top < region[1]:
                intermediate += [(region[0], bottom), (top, region[1])]
            else:
                intermediate += [region]
        result = intermediate
    return filter(lambda x: x[1]-x[0] > minLength, forwardResult), filter(lambda x: x[1]-x[0] > minLength, reverseResult)

def findPotentialGenes(genome, openLocations):
    result = []
    stop = None
    for region in openForwardLocations:
        for frame in xrange(3):
            for start in xrange(region[1]+frame, region[0]):
                if genome[start-3:start] in startCodons:
                    stop = start
                if genome[start-3:start] in startCodons and stop:
                    result.append((start-3, stop))
                    
def writePotentialGenes(genome, genes):
    output = open("intergenics.fas", "a")
    for location in genes:
        output.write("intergenic~" +
                     "-".join([location[0]+1, location[1]]) + "\n")
        proteins = translate(genome[location[0]:location[1]])
        for i in xrange(location[0], location[1], 50):
            output.write(genome[i:min(i+50, location[1])] + "\n")
    output.close()

def findIntergenics(query, genes, name, minLength, eValue):
    genome = utils.loadGenome(query)
    reverseComplementGenome = utils.reverseComplement(genome)
    openLocations = deleteLocations(genome.length, genes, minLength)

    potentialForwardGenes = findPotentialGenes(genome, openForwardLocations)
    potentialReverseGenes = map(lambda x: (len(genome) - x[0], len(genome) - x[1]), findPotentialGenes(reverseComplementGenome, openReverseLocations))

    os.remove("intergenics.fas")
    writePotentialGenes(genome, potentialForwardGenes)
    writePotentialGenes(reverseComplementGenome, potentialReverseGenes)

    result = utils.cachedBlast("intergenicBlasts/" + name + ".blastp.xml", blast, database, eValue, "intergenics.fas")
    os.remove("intergenics.fas")
    return result
