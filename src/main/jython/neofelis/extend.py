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

def getForwardExtensions(genome, locations, stops):
        result = {}
        for gene, location in locations.items():
                results[gene] = []
                lowerBound = max(filter(lambda x: x < location[0], stops))
                for i in xrange(location[0] - 3, lowerBound, -3):
                        if genome[i:i+3] in startCodons:
                                result[gene].append(i)
                        elif genome[i:i+3] in stopCodons:
                                break
        return result

def getReverseExtensions(genome, locations, stops):
        result = {}
        for gene, location in locations.items():
                results[gene] = []
                upperBound = max(filter(lambda x: x > location[1], stops))
                for i in xrange(location[1], upperBound, 3):
                        if reverseComplement(genome[i:i+3]) in startCodons:
                                result[gene].append(i)
                        elif reverseComplement(genome[i:i+3]) in stopCodons:
                                break
        return result

def writeForwardExtensions(genome, extensions):
        output = open("extensions.fas", "a")
        for gene, location in extension.items():
                output.write(gene + "~" +
                             "-".join([location[0]+1, location[1]]) + "\n")
                proteins = translate(genome[location[0]:location[1]])
                for i in xrange(location[0], location[1], 50):
                        output.write(genome[i:min(i+50, location[1])] + "\n")
        output.close()

def writeReverseExtensions(genome, extensions):
        output = open("extensions.fas", "a")
        for gene, location in extension.items():
                output.write(gene + "~" +
                             "-".join([location[1], location[0]+1]) + "\n")
                proteins = translate(reverseComplement(genome[location[0]:location[1]]))
                for i in xrange(location[0], location[1], 50):
                        output.write(genome[i:min(i+50, location[1])] + "\n")
        output.close()

def applyExtensions(genes, extendedGenes):
        def reduceFunction(x, y):
                if y[0].find(gene) == 0 and y[1].eValue < x.eValue: 
                        min(x, y[1], key = lambda x: x.eValue)
        result = {}
        for gene, geneData in genes.items():
                bestExtension = reduce(reduceFunction, extendedGenes.items(), GeneStruct())
                if bestExtension.eValue < geneData.eValue:
                        result[gene] = bestExtension
                else:
                        result[gene] = geneData
        return result

def extendGenes(query, genes, name, blast, database, eValue):
        genome = loadGenome(query)
        forwardLocations, reverseLocations = getGenelocations(genes)
        forwardStops = map(lambda x: x[1]-2, forwardLocations.values())
        reverseStops = map(lambda x: x[0]-2, reverseLocations.values())

        forwardExtensions = getForwardExtensions(genome, forwardLocations, forwardStops)
        reverseExtensions = getReverseExtensions(genome, reverseLocations, reverseStops)

        os.remove("extensions.fas")
        writeForwardExtensions(genome, forwardLocations, forwardExtensions)
        writeReverseExtensions(genome, reverseLocations, reverseExtensions)

        if not os.path.isfile("extendedBlasts/" + name + ".blastp.xml"):
                process = Runtime.getRuntime().exec(blast + "/bin/blastp" +
                                                    " -db " + database +
                                                    " -num_threads " + str(Runtime.getRuntime.availableProcessors()) +
                                                    " -evalue " + str(eValue) +
                                                    " -outfmt 5" +
                                                    " -query extensions.fas")
                input = PyFile(process.getInputStream())
                output = open("extendedBlasts/" + name + ".blastp.xml", "w")
                output.write(input.read)
                input.close()
                output.close()
        
        extendedGenes = parseBlast("extendedBlasts/" + name + ".blastp.xml")
        return applyExtensions(genes, extendedGenes)
