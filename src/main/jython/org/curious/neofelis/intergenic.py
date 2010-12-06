import copy
from java.lang import Runtime
from org.curious.neofelis import utils

def deleteLocations(original, deletions, minLength = 3):
    result = copy.deepcopy(original)
    for deletion in deletions:
        intermediate = []
        for region in result:
            if deletion[0] < region[1] and region[1] < deletion[1]:
                intermediate += [(region[0], deletion[0])]
            elif deletion[0] < region[0] and region[0] < deletion[1]:
                intermediate += [(deletion[1], region[1])]
            elif region[0] < deletion[0] and deletion[1] < region[1]:
                intermediate += [(region[0], deletion[0]), (deletion[1], region[1])]
            else:
                intermediate += [region]
        result = intermediate
    return filter(lambda x: x[1]-x[0] > minLength, result)

def findPotentialGenes(genome, openForwardLocations):
    result = []
    for region in openForwardLocations:
        for frame in xrange(3):
            for start in xrange(region[0]-frame, region[1]-2):
                if genome[start:start+3] in startCodons:
                    for stop in xrange(start+3, region[1]-2):
                        if genome[stop:stop+3] in stopCodons:
                            result += [(start, stop+3)]

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
    genome = loadGenome(query)
    reverseComplementGenome = reverseComplement(genome)
    forwardLocations, reverseLocations = getGeneLocations(genes)
    openForwardLocations = deleteLocations((0, genome.length), forwardLocations.values(), minLength)
    openReverseLocations = deleteLocations((0, genome.length), reverseLocations.values(), minLength)

    potentialForwardGenes = findPotentialGenes(genome, openForwardLocations)
    potentialReverseGenes = findPotentialGenes(reverseComplementGenome, openReverseLocations)

    os.remove("intergenics.fas")
    writePotentialGenes(genome, potentialForwardGenes)
    writePotentialGenes(reverseComplementGenome, potentialReverseGenes)

    if not os.path.isfile("intergenicBlasts/" + name + ".blastp.xml"):
        process = Runtime.getRuntime().exec(blast + "/bin/blastp" +
                                            " -db " + database +
                                            " -num_threads " + str(Runtime.getRuntime.availableProcessors()) +
                                            " -evalue " + str(eValue) +
                                            " -outfmt 5" +
                                            " -query intergenics.fas")
        input = PyFile(process.getInputStream())
        output = open("intergenicBlasts/" + name + ".blastp.xml", "w")
        output.write(input.read)
        input.close()
        output.close()

    intergenicGenes = parseBlast("intergenicBlasts/" + name + ".blastp.xml")
    return result