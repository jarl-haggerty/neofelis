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

def writePromoters(file, promoters):
    for promoter in promoters:
        if promoter.signal10.location[0] < promoter.signal10.location[1]:
            output.write("-10_signal\t" + "..".join(map(str, promoter.signal10.location)) + "\n")
        else:
            output.write("-10_signal\tcomplement(" + "..".join(map(str, promoter.signal10.location)) + ")\n")
        output.write("\t\t/note=Promoter Position:" + str(promoter.position) + "\tIDF:" + str(promoter.IDF) + "\n")
        output.write("\t\t/colour=255 0 255\n")
        if promoter.signal35.location[0] < promoter.signal35.location[1]:
            output.write("-35_signal\t" + "..".join(map(str, promoter.signal35.location)) + "\n")
        else:
            output.write("-35_signal\tcomplement(" + "..".join(map(str, promoter.signal35.location)) + ")\n")
        output.write("\t\t/note=Promoter Position:" + str(promoter.position) + "\tIDF:" + str(promoter.IDF) + "\n")
        output.write("\t\t/colour=255 0 255\n")

def writeTerminators(output, terminators):
    for terminator in terminators:
        if terminator.location[0] < terminator.location[1]:
            output.write("terminator\t" + "..".join(map(str, terminator.location)) + "\n")
        else:
            output.write("terminator\tcomplement(" + "..".join(map(str, terminator.location)) + ")\n")
        output.write("\t\t/note=\"confidence:" + str(terminator.confidence) + "\thp_score:" + terminator.hpScore + "\ttail_score:" + terminator.tailScore + "\t" + sequence + "\n")

def writeGenes(output, genes):
    for gene in genes:
        if gene.location[0] < gene.location[1]:
            output.write("     CDS             " + str(gene.location[0]) + ".." + str(gene.location[1]) + "\n")
        else:
            output.write("     CDS             complement(" + str(gene.location[0]) + ".." + str(gene.location[1]) + ")\n")
        output.write("                     /gene=\"" + gene.title + "\"\n")
        output.write("                     /note=\"" + gene.title + "\"\n")

def writeGenome(output, genome):
    output.write("\nORIGIN\n\n")
    for i in xrange(0, len(genome), 50):
        output.write(genome[i:min(i+50, len(genome))] + "\n")

def writeArtemisFile(fileName, genome, genes=[], promoters=[], terminators=[]):
    output = open(fileName, "w")
    writePromoters(output, promoters)
    writeTerminators(output, terminators)
    writeGenes(output, genes)
    writeGenome(output, genome)

    output.close()
