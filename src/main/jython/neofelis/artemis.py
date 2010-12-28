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

def writePromoters(file, promoters):
    for promoter in promoters:
        if promoter.signal10.location[0] < promoter.signal10.location[1]:
            output.write("-10_signal\t" + "..".join(map(str, promoter.signal10.location)) + "\n")
        else:
            output.write("-10_signal\tcomplement(" + "..".join(map(str, promoter.signal10.location)) + ")\n")
        output.write("\t\t/note=Promoter Position:" + str(promoter.position) + "\tIDF:" + str(promoter.LDF) + "\n")
        output.write("\t\t/colour=255 0 255\n")
        if promoter.signal35.location[0] < promoter.signal35.location[1]:
            output.write("-35_signal\t" + "..".join(map(str, promoter.signal35.location)) + "\n")
        else:
            output.write("-35_signal\tcomplement(" + "..".join(map(str, promoter.signal35.location)) + ")\n")
        output.write("\t\t/note=Promoter Position:" + str(promoter.position) + "\tIDF:" + str(promoter.LDF) + "\n")
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
        output.write("                     /note=\"" + gene.note + "\"\n")
        output.write("                     /colour=" + gene.color + "\n")

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
