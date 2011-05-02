"""
This module contains functions for writing artemis files.
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

def writePromoters(output, promoters):
    """
    output:    File object to write to.
    promoters: List of Promoter objects.

    Writes promoters to the file.
    """
    for promoter in promoters:
        if promoter.location[0] < promoter.location[1]:
            output.write("     promoter               " + "..".join(map(str, promoter.location)) + "\n")
        else:
            output.write("     promoter               complement(" + "..".join(map(str, promoter.location)) + ")\n")
        output.write("                            /note=Promoter Position:" + str(promoter.location[0]) + "\tscore:" + str(promoter.score) + "\n")
        output.write("                            /colour=255 0 255\n")

def writeTerminators(output, terminators):
    """
    output:      File object to write to.
    terminators: List of Terminator objects.

    Writes terminators to the file.
    """
    for terminator in terminators:
        if terminator.location[0] < terminator.location[1]:
            output.write("     terminator             " + "..".join(map(str, terminator.location)) + "\n")
        else:
            output.write("     terminator             complement(" + "..".join(map(str, terminator.location)) + ")\n")
        output.write("                            /note=\"confidence:" + str(terminator.confidence) + "\thp_score:" + str(terminator.hpScore) + "\ttail_score:" + str(terminator.tailScore) + "\"\n")

def writeTransferRNAs(output, transferRNAs):
    """
    output:      File object to write to.
    terminators: List of Terminator objects.

    Writes terminators to the file.
    """
    for transferRNA in transferRNAs:
        if transferRNA.location[0] < transferRNA.location[1]:
            output.write("     gene                   " + "..".join(map(str, transferRNA.location)) + "\n")
        else:
            output.write("     gene                   complement(" + "..".join(map(str, transferRNA.location)) + ")\n")
        output.write("                            /note=\"type:" + str(transferRNA.type) + "\tanti codon:" + str(transferRNA.antiCodon) + "\tcove_score:" + str(transferRNA.coveScore) + "\"\n")

def writeGenes(output, genes):
    """
    output: File object to write to.
    genes:  List of Iteration objects.

    Writes genes to the file.
    """
    for gene in genes:
        if gene.location[0] < gene.location[1]:
            output.write("     CDS             " + str(gene.location[0]) + ".." + str(gene.location[1]) + "\n")
        else:
            output.write("     CDS             complement(" + str(gene.location[0]) + ".." + str(gene.location[1]) + ")\n")
        output.write("                     /gene=\"" + gene.title + "\"\n")
        output.write("                     /note=\"" + gene.note + "\"\n")
        output.write("                     /colour=" + gene.color + "\n")

def writeGenome(output, genome):
    """
    output: File object to write to.
    genome: Genome as a string.

    Writes genome to the file.
    """
    output.write("\nORIGIN\n\n")
    for i in xrange(0, len(genome), 50):
        output.write(genome[i:min(i+50, len(genome))] + "\n")

def writeArtemisFile(fileName, genome, genes=[], promoters=[], terminators=[], transferRNAs=[]):
    """
    fileName:    Name of the artemis file to write.
    genome:      Genome as a string.
    genes:       List of Iteration objects.
    promoters:   List of Promoter objects.
    terminators: List of Terminator objects.

    Writes an artemis file consisting of the genome, genes, promoters, and terminators.
    """
    output = open(fileName, "w")
    writePromoters(output, promoters)
    writeTerminators(output, terminators)
    writeTransferRNAs(output, transferRNAs)
    writeGenes(output, genes)
    writeGenome(output, genome)

    output.close()
