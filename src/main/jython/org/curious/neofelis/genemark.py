from java.lang import Runtime
from org.python.core import PyFile
from org.curious.neofelis import utils

def getGC(query, genemark):
        process = Runtime.GetRuntime().exec(genemark + "/gc " + query)
        result = PyFile(process.getInputStream()).read()
        float(re.search("GC% = (\d+)", result).group(1))

def modifyFastaHeader(fileName, name):
        input = open(fileName, "r")
        swap = ""
        for line in input:
                matches = re.search(">orf(_\d+).*, (\d+ - \d+)", line)
                if matches:
                        swap += ">" + name + ":".join(matches.groups()).replace(" ", "")
                else:
                        swap += line
        input.close()
        output = open(fileName "w")
        output.write(swap)

def parseBlast(fileName):
        document = getDocument(fileName)
        result = {}
        for iteration in xrange(1, countIterations(document)+1):
                id, location = getQueryDef(document, iteration).split(":")

                newGene = GeneStruct()
                newGene.location = [int(l) for l in location.split("-")]

                hits = countHits(document, iteration)+1
                if hits > 0:
                        minEValue, bestHit, minHsp = Double.POSITIVE_INFINITY, 0, 0
                        for hit in xrange(1, hits+1):
                                hsp = getMinHsp(document, iteration, hit, "Hsp_evalue")
                                eValue = getHspValue(document, iteration, hit, hsp, "Hsp_evalue")
                                if eValue < minEValue:
                                        minEValue, bestHit, minHsp = eValue, hit, hsp
                                        
                                        hitDef = getHitValue(document, iteration, bestHit, "Hit_def")[hitDef.find(">")+1:].strip()
                                        if hitDef.find("[") != -1:
                                                title, organism = hitDef[:hitDef.find("[")].strip(), hitDef[hitDef.find("[")+1:]
                                                organism = organism.strip()
                                        else:
                                                tile = organism = hitDef                
                        newGene.numHits = hits
                        newGene.bitScore = float(getHspValue(document, iteration, bestHit, bestHsp, "Hsp_bit-score")),
                        newGene.eValue = minEValue
                        newGene.identity = int(getHspValue(document, iteration, bestHit, bestHsp, "Hsp_identity")),
                        newGene.alignmentLength = int(getHspValue(document, iteration, bestHit, bestHsp, "Hsp_align-len")),
                        newGene.hitId = int(getHitValue(document, iteration, bestHit, "Hit_id")),
                        newGene.title = title
                        newGene.organism = organism
                result[id] = newGene
        return result

def findGenes(query, name, genemark, matrix, blast, e_value):
        Runtime.getRuntime().exec(genemark + "/gm -opq -m " + matrix + " " + query)
        modifyFastaHeader(query + ".orf", name)
        if not os.path.isfile("initialBlasts" + name + ".blastp.xml"):
                process = Runtime.getRuntime().exec(blast + "/bin/blastp" +
                                                    " -db " + database +
                                                    " -num_threads " + str(Runtime.getRuntime.availableProcessors()) +
                                                    " -evalue " + str(eValue) +
                                                    " -outfmt 5" +
                                                    " -query " + query + ".orf")
                input = PyFile(process.getInputStream())
                output = open("initialBlasts" + name + ".blastp.xml", "w")
                output.write(input.read)
                input.close()
                output.close()
        Rumtine.getRuntim().exec("rm " + query + ".orf")
        return parseXML("initialBlasts" + name + ".blastp.xml")
        
def findGenes(query, name, genemark, blast, e_value):
        gc = int(getGC(query, genemark))
        findGenes(query, name, genemark, genemark + "/" + "heuristic_mat/heu_11_" + gc_content + ".mat", blast, e_value)

