startCodons = ("ATG", "GTG", "TTG")
stopCodons = ("TGA", "TAA", "TAG")

complementMap = {}
complementMap["A"] = "T"
complementMap["T"] = "A"
complementMap["G"] = "C"
complementMap["C"] = "G"

translationMap = {}
translationMap["TTT"] = "F"
translationMap["TTC"] = "F"
translationMap["TTA"] = "L"
translationMap["TTG"] = "L"
translationMap["CTT"] = "L"
translationMap["CTC"] = "L"
translationMap["CTA"] = "L"
translationMap["CTG"] = "L"
translationMap["TCT"] = "S"
translationMap["TCC"] = "S"
translationMap["TCA"] = "S"
translationMap["TCG"] = "S"
translationMap["TAT"] = "Y"
translationMap["TAC"] = "Y"
translationMap["TGA"] = "*"
translationMap["TAA"] = "*"
translationMap["TAG"] = "*"
translationMap["TGT"] = "C"
translationMap["TGC"] = "C"
translationMap["TGG"] = "W"
translationMap["CCA"] = "P"
translationMap["CCC"] = "P"
translationMap["CCG"] = "P"
translationMap["CCT"] = "P"
translationMap["CAC"] = "H"
translationMap["CAT"] = "H"
translationMap["CAA"] = "Q"
translationMap["CAG"] = "Q"
translationMap["CGA"] = "R"
translationMap["CGC"] = "R"
translationMap["CGG"] = "R"
translationMap["CGT"] = "R"
translationMap["AGA"] = "R"
translationMap["AGG"] = "R"
translationMap["ATT"] = "I"
translationMap["ATC"] = "I"
translationMap["ATA"] = "I"
translationMap["ATG"] = "M"
translationMap["ACA"] = "T"
translationMap["ACC"] = "T"
translationMap["ACG"] = "T"
translationMap["ACT"] = "T"
translationMap["AAT"] = "N"
translationMap["AAC"] = "N"
translationMap["AAA"] = "K"
translationMap["AAG"] = "K"
translationMap["AGT"] = "S"
translationMap["AGC"] = "S"
translationMap["GTA"] = "V"
translationMap["GTT"] = "V"
translationMap["GTG"] = "V"
translationMap["GTC"] = "V"
translationMap["GCA"] = "A"
translationMap["GCT"] = "A"
translationMap["GCG"] = "A"
translationMap["GCC"] = "A"
translationMap["GAT"] = "D"
translationMap["GAC"] = "D"
translationMap["GAA"] = "E"
translationMap["GAG"] = "E"
translationMap["GGA"] = "G"
translationMap["GGT"] = "G"
translationMap["GGG"] = "G"
translationMap["GGC"] = "G"

def translate(input):
    result = ""
    for i in xrange(len(input)-2):
        result += translationMap[input[i:i+3]]
    return result

def reverseComplement(input):
    result = map(lambda x: complementMap[x], input)
    result.reverse()
    return "".join(result)

xPath = XPatchFactory.newInstance().newXPath()
documentBuilder = DocumentBuilderFactory.newInstance().newDocumentBuilder()
iterationString = "/BlastOutput/BlastOutput_iterations/Iteration[%d]"
hitString = "/BlastOutput/BlastOutput_iterations/Iteration[%d]/Iteration_hits/Hit[%d]"
hspString = "/BlastOutput/BlastOutput_iterations/Iteration[%d]/Iteration_hits/Hit[%d]/Hit_hsps/Hsp[%d]"

def getDocument(fileName):
    return documentBuilder.parse(fileName)

def countIterations(document):
    return xPath.evaluate(iterationString[:-4], document, XPathConstants.NODESET).getLength()

def countIterationHits(document, iteration):
    return xPath.evaluate(hitString[:-4] % (iteration), document, XPathConstants.NODESET).getLength()

def countHsps(document, iteration, hit):
    return xPath.evaluate(hspString[:-4] % (iteration, hit), document, XPathConstants.NODESET).getLength()

def getIterationValue(document, iteration, value):
    return xPath.evaluate(iterationString % (iteration) + "/" + value, document)

def getHitValue(document, iteration, hit, value):
    return xPath.evaluate(hitString % (iteration, hit) + "/" + value, document)

def getHspValue(document, iteration, hit, hsp, value):
    return xPath.evaluate(hspString % (iterations, hit, hsp) + "/" + value, document)

def getMinHsp(document, iteration, hit, measure):
    minValue, minHsp = Double.POSITIVE_INFINITY, 0
    for hsp in xrange(countHsps(document, iteration, hit)):
        value = float(xPath.evaluate(hspString % (iteration, hit, hsp) + "/" + measure, document))
        if value < minValue:
            minValue, minHsp = value, hsp
    return minHsp

def getMaxHsp(document, iteration, hit, measure):
    maxValue, maxHsp = Double.NEGATIVE_INFINITY, 0
    for hsp in xrange(countHsps(document, iteration, hit)):
        value = float(xPath.evaluate(hspString % (iteration, hit, hsp) + "/" + measure, document))
        if value > maxValue:
            maxValue, maxHsp = value, hsp
    return maxHsp

class GeneStruct:
        location =        []
        numHits =         0
        bitScore =        0
        eValue =          Doube.POSITIVE_INFINITY
        identity =        0
        alignmentLength = 0
        hitId =           -1
        title =           "None"
        organism =        "None"

def loadGenome(fileName):
        input = open(fileName, "r")
        result = ""
        for line in input:
                match = re.match("([ACGT]+)\n", line)
                if match:
                        result += match.group(1)
        return result

def getGeneLocations(genes):
        forward = {}
        reverse = {}
        for k, v in genes.items():
                if v.location[0] < v.location[1]:
                        forward[k] = [v.location[0]-1, v.location[1]]
                else:
                        reverse[k] = [v.location[1]-1, v.location[0]]
        return forward, reverse
