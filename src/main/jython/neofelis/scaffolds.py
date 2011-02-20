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

import sys
import copy

class Scaffold():
    """
    This structure contains the information that describes a scaffold,
    where it starts, where is stops, and the genes contained in it.
    """
    def __init__(self, start, stop, genes):
        self.start = start
        self.stop = stop
        self.genes = genes
        
    def __str__(self):
        result = "<"
        result += "start = " + str(self.start) + ", "
        result += "stop = " + str(self.stop) + ", "
        result += "genes = " + str(self.genes) + ">"
        return result

def extractScaffolds(genes, scaffoldingDistance = 100):
    forwardScaffolds, reverseScaffolds = [], []
    for gene in genes:
        if gene.location[0] < gene.location[1]:
            newScaffold = Scaffold(gene.location[0], gene.location[1], [gene])
            scaffolds = forwardScaffolds
        else:
            newScaffold = Scaffold(gene.location[1], gene.location[0], [gene])
            scaffolds = reverseScaffolds
        running = True
        while running:
            running = False
            for scaffold in scaffolds:
                if abs(newScaffold.stop - scaffold.start) < scaffoldingDistance:
                    newScaffold.stop = scaffold.stop
                    newScaffold.genes += scaffold.genes
                    scaffolds.remove(scaffold)
                    running = True
                    break
                elif abs(newScaffold.start - scaffold.stop) < scaffoldingDistance:
                    newScaffold.start = scaffold.start
                    newScaffold.genes += scaffold.genes
                    scaffolds.remove(scaffold)
                    running = True
                    break
        scaffolds.append(newScaffold)
    map(lambda x: x.genes.sort(key = lambda y: (y.location[0] + y.location[1])/2), forwardScaffolds)
    map(lambda x: x.genes.sort(key = lambda y: (y.location[0] + y.location[1])/2), reverseScaffolds)
    return forwardScaffolds, reverseScaffolds

def overlap(intervalOne, intervalTwo):
    if intervalOne.start < intervalTwo.start and intervalTwo.start <= intervalOne.stop and intervalOne.stop < intervalTwo.stop:
        return intervalOne.stop - intervalTwo.start
    elif intervalTwo.start < intervalOne.start and intervalOne.start <= intervalTwo.stop and intervalTwo.stop < intervalOne.stop:
        return intervalTwo.stop - intervalOne.start
    elif intervalOne.start <= intervalTwo.start and intervalTwo.stop <= intervalOne.stop:
        return intervalTwo.stop - intervalTwo.start
    elif intervalTwo.start <= intervalOne.start and intervalOne.stop <= intervalTwo.stop:
        return intervalOne.stop - intervalOne.start
    else:
        return -1
    
def filterScaffolds(originalForwardScaffolds, originalReverseScaffolds):
    forwardScaffolds, reverseScaffolds = copy.deepcopy(originalForwardScaffolds), copy.deepcopy(originalReverseScaffolds)
    newForwardScaffolds, newReverseScaffolds = copy.copy(forwardScaffolds), copy.copy(reverseScaffolds)
    for forwardScaffold in forwardScaffolds:
        for reverseScaffold in copy.copy(newReverseScaffolds):
            forwardScaffoldRemoved = False
            while overlap(forwardScaffold, reverseScaffold) > 3:
                forwardHasGenemark = reduce(lambda x, y: x or not y.intergenic, forwardScaffold.genes, False)
                reverseHasGenemark = reduce(lambda x, y: x or not y.intergenic, reverseScaffold.genes, False)
                forwardCenter = (forwardScaffold.start + forwardScaffold.stop)/2
                reverseCenter = (reverseScaffold.start + reverseScaffold.stop)/2
                if forwardCenter < reverseCenter:
                    if forwardHasGenemark and not reverseHasGenemark:
                        removable = [reverseScaffold.genes[0]]
                    elif not forwardHasGenemark and reverseHasGenemark:
                        removable = [forwardScaffold.genes[-1]]
                    else:
                        removable = [forwardScaffold.genes[-1], reverseScaffold.genes[0]]
                else:
                    if forwardHasGenemark and not reverseHasGenemark:
                        removable = [reverseScaffold.genes[-1]]
                    elif not forwardHasGenemark and reverseHasGenemark:
                        removable = [forwardScaffold.genes[0]]
                    else:
                        removable = [reverseScaffold.genes[-1], forwardScaffold.genes[0]]
                removable = filter(lambda x: x.intergenic, removable)
                if removable:
                    toRemove = min(removable, key = lambda x: abs(x.location[1] - x.location[0]))
                    if toRemove.location[0] < toRemove.location[1]:
                        forwardScaffold.genes.remove(toRemove)
                        if not forwardScaffold.genes:
                            newForwardScaffolds.remove(forwardScaffold)
                            forwardScaffoldRemoved = True
                            break
                        elif forwardCenter < reverseCenter:
                            forwardScaffold.stop = forwardScaffold.genes[-1].location[1]
                        else:
                            forwardScaffold.start = forwardScaffold.genes[0].location[0]
                    else:
                        reverseScaffold.genes.remove(toRemove)
                        if not reverseScaffold.genes:
                            newReverseScaffolds.remove(reverseScaffold)
                            break
                        elif forwardCenter < reverseCenter:
                            reverseScaffold.start = reverseScaffold.genes[0].location[1]
                        else:
                            reverseScaffold.stop = reverseScaffold.genes[-1].location[0]
                elif forwardHasGenemark and reverseHasGenemark:
                    break
                elif forwardScaffold.stop - forwardScaffold.start < reverseScaffold.stop - reverseScaffold.start:
                    newForwardScaffolds.remove(forwardScaffold)
                    forwardScaffoldRemoved = True
                    break
                else:
                    newReverseScaffolds.remove(reverseScaffold)
                    break
            if forwardScaffoldRemoved:
                break
    return newForwardScaffolds, newReverseScaffolds

def refineScaffolds(genes, scaffoldingDistance):
    forwardScaffolds, reverseScaffolds = extractScaffolds(genes.values(), scaffoldingDistance)
    forwardFiltered, reverseFiltered = filterScaffolds(forwardScaffolds, reverseScaffolds)
    remainingGenes = reduce(lambda x, y: x + y, map(lambda x: x.genes, forwardFiltered), [])
    remainingGenes.extend(reduce(lambda x, y: x + y, map(lambda x: x.genes, reverseFiltered), []))
    return dict(map(lambda x: (x.query, x), remainingGenes))
