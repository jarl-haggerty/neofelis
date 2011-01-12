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
    def __init__(self, start, stop, genes, containsGenemarkGene):
        self.start = start
        self.stop = stop
        self.genes = initialGene
        self.containsGenemarkGene = containsGenemarkGene

    def __eq__(self, that):
        return self.start == that.start and \
               self.stop == that.stop and \
               self.containsGenemarkGene  == that.containsGenemarkGene

def extractScaffolds(genes, scaffoldingDistance = 100):
    forwardScaffolds, reverseScaffolds = [], []
    for gene in genes:
        if gene.location[0] < gene.location[1]:
            newScaffold = Scaffold(gene.location[0], gene.location[1], [gene], gene.note != "Intergenic")
            scaffolds = forwardScaffolds
        else:
            newScaffold = Scaffold(gene.location[1], gene.location[0], [gene], gene.note != "Intergenic")
            scaffolds = reverseScaffolds
        running = True
        while running:
            running = False
            for scaffold in scaffolds:
                if abs(newScaffold.stop - scaffold.start) < scaffoldingDistance:
                    newScaffold.stop = scaffold.stop
                    newScaffold.genes += scaffold.genes
                    newScaffold.containsGenemarkGene = newScaffold.containsGenemarkGene or scaffold.containsGenemarkGene
                    scaffolds.remove(scaffold)
                    running = True
                    break
                elif abs(newScaffold.start - scaffold.stop) < scaffoldingDistance:
                    newScaffold.start = scaffold.start
                    newScaffold.genes += scaffold.genes
                    newScaffold.containsGenemarkGene = newScaffold.containsGenemarkGene or scaffold.containsGenemarkGene
                    scaffolds.remove(scaffold)
                    running = True
                    break
        scaffolds.append(newScaffold)
    map(lambda x: x.genes.sort(x.genes, key = lambda x: (x.location[0] + x.location[1])/2), forwardScaffolds)
    map(lambda x: x.genes.sort(x.genes, key = lambda x: (x.location[0] + x.location[1])/2), reverseScaffolds)
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
    
def filterScaffolds(forwardScaffolds, reverseScaffolds):
    newForwardScaffolds, newReverseScaffolds = copy.copy(forwardScaffolds), copy.copy(reverseScaffolds)
    for forwardScaffold in forwardScaffolds:
        for reverseScaffold in reverseScaffolds:
            forwardScaffoldRemoved = False
            while overlap(forwardScaffold, reverseScaffold) > 3:
                removable = []
                forwardCenter, reverseCenter = (forwardScaffold.start + forwardScaffold.stop)/2, (reverseScaffold.start + reverseScaffold.stop)/2
                if forwardCenter < reverseCenter:
                    removable = filter(lambda x: x.note == "Intergenic", [forwardScaffold.genes[-1], reverseScaffold.genes[0]])
                else:
                    removable = filter(lambda x: x.note == "Intergenic", [reverseScaffold.genes[-1], forwardScaffold.genes[0]])
                if removable:
                    toRemove = min(removable, key = lambda x: abs(x.location[1] - x.location[0]))
                    if toRemove.location[0] < toRemove.location[1]:
                        forwardScaffold.genes.remove(toRemove)
                        if forwardCenter < reverseCenter:
                            forwardScaffold.stop = forwardScaffold.genes[-1].location[1]
                        else:
                            forwardScaffold.start = forwardScaffold.genes[0].location[0]
                    else:
                        reverseScaffold.genes.remove(toRemove)
                        if forwardCenter < reverseCenter:
                            reverseScaffold.stop = reverseScaffold.genes[0].location[0]
                        else:
                            reverseScaffold.start = reverseScaffold.genes[-1].location[1]
                elif forwardScaffold.stop - forwardScaffold.start < reverseScaffold.stop - reverseScaffold.start:
                    newForwardScaffolds.remove(forwardScaffold)
                    forwardScaffoldRemoved = True
                    break
                elif reverseScaffold in newReverseScaffolds:
                    newReverseScaffolds.remove(reverseScaffold)
                    break
            if forwardScaffoldRemoved:
                break
    return newForwardScaffolds, newReverseScaffolds

def maskedGenes(genes, masks):
    result = []
    for gene in genes:
        for mask in masks:
            center = (gene.location[0] + gene.location[1])/2
            if mask.start < center and center < mask.stop:
                result.append(gene)
    return result

def refineScaffolds(genes, scaffoldingDistance):
    forwardScaffolds, reverseScaffolds = extractScaffolds(genes.values(), scaffoldingDistance)
    forwardFiltered, reverseFiltered = filterScaffolds(forwardScaffolds, reverseScaffolds)
    masked = maskedGenes(filter(lambda x: x.location[0] < x.location[1], genes.values()), forwardFiltered)
    masked.extend(maskedGenes(filter(lambda x: x.location[0] > x.location[1], genes.values()), reverseFiltered))
    return dict(map(lambda x: (x.query, x), masked))
