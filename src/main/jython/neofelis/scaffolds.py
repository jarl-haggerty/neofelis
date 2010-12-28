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

def extractScaffolds(genes, scaffoldingDistance = 100):
    forwardScaffolds, reverseScaffolds = [], []
    for gene in genes:
        if gene.location[0] < gene.location[1]:
            newScaffold = [gene.location[0], gene.location[1]]
            scaffolds = forwardScaffolds
        else:
            newScaffold = [gene.location[1], gene.location[0]]
            scaffolds = reverseScaffolds
        running = True
        while running:
            running = False
            for scaffold in scaffolds:
                if abs(newScaffold[1] - scaffold[0]) < scaffoldingDistance:
                    newScaffold[1] = scaffold[1]
                    scaffolds.remove(scaffold)
                    running = True
                    break
                elif abs(newScaffold[0] - scaffold[1]) < scaffoldingDistance:
                    newScaffold[0] = scaffold[0]
                    scaffolds.remove(scaffold)
                    running = True
                    break
        scaffolds.append(newScaffold)
    return forwardScaffolds, reverseScaffolds

def overlap(intervalOne, intervalTwo):
    if intervalOne[0] < intervalTwo[0] and intervalTwo[0] <= intervalOne[1] and intervalOne[1] < intervalTwo[1]:
        return intervalOne[1] - intervalTwo[0]
    elif intervalTwo[0] < intervalOne[0] and intervalOne[0] <= intervalTwo[1] and intervalTwo[1] < intervalOne[1]:
        return intervalTwo[1] - intervalOne[0]
    elif intervalOne[0] <= intervalTwo[0] and intervalTwo[1] <= intervalOne[1]:
        return intervalTwo[1] - intervalTwo[0]
    elif intervalTwo[0] <= intervalOne[0] and intervalOne[1] <= intervalTwo[1]:
        return intervalOne[1] - intervalOne[0]
    else:
        return -1
    
def filterScaffolds(forwardScaffolds, reverseScaffolds):
    newForwardScaffolds, newReverseScaffolds = copy.deepcopy(forwardScaffolds), copy.deepcopy(reverseScaffolds)
    for forwardScaffold in forwardScaffolds:
        for reverseScaffold in reverseScaffolds:
            if overlap(forwardScaffold, reverseScaffold) > 3:
                if forwardScaffold[1] - forwardScaffold[0] < reverseScaffold[1] - reverseScaffold[0]:
                    newForwardScaffolds.remove(forwardScaffold)
                    break
                elif reverseScaffold in newReverseScaffolds:
                    newReverseScaffolds.remove(reverseScaffold)
    return newForwardScaffolds, newReverseScaffolds

def maskedGenes(genes, masks):
    result = []
    for gene in genes:
        for mask in masks:
            center = (gene.location[0] + gene.location[1])/2
            if mask[0] < center and center < mask[1]:
                result.append(gene)
    return result

def refineScaffolds(genes, scaffoldingDistance):
    forwardScaffolds, reverseScaffolds = extractScaffolds(genes.values(), scaffoldingDistance)
    forwardFiltered, reverseFiltered = filterScaffolds(forwardScaffolds, reverseScaffolds)
    masked = maskedGenes(filter(lambda x: x.location[0] < x.location[1], genes.values()), forwardFiltered)
    masked.extend(maskedGenes(filter(lambda x: x.location[0] > x.location[1], genes.values()), reverseFiltered))
    return dict(map(lambda x: (x.query, x), masked))
