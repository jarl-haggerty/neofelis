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

def extractScaffolds(genes, scaffoldingDistance = 100):
    forwardScaffolds, reverseScaffolds = [], []
    for gene in genes:
        if gene.location[0] < gene.location[1]:
            newScaffold = gene.location[0], gene.location[1]
            scaffolds = forwardScaffolds
        else:
            newScaffold = gene.location[1], gene.location[0]
            scaffolds = reverseScaffolds
        running = True
        while running:
            running = False
            for scaffold in scaffolds:
                if abs(newScaffold[1] - scaffold[0]) < scaffoldingDistance:
                    newScaffold[1] = scaffold[1]
                    running = True
                    break
                elif abs(scaffold[1] - newScaffold[0]) < scaffoldingDistance:
                    newScaffold[0] = scaffold[0]
                    running = True
                    break
        scaffolds.append(newScaffold)
    return forwardScaffolds, reverseScaffolds

def overlap(intervalOne, intervalTwo):
    if intervalOne[0] < intervalTwo[0] and intervalTwo[0] <= intervalOne[1] and intervalOne[1] < intervalTwo[1]:
        return intervalOne[1] - intervalTwo[0]
    elif intervalTwo[0] < intervalOne[0] and intervalOne[0] <= intervalTwo[1] and intervalTwo[1] < intervalOne[1]:
        return intervalTwo[1] - intervalOne[0]
    elif intervalOne[0] < intervalTwo[0] and intervalTwo[1] < intervalOne[1]:
        return intervalTwo[1] - intervalTwo[0]
    else:
        return intervalOne[1] - intervalOne[0]
    
def filterScaffolds(forwardScaffolds, reverseScaffolds):
    newForwardScaffolds, newReverseScaffolds = copy.deepcopy(forwardScaffolds), copy.deepcopy(reverseScaffolds)
    for forwardScaffold in fowardScaffolds:
        for reverseScaffold in reverseScaffolds:
            if overlap(forwardScaffold, reverseScaffold) > 3:
                if forwardScaffold[1] - forwardScaffold[0] < reverseScaffold[1] - reverseScaffold[0]:
                    newForwardScaffolds.remove(forwardScaffold)
                    break
                else:
                    newReverseScaffolds.remove(reverseScaffold)
    return newForwardScaffolds, newReverseScaffolds

def refineScaffolds(genes, scaffoldingDistance):
    forwardScaffolds, reverseScaffolds = extractScaffolds(genes.values(), scaffoldingDistance)
    return filterScaffolds(forwardScaffolds, reverseScaffolds)
    
