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

import itertools
import copy

def removeSignals(genes, signals):
    remainingSignals = copy.copy(signals)
    for signal in signals:
        for gene in genes:
            center = (signal.location[0] + signal.location[1])/2
            if min(gene.location) < center and center < max(gene.location):
                remainingSignals.remove(signal)
                break
    return remainingSignals

def filterSignals(genes, signals):
    forwardGenes = filter(lambda x: x.location[0] < x.location[1], genes)
    forwardSignals = filter(lambda x: x.location[0] < x.location[1], signals)
    forwardResult = removeSignals(forwardGenes, forwardSignals)

    reverseGenes = filter(lambda x: x.location[1] < x.location[0], genes)
    reverseSignals = filter(lambda x: x.location[1] < x.location[0], signals)
    reverseResult = removeSignals(reverseGenes, reverseSignals)
    
    return forwardResult + reverseResult
