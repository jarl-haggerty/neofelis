"""
Module for cleaning excess signals.
"""

import itertools
import copy

def removeSignals(genes, signals):
    """
    genes:   List of Iteration objects.
    signals: List of signals(any object with a location tuple) to filter

    return:  A filtered list of the singals which does not contain any signal which has a center of
             mass inside a gene

    A helper function for filter signals that expects genes and signals on the same direction.
    """
    remainingSignals = copy.copy(signals)
    
    for signal in signals:
        nearStart = False
        for gene in genes:
            center = (signal.location[0] + signal.location[1])/2
            nearStart = nearStart or abs(gene.location[0] - center) < 100
            if min(gene.location) < center and center < max(gene.location):
                remainingSignals.remove(signal)
                break
        if not nearStart and signal in remainingSignals:
            remainingSignals.remove(signal)
            
    return remainingSignals

def filterSignals(genes, signals):
    """
    genes:   List of Iteration objects.
    signals: List of signals(any object with a location tuple) to filter

    return:  A filtered list of the singals which does not contain any signal which has a center of
             mass inside a gene
    """
    forwardGenes = filter(lambda x: x.location[0] < x.location[1], genes)
    forwardSignals = filter(lambda x: x.location[0] < x.location[1], signals)
    forwardResult = removeSignals(forwardGenes, forwardSignals)

    reverseGenes = filter(lambda x: x.location[1] < x.location[0], genes)
    reverseSignals = filter(lambda x: x.location[1] < x.location[0], signals)
    reverseResult = removeSignals(reverseGenes, reverseSignals)
    
    return forwardResult + reverseResult
