"""
Function for finding significant RNA sequences.
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

import os
import re
import subprocess

class TransferRNA():
    """
    Class representing a transfer RNA.
    """
    def __init__(self, start, stop, type, antiCodon, coveScore):
        self.location, self.type, self.antiCodon, self.coveScore = [start, stop], type, antiCodon, coveScore

def findtRNAs(tRNAscanLocation, query):
    """
    tRNAscanLocation: Directory that tRNAscan resides in.
    query:            Name of the fasta file to scan.

    return: List of TransferRNA objects.

    Uses tRNAscan to find transfer RNAs in a fasta file.
    """
    tRNAscanProcess = subprocess.Popen([tRNAscanLocation + "/tRNAscan-SE", "-P", os.getcwd() + "/" + query], stdout = subprocess.PIPE, stderr = subprocess.PIPE, env = {"PATH" : os.getenv("PATH") + ":" + tRNAscanLocation}, cwd = tRNAscanLocation)
    result = ""
    while tRNAscanProcess.poll() == None:
        result += tRNAscanProcess.stdout.read()
        tRNAscanProcess.stderr.read()
    nextRead = tRNAscanProcess.stdout.read()
    while nextRead:
        result += nextRead
        nextRead = tRNAscanProcess.stdout.read()

    transferRNAs = []
    for line in result.split("\n"):
        match = re.match(".+\s+\d+\s+(\d+)\s+(\d+)\s+(\w+)\s+([ACTG?]+)\s+\d+\s+\d+\s+(\d*\.\d*)", line)
        if match:
            transferRNAs.append(TransferRNA(int(match.group(1)), int(match.group(2)), match.group(3), match.group(4), float(match.group(5))))
    return transferRNAs
    
