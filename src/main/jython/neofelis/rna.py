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

import subprocess

class TransferRNA():
    def __init__(self, start, stop):
        self.start, self.stop = start, stop

def findtRNAs(tRNAscanLocation, query):
    tRNAscanProcess = subprocess.Popen([tRNAscanLocation + "/tRNAscan-SE", "-P", query], stdout = subprocess.PIPE)
    result = ""
    while tRNAscanProcess.poll() == None:
        result += tRNAscanProcess.stdout.read()
    nextRead = tRNAscanProcess.stdout.read()
    while nextRead:
        result += nextRead
        nextRead = tRNAscanProcess.stdout.read()

    transferRNAs = []
    for line in result.split("\n"):
        match = re.match(".+\s+\d+\s+(\d+)\s+(\d+)\s+\w+\s+[ACTG]+\s+\d+\s+\d+\s+\d*\.\d*\n", line)
        if match:
            transferRNAs += TransferRNA(int(match.groups(0)), int(match.groups(1)))
    return transferRNAs
    
