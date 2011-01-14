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
import sys
from neofelis import genemark
from neofelis import extend
from neofelis import intergenic
from neofelis import promoters
from neofelis import terminators
from neofelis import artemis
from neofelis import utils
from neofelis import scaffolds
from neofelis import signals
from javax.swing import JFrame
from javax.swing import JPanel
from javax.swing import JLabel
from javax.swing import JButton
from javax.swing import JProgressBar
from javax.swing import AbstractAction
from java.awt import GridBagLayout
from java.awt import GridBagConstraints

frame = None
jobCount = 0
currentJob = ""
message = 0
messages = ["Searching for genes via genemark",
            "Extending genes found via genemark",
            "Searching for intergenic genes",
            "Removing overlapping genes",
            "Using BPROM to find promoters",
            "Using transterm to find terminators",
            "Removing transcription signals which conflict with genes"]

class DoneAction(AbstractAction):
  def __init__(self):
    AbstractAction.__init__(self, "Done")

  def actionPerformed(self, event):
    global frame
    sys.exit(0)

def initializeDisplay(queries, swing):
  global numJobs, frame, currentLabel, currentProgress, globalLabel, globalProgress, doneButton, messages
  
  numJobs = len(queries)
  if swing:
    frame = JFrame("Neofelis")
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE)
    contentPane = JPanel(GridBagLayout())
    frame.setContentPane(contentPane)
    globalLabel = JLabel(max(queries, key = len))
    globalProgress = JProgressBar(0, numJobs)
    currentLabel = JLabel(max(messages, key = len))
    currentProgress = JProgressBar(0, len(messages))
    doneButton = JButton(DoneAction())
    doneButton.setEnabled(False)

    constraints = GridBagConstraints()

    constraints.gridx, constraints.gridy = 0, 0
    constraints.gridwidth, constraints.gridheight = 1, 1
    constraints.weightx = 1
    constraints.fill = GridBagConstraints.HORIZONTAL
    contentPane.add(globalLabel, constraints)
    constraints.gridy = 1
    contentPane.add(globalProgress, constraints)
    constraints.gridy = 2
    contentPane.add(currentLabel, constraints)
    constraints.gridy = 3
    contentPane.add(currentProgress, constraints)
    constraints.gridy = 4
    constraints.weightx = 0
    constraints.fill = GridBagConstraints.NONE
    constraints.anchor = GridBagConstraints.LINE_END
    contentPane.add(doneButton, constraints)
    
    frame.pack()
    frame.setResizable(False)
    globalLabel.setText(" ")
    globalLabel.setText(" ")
    frame.setLocationRelativeTo(None)
    frame.setVisible(True)

def updateProgress(job):
  global numJobs, currentJob, jobCount, message, messages, frame, currentLabel, currentProgress, globalLabel, globalProgress
  if frame:
    if job != currentJob:
      currentProgress.setValue(currentProgress.getMaximum())
      globalLabel.setText(job)
      globalProgress.setValue(jobCount)
      jobCount += 1
      currentJob = job
      message = -1
    message += 1
    currentProgress.setValue(message)
    currentLabel.setText(messages[message])
  else:
    if job != currentJob:
      print "Processing %s, %.2f%% done" % (job, 100.0*jobCount/numJobs)
      jobCount += 1
      currentJob = job
      message = -1
    message += 1
    print "    %s, %.2f%% done" % (messages[message], 100.0*message/len(messages))

def finished():
  global frame, currentLabel, currentProgress, globalLabel, globalProgress, doneButton
  if frame:
    globalLabel.setText("Finished")
    globalProgress.setValue(globalProgress.getMaximum())
    currentLabel.setText(" ")
    currentProgress.setValue(currentProgress.getMaximum())
    doneButton.setEnabled(True)
    while frame.isVisible():
      pass
  else:
    print "Processing 100.00% done"

def run(blastLocation, genemarkLocation, transtermLocation, database, eValue, matrix, minLength, scaffoldingDistance, remote, ldfCutoff, queries, swing = False):
  initializeDisplay(queries, swing)
    
  for query in queries:
    name = os.path.splitext(query)[0]
    queryDirectory, name = os.path.split(name)
    genome = utils.loadGenome(query)
    queryFile = open("query.fas", "w")
    queryFile.write(">" + name + "\n")
    for i in range(0, len(genome), 50):
      queryFile.write(genome[i:min(i+50, len(genome))] + "\n")
    queryFile.close()

    updateProgress(query)
    initialGenes = genemark.findGenes("query.fas", name, blastLocation, database, eValue, genemarkLocation, matrix, remote)
    artemis.writeArtemisFile(name + "genemark.art", genome, initialGenes.values())
    
    updateProgress(query)
    extendedGenes = extend.extendGenes("query.fas", initialGenes, name, blastLocation, database, eValue, remote)
    artemis.writeArtemisFile(name + "extended.art", genome, extendedGenes.values())
    
    updateProgress(query)
    intergenicGenes = intergenic.findIntergenics("query.fas", extendedGenes, name, minLength, blastLocation, database, eValue, remote)
    genes = {}
    for k, v in extendedGenes.items() + intergenicGenes.items():
      genes[k] = v
    artemis.writeArtemisFile(name + "intergenic.art", genome, genes.values())
    
    updateProgress(query)
    scaffolded = scaffolds.refineScaffolds(genes, scaffoldingDistance)
    artemis.writeArtemisFile(name + "scaffolds.art", genome, scaffolded.values())
    artemis.writeArtemisFile(os.path.splitext(query)[0] + ".art", genome, scaffolded.values())
 
    updateProgress(query)
    initialPromoters = promoters.findPromoters("query.fas", name, ldfCutoff)
    artemis.writeArtemisFile(name + "promoters.art", genome, scaffolded.values(), initialPromoters)
    
    updateProgress(query)
    initialTerminators = terminators.findTerminators("query.fas", name, genes.values(), transtermLocation)
    artemis.writeArtemisFile(name + "promoters-and-terminators.art", genome, scaffolded.values(), initialPromoters, initialTerminators)

    updateProgress(query)
    filteredSignals = signals.filterSignals(genes.values(), initialPromoters + initialTerminators)
    filteredPromoters = filter(lambda x: isinstance(x, promoters.Promoter), filteredSignals)
    filteredTerminators = filter(lambda x: isinstance(x, terminators.Terminator), filteredSignals)

    artemis.writeArtemisFile(name + "final.art", genome, scaffolded.values(), filteredPromoters, filteredTerminators)
    
  finished()
