"""
The central component of the Neofelis pipeline.
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
import sys
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from neofelis import genemark
from neofelis import extend
from neofelis import intergenic
from neofelis import promoters
from neofelis import terminators
from neofelis import artemis
from neofelis import utils
from neofelis import scaffolds
from neofelis import signals
from neofelis import report
from javax.swing import JFrame
from javax.swing import JPanel
from javax.swing import JLabel
from javax.swing import JButton
from javax.swing import JProgressBar
from javax.swing import AbstractAction
from java.awt import GridBagLayout
from java.awt import GridBagConstraints

#If a swing interface is asked for this will be the JFrame.
frame = None
#Keeps track of the number of queries processed.
jobCount = 0
#Keeps track of the query currently being processed.
currentJob = ""
#Keeps track of the massage to be displayed.
message = 0
#Messages to be displayed at each stage in the processing of a single query.
messages = ["Searching for genes via genemark",
            "Extending genes found via genemark",
            "Searching for intergenic genes",
            "Removing overlapping genes",
            "Using BPROM to find promoters",
            "Using transterm to find terminators",
            "Removing transcription signals which conflict with genes"]

class DoneAction(AbstractAction):
  """
  Action for finishing the pipeline.  Simply exits the program.
  """
  def __init__(self):
    AbstractAction.__init__(self, "Done")

  def actionPerformed(self, event):
    global frame
    sys.exit(0)

def initializeDisplay(queries, swing):
  """
  queries: A list of the fasts files to be processed.
  swing:   If true then updates about progress will be displayed in a swing window, otherwise they will be written to stdout.
  
  Initializes the interface for telling the user about progress in the pipeline.  Queries is used to count the
  number of queries the pipeline will process and to size the swing display(if it is used) so that text
  isn't cutoff at the edge of the window.  The swing display is setup if swing is true.
  """
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
    currentLabel.setText(" ")
    frame.setLocationRelativeTo(None)
    frame.setVisible(True)

def updateProgress(job):
  """
  query: Name of the query currently being processed.
  
  This function use used for updating the progress shown in the interface.  If job is not equal to currentJob then
  global progress is incremented and shown and the currentProgress is reset and shown.  If job is equal to currentJob
  then the globalProgress does not change and currentProgress is incremented.
  """
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
  """
  This function is to be called at the end of the pipeline.  Informs the user that the pipeline is finished
  and if a swing interface is being used the Done button is enabled.
  """
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

def run(blastLocation, genemarkLocation, transtermLocation, database, eValue, matrix, minLength, scaffoldingDistance, ldfCutoff, queries, swing = False, email = ""):
  """
  blastLocation:       Directory blast was installed in.
  genemarkLocation:    Directory genemark was installed in.
  transtermLocation:   Directory transterm was installed in.
  database:            Name of the blast database to use.
  eValue:              The e value used whenever a blast search is done.
  matrix:              The matrix to use when running genemark.  If None then genemark is run heuristically.
  minLength:           Minimum length of any genes included in the resulting annotation.
  scaffoldingDistance: The maximum length allowed between genes when contiguous regions of genes are being identified
  ldfCutoff:           Minimum LDF allowed for any promoters included in the resulting annotation
  queries:             A list of faster files to process.
  swing:               If true a swing window will be used to updated the user about the pipeline's progress.
  email:               If this is a non-empty string an email will be sent to the address in the string when the pipeline is done.  The local machine will be used as
                       an SMTP server and this will not work if it isn't.
  
  The main pipeline function.  For every query genemark is used to predict genes, these genes are then extended to any preferable starts.  Then the pipeline searches
  for any intergenic genes(genes between those found by genemark) and these are combined with the extended genemark genes.  Then the genes are pruned to remove
  any undesirable genes found in the intergenic stage.  BPROM and Transterm are used to find promoters and terminators, which are then pruned to remove any
  signals which are inside or too far away from any genes.  Finally, all the remaining genes, promoters, and terminators ar written to an artemis file in the directory
  of the query with the same name but with a .art extension, and .dat and .xls files will be generating describing the blast results of the final genes.
  """
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
    initialGenes = genemark.findGenes("query.fas", name, blastLocation, database, eValue, genemarkLocation, matrix)
    #artemis.writeArtemisFile(os.path.splitext(query)[0] + ".genemark.art", genome, initialGenes.values())
    
    updateProgress(query)
    extendedGenes = extend.extendGenes("query.fas", initialGenes, name, blastLocation, database, eValue)
    #artemis.writeArtemisFile(os.path.splitext(query)[0] + ".extended.art", genome, extendedGenes.values())
    
    updateProgress(query)
    intergenicGenes = intergenic.findIntergenics("query.fas", extendedGenes, name, minLength, blastLocation, database, eValue)
    #artemis.writeArtemisFile(os.path.splitext(query)[0] + ".intergenic.art", genome, intergenicGenes.values())
    genes = {}
    for k, v in extendedGenes.items() + intergenicGenes.items():
      genes[k] = v
    
    updateProgress(query)
    scaffolded = scaffolds.refineScaffolds(genes, scaffoldingDistance)
 
    updateProgress(query)
    initialPromoters = promoters.findPromoters("query.fas", name, ldfCutoff)
    
    updateProgress(query)
    initialTerminators = terminators.findTerminators("query.fas", name, genes.values(), transtermLocation)

    updateProgress(query)
    filteredSignals = signals.filterSignals(scaffolded.values(), initialPromoters + initialTerminators)
    filteredPromoters = filter(lambda x: isinstance(x, promoters.Promoter), filteredSignals)
    filteredTerminators = filter(lambda x: isinstance(x, terminators.Terminator), filteredSignals)

    artemis.writeArtemisFile(os.path.splitext(query)[0] + ".art", genome, scaffolded.values(), filteredPromoters, filteredTerminators)
    
    report.report(name, scaffolded, os.path.splitext(query)[0])

  if email:
    message = MIMEText("Your genome has been annotated.")
    message["Subject"] = "Annotation complete"
    message["From"] = "Neofelis"
    message["To"] = email
    

    smtp = smtplib.SMTP("tmpl.arizona.edu", 587)
    smtp.ehlo()
    smtp.starttls()
    smtp.ehlo
    smtp.login("tmpl", "gobanana23")
    smtp.sendmail("Neofelis", [email], message.as_string())
    smtp.close()
    
  finished()
