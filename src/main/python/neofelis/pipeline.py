"""
The central component of the Neofelis pipeline.
"""

import os
import sys
import subprocess
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
from neofelis import rna
from javax.swing import JFrame
from javax.swing import JPanel
from javax.swing import JLabel
from javax.swing import JButton
from javax.swing import JProgressBar
from javax.swing import AbstractAction
from java.awt.event import WindowAdapter
from java.awt import GridBagLayout
from java.awt import GridBagConstraints

class PipelineException(Exception):
  """
  Used to terminate the pipeline early.
  """
  def __init__(self):
    Exception.__init__(self, "Stop the Pipeline")

class DoneAction(AbstractAction):
  """
  Action for finishing the pipeline.  Simply exits the program.
  """
  def __init__(self, frame):
    AbstractAction.__init__(self, "Done")
    self.frame = frame

  def actionPerformed(self, event):
    self.frame.dispose()

class PipelineWindowAdapter(WindowAdapter):
  """
  Window adapter for the swing interface.
  """
  def __init__(self, pipeline):
    self.pipeline = pipeline

  def windowClosing(self, event):
    self.pipeline.frame.dispose()
    
  def windowClosed(self, event):
    self.pipeline.exception = PipelineException()

class Pipeline():
  def __init__(self):
    #If a swing interface is asked for this will be the JFrame.
    self.frame = None
    #Keeps track of the number of queries processed.
    self.jobCount = 0
    #Keeps track of the query currently being processed.
    self.currentJob = ""
    #Keeps track of the massage to be displayed.
    self.message = 0
    #Messages to be displayed at each stage in the processing of a single query.
    self.messages = ["Searching for genes via genemark",
                     "Extending genes found via genemark",
                     "Searching for intergenic genes",
                     "Removing overlapping genes",
                     "Searching for promoters",
                     "Using transterm to find terminators",
                     "Removing transcription signals which conflict with genes",
                     "Using tRNAscan to find transfer RNAs",
                     "Writing Artemis file",
                     "Writing summary .xml, .html, and .xls files"]
    self.exception = None

  def initializeDisplay(self, queries, swing):
    """
    queries: A list of the fasts files to be processed.
    swing:   If true then updates about progress will be displayed in a swing window, otherwise they will be written to stdout.
    
    Initializes the interface for telling the user about progress in the pipeline.  Queries is used to count the
    number of queries the pipeline will process and to size the swing display(if it is used) so that text
    isn't cutoff at the edge of the window.  The swing display is setup if swing is true.
    """
  
    self.numJobs = len(queries)
    if swing:
      self.frame = JFrame("Neofelis")
      self.frame.addWindowListener(PipelineWindowAdapter(self))
      contentPane = JPanel(GridBagLayout())
      self.frame.setContentPane(contentPane)
      self.globalLabel = JLabel(max(queries, key = len))
      self.globalProgress = JProgressBar(0, self.numJobs)
      self.currentLabel = JLabel(max(self.messages, key = len))
      self.currentProgress = JProgressBar(0, len(self.messages))
      self.doneButton = JButton(DoneAction(self.frame))
      self.doneButton.setEnabled(False)

      constraints = GridBagConstraints()
      
      constraints.gridx, constraints.gridy = 0, 0
      constraints.gridwidth, constraints.gridheight = 1, 1
      constraints.weightx = 1
      constraints.fill = GridBagConstraints.HORIZONTAL
      contentPane.add(self.globalLabel, constraints)
      constraints.gridy = 1
      contentPane.add(self.globalProgress, constraints)
      constraints.gridy = 2
      contentPane.add(self.currentLabel, constraints)
      constraints.gridy = 3
      contentPane.add(self.currentProgress, constraints)
      constraints.gridy = 4
      constraints.weightx = 0
      constraints.fill = GridBagConstraints.NONE
      constraints.anchor = GridBagConstraints.LINE_END
      contentPane.add(self.doneButton, constraints)
    
      self.frame.pack()
      self.frame.setResizable(False)
      self.globalLabel.setText(" ")
      self.currentLabel.setText(" ")
      self.frame.setLocationRelativeTo(None)
      self.frame.setVisible(True)

  def updateProgress(self, job):
    """
    query: Name of the query currently being processed.
    
    This function use used for updating the progress shown in the interface.  If job is not equal to currentJob then
    global progress is incremented and shown and the currentProgress is reset and shown.  If job is equal to currentJob
    then the globalProgress does not change and currentProgress is increased.
    """
    if self.exception:
      raise self.exception
    
    if self.frame:
      if job != self.currentJob:
        self.currentProgress.setValue(self.currentProgress.getMaximum())
        self.globalLabel.setText(job)
        self.globalProgress.setValue(self.jobCount)
        print "Processing %s, %.2f%% done" % (job, 100.0*self.jobCount/self.numJobs)
        self.jobCount += 1
        self.currentJob = job
        self.message = -1
      self.message += 1
      print "    %s, %.2f%% done" % (self.messages[self.message], 100.0*self.message/len(self.messages))
      self.currentProgress.setValue(self.message)
      self.currentLabel.setText(self.messages[self.message])
    else:
      if job != self.currentJob:
        print "Processing %s, %.2f%% done" % (job, 100.0*self.jobCount/self.numJobs)
        self.jobCount += 1
        self.currentJob = job
        self.message = -1
      self.message += 1
      print "    %s, %.2f%% done" % (self.messages[self.message], 100.0*self.message/len(self.messages))

  def finished(self):
    """
    This function is to be called at the end of the pipeline.  Informs the user that the pipeline is finished
    and if a swing interface is being used the Done button is enabled.
    """
    print "Processing 100.00% done"
    if self.frame:
      self.globalLabel.setText("Finished")
      self.globalProgress.setValue(self.globalProgress.getMaximum())
      self.currentLabel.setText(" ")
      self.currentProgress.setValue(self.currentProgress.getMaximum())
      self.doneButton.setEnabled(True)
      while self.frame.isVisible():
        pass

  def run(self, blastLocation, genemarkLocation, transtermLocation, tRNAscanLocation, database, eValue, matrix, minLength, scaffoldingDistance, promoterScoreCutoff, queries, swing = False, email = ""):
    """
    blastLocation:       Directory blast was installed in.
    genemarkLocation:    Directory genemark was installed in.
    transtermLocation:   Directory transterm was installed in.
    tRNAscanLocation:    Directory tRNAscan was installed in.
    database:            Name of the blast database to use.
    eValue:              The e value used whenever a blast search is done.
    matrix:              The matrix to use when running genemark.  If None then genemark is run heuristically.
    minLength:           Minimum length of any genes included in the resulting annotation.
    scaffoldingDistance: The maximum length allowed between genes when contiguous regions of genes are being identified
    promoterScoreCutoff: Minimum score allowed for any promoters included in the resulting annotation
    queries:             A list of faster files to process.
    swing:               If true a swing window will be used to updated the user about the pipeline's progress.
    email:               If this is a non-empty string an email will be sent to the address in the string when the pipeline is done.  This will be attempted with the sendmail command on the local computer.
    
    The main pipeline function.  For every query genemark is used to predict genes, these genes are then extended to any preferable starts.  Then the pipeline searches
    for any intergenic genes(genes between those found by genemark) and these are combined with the extended genemark genes.  Then the genes are pruned to remove
    any undesirable genes found in the intergenic stage.  BPROM and Transterm are used to find promoters and terminators, which are then pruned to remove any
    signals which are inside or too far away from any genes.  Next, tRNAscan is used to find any transfer RNAs in the genome.  Finally, all the remaining genes,
    promoters, and terminators are written to an artemis file in the directory of the query with the same name but with a .art extension, and .xml, .html, and
    .xls files will be generating describing the blast results of the final genes.
    """
    self.initializeDisplay(queries, swing)

    try:
      for query in queries:
        name = os.path.splitext(query)[0]
        queryDirectory, name = os.path.split(name)
        
        genome = utils.loadGenome(query)
        swapFileName = "query" + str(id(self)) + ".fas"
        queryFile = open(swapFileName, "w")
        queryFile.write(">" + name + "\n")
        for i in range(0, len(genome), 50):
          queryFile.write(genome[i:min(i+50, len(genome))] + "\n")
        queryFile.close()

        self.updateProgress(query)
        initialGenes = genemark.findGenes(swapFileName, name, blastLocation, database, eValue, genemarkLocation, matrix, self)
      
        self.updateProgress(query)
        extendedGenes = extend.extendGenes(swapFileName, initialGenes, name, blastLocation, database, eValue, self)
    
        self.updateProgress(query)
        intergenicGenes = intergenic.findIntergenics(swapFileName, extendedGenes, name, minLength, blastLocation, database, eValue, self)

        genes = {}
        for k, v in extendedGenes.items() + intergenicGenes.items():
          genes[k] = v
        
        self.updateProgress(query)
        scaffolded = scaffolds.refineScaffolds(genes, scaffoldingDistance)
 
        self.updateProgress(query)
        initialPromoters = promoters.findPromoters(swapFileName, name, promoterScoreCutoff, self.frame)
    
        self.updateProgress(query)
        initialTerminators = terminators.findTerminators(swapFileName, name, genes.values(), transtermLocation)
      
        self.updateProgress(query)
        filteredSignals = signals.filterSignals(scaffolded.values(), initialPromoters + initialTerminators)
        filteredPromoters = filter(lambda x: isinstance(x, promoters.Promoter), filteredSignals)
        filteredTerminators = filter(lambda x: isinstance(x, terminators.Terminator), filteredSignals)

        self.updateProgress(query)
        transferRNAs = rna.findtRNAs(tRNAscanLocation, swapFileName)

        os.remove(swapFileName)

        self.updateProgress(query)
        artemis.writeArtemisFile(os.path.splitext(query)[0] + ".art", genome, scaffolded.values(), filteredPromoters, filteredTerminators, transferRNAs)

        self.updateProgress(query)
        report.report(name, scaffolded, os.path.splitext(query)[0])

      if email:
        if not os.path.isfile("EMAIL_MESSAGE"):
          message = open("EMAIL_MESSAGE", "w")
          message.write("Subject: Annotation Complete\nYour genome has been annotated.\n")
          message.close()
        
        sent = False
        while not sent:
          message = open("EMAIL_MESSAGE", "r")
          sendmailProcess = subprocess.Popen(["/usr/sbin/sendmail", "-F", "Neofelis", "-f", "neofelis@tmpl.arizona.edu", email],
                                             stdin = message,
                                             stdout = subprocess.PIPE)
          result = ""
          nextRead = sendmailProcess.stdout.read()
          while nextRead:
            result += nextRead
            nextRead = sendmailProcess.stdout.read()
          sent = not result.strip()
          message.close()
    
      self.finished()
    except PipelineException:
      return
