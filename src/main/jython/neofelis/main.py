"""
This Module is the entry point for Neofelis when used on the desktop.  It parses any command line arguments and if any required arguments
are missing a window is displayed to collect the remaining arguments.  If a directory is specified as the query then that directory and all subdirectories
will be searched for fasta files with a single genome, these files will then be used as queries.  All the arguments are then passed onto the pipeline.
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
import sys
import socket
import threading
from getopt import getopt
from neofelis import pipeline
from neofelis import utils
from javax.swing import JFrame
from javax.swing import JPanel
from javax.swing import JFileChooser
from javax.swing import JButton
from javax.swing import JTextField
from javax.swing import AbstractAction
from javax.swing import JLabel
from javax.swing import JEditorPane
from javax.swing import JScrollPane
from java.awt import GridBagLayout
from java.awt import GridBagConstraints
from java.awt import Toolkit
from java.lang import ClassLoader
from java.lang import System
from java.lang import Class
from java.lang import Runtime
from java.lang import String

class NeofelisThread(threading.Thread):
  def __init__(self, arguments):
    threading.Thread.__init__(self)
    self.arguments = arguments
    self.name = "NeofelisThread"
    
  def run(self):
    self.program = Main()
    print "NeofelisThread", self.arguments
    self.program.run(self.arguments)

  def terminate(self):
    self.program.pipeline.exception = pipeline.PipelineException()

class Main():
  def __init__(self):
    self.shutDown = False
    
  def getArguments(self):
    """
    This function brings up a window to retreive any required arguments.  This function brings up a window with fields for each argument, filled with any arguments already given.
    While this window is visible the program will wait, once it is no longer visible all the arguments will be filled with the entries in the fields.
    """

    class LocationAction(AbstractAction):
      """
      Action for selecting the database.  Brings up a file selection dialog and fills the text field for the database with the selection.
      """
      def __init__(self, field):
        AbstractAction.__init__(self, "...")
        self.field = field

      def actionPerformed(self, event):
        fileChooser = JFileChooser()
        fileChooser.fileSelectionMode = JFileChooser.FILES_AND_DIRECTORIES
        if fileChooser.showOpenDialog(None) == JFileChooser.APPROVE_OPTION:
          field.text = fileChooser.selectedFile.absolutePath
    

    class HelpAction(AbstractAction):
      """
      Displays a help dialog.
      """
      def __init__(self):
        AbstractAction.__init__(self, "Help")

      def actionPerformed(self, event):
        browsers = ["google-chrome", "firefox", "opera", "epiphany", "konqueror", "conkeror", "midori", "kazehakase", "mozilla"]
        osName = System.getProperty("os.name")
        helpHTML = ClassLoader.getSystemResource("help.html").toString()
        if osName.find("Mac OS") == 0:
          Class.forName("com.apple.eio.FileManager").getDeclaredMethod( "openURL", [String().getClass()]).invoke(None, [helpHTML])
        elif osName.find("Windows") == 0:
          Runtime.getRuntime().exec( "rundll32 url.dll,FileProtocolHandler " + helpHTML)
        else:
          browser = None
          for b in browsers:
            if browser == None and Runtime.getRuntime().exec(["which", b]).getInputStream().read() != -1:
              browser = b
              Runtime.getRuntime().exec([browser, helpHTML])

    class OKAction(AbstractAction):
      """
      Action for starting the pipeline.  This action will simply make the window invisible.
      """
      def __init__(self):
        AbstractAction.__init__(self, "Ok")

      def actionPerformed(self, event):
        frame.setVisible(False)

    class CancelAction(AbstractAction):
      """
      Action for canceling the pipeline.  Exits the program.
      """
      def __init__(self):
        AbstractAction.__init__(self, "Cancel")

      def actionPerformed(self, event):
        sys.exit(0)

    frame = JFrame("Neofelis")
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE)
    constraints = GridBagConstraints()
    contentPane = JPanel(GridBagLayout())
    frame.setContentPane(contentPane)
    blastField = JTextField(self.blastLocation)
    genemarkField = JTextField(self.genemarkLocation)
    transtermField = JTextField(self.transtermLocation)
    tRNAscanField = JTextField(self.tRNAscanLocation)
    databaseLocationField = JTextField(os.path.split(self.database)[0])
    databaseField = JTextField(os.path.split(self.database)[1])
    matrixField = JTextField(str(self.matrix))
    eValueField = JTextField(str(self.eValue))
    minLengthField = JTextField(str(self.minLength))
    scaffoldingDistanceField = JTextField(str(self.scaffoldingDistance))
    promoterScoreField = JTextField(str(self.promoterScoreCutoff))
    queryField = JTextField(self.sources[0])

    constraints.gridx = 0
    constraints.gridy = 0
    constraints.gridwidth = 1
    constraints.gridheight = 1
    constraints.fill = GridBagConstraints.HORIZONTAL
    constraints.weightx = 0
    constraints.weighty = 0
    contentPane.add(JLabel("Blast Location"), constraints)
    constraints.gridy = 1
    contentPane.add(JLabel("Genemark Location"), constraints)
    constraints.gridy = 2
    contentPane.add(JLabel("Transterm Location"), constraints)
    constraints.gridy = 3
    contentPane.add(JLabel("tRNAscan Location"), constraints)
    constraints.gridy = 4
    contentPane.add(JLabel("Databases Location"), constraints)
    constraints.gridy = 5
    contentPane.add(JLabel("Database"), constraints)
    constraints.gridy = 6
    contentPane.add(JLabel("Matrix(Leave blank to use heuristic matrix)"), constraints)
    constraints.gridy = 7
    contentPane.add(JLabel("E Value"), constraints)
    constraints.gridy = 8
    contentPane.add(JLabel("Minimum Intergenic Length"), constraints)
    constraints.gridy = 9
    contentPane.add(JLabel("Scaffold Distance"), constraints)
    constraints.gridy = 0
    contentPane.add(JLabel("Promoter Score Cutoff"), constraints)
    constraints.gridy = 11
    contentPane.add(JLabel("Query"), constraints)
    constraints.gridx = 1
    constraints.gridy = 0
    constraints.weightx = 1
    contentPane.add(blastField, constraints)
    constraints.gridy = 1
    contentPane.add(genemarkField, constraints)
    constraints.gridy = 2
    contentPane.add(transtermField, constraints)
    constraints.gridy = 3
    contentPane.add(tRNAscanField, constraints)
    constraints.gridy = 4
    contentPane.add(databaseLocationField, constraints)
    constraints.gridy = 5
    contentPane.add(databaseField, constraints)
    constraints.gridy = 6
    contentPane.add(matrixField, constraints)
    constraints.gridy = 7
    contentPane.add(eValueField, constraints)
    constraints.gridy = 8
    contentPane.add(minLengthField, constraints)
    constraints.gridy = 9
    contentPane.add(scaffoldingDistanceField, constraints)
    constraints.gridy = 10
    contentPane.add(promoterScoreField, constraints)
    constraints.gridy = 11
    contentPane.add(queryField, constraints)
    constraints.gridx = 2
    constraints.gridy = 0
    constraints.weightx = 0
    constraints.fill = GridBagConstraints.NONE
    constraints.anchor = GridBagConstraints.LINE_END
    contentPane.add(JButton(LocationAction(blastField)), constraints)
    constraints.gridy = 1
    contentPane.add(JButton(LocationAction(genemarkField)), constraints)
    constraints.gridy = 2
    contentPane.add(JButton(LocationAction(transtermField)), constraints)
    constraints.gridy = 3
    contentPane.add(JButton(LocationAction(tRNAscanField)), constraints)
    constraints.gridy = 4
    contentPane.add(JButton(LocationAction(databaseLocationField)), constraints)
    constraints.gridy = 11
    contentPane.add(JButton(LocationAction(queryField)), constraints)

    constraints.gridx = 0
    constraints.gridy = 12
    constraints.anchor = GridBagConstraints.LINE_START
    contentPane.add(JButton(HelpAction()), constraints)
    constraints.gridx = 1
    constraints.anchor = GridBagConstraints.LINE_END
    contentPane.add(JButton(OKAction()), constraints)
    constraints.gridx = 2
    contentPane.add(JButton(CancelAction()), constraints)

    frame.pack()
    frame.setLocationRelativeTo(None)
    frame.setVisible(True)

    while frame.isVisible():
      pass

    self.blastLocation = blastField.getText()
    self.genemarkLocation = genemarkField.getText()
    self.transtermLocation = transtermField.getText()
    self.database = databaseLocationField.getText() + "/" + databaseField.getText()
    self.matrix = matrixField.getText()
    self.eValue = float(eValueField.getText())
    self.minLength = int(minLengthField.getText())
    self.scaffoldingDistance = int(scaffoldingDistanceField.getText())
    self.promoterScoreCutoff = float(promoterScoreField.getText())
    self.sources = [queryField.getText()]


  def run(self, arguments):
    documentation = """
-m --matrix                Matrix with which to run genemark
-d --database              Database to use when running blast
-g --genemark              Location of Genemark
-b --blast                 Location of Blast+
-e --e-value               Minimal evalue for any genes detected
-l --min-length            Minimum length of any genes discovered
-t --transterm             Location of Transterm
-d --promoter-score-cutoff Minimum promoter score of any promoters selected from a promoter search
-c --scaffolding-distance  Distance to allow between genes when determining scaffolds
-q --query                 Genome or directory of genomes to run pipeline on
-h --help                  Print help documentation
-s --swing                 Use a swing interface
-v --server                Neofelis will be set to read lines of command line arguments from the socket 1122 on localhost.  For each line read a new thread will be spawned to process the query.  Only one connection will be accepted per process.
-u --shutdown              If Neofelis is already running as a server then this will shut it down.
-n --no-swing              If any required arguments are missing then the program will exit instead of using a Swing interface to get the missing arguments
-a --email                 Email address that results will be sent to.
"""
    try:
      opts, args = getopt(arguments, "m:d:g:b:e:l:t:c:q:hsna:vw:r:y:", ["matrix=", "database=", "genemark=", "blast=", "e-value=", "min-length=", "transterm=", "ldf-cutoff=", "scaffolding-distance=", "query=", "help", "swing", "no-swing", "email=", "server", "smtp-server=", "smtp-user=", "smtp-password="])
    except GetoptError:
      print documentation
      sys.exit(0)
        
    self.matrix = ""
    self.database = ""
    self.genemarkLocation = ""  
    self.blastLocation = ""
    self.eValue = 0.1
    self.minLength = 100  
    self.transtermLocation = ""
    self.tRNAscanLocation = ""
    self.promoterScoreCutoff = 0
    self.scaffoldingDistance = 100
    self.sources = [""]
    self.swingInterface = False
    self.noSwing = False
    self.email = ""
    self.remote = False
    self.server = False
    self.smtpServer = ""
    self.smtpUser = ""
    self.smtpPassword = ""
    
    for opt, arg in opts:
      if opt in ("-q", "--query"):
        self.sources = [arg]
      elif opt in ("-g", "--genemark"):
        self.genemarkLocation = arg
      elif opt in ("-b", "--blast"):
        self.blastLocation = arg
      elif opt in ("-d", "--database"):
        self.database = arg
      elif opt in ("-m", "--matrix"):
        self.matrix = arg
      elif opt in ("-e", "--e-value"):
        self.eValue = float(arg)
      elif opt in ("-l", "--min-length"):
        self.minLength = int(arg)
      elif opt in ("-d", "--promoter-score-cutoff"):
        self.promoterScoreCutoff = float(arg)
      elif opt in ("-t", "--transterm"):
        self.transtermLocation = arg
      elif opt in ("-c", "--scaffolding-distance"):
        self.scaffoldingDistance = int(arg)
      elif opt in ("-s", "--swing"):
        self.swingInterface = True
      elif opt in ("-n", "--no-swing"):
        self.noSwing = True
      elif opt in ("-a", "--email"):
        self.email = arg
      elif opt in ("-v", "--server"):
        self.server = True
      elif opt in ("-w", "--smtp-server"):
        self.smtpServer = arg
      elif opt in ("-r", "--smtp-user"):
        self.smtpUser = arg
      elif opt in ("-y", "--smtp-password"):
        self.smtpPassword = arg
      elif opt in ("-h", "--help"):
        print documentation
        sys.exit(0)
        
    if self.server:
      print self.server
      s = socket.socket()
      s.bind(("localhost", 1122))
      s.listen(5)
      running = True
      try:
        while running:
          connection, address = s.accept()
          recieved, abort = "", False
          print "listening"
          while not "\n" in recieved:
            data = connection.recv(4096)
            if data:
              recieved += data
            elif data == None:
              abort = True
              break
          connection.close()
          print recieved, abort
          if abort:
            continue
          arguments = re.split(r"\s+", recieved[:recieved.find("\n")])
          print arguments, recieved, running
          if "--server" in arguments or "-v" in arguments:
            print "ERROR: Can't run a neofelis server within a neofelis server, you just had to try that didn't you?"
          elif "--shutdown" in arguments or "-u" in arguments:
            running = False
          else:
            print "Threading"
            NeofelisThread(arguments).start()
          print "Continue", running
      except Exception, e:
        print e
      finally:
        for thread in threading.enumerate():
          if isinstance(thread, NeofelisThread):
            thread.terminate()
        s.close()
        return

    print "Getting", self.blastLocation, self.database, self.genemarkLocation, self.transtermLocation, self.sources, self.noSwing
    if not self.blastLocation or not self.database or not self.genemarkLocation or not self.transtermLocation or self.sources == [""]:
      if self.noSwing:
        sys.exit(1)
      else:
        self.getArguments()
        self.swingInterface = True

    self.queries = []
    while self.sources:
      source = self.sources.pop()
      if os.path.isdir(source):
        newSources = map(lambda x: os.path.join(source, x), os.listdir(source))
        self.sources.extend(newSources)
      elif utils.isGenome(source):
        self.queries.append(source)
        
    print self.blastLocation, self.genemarkLocation, self.transtermLocation, self.database, self.eValue, self.matrix, self.minLength, self.scaffoldingDistance, self.promoterScoreCutoff, self.queries, self.swingInterface, self.email
    self.pipeline = pipeline.Pipeline()
    self.pipeline.run(self.blastLocation, self.genemarkLocation, self.transtermLocation, self.tRNAscanLocation, self.database, self.eValue, self.matrix, self.minLength, self.scaffoldingDistance, self.promoterScoreCutoff, self.queries, self.swingInterface, self.email, self.smtpServer, self.smtpUser, self.smtpPassword)

if __name__ == "__main__":
  Main().run(sys.argv)
