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
import sys
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
from java.awt import GridBagLayout
from java.awt import GridBagConstraints

def getArguments():
  """
  This function brings up a window to retreive any required arguments.  This function brings up a window with fields for each argument, filled with any arguments already given.
  While this window is visible the program will wait, once it is no longer visible all the arguments will be filled with the entries in the fields.
  """
  global blastLocation, genemarkLocation, transtermLocation, database, matrix, eValue, minLength, scaffoldingDistance, ldfCutoff, sources, email

  class BlastAction(AbstractAction):
    """
    Action for selecting the location of Blast+.  Brings up a file selection dialog and fills the text field for blast with the selection.
    """
    def __init__(self):
      AbstractAction.__init__(self, "...")

    def actionPerformed(self, event):
      fileChooser = JFileChooser()
      fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY)
      if fileChooser.showOpenDialog(None) == JFileChooser.APPROVE_OPTION:
        blastField.setText(fileChooser.getSelectedFile().getAbsolutePath())

  class GenemarkAction(AbstractAction):
    """
    Action for selecting the location of Genemark.  Brings up a file selection dialog and fills the text field for genemark with the selection.
    """
    def __init__(self):
      AbstractAction.__init__(self, "...")

    def actionPerformed(self, event):
      fileChooser = JFileChooser()
      fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY)
      if fileChooser.showOpenDialog(None) == JFileChooser.APPROVE_OPTION:
        genemarkField.setText(fileChooser.getSelectedFile().getAbsolutePath())

  class TranstermAction(AbstractAction):
    """
    Action for selecting the location of Transterm.  Brings up a file selection dialog and fills the text field for transterm with the selection.
    """
    def __init__(self):
      AbstractAction.__init__(self, "...")

    def actionPerformed(self, event):
      fileChooser = JFileChooser()
      fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY)
      if fileChooser.showOpenDialog(None) == JFileChooser.APPROVE_OPTION:
        transtermField.setText(fileChooser.getSelectedFile().getAbsolutePath())

  class QueryAction(AbstractAction):
    """
    Action for selecting the query file or directory.  Brings up a file selection dialog and fills the text field for the query with the selection.
    """
    def __init__(self):
      AbstractAction.__init__(self, "...")

    def actionPerformed(self, event):
      fileChooser = JFileChooser()
      fileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES)
      if fileChooser.showOpenDialog(None) == JFileChooser.APPROVE_OPTION:
        queryField.setText(fileChooser.getSelectedFile().getAbsolutePath())

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
  blastField = JTextField(blastLocation)
  genemarkField = JTextField(genemarkLocation)
  transtermField = JTextField(transtermLocation)
  databaseField = JTextField(database)
  matrixField = JTextField(str(matrix))
  eValueField = JTextField(str(eValue))
  minLengthField = JTextField(str(minLength))
  scaffoldingDistanceField = JTextField(str(scaffoldingDistance))
  ldfField = JTextField(str(ldfCutoff))
  queryField = JTextField(sources[0])
  #emailField = JTextField(email)

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
  contentPane.add(JLabel("Database"), constraints)
  constraints.gridy = 4
  contentPane.add(JLabel("Matrix(Leave blank to use heuristic matrix)"), constraints)
  constraints.gridy = 5
  contentPane.add(JLabel("E Value"), constraints)
  constraints.gridy = 6
  contentPane.add(JLabel("Minimum Intergenic Length"), constraints)
  constraints.gridy = 7
  contentPane.add(JLabel("Scaffold Distance"), constraints)
  constraints.gridy = 8
  contentPane.add(JLabel("LDF Cutoff"), constraints)
  constraints.gridy = 9
  contentPane.add(JLabel("Query"), constraints)
  #constraints.gridy = 10
  #contentPane.add(JLabel("Email"), constraints)
  constraints.gridx = 1
  constraints.gridy = 0
  constraints.weightx = 1
  contentPane.add(blastField, constraints)
  constraints.gridy = 1
  contentPane.add(genemarkField, constraints)
  constraints.gridy = 2
  contentPane.add(transtermField, constraints)
  constraints.gridy = 3
  contentPane.add(databaseField, constraints)
  constraints.gridy = 4
  contentPane.add(matrixField, constraints)
  constraints.gridy = 5
  contentPane.add(eValueField, constraints)
  constraints.gridy = 6
  contentPane.add(minLengthField, constraints)
  constraints.gridy = 7
  contentPane.add(scaffoldingDistanceField, constraints)
  constraints.gridy = 8
  contentPane.add(ldfField, constraints)
  constraints.gridy = 9
  contentPane.add(queryField, constraints)
  #constraints.gridy = 10
  #contentPane.add(emailField, constraints)
  constraints.gridx = 2
  constraints.gridy = 0
  constraints.weightx = 0
  constraints.fill = GridBagConstraints.NONE
  constraints.anchor = GridBagConstraints.LINE_END
  contentPane.add(JButton(BlastAction()), constraints)
  constraints.gridy = 1
  contentPane.add(JButton(GenemarkAction()), constraints)
  constraints.gridy = 2
  contentPane.add(JButton(TranstermAction()), constraints)
  constraints.gridy = 9
  contentPane.add(JButton(QueryAction()), constraints)

  constraints.gridx = 1
  constraints.gridy = 11
  contentPane.add(JButton(OKAction()), constraints)
  constraints.gridx = 2
  contentPane.add(JButton(CancelAction()), constraints)

  frame.pack()
  frame.setLocationRelativeTo(None)
  frame.setVisible(True)

  while frame.isVisible():
    pass

  blastLocation = blastField.getText()
  genemarkLocation = genemarkField.getText()
  transtermLocation = transtermField.getText()
  database = databaseField.getText()
  matrix = matrixField.getText()
  eValue = float(eValueField.getText())
  minLength = int(minLengthField.getText())
  scaffoldingDistance = int(scaffoldingDistanceField.getText())
  ldfCutoff = float(ldfField.getText())
  sources = [queryField.getText()]

def main(arguments):
  global blastLocation, genemarkLocation, transtermLocation, database, matrix, eValue, minLength, scaffoldingDistance, ldfCutoff, sources, email
  
  documentation = """
-m --matrix               Matrix with which to run genemark
-d --database             Database to use when running blast
-g --genemark             Location of Genemark
-b --blast                Location of Blast+
-e --e-value              Minimal evalue for any genes detected
-l --min-length           Minimum length of any genes discovered
-t --transterm            Location of Transterm
-d --ldf-cutoff           Minimum LDF value of any promoters selected from a BPROM search
-c --scaffolding-distance Distance to allow between genes when determining scaffolds
-q --query                Genome or directory of genomes to run pipeline on
-h --help                 Print help documentation
-s --swing                Use a swing interface
-p --pipe                 Neofelis will be set to read lines of command line arguments from the pipe, "neofelis_pipe".  For each line read a new thread will be spawned to process the query.
-n --no-swing             If any required arguments are missing then the program will exit instead of using a Swing interface to get the missing arguments
-a --email                Email address that results will be sent to.
"""
  try:
    opts, args = getopt(arguments, "m:d:g:b:e:l:t:c:q:hsna:", ["matrix=", "database=", "genemark=", "blast=", "e-value=", "min-length=", "transterm=", "ldf-cutoff=", "scaffolding-distance=", "query=", "help", "swing", "no-swing", "email="])
  except GetoptError:
    print documentation
    sys.exit(0)

  matrix = ""
  database = ""
  genemarkLocation = ""  
  blastLocation = ""
  eValue = 0.1
  minLength = 100  
  transtermLocation = ""
  ldfCutoff = 0
  scaffoldingDistance = 100
  sources = [""]
  swingInterface = False
  noSwing = False
  email = ""
  remote = False
  pipe = False
	
  for opt, arg in opts:
    if opt in ("-q", "--query"):
      sources = [arg]
    elif opt in ("-g", "--genemark"):
      genemarkLocation = arg
    elif opt in ("-b", "--blast"):
      blastLocation = arg
    elif opt in ("-d", "--database"):
      database = arg
    elif opt in ("-m", "--matrix"):
      matrix = arg
    elif opt in ("-e", "--e-value"):
      eValue = float(arg)
    elif opt in ("-l", "--min-length"):
      minLength = int(arg)
    elif opt in ("-l", "--ldf-cutoff"):
      ldfCutoff = float(arg)
    elif opt in ("-t", "--transterm"):
      transtermLocation = arg
    elif opt in ("-c", "--scaffolding-distance"):
      scaffoldingDistance = int(arg)
    elif opt in ("-s", "--swing"):
      swingInterface = True
    elif opt in ("-n", "--no-swing"):
      noSwing = True
    elif opt in ("-a", "--email"):
      email = arg
    elif opt in ("-p", "--pipe"):
      pipe = True
    elif opt in ("-h", "--help"):
      print documentation
      sys.exit(0)

  if pipe:
    os.mkfifo("neofelis_pipe", "r")
    pipe = open("neofelis_pipe", "r")
    while True:
      line = pipe.readLine()
      if line:
        threading.Thread(target = main, args = re.split(r"\s+", line))
        
  if not blastLocation or not database or not genemarkLocation or not transtermLocation:
    if noSwing:
      sys.exit(1)
    else:
      getArguments()
      swingInterface = True

  queries = []
  while sources:
    source = sources.pop()
    if os.path.isdir(source):
      newSources = map(lambda x: os.path.join(source, x), os.listdir(source))
      sources.extend(newSources)
    elif utils.isGenome(source):
      queries.append(source)

  pipeline.run(blastLocation, genemarkLocation, transtermLocation, database, eValue, matrix, minLength, scaffoldingDistance, ldfCutoff, queries, swingInterface, email)

if __name__ == "__main__":
  main(sys.argv)
