"""
This Module is the main file for the Neofelis Genome Annotator
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

def getParameters():
  global blastLocation, genemarkLocation, database, matrix, eValue, minLength, scaffoldingDistance, sources

  class BlastAction(AbstractAction):
    def __init__(self):
      AbstractAction.__init__(self, "...")

    def actionPerformed(self, event):
      fileChooser = JFileChooser()
      fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY)
      if fileChooser.showOpenDialog(None) == JFileChooser.APPROVE_OPTION:
        blastField.setText(fileChooser.getSelectedFile().getAbsolutePath())

  class GenemarkAction(AbstractAction):
    def __init__(self):
      AbstractAction.__init__(self, "...")

    def actionPerformed(self, event):
      fileChooser = JFileChooser()
      fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY)
      if fileChooser.showOpenDialog(None) == JFileChooser.APPROVE_OPTION:
        genemarkField.setText(fileChooser.getSelectedFile().getAbsolutePath())

  class QueryAction(AbstractAction):
    def __init__(self):
      AbstractAction.__init__(self, "...")

    def actionPerformed(self, event):
      fileChooser = JFileChooser()
      fileChooser.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES)
      if fileChooser.showOpenDialog(None) == JFileChooser.APPROVE_OPTION:
        queryField.setText(fileChooser.getSelectedFile().getAbsolutePath())

  class OKAction(AbstractAction):
    def __init__(self):
      AbstractAction.__init__(self, "Ok")

    def actionPerformed(self, event):
      frame.setVisible(False)

  class CancelAction(AbstractAction):
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
  databaseField = JTextField(database)
  matrixField = JTextField(str(matrix))
  eValueField = JTextField(str(eValue))
  minLengthField = JTextField(str(minLength))
  scaffoldingDistanceField = JTextField(str(scaffoldingDistance))
  queryField = JTextField(sources[0])

  constraints.gridx = 0
  constraints.gridy = 0
  constraints.gridwidth = 1
  constraints.gridheight = 1
  constraints.fill = GridBagConstraints.HORIZONTAL
  constraints.weightx = 0
  constraints.weighty = 0
  contentPane.add(JLabel("Blast Location"), constraints)
  constraints.gridy = 1
  contentPane.add(JLabel("Genemark Label"), constraints)
  constraints.gridy = 2
  contentPane.add(JLabel("Database"), constraints)
  constraints.gridy = 3
  contentPane.add(JLabel("Matrix(Leave blank to use heuristic matrix)"), constraints)
  constraints.gridy = 4
  contentPane.add(JLabel("E Value"), constraints)
  constraints.gridy = 5
  contentPane.add(JLabel("Minimum Intergenic Length"), constraints)
  constraints.gridy = 6
  contentPane.add(JLabel("Scaffold Distance"), constraints)
  constraints.gridy = 7
  contentPane.add(JLabel("Query"), constraints)
  constraints.gridx = 1
  constraints.gridy = 0
  constraints.weightx = 1
  contentPane.add(blastField, constraints)
  constraints.gridy = 1
  contentPane.add(genemarkField, constraints)
  constraints.gridy = 2
  contentPane.add(databaseField, constraints)
  constraints.gridy = 3
  contentPane.add(matrixField, constraints)
  constraints.gridy = 4
  contentPane.add(eValueField, constraints)
  constraints.gridy = 5
  contentPane.add(minLengthField, constraints)
  constraints.gridy = 6
  contentPane.add(scaffoldingDistanceField, constraints)
  constraints.gridy = 7
  contentPane.add(queryField, constraints)
  constraints.gridx = 2
  constraints.gridy = 0
  constraints.weightx = 0
  constraints.fill = GridBagConstraints.NONE
  constraints.anchor = GridBagConstraints.LINE_END
  contentPane.add(JButton(BlastAction()), constraints)
  constraints.gridy = 1
  contentPane.add(JButton(GenemarkAction()), constraints)
  constraints.gridy = 7
  contentPane.add(JButton(QueryAction()), constraints)

  constraints.gridx = 1
  constraints.gridy = 8
  contentPane.add(JButton(OKAction()), constraints)
  constraints.gridx = 2
  contentPane.add(JButton(CancelAction()), constraints)

  frame.pack()
  frame.setResizable(False)
  frame.setLocationRelativeTo(None)
  frame.setVisible(True)

  while frame.isVisible():
    pass

  blastLocation = blastField.getText()
  genemarkLocation = genemarkField.getText()
  database = databaseField.getText()
  matrix = matrixField.getText()
  eValue = float(eValueField.getText())
  minLength = int(minLengthField.getText())
  scaffoldingDistance = int(scaffoldingDistanceField.getText())

if __name__ == "__main__":
  documentation = """
-m --matrix               :Matrix with which to run genemark
-d --database             :Database to use when running blast
-g --genemark             :Location of genemark
-b --blast                :Location of blast
-l --min-length           :Minimum length of any genes discovered
-e --eValue               :Minimal evalue for any genes detected
-r --remote               :Run blast with remote NCBI servers.
-c --scaffolding-distance :Distance to allow between genes when determining scaffolds
-h --help                 :Print help documentation
-q --query                :Genome or directory of genomes to run pipeline on
-s --swing                :Use a swing interface
"""
  try:
    opts, args = getopt(sys.argv, "m:d:g:b:e:rq:hs", ["matrix=", "database=", "genemark=", "blast=", "eValue=", "remote", "query=", "help", "swing"])
  except GetoptError:
    print documentation
    sys.exit(0)
  
  blastLocation = ""
  database = ""
  genemarkLocation = ""
  eValue = 0.1
  matrix = ""
  minLength = 100
  scaffoldingDistance = 100
  sources = [""]
  remote = False
  swingInterface = False
	
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
    elif opt in ("-e", "--eValue"):
      eValue = float(arg)
    elif opt in ("-l", "--min-length"):
      minLength = int(arg)
    elif opt in ("-r", "--remote"):
      remote = True
    elif opt in ("-c", "--scaffolding-distance"):
      scaffoldingDistance = int(arg)
    elif opt in ("-s", "--swing"):
      swingInterface = True
    elif opt in ("-h", "--help"):
      print documentation
      sys.exit(0)

  if not blastLocation or not database or not genemarkLocation:
    getParameters()
    swingInterface = True

  queries = []
  while sources:
    source = sources.pop()
    if os.path.isdir(source):
      newSources = filter(lambda x: utils.isGenome(os.path.join(source, x)), os.listdir(source))
      newSources = map(lambda x: os.path.join(source, x), newSources)
      sources.extend(newSources)
    else:
      queries.append(source)

  print queries

  pipeline.run(blastLocation, genemarkLocation, database, eValue, matrix, minLength, scaffoldingDistance, remote, queries, swingInterface)

  """
		#remove all signals inside orfs
		for f in os.listdir("artemis_complete/"):
			for c in cut_offs:
				print name + "_" + c
				print f
				if name + "_" + c in f:
					signal_command_now = signal_command_now + "artemis_complete/" + f + " "
					
		print signal_command_now
					
		os.system(signal_command_now)
		
		#build fasta file of all orfs
		fasta_command_now = "python make_fasta.py " + q
		print fasta_command_now
		os.system(fasta_command_now)
"""	
