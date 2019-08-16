import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
from auto3dgm_nazar.mesh.meshexport import MeshExport
import auto3dgm_nazar
from auto3dgm_nazar.dataset.datasetfactory import DatasetFactory
import numpy as np
#
# Auto3dgm
#

class Auto3dgm(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Auto3dgm" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Examples"]
    self.parent.dependencies = []
    self.parent.contributors = ["John Doe (AnyWare Corp.)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
This is an example of scripted loadable module bundled in an extension.
It performs a simple thresholding on the input volume and optionally captures a screenshot.
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
and Steve Pieper, Isomics, Inc. and was partially funded by NIH grant 3P41RR013218-12S1.
""" # replace with organization, grant and thanks.

    # module level data

#
# Auto3dgmModuleData
#
  """May be a temporary stand-in class, this stores module-level results and data
  not stored in MRML scene"""
class Auto3dgmData():
  def __init__(self):
    self.dataset = None
    """ To be discussed?
    1) Store all the datasets into dataset collection
    2) Add phase1 and phase2 sampled points. The idea here is that
    once we have subsampled the mesh, we dont want to access thesbased    on the slider value (that the user can change, resulting in a crash if not subsampled again). Instead: Record the keys everytime subsampling is done
    """
    self.datasetCollection = None
    self.phase1SampledPoints = None
    self.phase2SampledPoints = None
    self.aligned_meshes = []
#
# Auto3dgmWidget
#

class Auto3dgmWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)

    #Widget variables
    self.meshFolder = None
    self.outputFolder = None
    self.visualizationMeshFolder = None
    self.outputFolderPrepared = False
    self.Auto3dgmData = Auto3dgmData()

    # Instantiate and connect widgets ...
    tabsWidget = qt.QTabWidget()
    setupTab = qt.QWidget()
    setupTabLayout = qt.QFormLayout(setupTab)
    
    runTab = qt.QWidget()
    runTabLayout = qt.QFormLayout(runTab)
    outTab = qt.QWidget()
    outTabLayout = qt.QFormLayout(outTab)

    tabsWidget.addTab(setupTab, "Setup")
    tabsWidget.addTab(runTab, "Run")
    tabsWidget.addTab(outTab, "Visualize/Output")
    self.layout.addWidget(tabsWidget)

    self.setupSetupTab(setupTabLayout)
    self.setupRunTab(runTabLayout)
    self.setupOutTab(outTabLayout)
  
  ### SETUP TAB WIDGETS AND BEHAVIORS ###

  def setupSetupTab(self, setupTabLayout):
    inputfolderWidget=ctk.ctkCollapsibleButton()
    inputfolderLayout=qt.QFormLayout(inputfolderWidget)
    inputfolderWidget.text = "Input and output folder"
    setupTabLayout.addRow(inputfolderWidget)

    self.meshInputText, volumeInLabel, self.inputFolderButton=self.textIn('Input folder','Choose input folder', '')

    self.inputFolderButton.connect('clicked(bool)', self.selectMeshFolder)
    self.loadButton=qt.QPushButton("Load Data")
    self.loadButton.enabled=False
    self.loadButton.connect('clicked(bool)', self.onLoad)
    self.LMText, volumeInLabel, self.LMbutton=self.textIn('Input Directory','..', '')

    self.meshOutputText, volumeOutLabel, self.outputFolderButton=self.textIn('Output folder','Choose output folder', '')
    self.outputFolderButton.connect('clicked(bool)', self.selectOutputFolder)

    inputfolderLayout.addRow(volumeInLabel)#,1,1)    
    inputfolderLayout.addRow(self.meshInputText)#,1,2)
    inputfolderLayout.addRow(self.inputFolderButton)#,1,3)
    inputfolderLayout.addRow(self.loadButton)
    inputfolderLayout.addRow(self.meshOutputText)
    inputfolderLayout.addRow(self.outputFolderButton)
    self.LMbutton.connect('clicked(bool)', self.selectMeshFolder)

    self.parameterWidget = ctk.ctkCollapsibleButton()
    self.parameterLayout = qt.QFormLayout(self.parameterWidget)
    self.parameterWidget.text = "Parameters"
    setupTabLayout.addRow(self.parameterWidget)

    self.maxIterSliderWidget = ctk.ctkSliderWidget()
    self.maxIterSliderWidget.setDecimals(0)
    self.maxIterSliderWidget.singleStep = 1
    self.maxIterSliderWidget.minimum = 1
    self.maxIterSliderWidget.maximum = 5000
    self.maxIterSliderWidget.value = 1000
    self.maxIterSliderWidget.setToolTip("Maximum possible number of iterations for pairwise alignment optimization.")
    self.parameterLayout.addRow("Maximum iterations", self.maxIterSliderWidget)

    self.reflectionCheckBox = qt.QCheckBox()
    self.reflectionCheckBox.checked = 0
    self.reflectionCheckBox.setToolTip("Whether meshes can be reflected/mirrored to achieve more optimal alignments.")
    self.parameterLayout.addRow("Allow reflection", self.reflectionCheckBox)

    self.subsampleComboBox = qt.QComboBox()
    self.subsampleComboBox.addItem("FPS (Furthest Point Sampling)")
    self.subsampleComboBox.addItem("GPL (Gaussian Process Landmarks)")
    self.subsampleComboBox.addItem("FPS/GPL Hybrid")
    self.parameterLayout.addRow("Subsampling", self.subsampleComboBox)

    self.fpsSeed = qt.QSpinBox()
    self.fpsSeed.setSpecialValueText('-')
    self.parameterLayout.addRow("Optional FPS Seed", self.fpsSeed)

    self.hybridPoints = qt.QSpinBox()
    self.hybridPoints.setSpecialValueText('-')
    self.hybridPoints.setMinimum(1)
    self.hybridPoints.setMaximum(1000)
    self.parameterLayout.addRow("Hybrid GPL Points", self.hybridPoints)

    self.phaseChoiceComboBox = qt.QComboBox()
    self.phaseChoiceComboBox.addItem("1 (Single Alignment Pass)")
    self.phaseChoiceComboBox.addItem("2 (Double Alignment Pass)")
    self.phaseChoiceComboBox.setCurrentIndex(1)
    self.parameterLayout.addRow("Analysis phases", self.phaseChoiceComboBox)

    self.phase1PointNumber = qt.QSpinBox()
    self.phase1PointNumber.setMinimum(1)
    self.phase1PointNumber.setMaximum(100000)
    self.phase1PointNumber.setValue(200)
    self.parameterLayout.addRow("Phase 1 Points", self.phase1PointNumber)

    self.phase2PointNumber = qt.QSpinBox()
    self.phase2PointNumber.setMinimum(1)
    self.phase2PointNumber.setMaximum(1000000)
    self.phase2PointNumber.setValue(1000)
    self.parameterLayout.addRow("Phase 2 Points", self.phase2PointNumber)

    self.processingComboBox = qt.QComboBox()
    self.processingComboBox.addItem("Local Single CPU Core")
    self.processingComboBox.addItem("Local Multiple CPU Cores")
    self.processingComboBox.addItem("Cluster/Grid")
    self.parameterLayout.addRow("Processing", self.processingComboBox)

  def selectMeshFolder(self):
      self.meshFolder=qt.QFileDialog().getExistingDirectory()
      self.meshInputText.setText(self.meshFolder)
      try:
        self.loadButton.enabled = bool(self.meshFolder)
      except AttributeError:
        self.loadButton.enable = False

  def selectOutputFolder(self):
    self.outputFolder=qt.QFileDialog().getExistingDirectory()
    self.meshOutputText.setText(self.outputFolder)

  def onLoad(self):
    #print("Mocking a call to the Logic service AL001.1  with directory" + str(self.mesh_folder))
    #self.Auto3dgmData.datasetCollection=Auto3dgmLogic.createDatasetCollection(Auto3dgmLogic.createDataset(self.meshFolder),"original")
    self.Auto3dgmData.datasetCollection=Auto3dgmLogic.createDataset(self.meshFolder)
    print(self.Auto3dgmData.datasetCollection)
    try:
      self.subStepButton.enabled = bool(self.meshFolder)
    except AttributeError:
      self.subStepButton.enable = False


  def prepareOutputFolder(self):
    if (not os.path.exists(self.outputFolder+"/lowres/")):
      os.mkdir(self.outputFolder+"/lowres/")
    if (not os.path.exists(self.outputFolder+"/highres/")):
      os.mkdir(self.outputFolder+"/highres/")
    if (not os.path.exists(self.outputFolder+"/aligned/")):
      os.mkdir(self.outputFolder+"/aligned/")
    self.outputFolderPrepared=True

  ### RUN TAB WIDGETS AND BEHAVIORS ###

  def setupRunTab(self, runTabLayout):
    self.singleStepGroupBox = qt.QGroupBox("Run individual steps")
    self.singleStepGroupBoxLayout = qt.QVBoxLayout()
    self.singleStepGroupBoxLayout.setSpacing(5)
    self.singleStepGroupBox.setLayout(self.singleStepGroupBoxLayout)
    runTabLayout.addRow(self.singleStepGroupBox)

    self.subStepButton = qt.QPushButton("Subsample")
    self.subStepButton.toolTip = "Run subsample step of analysis, creating collections of subsampled points per mesh."
    self.subStepButton.connect('clicked(bool)', self.subStepButtonOnLoad)
    self.subStepButton.enabled=False
    self.singleStepGroupBoxLayout.addWidget(self.subStepButton)

    self.phase1StepButton = qt.QPushButton("Phase 1")
    self.phase1StepButton.toolTip = "Run the first analysis phase, aligning meshes based on low resolution subsampled points."
    self.phase1StepButton.connect('clicked(bool)', self.phase1StepButtonOnLoad)
    self.singleStepGroupBoxLayout.addWidget(self.phase1StepButton)

    self.phase2StepButton = qt.QPushButton("Phase 2")
    self.phase2StepButton.toolTip = "Run the second analysis phase, aligning meshes based on high resolution subsampled points."
    self.phase2StepButton.connect('clicked(bool)', self.phase2StepButtonOnLoad)
    self.singleStepGroupBoxLayout.addWidget(self.phase2StepButton)

    self.allStepsButton = qt.QPushButton("Run all steps")
    self.allStepsButton.toolTip = "Run all possible analysis steps and phases."
    self.allStepsButton.connect('clicked(bool)', self.allStepsButtonOnLoad)
    runTabLayout.addRow(self.allStepsButton)

    runTabLayout.setVerticalSpacing(15)

  def subStepButtonOnLoad(self):
    print("Mocking a call to the Logic service AL002.1")
    # Store the keys
    self.Auto3dgmData.phase1SampledPoints = self.phase1PointNumber.value
    self.Auto3dgmData.phase2SampledPoints = self.phase2PointNumber.value
    self.Auto3dgmData.fpsSeed=self.fpsSeed.value
    Auto3dgmLogic.subsample(self.Auto3dgmData,list_of_pts = [self.phase1PointNumber.value,self.phase2PointNumber.value], meshes=self.Auto3dgmData.datasetCollection.datasets[0])
    print("Dataset collection updated")
    print(self.Auto3dgmData.datasetCollection.datasets)

  def phase1StepButtonOnLoad(self):
    #corr = Auto3dgmLogic.correspondence(self, phase = 1)
    self.Auto3dgmData.datasetCollection.add_analysis_set(Auto3dgmLogic.correspondence(self.Auto3dgmData, phase=1),"Phase 1")
    print("Mocking a call to the Logic service AL002.2")

  def phase2StepButtonOnLoad(self):
    #corr = Auto3dgmLogic.correspondence(self, phase=2)
    self.Auto3dgmData.datasetCollection.add_analysis_set(Auto3dgmLogic.correspondence(self.Auto3dgmData, phase=2),"Phase 2")
    print("Mocking a call to the Logic service AL002.2")

  def allStepsButtonOnLoad(self):
    print("Mocking a call to the Logic service AL000.0")
    self.Auto3dgmData.phase1SampledPoints = self.phase1PointNumber.value
    self.Auto3dgmData.phase2SampledPoints = self.phase2PointNumber.value
    #Auto3dgmLogic.subsample(self,list_of_pts = [self.phase1PointNumber.value,self.phase2PointNumber.value], meshes=self.Auto3dgmData.datasetCollection.datasets[0])
    #self.Auto3dgmData.datasetCollection.add_analysis_set(Auto3dgmLogic.correspondence(self, phase = 1),"Phase 1")
    #self.Auto3dgmData.datasetCollection.add_analysis_set(Auto3dgmLogic.correspondence(self, phase = 2),"Phase 2")
    Auto3dgmLogic.runAll(self.Auto3dgmData)

  ### OUTPUT TAB WIDGETS AND BEHAVIORS

  def setupOutTab(self, outTabLayout):
    visualizationinputfolderWidget=ctk.ctkCollapsibleButton()
    visualizationinputfolderLayout=qt.QFormLayout(visualizationinputfolderWidget)
    visualizationinputfolderWidget.text = ""
    outTabLayout.addRow(visualizationinputfolderWidget)

    self.visualizationmeshInputText, visualizationvolumeInLabel, self.visualizationinputFolderButton=self.textIn('Input folder','Choose input folder', '')

    visualizationinputfolderLayout.addRow(visualizationvolumeInLabel)#,1,1)    
    visualizationinputfolderLayout.addRow(self.visualizationmeshInputText)#,1,2)
    visualizationinputfolderLayout.addRow(self.visualizationinputFolderButton)#,1,3)

    self.visualizationinputFolderButton.connect('clicked(bool)', self.visualizationSelectMeshFolder)

    self.visGroupBox = qt.QGroupBox("Visualize results")
    self.visGroupBoxLayout = qt.QVBoxLayout()
    self.visGroupBoxLayout.setSpacing(5)
    self.visGroupBox.setLayout(self.visGroupBoxLayout)
    outTabLayout.addRow(self.visGroupBox)

    self.visSubButton = qt.QPushButton("Subsample")
    self.visSubButton.toolTip = "Visualize collections of subsampled points per mesh."
    self.visSubButton.connect('clicked(bool)', self.visSubButtonOnLoad)
    self.visGroupBoxLayout.addWidget(self.visSubButton) 

    self.visPhase1Button = qt.QPushButton("Phase 1")
    self.visPhase1Button.toolTip = "Visualize aligned meshes with alignment based on low resolution subsampled points."
    self.visPhase1Button.connect('clicked(bool)', self.visPhase1ButtonOnLoad)
    self.visGroupBoxLayout.addWidget(self.visPhase1Button)

    self.visPhase2Button = qt.QPushButton("Phase 2")
    self.visPhase2Button.toolTip = "Visualize aligned meshes with alignment based on high resolution subsampled points."
    self.visPhase2Button.connect('clicked(bool)', self.visPhase2ButtonOnLoad)
    self.visGroupBoxLayout.addWidget(self.visPhase2Button)

    self.outGroupBox = qt.QGroupBox("Output results")
    self.outGroupBoxLayout = qt.QVBoxLayout()
    self.outGroupBoxLayout.setSpacing(5)
    self.outGroupBox.setLayout(self.outGroupBoxLayout)
    outTabLayout.addRow(self.outGroupBox)

    self.outPhase1Button = qt.QPushButton("Phase 1 results to GPA")
    self.outPhase1Button.toolTip = "Output aligned mesh results (based on low resolution subsampled points) to GPA toolkit extension for PCA and other analysis."
    self.outPhase1Button.connect('clicked(bool)', self.outPhase1ButtonOnLoad)
    self.outGroupBoxLayout.addWidget(self.outPhase1Button)

    self.outPhase2Button = qt.QPushButton("Phase 2 results to GPA")
    self.outPhase2Button.toolTip = "Output aligned mesh results (based on high resolution subsampled points)to GPA toolkit extension for PCA and other analysis."
    self.outPhase2Button.connect('clicked(bool)', self.outPhase2ButtonOnLoad)
    self.outGroupBoxLayout.addWidget(self.outPhase2Button)

    self.outVisButton = qt.QPushButton("Visualization(s)")
    self.outVisButton.toolTip = "Output all visualizations produced."
    self.outVisButton.connect('clicked(bool)', self.outVisButtonOnLoad)
    self.outGroupBoxLayout.addWidget(self.outVisButton)

    self.importAlignedButton = qt.QPushButton("Import aligned meshes")
    self.importAlignedButton.toolTip = "Imports aligned meshes to be plugged to the webviewer"
    self.importAlignedButton.connect('clicked(bool)', self.onImportAligned)
    self.outGroupBoxLayout.addWidget(self.importAlignedButton)



    outTabLayout.setVerticalSpacing(15)

  def visSubButtonOnLoad(self):
    print("Mocking a call to the Logic service AL003.3, maybe others")

  def visPhase1ButtonOnLoad(self):
    print("Mocking a call to the Logic service AL003.3, AL003.4")

  def visPhase2ButtonOnLoad(self):
    print("Mocking a call to the Logic service AL003.3, AL003.4")

  def outPhase1ButtonOnLoad(self):
    print("Mocking a call to the Logic service AL003.1")
    if self.outputFolderPrepared == False:
      self.prepareOutputFolder()
    for mesh in self.Auto3dgmData.datasetCollection.datasets[self.Auto3dgmData.phase1SampledPoints][self.Auto3dgmData.phase1SampledPoints]:
        print(self.Auto3dgmData.datasetCollection.datasets[self.Auto3dgmData.phase1SampledPoints])
        print(mesh)
        cfilename=self.outputFolder+"/lowres/"+mesh.name
        Auto3dgmLogic.saveNumpyArrayToCsv(mesh.vertices,cfilename)

  def outPhase2ButtonOnLoad(self):
    print("Mocking a call to the Logic service AL003.1")
    if self.outputFolderPrepared==False:
      self.prepareOutputFolder()
    for mesh in self.Auto3dgmData.datasetCollection.datasets[self.Auto3dgmData.phase2SampledPoints][self.Auto3dgmData.phase2SampledPoints]:
        cfilename=self.outputFolder+"/highres/"+mesh.name
        Auto3dgmLogic.saveNumpyArrayToCsv(mesh.vertices,cfilename)

  def outVisButtonOnLoad(self):
    print("Mocking a call to the Logic service unnamed visualization output service")

  def textIn(self,label, dispText, toolTip):
    """ a function to set up the appearnce of a QlineEdit widget.
    the widget is returned.
    """
    # set up text line
    textInLine=qt.QLineEdit();
    textInLine.setText(dispText)
    textInLine.toolTip = toolTip
    # set up label
    lineLabel=qt.QLabel()
    lineLabel.setText(label)

    # make clickable button
    button=qt.QPushButton(dispText)
    return textInLine, lineLabel, button

  def visualizationSelectMeshFolder(self):
    self.visualizationMeshFolder=qt.QFileDialog().getExistingDirectory()
    self.visualizationmeshInputText.setText(self.visualizationMeshFolder)
    slicer.app.coreIOManager().loadFile(self.visualizationMeshFolder+"/MSTmatrix.csv")
    slicer.app.coreIOManager().loadFile(self.visualizationMeshFolder+"/distancematrix.csv")
    print(self.visualizationMeshFolder)

  def cleanup(self):
    pass

  def onImportAligned(self):
    Auto3dgmLogic.alignOriginalMeshes(self.Auto3dgmData)
    Auto3dgmLogic.saveAlignedMeshesForViewer(self.Auto3dgmData, self.outputFolder)
    print("Meshes exported")

#
# Auto3dgmLogic
#

class Auto3dgmLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """
  def runAll(Auto3dgmData):
    Auto3dgmLogic.subsample(Auto3dgmData,[Auto3dgmData.phase1SampledPoints,Auto3dgmData.phase2SampledPoints],Auto3dgmData.datasetCollection.datasets[0])
    print("Subsampling complete.")
    Auto3dgmData.datasetCollection.add_analysis_set(Auto3dgmLogic.correspondence(Auto3dgmData, phase=1),"Phase 1")
    print("Phase 1 complete.")
    Auto3dgmData.datasetCollection.add_analysis_set(Auto3dgmLogic.correspondence(Auto3dgmData, phase=2),"Phase 2")
    print("Phase 2 complete.")

  # Logic service function AL001.001 Create dataset
  def createDataset(inputdirectory):
    dataset = DatasetFactory.ds_from_dir(inputdirectory)
    return dataset

  # Logic service function AL002.1 Subsample

  # In: List of points, possibly just one
  # list of meshes
  def subsample(Auto3dgmData,list_of_pts, meshes):
    print(list_of_pts)
    for mesh in meshes:
        print(len(mesh.vertices))
    ss = auto3dgm_nazar.mesh.subsample.Subsample(pointNumber=list_of_pts, meshes=meshes, seed={})
    for point in list_of_pts:
      names = []
      meshes = []
      for key in ss.ret[point]['output']['output']:
        mesh = ss.ret[point]['output']['output'][key]
        meshes.append(mesh)
      dataset = {}
      dataset[point] = meshes
      Auto3dgmData.datasetCollection.add_dataset(dataset,point)
    return(Auto3dgmData)

  def createDatasetCollection(dataset,name):
    datasetCollection=auto3dgm_nazar.dataset.datasetcollection.DatasetCollection(datasets = [dataset],dataset_names = [name])
    return datasetCollection

  def correspondence(Auto3dgmData, phase = 1):
    if phase == 1:
      npoints = Auto3dgmData.phase1SampledPoints
    else:
      npoints = Auto3dgmData.phase2SampledPoints
    meshes = Auto3dgmData.datasetCollection.datasets[npoints][npoints]
    corr = auto3dgm_nazar.analysis.correspondence.Correspondence(meshes=meshes)
    print("Correspondence compute for Phase " + str(phase))
    return(corr)

  def saveNumpyArrayToCsv(array,filename):
    print(array)
    print(filename)
    np.savetxt(filename+".csv",array,delimiter = ",",fmt = "%s")
    print(str(array) + " saved to file " + str(filename))

  def alignOriginalMeshes(Auto3dgmData, phase = 2):
    if 'Phase 2' in Auto3dgmData.datasetCollection.analysis_sets:
      corr = Auto3dgmData.datasetCollection.analysis_sets['Phase 2']
    elif 'Phase 1' in Auto3dgmData.datasetCollection.analysis_sets:
      corr = Auto3dgmData.datasetCollection.analysis_sets['Phase 1']
      print("Phase 2 results do not exist, computing with Phase 1")
    else:
      print("No alignment has been computed")
      return(0)
    meshes = Auto3dgmData.datasetCollection.datasets[0]
    Auto3dgmData.aligned_meshes = []
    for t in range(len(meshes)):
      R = corr.globalized_alignment['r'][t]
      aligned_mesh = meshes[t]
      aligned_mesh.rotate(arr = R)
      Auto3dgmData.aligned_meshes.append(aligned_mesh)
    return(Auto3dgmData)

  def saveAlignedMeshesForViewer(Auto3dgmData,outputFolder):
    outputdir=outputFolder+'/aligned/'
    for mesh in Auto3dgmData.aligned_meshes:
      print(outputdir)
      print(mesh.name)
      MeshExport.writeToFile(outputdir,mesh,format='obj')

class Auto3dgmTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_Auto3dgm1()

  def test_Auto3dgm1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import SampleData
    SampleData.downloadFromURL(
      nodeNames='FA',
      fileNames='FA.nrrd',
      uris='http://slicer.kitware.com/midas3/download?items=5767')
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = Auto3dgmLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
