#!/usr/bin/env python
# coding: utf-8

# # Interactive Unit Test

import ROOT
from bokeh.io import output_notebook
from RootInteractive.Tools.aliTreePlayer import *
from RootInteractive.InteractiveDrawing.bokeh.bokehDrawSA import *
from bokeh.io import curdoc
import numpy as np
import os
import sys
from ROOT import TFile, gSystem
#import plotly.express as px
#output_notebook()
from IPython.display import Image
from IPython.core.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import seaborn as sns
import matplotlib.pyplot as plt
from RootInteractive.Tools.compressArray import arrayCompressionRelative16,arrayCompressionRelative8

import pprint
pp = pprint.PrettyPrinter(indent=4)

from RootInteractive.InteractiveDrawing.bokeh.bokehInteractiveParameters import *
#pp.pprint(figureParameters["legend"])

print("Reading tree...")

# Load some macros (to reed tree and to set some aliases)
ROOT.gROOT.LoadMacro("readTree.C")
ROOT.gROOT.LoadMacro("kinematicFunctions.C")

readTreeStatus = ROOT.readTree();
if readTreeStatus != 0:
    sys.exit("readTree() has been failed with error code " + str(readTreeStatus) + ". Exit.")

#ROOT.treeMC.SetAlias("nTracks", "nTracks")
ROOT.treeMC.SetAlias("initialPT","PT(initialPX, initialPY)")
ROOT.treeMC.SetAlias("finalPT","PT(finalPX, finalPY)")
ROOT.treeMC.SetAlias("initialP","P(initialPX, initialPY, initialPZ)")
ROOT.treeMC.SetAlias("finalP","P(finalPX, finalPY, finalPZ)")
ROOT.treeMC.SetAlias("initialE","E(initialPX, initialPY, initialPZ, mass)")
ROOT.treeMC.SetAlias("finalE","E(finalPX, finalPY, finalPZ, mass)")
ROOT.treeMC.SetAlias("initialEta","Eta(initialPX, initialPY, initialPZ)")
ROOT.treeMC.SetAlias("finalEta","Eta(finalPX, finalPY, finalPZ)")
ROOT.treeMC.SetAlias("initialPhi","Phi(initialPX, initialPY)")
ROOT.treeMC.SetAlias("finalPhi","Phi(finalPX, finalPY)")
ROOT.treeMC.SetAlias("initialTheta","Theta(initialPX, initialPY, initialPZ)")
ROOT.treeMC.SetAlias("finalTheta","Theta(finalPX, finalPY, finalPZ)")

nEntriesMC=ROOT.treeMC.GetEntries()

varListMC = [
    'nTracks', 'trackID', 'parentID', 'pdg', 'charge', 'mass', 'initialX', 'initialY', 'initialZ',
    'finalX', 'finalY', 'finalZ', 'trackLength', 'initialPX', 'initialPY', 'initialPZ',
    'finalPX', 'finalPY', 'finalPZ', 'initialGlobalTime', 'finalGlobalTime', 
    'initialPT', 'finalPT', 'initialP', 'finalP', 'initialE', 'finalE',
    'initialEta', 'finalEta', 'initialPhi', 'finalPhi', 'initialTheta', 'finalTheta'
]

columnMask=[["particle__particle_f","p_"]] #????????
particlesMC = tree2Panda(ROOT.treeMC, varListMC, "1", exclude=[], columnMask=columnMask,nEntries=nEntriesMC) #?????

pd.set_option('display.max_columns', None)

nEntriesKFAR=ROOT.treeKFAR.GetEntries()

varListKFAR = [
    'nTracks', 'numOfMCEvent', 'isAllCovMatOK', 'constructMethod',
    'ID', 'IDMC', 'pdgID', 'isCovMatOK', 'massConstraint', 'topoConstraint',
    'x', 'y', 'z', 'pX', 'pY', 'pZ', 'e', 's', 'charge', 'chi2',
    'ndf', 'massHypo', 'isMother', 'isDaughter', 'p', 'pT', 'eta',
    'phi', 'mass', 'decayLength', 'decayLengthXY', 'lifeTime',
    'r', 'rapidity', 'theta', 'xErr', 'yErr', 'zErr', 'pXErr',
    'pYErr', 'pZErr', 'eErr', 'sErr', 'pErr', 'pTErr', 'etaErr',
    'phiErr', 'massErr', 'decayLengthErr', 'decayLengthXYErr',
    'lifeTimeErr', 'rErr',
]

columnMask=[["particle__particle_f","p_"]] #????????
particlesKFAR = tree2Panda(ROOT.treeKFAR, varListKFAR, "1", exclude=[], columnMask=columnMask,nEntries=nEntriesKFAR) #?????

#import pandas as pd
particles = particlesMC.join(particlesKFAR,lsuffix="_MC", rsuffix="_KFAR")
initMetadata(particles)


particles['Xres'] = particles['x'] - np.where(particles['isMother']==1.0, particles['finalX'], particles['initialX'])
particles['Yres'] = particles['y'] - np.where(particles['isMother']==1.0, particles['finalY'], particles['initialY'])
particles['Zres'] = particles['z'] - np.where(particles['isMother']==1.0, particles['finalZ'], particles['initialZ'])

particles['PXres'] = particles['pX'] - np.where(particles['isMother']==1.0, particles['finalPX'], particles['initialPX'])
particles['PYres'] = particles['pY'] - np.where(particles['isMother']==1.0, particles['finalPY'], particles['initialPY'])
particles['PZres'] = particles['pZ'] - np.where(particles['isMother']==1.0, particles['finalPZ'], particles['initialPZ'])

particles['Massres'] = particles['mass_KFAR'] - particles['mass_MC']
particles['PTres'] = particles['pT'] - np.where(particles['isMother']==1.0, particles['initialPT'], particles['finalPT'])
particles['Eres'] = particles['e'] - np.where(particles['isMother']==1.0, particles['finalE'], particles['initialE'])

#-----------

# Pulls

particles['Xpull'] = particles['Xres'] / particles['xErr']
particles['Ypull'] = particles['Yres'] / particles['yErr']
particles['Zpull'] = particles['Zres'] / particles['zErr']

particles['PXpull'] = particles['PXres'] / particles['pXErr']
particles['PYpull'] = particles['PYres'] / particles['pYErr']
particles['PZpull'] = particles['PZres'] / particles['pZErr']

particles['Masspull'] = particles['Massres'] / particles['massErr']
particles['PTpull'] = particles['PTres'] / particles['pTErr']
particles['Epull'] = particles['Eres'] / particles['eErr']

#----------------

particles["OneOverPT"] = 1/np.where(particles['isMother']==1.0, particles['initialPT'], particles['finalPT'])

particles["indexMI"]=particles.index


#AxisTitles
#particles.meta.metaData = {'Xres.AxisTitle': "Xres (cm)", }


# Create Rootinteractive interface

print("Creating RootInteractive interface...")
# statistic parameterization
figureParameters["StatParam"]={}
figureParameters["StatParam"]["parameterArray"]=[
        {"name": "nPointsRender", "range":[10, 5000], "value": 350},
]

figureParameters["StatParam"]["widgets"]=[
    ['slider', ['nPointsRender'], {"name": "nPointsRender"}],
#    ['range', ['indexMI'], {"index": True}],
]
figureParameters["StatParam"]["widgetLayout"]=[
    ["nPointsRender"]
]

parameterArray = [
    {"name": "paramX", "value":10, "range": [-20, 20]},
    {"name": "C_cut", "value": 1, "range": [0, 1]},
    {"name": "nPoints", "range":[0, 1200], "value": 1000},
]
parameterArray.extend(figureParameters["legend"]['parameterArray'])   # add legend configuration as in the default 
parameterArray.extend(figureParameters["markers"]['parameterArray'])   # add legend configuration as in the default 
parameterArray.extend(figureParameters["StatParam"]['parameterArray'])   # add legend configuration as in the default 
#
widgetParams=[
              ['range', ['finalPT']], #0
              ['range', ['OneOverPT']], #1
              ['range', ['finalEta']], #2
              ['range', ['finalPhi']], #3
              #
              ['multiSelect', ['pdg']], #4
              ['multiSelect', ['isMother']], #5
              ['multiSelect', ['isAllCovMatOK']], #6
              ['multiSelect', ['constructMethod']], #7
              ['multiSelect', ['massConstraint']], #8
              ['multiSelect', ['topoConstraint']], #9
]

widgetParams.extend(figureParameters["legend"]["widgets"])
widgetParams.extend(figureParameters["markers"]["widgets"])
widgetParams.extend(figureParameters["StatParam"]["widgets"])

widgetLayoutDesc0=[[0,1,2,3],[4,5,6,7,8,9], {'sizing_mode': 'scale_width'}] 
#
widgetLayoutDesc={
    "Data selection":widgetLayoutDesc0,
    "Legend": figureParameters['legend']['widgetLayout'],
    "Markers":["markerSize"],
    "Stat":figureParameters['StatParam']['widgetLayout'],
}

#figureGlobalOption={"size":"size"}
figureGlobalOption={}
figureGlobalOption=figureParameters["legend"]["figureOptions"]
figureGlobalOption["size"]="markerSize"

histoArray=[
    # QA histos
    {"name": "X_KFAR", "variables": ["x"], "nbins":100,"range": [-8, 8],  },#
     {"name": "Y_KFAR", "variables": ["y"], "nbins":100,"range": [-8, 8],  },#
     {"name": "Z_KFAR", "variables": ["z"], "nbins":100,"range": [-30, 30],  },#
    
    {"name": "PX_KFAR", "variables": ["pX"], "nbins":100,"range": [-8, 8],  },#
     {"name": "PY_KFAR", "variables": ["pY"], "nbins":100,"range": [-8, 8],  },#
     {"name": "PZ_KFAR", "variables": ["pZ"], "nbins":100,"range": [-8, 8],  },#
    
    {"name": "Mass_KFAR", "variables": ["mass_KFAR"], "nbins":100,"range": [0, 5],  },#
     {"name": "PT_KFAR", "variables": ["pT"], "nbins":100,"range": [0, 10],  },#
     {"name": "E_KFAR", "variables": ["e"], "nbins":100,"range": [0, 15],  },#
    
    # Residuals
     {"name": "Xres_hist", "variables": ["Xres"], "nbins":100,"range": [-0.1, 0.1],  },#
     {"name": "Yres_hist", "variables": ["Yres"], "nbins":100,"range": [-0.1, 0.1],  },#
     {"name": "Zres_hist", "variables": ["Zres"], "nbins":100,"range": [-0.1, 0.1],  },#
    
    {"name": "PXres_hist", "variables": ["PXres"], "nbins":100,"range": [-0.1, 0.1],  },#
     {"name": "PYres_hist", "variables": ["PYres"], "nbins":100,"range": [-0.1, 0.1],  },#
     {"name": "PZres_hist", "variables": ["PZres"], "nbins":100,"range": [-0.1, 0.1],  },#
    
    {"name": "Massres_hist", "variables": ["Massres"], "nbins":100,"range": [-0.1, 0.1],  },#
     {"name": "PTres_hist", "variables": ["PTres"], "nbins":100,"range": [-0.1, 0.1],  },#
     {"name": "Eres_hist", "variables": ["Eres"], "nbins":100,"range": [-0.1, 0.1],  },#
    
    # Pulls
     {"name": "Xpull_hist", "variables": ["Xpull"], "nbins":100,"range": [-8, 8],  },#
     {"name": "Ypull_hist", "variables": ["Ypull"], "nbins":100,"range": [-8, 8],  },#
     {"name": "Zpull_hist", "variables": ["Zpull"], "nbins":100,"range": [-8, 8],  },#
    
    {"name": "PXpull_hist", "variables": ["PXpull"], "nbins":100,"range": [-8, 8],  },#
     {"name": "PYpull_hist", "variables": ["PYpull"], "nbins":100,"range": [-8, 8],  },#
     {"name": "PZpull_hist", "variables": ["PZpull"], "nbins":100,"range": [-8, 8],  },#
    
    {"name": "Masspull_hist", "variables": ["Masspull"], "nbins":100,"range": [-5, 5],  },#
     {"name": "PTpull_hist", "variables": ["PTpull"], "nbins":100,"range": [-8, 8],  },#
     {"name": "Epull_hist", "variables": ["Epull"], "nbins":100,"range": [-8, 8],  },#
    
]

jsFunctionArray=[]
aliasArray=[]


tooltips = [("initialPT", "(@initialPT)")]

figureArray = [
    # residuals
    [['initialPT'], ['Xres'], { "colorZvar": "pdg","xAxisTitle":"pT (GeV)", "yAxisTitle":"Xres (cm)",}], #0
    [['initialPT'], ['Yres'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"Yres (cm)",}], #1
    [['initialPT'], ['Zres'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"Zres (cm)",}], #2
    
    [['initialPT'], ['PXres'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"PXres (GeV)",}], #3
    [['initialPT'], ['PYres'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"PYres (GeV)",}], #4
    [['initialPT'], ['PZres'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"PZres (GeV)",}], #5
    
    [['initialPT'], ['Massres'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"Massres (GeV)",}], #6
    [['initialPT'], ['PTres'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"PTres (GeV)",}], #7
    [['initialPT'], ['Eres'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"Eres (GeV)",}], #8
    
    # pulls
    [['initialPT'], ['Xpull'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"Xpull",}], #9
    [['initialPT'], ['Ypull'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"Ypull",}], #10
    [['initialPT'], ['Zpull'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"Zpull",}], #11
    
    [['initialPT'], ['PXpull'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"PXpull",}], #12
    [['initialPT'], ['PYpull'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"PYpull",}], #13
    [['initialPT'], ['PZpull'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"PZpull",}], #14
    
    [['initialPT'], ['Masspull'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"Masspull",}], #15
    [['initialPT'], ['PTpull'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"PTpull",}], #16
    [['initialPT'], ['Epull'], { "colorZvar": "pdg", "xAxisTitle":"pT (GeV)", "yAxisTitle":"Epull",}], #17
    
    # QA KFAR histos
    [['X_KFAR'], ['X_KFAR'], {"show_histogram_error": True}], #18
    [['Y_KFAR'], ['Y_KFAR'], {"show_histogram_error": True}], #19
    [['Z_KFAR'], ['Z_KFAR'], {"show_histogram_error": True}], #20
    
    [['PX_KFAR'], ['PX_KFAR'], {"show_histogram_error": True}], #21
    [['PY_KFAR'], ['PY_KFAR'], {"show_histogram_error": True}], #22
    [['PZ_KFAR'], ['PZ_KFAR'], {"show_histogram_error": True}], #23
    
    [['Mass_KFAR'], ['Mass_KFAR'], {"show_histogram_error": True}], #24
    [['PT_KFAR'], ['PT_KFAR'], {"show_histogram_error": True}], #25
    [['E_KFAR'], ['E_KFAR'], {"show_histogram_error": True}], #26
    
#     # Residuals histos
     [['Xres'], ['Xres_hist'], {"show_histogram_error": True}], #27
    [['Yres'], ['Yres_hist'], {"show_histogram_error": True}], #28
    [['Zres'], ['Zres_hist'], {"show_histogram_error": True}], #29
    
    [['PXres'], ['PXres_hist'], {"show_histogram_error": True}], #30
    [['PYres'], ['PYres_hist'], {"show_histogram_error": True}], #31
    [['PZres'], ['PZres_hist'], {"show_histogram_error": True}], #32
    
    [['Massres'], ['Massres_hist'], {"show_histogram_error": True}], #33
    [['PTres'], ['PTres_hist'], {"show_histogram_error": True}], #34
    [['Eres'], ['Eres_hist'], {"show_histogram_error": True}], #35
    
     # Pulls histos
    [['Xpull'], ['Xpull_hist'], {"show_histogram_error": True}], #36
    [['Ypull'], ['Ypull_hist'], {"show_histogram_error": True}], #37
    [['Zpull'], ['Zpull_hist'], {"show_histogram_error": True}], #38
    
    [['PXpull'], ['PXpull_hist'], {"show_histogram_error": True}], #39
    [['PYpull'], ['PYpull_hist'], {"show_histogram_error": True}], #40
    [['PZpull'], ['PZpull_hist'], {"show_histogram_error": True}], #41
    
    [['Masspull'], ['Masspull_hist'], {"show_histogram_error": True}], #42
    [['PTpull'], ['PTpull_hist'], {"show_histogram_error": True}], #43
    [['Epull'], ['Epull_hist'], {"show_histogram_error": True}], #44
    
    # Histo tables
    ["tableHisto", {"rowwise": False}], #45
    

    
]

figureLayoutDesc={
    "Residuals vs PT": [
        [0, 1, 2,   {  'plot_height': 250}],
        [3, 4, 5,   {  'plot_height': 250}],
        [6, 7, 8,   {  'plot_height': 250}],
        {'plot_height': 200, 'sizing_mode': 'scale_width'} 
    ],
    "Pulls vs PT": [
        [9, 10, 11,   {  'plot_height': 250}],
        [12, 13, 14,   {  'plot_height': 250}],
        [15, 16, 17,   {  'plot_height': 250}],
        {'plot_height': 200, 'sizing_mode': 'scale_width'} 
    ],
    "QA histograms": [
        [18,19,20,   {  'plot_height': 250}],
        [21,22,23,   {  'plot_height': 250}],
        [24,25,26,   {  'plot_height': 250}],        
        {'plot_height': 200, 'sizing_mode': 'scale_width'} 
    ],
     "Residuals histograms": [
        [27,28,29,   {  'plot_height': 250}],
        [30,31,32,   {  'plot_height': 250}],
        [33,34,35,   {  'plot_height': 250}],        
        {'plot_height': 200, 'sizing_mode': 'scale_width'} 
    ],
    "Pulls histograms": [
        [36,37,38,   {  'plot_height': 250}],
        [39,40,41,   {  'plot_height': 250}],
        [42,43,44,   {  'plot_height': 250}],        
        {'plot_height': 200, 'sizing_mode': 'scale_width'} 
    ],
    "Histogram statistics": [
     [45],
        {'plot_height': 200, 'sizing_mode': 'scale_width'} 
    ]
    
}


#%pdb

outputFileName = "InteractiveTest.html"
output_file(outputFileName)
fig=bokehDrawSA.fromArray(particles, "finalPT>-1", figureArray, widgetParams, layout=figureLayoutDesc, tooltips=tooltips, parameterArray=parameterArray,
                          widgetLayout=widgetLayoutDesc, sizing_mode="scale_width", nPointRender="nPointsRender", jsFunctionArray=jsFunctionArray,
                           aliasArray=aliasArray,
                          histogramArray=histoArray,
                          arrayCompression=arrayCompressionRelative16)

#testBokehClientHistogramRowwiseTable()



print("Interactive test is created and saved in "+outputFileName)


