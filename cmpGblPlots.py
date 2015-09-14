#!/usr/bin/python

import sys, os, re
import argparse
from ROOT import TFile, TH1F, TF1, TCanvas, gPad, gDirectory, TGraph, TGraphErrors
sys.path.append('pythonutils')
import utils
import compareRootHists 
import plotutils

def getArgs():
    parser = argparse.ArgumentParser(description='Run comparison')
    parser.add_argument('--file','-f',nargs='+', required=True,help='Input files.')
    parser.add_argument('--tag','-t',default='',help='Tag to output files.')
    args = parser.parse_args();
    print args
    return args


def getFileName(f):
    print  os.path.splitext(os.path.basename(f.GetName()))[0]
    return os.path.splitext(os.path.basename(f.GetName()))[0]
    


def getSensor(name):
    m = re.match('.*_module_(.*)_sensor0.*',name)
    if m!=None:
        return m.group(1)
    else:
        print 'error get sensor from ', name


def compareSensorHists(files,histName,tag='',half=None,singleCanvas=None):
    print 'compareSensorHists for ', len(files), ' for histName ', histName, ' tag ', tag
    sensorHistNames = []
    sensors = utils.getSensorNames()
    for s in sensors:
        if utils.getHalf(s) == 't' and half != 'top':
            continue
        if utils.getHalf(s) == 'b' and half != 'bottom':
            continue
        addOn = ''
        if utils.getHalf(s) == 't':
            addOn = '_top'
        else:
            addOn = '_bot'
        sensorHistNames.append(histName + s + addOn)
    print sensorHistNames
    print 'open TFiles'
    tFiles = []
    for f in files:
        tFiles.append(TFile(f))
    print 'opened ', len(tFiles)
    c = None
    print 'singleCanvas ', singleCanvas
    if singleCanvas != None:
        cName = 'c_' + histName + '_' + half + '_' + tag
        c = TCanvas(cName, cName, 10, 10, 690*2, 390*2)
        c.Divide(2,12)
    currentPad = None

    graphBinNames = {}
    graphMean = []
    graphRMS = []
    for tf in tFiles:
        print 'tf ', tf.GetName()
        grM = TGraphErrors()
        grM.SetName('grMean_' + tf.GetName())
        grR = TGraphErrors()
        grR.SetName('grRMS_' + tf.GetName())
        graphMean.append( grM )
        graphRMS.append( grR )
    
    for hName in sensorHistNames:
        histos = []
        print 'get ', hName
        for tf in tFiles:
            print 'tf ', tf.GetName()
            h = tf.Get(hName)
            histos.append(h)
            
        print 'compare ', len(histos),' root histos'
        if c != None:
            i = utils.getCanvasIdxTwoCols(hName)
            currentPad = c.cd(i)
        print currentPad
        compareRootHists.compareHists(histos,legends=None,normalize=None,fitName='gaus',t=tag,pad=currentPad)
        print 'make graphs of mean and RMS'
        for ih in range(len(histos)):
            h = histos[ih]
            f = h.GetFunction('gaus')
            if f != None:
                mean = f.GetParameter(1)
                meanError = f.GetParError(1)
                rms = f.GetParameter(2)
                rmsError = f.GetParError(2)
                ipoint = graphMean[ih].GetN()
                graphMean[ih].SetPoint(ipoint, ipoint, mean)
                graphRMS[ih].SetPoint(ipoint, ipoint, rms)
                graphMean[ih].SetPointError(ipoint, 0., meanError)
                graphRMS[ih].SetPointError(ipoint, 0., rmsError)
                graphBinNames[ ipoint ] = utils.getShortSensorName( hName )
                #print 'mean ', mean, ' RMS ', rms
        
        print 'done comparing ', len(histos),' root histos'
    if c != None:
        c.SaveAs(cName + '.png')

    cGrName = 'c_summary' + histName + '_' + half + '_' + tag 
    cGr = TCanvas(cGrName, cGrName, 10, 10, 690*2, 390*2)
    cGr.Divide(1,2)
    grMaxValMean = -1000000.
    grMaxValRMS = -10000000.
    grMinValMean = 1000000.
    grMinValRMS = 10000000.
    
    for igr in range(len(graphMean)):
        grM = graphMean[ igr ]
        grR = graphRMS [ igr ]
        limVals = plotutils.getGraphMaxMinVal(grM)
        print limVals
        if limVals[1] > grMaxValMean:
            grMaxValMean = limVals[1]
        if limVals[0] < grMinValMean:
            grMinValMean = limVals[0]
        limVals = plotutils.getGraphMaxMinVal(grR)
        if limVals[1] > grMaxValRMS:
            grMaxValRMS = limVals[1]
        if limVals[0] < grMinValMean:
            grMinValRMS = limVals[0]

    print grMaxValMean, grMinValMean
    
    for igr in range(len(graphMean)):
        grM = graphMean[ igr ]
        grR = graphRMS [ igr ]
        #print 'igr ', igr, ' name ', grM.GetName() 
        plotutils.setGraphStyle(grM,igr+1)
        plotutils.setGraphStyle(grR,igr+1)
        plotutils.setBinLabelsDict(grM,graphBinNames)
        plotutils.setBinLabelsDict(grR,graphBinNames)
        if igr == 0:
            cGr.cd(1)
            gPad.SetBottomMargin(0.3)
            grM.GetHistogram().SetMaximum(grMaxValMean*1.2)
            grM.GetHistogram().SetMinimum(grMinValMean*1.2)
            grM.Draw('AXPL')
            cGr.cd(2)
            gPad.SetBottomMargin(0.3)
            grR.SetMaximum(grMaxValRMS)
            grR.SetMinimum(grMinValRMS)
            grR.Draw('AXPL')
        else:
            cGr.cd(1)
            grM.Draw('PL,same')
            cGr.cd(2)
            grR.Draw('PL,same')
    cGr.SaveAs(cGr.GetName() + '.png')

    #ans = raw_input('press anything')
    
        
    for tf in tFiles:        
        print 'closing ', tf.GetName()
        tf.Close()
    return

    

def compareGblHists(files,tag):

    print 'compareGblHists for ', len(files), ' files'
    compareSensorHists(files,'h_res_gbl_',tag,'top',TCanvas())
    ans = raw_input('press anything')
    compareSensorHists(files,'h_res_gbl_',tag,'bottom',TCanvas())
    ans = raw_input('press anything')


    compareSensorHists(files,'h_corrdiff_lambda_',tag,'top',TCanvas())
    ans = raw_input('press anything')
    compareSensorHists(files,'h_corrdiff_lambda_',tag,'bottom',TCanvas())
    ans = raw_input('press anything')

    compareSensorHists(files,'h_corrdiff_phi_',tag,'top',TCanvas())
    ans = raw_input('press anything')
    compareSensorHists(files,'h_corrdiff_phi_',tag,'bottom',TCanvas())
    ans = raw_input('press anything')



def main(args):

    if len(args.file)>1:

        compareGblHists(args.file,args.tag)
    



if __name__ == '__main__':

    args = getArgs()

    main(args)
