#!/usr/bin/python

import sys, os, re
import argparse
from ROOT import TFile, TH1F, TF1, TCanvas, gPad, gDirectory, TGraph, TGraphErrors, gROOT
sys.path.append('pythonutils')
import utils
import compareRootHists 
import plotutils
import hps_utils

def getArgs():
    parser = argparse.ArgumentParser(description='Run comparison')
    parser.add_argument('--file','-f',nargs='+', required=True,help='Input files.')
    parser.add_argument('--tag','-t',default='',help='Tag to output files.')
    parser.add_argument('--beamspot','-b',action='store_true',help='Beamspot is included.')
    parser.add_argument('--legend','-l',help='Legend')
    parser.add_argument('--batch', action='store_true',help='ROOT batch mode.')
    args = parser.parse_args();
    print args
    return args


def getFileName(f):
    print  os.path.splitext(os.path.basename(f.GetName()))[0]
    return os.path.splitext(os.path.basename(f.GetName()))[0]
    


def compareSensorHists(files,histName,tag='',half=None,singleCanvas=None,beamspot=False,legends=None):
    print 'compareSensorHists for ', len(files), ' for histName ', histName, ' tag ', tag, ' beamspot ', beamspot
    sensorHistNames = []
    sensors = hps_utils.getSensorNames(beamspot)
    for s in sensors:
        if hps_utils.getLayer > 3 and hps_utils.getHoleSlot(s) == 'slot':
            continue
        #if hps_utils.getHalf(s) != 'b' or hps_utils.getLayer(s) !=2 or hps_utils.getAxialStereo(s) != 'axial':
        #    continue
        #if hps_utils.getHalf(s) != 't' or hps_utils.getLayer(s) !=6 or hps_utils.getHoleSlot(s) != 'slot' and hps_utils.getAxialStereo(s) != 'slot':
        #    continue
        #if beamspot and hps_utils.getLayer(s)==0:
        #    if hps_utils.getHalf(s) != 'b':
        #        print 'this beamspot sensor is weird! ', s
        #        sys.exit(1)
        if hps_utils.getHalf(s) == 't' and half != 'top':
            continue
        if hps_utils.getHalf(s) == 'b' and half != 'bottom':
            continue
        addOn = ''
        if hps_utils.getHalf(s) == 't':
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
            i = hps_utils.getCanvasIdxTwoCols(hName,beamspot)
            currentPad = c.cd(i)
        print currentPad
        compareRootHists.compareHists(histos,legends=legends,normalize=True,fitName={'name':'gaus','rms':2},t=tag,pad=currentPad,myTextSize=0.1)
        print 'make graphs of mean and RMS'
        for ih in range(len(histos)):
            h = histos[ih]
            f = h.GetFunction('fg_' + h.GetName())
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
                graphBinNames[ ipoint ] = hps_utils.getshortsensorname( hName )
                print 'mean ', mean, ' RMS ', rms, ' fg ', f.GetName()
            else:
                print 'No fg in histo ', ih, ' name ', h.GetName()
        
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

    if legends:
        styles = [ 'PL' for x in range( len( graphMean ) ) ] 
        l = plotutils.getLegendList(0.13,0.75,0.25,0.85,graphMean,legends,styles)
        l.Draw()
    

    cGr.SaveAs(cGr.GetName() + '.png')

    #ans = raw_input('press anything')
    
        
    for tf in tFiles:        
        print 'closing ', tf.GetName()
        tf.Close()
    return




def compareHists(files,histName,tag='',fitFunctionName=None, legends=None):
    print 'compareHists for ', len(files), ' for histName ', histName, ' tag ', tag
    print 'open TFiles'
    tFiles = []
    for f in files:
        tFiles.append(TFile(f))
    print 'opened ', len(tFiles)

    histos = []
    for tf in tFiles:
        print 'tf ', tf.GetName()
        h = tf.Get(histName)
        histos.append(h)
            
    print 'compare ', len(histos),' root histos'

    compareRootHists.compareHists(histos,legends=legends,normalize=True,fitName=fitFunctionName,t=tag,pad=None,myTextSize=0.05)

    for tf in tFiles:        
        print 'closing ', tf.GetName()
        tf.Close()
    return

def compareTwoHists(file1,histName1,file2,histName2,tag='',fitFunctionName=None,legs=None):
    print 'open TFiles'
    tFiles = [ TFile(file1), TFile(file2)]
    print 'opened ', len(tFiles), ' TFiles'
    histos = [tFiles[0].Get(histName1), tFiles[1].Get(histName2)]
    print 'compare ', len(histos),' root histos '

    compareRootHists.compareHists(histos,legends=legs,normalize=True,fitName=fitFunctionName,t=tag,pad=None,myTextSize=0.05)

    for tf in tFiles:        
        print 'closing ', tf.GetName()
        tf.Close()
    return



    

def compareGblHists(files,tag,beamspot=False,legends=None):

    print 'compareGblHists for ', len(files), ' files'

    compareSensorHists(files,'h_res_gbl_',tag,'bottom',None,beamspot,legends=legends)

    compareSensorHists(files,'h_res_gbl_',tag,'top',None,beamspot,legends=legends)
    compareSensorHists(files,'h_res_gbl_',tag,'bottom',None,beamspot,legends=legends)
    compareSensorHists(files,'h_corrdiff_lambda_',tag,'top',None,beamspot,legends=legends)
    compareSensorHists(files,'h_corrdiff_lambda_',tag,'bottom',None,beamspot,legends=legends)
    compareSensorHists(files,'h_corrdiff_phi_',tag,'top',None,beamspot,legends=legends)
    compareSensorHists(files,'h_corrdiff_phi_',tag,'bottom',None,beamspot,legends=legends)

    compareSensorHists(files,'h_res_gbl_',tag,'top',TCanvas(),beamspot,legends=legends)
    compareSensorHists(files,'h_res_gbl_',tag,'bottom',TCanvas(),beamspot,legends=legends)
    compareSensorHists(files,'h_corrdiff_lambda_',tag,'top',TCanvas(),beamspot,legends=legends)
    compareSensorHists(files,'h_corrdiff_lambda_',tag,'bottom',TCanvas(),beamspot,legends=legends)
    compareSensorHists(files,'h_corrdiff_phi_',tag,'top',TCanvas(),beamspot,legends=legends)
    compareSensorHists(files,'h_corrdiff_phi_',tag,'bottom',TCanvas(),beamspot,legends=legends)





def main(args):

    if len(args.file)>1:
        l = None
        if args.legend != None: l = args.legend.split()

        compareGblHists(args.file,args.tag,args.beamspot,l)

        compareHists(args.file,'h_d0_initial_top',args.tag,{'name':'gaus','peakbin':1,'rms':2},l)
        compareHists(args.file,'h_d0_initial_bot',args.tag,{'name':'gaus','peakbin':1,'rms':2},l)
        compareHists(args.file,'h_d0_gbl_top',args.tag,{'name':'gaus','peakbin':1,'rms':2},l)
        compareHists(args.file,'h_d0_gbl_bot',args.tag,{'name':'gaus','peakbin':1,'rms':2},l)


        compareHists(args.file,'h_z0_initial_top',args.tag,{'name':'gaus','peakbin':1,'rms':2},l)
        compareHists(args.file,'h_z0_initial_bot',args.tag,{'name':'gaus','peakbin':1,'rms':2},l)
        compareHists(args.file,'h_z0_gbl_top',args.tag,{'name':'gaus','peakbin':1,'rms':2},l)
        compareHists(args.file,'h_z0_gbl_bot',args.tag,{'name':'gaus','peakbin':1,'rms':2},l)

        compareHists(args.file,'h_p_top',args.tag,{'name':'gaus','peakbin':1,'xmin':0.95,'xmax':1.14},l)
        compareHists(args.file,'h_p_bot',args.tag,{'name':'gaus','peakbin':1,'xmin':0.95,'xmax':1.14},l)
        compareHists(args.file,'h_p_gbl_top',args.tag,{'name':'gaus','peakbin':1,'xmin':0.95,'xmax':1.14},l)
        compareHists(args.file,'h_p_gbl_bot',args.tag,{'name':'gaus','peakbin':1,'xmin':0.95,'xmax':1.14},l)

        compareHists(args.file,'h_chi2_initial_top',args.tag,None,l)
        compareHists(args.file,'h_chi2_initial_bot',args.tag,None,l)
        compareHists(args.file,'h_chi2_top',args.tag,None,l)
        compareHists(args.file,'h_chi2_bot',args.tag,None,l)

        compareHists(args.file,'h_chi2prob_initial_top',args.tag,None,l)
        compareHists(args.file,'h_chi2prob_initial_bot',args.tag,None,l)
        compareHists(args.file,'h_chi2prob_top',args.tag,None,l)
        compareHists(args.file,'h_chi2prob_bot',args.tag,None,l)


        ifile = 0
        for f in args.file:
            compareTwoHists(f,'h_d0_initial_top', f,'h_d0_initial_bot',args.tag + '-top_vs_bottom-file-' + str(ifile),{'name':'gaus','rms':2},['top','bot'])
            compareTwoHists(f,'h_z0_initial_top', f,'h_z0_initial_bot',args.tag + '-top_vs_bottom-file-' + str(ifile),{'name':'gaus','rms':2},['top','bot'])
            compareTwoHists(f,'h_d0_gbl_top', f,'h_d0_gbl_bot',args.tag + '-top_vs_bottom-file-' + str(ifile),{'name':'gaus','rms':2},['top','bot'])
            compareTwoHists(f,'h_z0_gbl_top', f,'h_z0_gbl_bot',args.tag + '-top_vs_bottom-file-' + str(ifile),{'name':'gaus','rms':2},['top','bot'])
            compareTwoHists(f,'h_p_top', f,'h_p_bot',args.tag + '-top_vs_bottom-file-' + str(ifile),{'name':'gaus','peakbin':1,'xmin':0.95,'xmax':1.14},['top','bot'])
            compareTwoHists(f,'h_p_gbl_top', f,'h_p_gbl_bot',args.tag + '-top_vs_bottom-file-' + str(ifile),{'name':'gaus','peakbin':1,'xmin':0.95,'xmax':1.14},['top','bot'])
            compareTwoHists(f,'h_chi2prob_initial_top', f,'h_chi2prob_initial_bot',args.tag + '-top_vs_bottom-file-' + str(ifile),None,['top','bot'])
            compareTwoHists(f,'h_chi2_initial_top', f,'h_chi2_initial_bot',args.tag + '-top_vs_bottom-file-' + str(ifile),None,['top','bot'])
            compareTwoHists(f,'h_chi2_top', f,'h_chi2_bot',args.tag + '-top_vs_bottom-file-' + str(ifile),None,['top','bot'])
            compareTwoHists(f,'h_chi2prob_top', f,'h_chi2prob_bot',args.tag + '-top_vs_bottom-file-' + str(ifile),None,['top','bot'])
            ifile = ifile + 1

        


if __name__ == '__main__':

    args = getArgs()

    gROOT.SetBatch( args.batch )

    main(args)
