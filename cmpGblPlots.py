#!/usr/bin/python

import sys, os, re
import argparse
from ROOT import TFile, TH1F, TF1, TCanvas, gDirectory, TGraph
#sys.path.append('pythonutils')
import utils
import compareRootHists 

def getArgs():
    parser = argparse.ArgumentParser(description='Run comparison')
    parser.add_argument('--file','-f',nargs='+',help='Input files.')
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
        print 'done comparing ', len(histos),' root histos'
    if c != None:
        c.SaveAs(cName + '.png')
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
