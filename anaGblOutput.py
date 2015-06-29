#!/usr/bin/python

import sys, os, re
import argparse
from ROOT import TFile, TH1F, TF1, TCanvas, gDirectory, TGraph
sys.path.append('pythonutils')
from plotutils import myText
from plotutils import  setBinLabels
import compareRootHists 

def getArgs():
    parser = argparse.ArgumentParser(description='Run comparison')
    parser.add_argument('--file','-f',required=True,help='Input file.')
    parser.add_argument('--tag','-t',default='',help='Tag to output files.')
    args = parser.parse_args();
    print args
    return args


def findFEEPeakBin(h):
    m = -1.0
    mb = -1
    for b in range(1,h.GetNbinsX()+1):
        x = h.GetBinCenter(b)
        if x>0.85:
            if h.GetBinContent(b) > m:
                m = h.GetBinContent(b)
                mb = b
    return  mb

def findPeakBin(h):
    m = -1.0
    mb = -1
    for b in range(1,h.GetNbinsX()+1):
        x = h.GetBinCenter(b)
        if h.GetBinContent(b) > m:
            m = h.GetBinContent(b)
            mb = b
    return  mb

                

def setStyle(h):
    h.SetFillColor(2)
    h.SetFillStyle(3003)
    

def getFileName(f):
    print  os.path.splitext(os.path.basename(f.GetName()))[0]
    return os.path.splitext(os.path.basename(f.GetName()))[0]
    

def fitMom(f,name):    
    c = TCanvas('c','c',10,10,700,500)
    h_p = f.Get(name)
    setStyle(h_p)
    peakbin = findFEEPeakBin(h_p)
    l = h_p.GetBinCenter(peakbin)-0.11
    h = h_p.GetBinCenter(peakbin)+0.2
    f_p = TF1('f'+'_'+name,'gaus',l,h)
    h_p.Fit(f_p,'R')
    h_p.Draw()
    myText(0.5,0.85,'<m>=%.2f #sigma=%.2f'%(f_p.GetParameter(1),f_p.GetParameter(2)), 0.05, 2)
    c.SaveAs(name + '-' + getFileName(f) + args.tag + '.png')
    ans = raw_input('continue?')



def getSensor(name):
    m = re.match('.*_module_(.*)_sensor0.*',name)
    if m!=None:
        return m.group(1)
    else:
        print 'error get sensor from ', name

def plotResiduals(f,name,half):    
    f.cd()
    
    c = TCanvas('c','c',10,10,700,500)
    c.Print('res-' +name+'-'+half + '-' + getFileName(f) + args.tag + '.ps[')
    histos = compareRootHists.getHistograms(gDirectory)
    print 'found ', len(histos), ' histograms'
    res_map = {}
    for h_p in histos:
        if re.match('^' + name + '.*' + half + '$', h_p.GetName() ) and not 'large' in h_p.GetName() and not 'truth' in h_p.GetName() :
            print 'found ' , h_p.GetName()
            sensorName = getSensor(h_p.GetName())
            print 'sensor ', sensorName
            setStyle(h_p)
            peakbin = findPeakBin(h_p)
            l = h_p.GetBinCenter(peakbin) - 5*h_p.GetRMS()
            h = h_p.GetBinCenter(peakbin) + 5*h_p.GetRMS()
            print 'peakbin ', peakbin, ' c ' , h_p.GetBinCenter(peakbin), ' l ', l, ' h', h
            f_p = TF1('f'+'_'+name+'-'+half+'_'+h_p.GetName(),'gaus',l,h)
            h_p.Fit(f_p,'R')
            c1 = TCanvas('c1-'+h_p.GetName(),'c1-'+h_p.GetName(),10,10,700,500)
    
            h_p.Draw()
            myText(0.5,0.85,'<m>=%.2f #sigma=%.2f'%(f_p.GetParameter(1),f_p.GetParameter(2)), 0.05, 2)
            res_map[sensorName] = f_p

            c1.Print('res-' +name+'-'+half + '-' + getFileName(f) + args.tag + '.ps')
            #c.SaveAs(name + '-' + getFileName(f) + args.tag + '.png')
            #ans = raw_input('continue?')

    grMean = TGraph()
    grSigma = TGraph()
    i=0
    binLabels = []
    res_map_sort = sorted(res_map.keys())
    print res_map
    print 'found ', len(res_map), ' res plots' 
    for sensor in res_map_sort:
        fit = res_map[sensor]
        grMean.SetPoint(i,i,fit.GetParameter(1))
        grSigma.SetPoint(i,i,fit.GetParameter(2))
        binLabels.append(sensor)
        i=i+1
    setBinLabels(grMean,binLabels)
    setBinLabels(grSigma,binLabels)
    grMean.SetTitle(name+'-'+half + '-' + getFileName(f) + args.tag + ';;'+name+'-'+half+' mean')
    grMean.SetMarkerStyle(20)
    c1 = TCanvas('c1-res','c1-res',10,10,700,500)
    c1.SetBottomMargin(0.4)
    grMean.Draw('ALP')
    c1.Print('res-' +name+'-'+half + '-' + getFileName(f) + args.tag + '.ps')

    grSigma.SetTitle(name+'-'+half + '-' + getFileName(f) + args.tag + ';;'+name+'-'+half+' #sigma')
    grSigma.SetMarkerStyle(20)
    c1 = TCanvas('c1-ress','c1-ress',10,10,700,500)
    c1.SetBottomMargin(0.4)
    grSigma.Draw('ALP')
    c1.Print('res-' +name+'-'+half + '-' + getFileName(f) + args.tag + '.ps')

    ans = raw_input('continue?')


    c.Print('res-' +name+'-'+half + '-' + getFileName(f) + args.tag + '.ps]')



def main(args):

    f = TFile(args.file)

    #fitMom(f,'h_p')
    #fitMom(f,'h_p_gbl')
    #fitMom(f,'h_p_top')
    #fitMom(f,'h_p_gbl_top')
    #fitMom(f,'h_p_bot')
    #fitMom(f,'h_p_gbl_bot')

    plotResiduals(f,'h_res','top')
    plotResiduals(f,'h_res','bot')
    plotResiduals(f,'h_res_gbl','top')
    plotResiduals(f,'h_res_gbl','bot')

    ans = raw_input('continue?')
    
    f.Close()



if __name__ == '__main__':

    args = getArgs()

    main(args)
