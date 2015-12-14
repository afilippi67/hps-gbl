#!/usr/bin/python

import sys, os, re
import argparse
from ROOT import TFile, TH1F, TF1, TCanvas, gDirectory, TGraph, TGraphErrors, gROOT
sys.path.append('pythonutils')
import plotutils
import hps_utils

def getArgs():
    parser = argparse.ArgumentParser(description='Run analysis')
    parser.add_argument('--file','-f',nargs=1, required=True,help='Input files.')
    parser.add_argument('--tag','-t',default='',help='Tag to output files.')
    parser.add_argument('-b','--batch', action='store_true',help='batch mode')
    parser.add_argument('-r','--regexp', help='regexp to down select')
    parser.add_argument('-u','--uflip', action='store_true', help='Flip sign of u-residuals for stereo sensors.')
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
    plotutils.myText(0.5,0.85,'<m>=%.2f #sigma=%.2f'%(f_p.GetParameter(1),f_p.GetParameter(2)), 0.05, 2)
    c.SaveAs(name + '-' + getFileName(f) + args.tag + '.png')
    ans = raw_input('continue?')





def plotResiduals(f,name,half):    
    f.cd()
    
    c = TCanvas('c','c',10,10,700,500)
    c.Print('res-' +name+'-'+half + '-' + getFileName(f) + args.tag + '.ps[')
    histos = plotutils.getHistograms(gDirectory)
    print 'found ', len(histos), ' histograms'
    res_map = {}
    for h_p in histos:
        if re.match('^' + name + '.*' + half + '$', h_p.GetName() ) and not 'large' in h_p.GetName() and not 'truth' in h_p.GetName() :
            print 'found ' , h_p.GetName()
            sensorName = utils.getSensor(h_p.GetName())
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
            plotutils.myText(0.5,0.85,'<m>=%.2f #sigma=%.2f'%(f_p.GetParameter(1),f_p.GetParameter(2)), 0.05, 2)
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
    plotutils.setBinLabels(grMean,binLabels)
    plotutils.setBinLabels(grSigma,binLabels)
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







def plotRotations(t_file):    

    names = hps_utils.getSensorNames()

    print ' got ', len(names), ' sensor names'

    #hname = 'h_res_gbl_vs_u_'
    hname = 'h_res_gbl_vs_vpred_'

    
    c = TCanvas('c','c',10,10,700,500)

    grSlopes = TGraphErrors()
    grSlopes.SetName('slopes_' + hname)
    grSlopesBinLabels = []
    for sensorname in names:

        name = hname + sensorname

        if args.regexp != None:
            m = re.match(args.regexp, name)
            if m == None:
                print 'skip this histogram \"', name, '\"'
                continue
        
        h = t_file.Get(name)
        if h == None:
            print 'no histogram \"', name, '\"'
            sys.exit(1)

        print 'process \"', name , '\"'
        c.Clear()
        h.Draw('colz')
        c.SaveAs(name + '-' + args.tag + '.png')
        #ans = raw_input('continue?')

        grMean = TGraphErrors()
        grMean.SetName( h.GetName() + '_prjYmean' )
        for b in range(1, h.GetNbinsX()+1):
            h_p = h.ProjectionY(h.GetName() + '_prjy' + str(b), b, b, 'E')
            f = None
            if h_p.GetEntries() > 77.:
                peakbin = findPeakBin(h_p)
                x =  h_p.GetBinCenter(peakbin)
                minv = x - 1.5*h_p.GetRMS()
                maxv = x + 1.5*h_p.GetRMS()
                print 'vpred ', h.GetXaxis().GetBinCenter(b) , ' peakbin ', peakbin, ' c ' , h_p.GetBinCenter(peakbin), ' minv ', minv, ' maxv ', maxv
                f = TF1('f'+'_'+name,'gaus',minv,maxv)
                h_p.Fit(f,'RQ')
                ipoint = grMean.GetN()
                grMean.SetPoint(ipoint, h.GetXaxis().GetBinCenter(b), f.GetParameter(1))
                grMean.SetPointError(ipoint, 0., f.GetParError(1))
            c.Clear()
            if f != None:
                plotutils.myText(0.5,0.85,'<m>=%.2f #sigma=%.2f'%(f.GetParameter(1),f.GetParameter(2)), 0.05, 2)
            h_p.Draw()
            c.SaveAs(h_p.GetName() + '-' + args.tag + '.png')
            #ans = raw_input('continue?')

        grMean.SetTitle(name + ';;mean')
        grMean.SetMarkerStyle(20)
        c.Clear()
        fpol1 = TF1('f' + grMean.GetName(),'pol1')
        grMean.Fit(fpol1)
        grMean.Draw('ALP')
        plotutils.myText(0.35,0.8,'slope=%.2e m=%.2e'%(fpol1.GetParameter(1),fpol1.GetParameter(0)), 0.05, 2)
        ipoint = grSlopes.GetN()
        if args.uflip and hps_utils.getAxialStereo(sensorname) == 'stereo':
            print 'flip ', sensorname, ' ', fpol1.GetParameter(1) , ' -> ',  -1.*fpol1.GetParameter(1)
            grSlopes.SetPoint( ipoint, ipoint, -1.*fpol1.GetParameter(1) )
        else:
            print 'NO flip ', sensorname, ' ', fpol1.GetParameter(1) 
            grSlopes.SetPoint( ipoint, ipoint, fpol1.GetParameter(1) )
        grSlopes.SetPointError( ipoint, 0., fpol1.GetParError(1) )
        grSlopesBinLabels.append( hps_utils.getshortsensorname( sensorname ) )
        c.SaveAs(name + '-mean-' + args.tag + '.png')
        #ans = raw_input('continue?')

    c.Clear()
    grSlopes.Draw('ALP')
    plotutils.setBinLabels( grSlopes, grSlopesBinLabels )
    c.SetBottomMargin(0.2)
    if args.uflip:
        plotutils.myText(0.35,0.8,'u-flipped', 0.05, 2)
    else:
        plotutils.myText(0.35,0.8,'NOT u-flipped', 0.05, 2)
    if args.uflip:
        c.SaveAs(grSlopes.GetName() + '-uflipped-' + args.tag + '.png')
    else:
        c.SaveAs(grSlopes.GetName() + '-' + args.tag + '.png')
    ans = raw_input('continue?')
    






def main(args):

    if len(args.file)==1:

        f = TFile(args.file[0])
        #fitMom(f,'h_p')
        #fitMom(f,'h_p_gbl')
        #fitMom(f,'h_p_top')
        #fitMom(f,'h_p_gbl_top')
        #fitMom(f,'h_p_bot')
        #fitMom(f,'h_p_gbl_bot')
        #plotResiduals(f,'h_res','top')
        #plotResiduals(f,'h_res','bot')
        #plotResiduals(f,'h_res_gbl','top')
        #plotResiduals(f,'h_res_gbl','bot')
        plotRotations(f)
        ans = raw_input('continue?')
        f.Close()





if __name__ == '__main__':

    args = getArgs()

    gROOT.SetBatch(args.batch)

    main(args)
