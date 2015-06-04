#!/usr/bin/python

import sys
import argparse
from ROOT import TFile, gDirectory, TIter, TCanvas

def getHistograms(direc):
    histos = []
    print direc.GetList().GetSize()
    iter = TIter(direc.GetListOfKeys())
    print iter.GetCollection().GetSize(), ' items in collection'
    while True:
        key = iter.Next()
        if not key:
            break
        obj = key.ReadObj()
        if obj.InheritsFrom('TH1') and not obj.InheritsFrom('TH2') and not obj.InheritsFrom('TH3'):
            print 'Getting TH1 %s' % obj.GetName()
            histos.append(obj)
    print 'Got ', len(histos), ' TH1s'
    return histos

def main(args):
    print 'hej'

    f1 = TFile(args.files[0])
    histos1 = getHistograms(gDirectory)

    f2 = TFile(args.files[1])
    histos2 = getHistograms(gDirectory)

    n = 0
    nm = 0
    for h1 in histos1:
        if args.skip != None:
            if args.skip in h1.GetName():
                continue
        if args.only != None:
            if args.only not in h1.GetName():
                continue
        for h2 in histos2:
            #if h1.GetName()==h2.GetName():            
            #if (h1.GetName()+'_bot')==h2.GetName():            
            if (h1.GetName()+'_top')==h2.GetName():            
                if args.skipempty:
                    if h1.GetEntries()==0.0 or h2.GetEntries()==0.0 or h1.Integral()==0.0 or h2.Integral()==0.0:
                        continue
                print 'match ', h1.GetName(), ' ',  h1.GetEntries(), ' ', h2.GetEntries()
                c = TCanvas('c_'+h1.GetName(),'c_'+h1.GetName(),10,10,700,500)
                h1.SetLineColor(4)
                h1.SetFillColor(4)
                h1.SetFillStyle(3004)
                #h1.DrawNormalized("hist",1.0/h1.Integral(-1,999999))
                h1.DrawNormalized("hist",1.0/h1.Integral())
                h2.SetFillStyle(3005)
                h2.SetFillColor(2)
                h2.SetLineColor(2)
                #h2.Scale(1.0/h2.GetEntries())
                #h2.Scale(h1.GetEntries())
                h2.DrawNormalized("same,hist",1.0/h2.Integral())
                c.SaveAs(c.GetName()+'.png')
                n=n+1
                if args.pause:
                    ans = raw_input('press anywhere to continue')
            else:
                #print 'NO match ', h1.GetName()
                nm=nm+1
    print 'compared ', n, ' missed ', nm
    

    f1.Close()
    f2.Close()
    return 0



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run comparison.')
    parser.add_argument('--files',nargs=2, help='Input files.')
    parser.add_argument('--skipempty',action='store_true', help='Input files.')
    parser.add_argument('--skip', help='Skip histogram containing this text.')
    parser.add_argument('--only', help='Use only histograms containing this text..')
    parser.add_argument('--pause', action='store_true',help='Pause b/w histograms')
    args = parser.parse_args();
    print args
    
    main(args)
