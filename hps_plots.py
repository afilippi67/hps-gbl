import re, sys, subprocess
from ROOT import TH1F, TH2F, TGraph, TGraphErrors, TCanvas, TLegend, TLatex, gStyle, gDirectory, TIter, TFile, gPad
from ROOT import Double as ROOTDouble
from math import sqrt
sys.path.append('pythonutils')
from plotutils import myText
from utils import getLayer, getHalf, getAxialStereo,getHoleSlot,setGraphXLabels,getCanvasIdxTwoCols



class plotter:
    def __init__(self,tag,picExt,testRunFlag,isTop=False,isBot=True):
        self.picExt = picExt
        self.testRun = testRunFlag
        self.isTop = isTop
        self.isBot = isBot
        if not self.isTop and not self.isBot:
            print 'Cannot create a plotter with no half!?'
            sys.exit(1)
        if self.isTop and self.isBot:
            self.halftag = ''
        elif self.isTop:
            self.halftag = '_top'
        else:
            self.halftag = '_bot'
        self.doDetailed = False
        self.tag = tag
        self.makePlots()
    
    def getTag(self):
        str = ''
        if self.tag != '':
            str += '-'+self.tag
        if self.isTop:
            str += '-top'
        if self.isBot:
            str += '-bot'
        return str
    
            

    def makePlots(self):
        self.h_chi2_gbl_truth = TH1F('h_chi2_gbl_truth'+self.halftag,';GBL truth #chi^{2};Entries',50,0.,50.)
        self.h_chi2ndf_gbl_truth = TH1F('h_chi2ndf_gbl_truth'+self.halftag,';GBL truth #chi^{2}/ndf;Entries',50,0.,20.)
        self.h_chi2prob_gbl_truth = TH1F('h_chi2prob_gbl_truth'+self.halftag,';GBL truth #chi^{2} prob.;Entries',50,0.,1.)
        self.h_chi2 = TH1F('h_chi2'+self.halftag,';GBL #chi^{2};Entries',50,0.,50.)
        self.h_chi2ndf = TH1F('h_chi2ndf'+self.halftag,';GBL #chi^{2}/ndf;Entries',50,0.,8.)
        self.h_chi2prob = TH1F('h_chi2prob'+self.halftag,';GBL #chi^{2} prob;Entries',50,0.,1.)
        self.h_chi2_initial = TH1F('h_chi2_initial'+self.halftag,';Track #chi^{2};Entries',50,0.,50.)
        self.h_chi2ndf_initial = TH1F('h_chi2ndf_initial'+self.halftag,';Track #chi^{2}/ndf;Entries',50,0.,8.)
        self.h_chi2prob_initial = TH1F('h_chi2prob_initial'+self.halftag,';Track #chi^{2} prob;Entries',50,0.,1.)
        self.h_chi2_initial_truth = TH1F('h_chi2_initial_truth'+self.halftag,';Track #chi^{2};Entries',50,0.,50.)
        self.h_chi2ndf_initial_truth = TH1F('h_chi2ndf_initial_truth'+self.halftag,';Track #chi^{2}/ndf;Entries',50,0.,8.)
        self.h_chi2prob_initial_truth = TH1F('h_chi2prob_initial_truth'+self.halftag,';Track #chi^{2} prob.;Entries',50,0.,1.)
        self.h_p = TH1F('h_p'+self.halftag,';Track momentum;Entries',50,0.,1.5)
        self.h_p_gbl = TH1F('h_p_gbl'+self.halftag,';Track momentum;Entries',50,0.,1.5)
        self.h_p_truth = TH1F('h_p_truth'+self.halftag,';Track momentum;Entries',50,0.,1.5)
        self.h_p_truth_res = TH1F('h_p_truth_res'+self.halftag,';Track - Truth momentum;Entries',50,-0.3,0.3)
        self.h_p_truth_res_vs_p = TH2F('h_p_truth_res_vs_p'+self.halftag,';Truth momentum;Track - Truth momentum',6,0.4,1.6,50,-0.3,0.3)
        self.h_p_truth_res_gbl = TH1F('h_p_truth_res_gbl'+self.halftag,';Track - Truth momentum;Entries',50,-0.3,0.3)
        self.h_p_truth_res_gbl_vs_p = TH2F('h_p_truth_res_gbl_vs_p'+self.halftag,';Truth momentum;Track - Truth momentum',6,0.4,1.6,50,-0.3,0.3)
        self.h_qOverP_truth_res_gbl = TH1F('h_qOverP_truth_res_gbl'+self.halftag,';q/p - Truth q/p;Entries',50,-0.3,0.3)
        self.h_qOverP_truth_res = TH1F('h_qOverP_truth_res'+self.halftag,';q/p - Truth q/p;Entries',50,-0.3,0.3)
        self.h_qOverP_corr = TH1F('h_qOverP_corr'+self.halftag,';q/p correction;Entries',50,-6.0e-1,6.0e-1)
        self.h_qOverP = TH1F('h_qOverP'+self.halftag,';q/p;Entries',50,-5,5)
        self.h_qOverP_gbl = TH1F('h_qOverP_gbl'+self.halftag,';q/p;Entries',50,-5,5)
        self.h_vtx_xT_corr = TH1F('h_vtx_xT_corr'+self.halftag,';x_{T} vtx corr',50,-5.0e-1,5.0e-1)
        self.h_vtx_yT_corr = TH1F('h_vtx_yT_corr'+self.halftag,';y_{T} vtx corr',50,-5.0e-1,5.0e-1)
        self.h_d0_corr = TH1F('h_d0_corr'+self.halftag,';d_{0} corr [mm]',50,-5.0e-1,5.0e-1)
        self.h_z0_corr = TH1F('h_z0_corr'+self.halftag,';z_{0} corr [mm]',50,-5.0e-1,5.0e-1)
        if self.testRun:
            self.h_d0_initial = TH1F('h_d0_initial'+self.halftag,';d_{0} [mm]',50,-2.0e1,6.0e1)
            self.h_d0_gbl = TH1F('h_d0_gbl'+self.halftag,';d_{0} [mm]',50,-2.0e1,6.0e1)
            self.h_clPar_initial_xT = TH1F('h_clPar_initial_xT'+self.halftag,';x_{T} (mm)',50,-10.,40.)
            if self.isTop: 
                self.h_z0_initial = TH1F('h_z0_initial'+self.halftag,';z_{0} [mm]',50,1.0e1,4.0e1)
                self.h_z0_gbl = TH1F('h_z0_gbl'+self.halftag,';z_{0} [mm]',50,1.0e1,4.0e1)
                self.h_clPar_initial_yT = TH1F('h_clPar_initial_yT'+self.halftag,';y_{T} (mm)',50,10.,40.)
                self.h_clPar_initial_lambda = TH1F('h_clPar_initial_lambda'+self.halftag,';#lambda (rad)',50,0,0.08)
            elif self.isBot: 
                self.h_z0_initial = TH1F('h_z0_initial'+self.halftag,';z_{0} [mm]',50,-4.0e1,-1.0e1)
                self.h_z0_gbl = TH1F('h_z0_gbl'+self.halftag,';z_{0} [mm]',50,-4.0e1,-1.0e1)
                self.h_clPar_initial_yT = TH1F('h_clPar_initial_yT'+self.halftag,';y_{T} (mm)',50,-40.,-10.)
                self.h_clPar_initial_lambda = TH1F('h_clPar_initial_lambda'+self.halftag,';#lambda (rad)',50,-0.08,0.)
            else: 
                self.h_z0_initial = TH1F('h_z0_initial'+self.halftag,';z_{0} [mm]',50,-4.0e1,4.0e1)
                self.h_z0_gbl = TH1F('h_z0_gbl'+self.halftag,';z_{0} [mm]',50,-4.0e1,4.0e1)
                self.h_clPar_initial_yT = TH1F('h_clPar_initial_yT'+self.halftag,';y_{T} (mm)',50,-40.,40.)
                self.h_clPar_initial_lambda = TH1F('h_clPar_initial_lambda'+self.halftag,';#lambda (rad)',50,-0.1,0.1)
        else:
            self.h_d0_initial = TH1F('h_d0_initial'+self.halftag,';d_{0} [mm]',50,-3.0e0,3.0e0)
            self.h_z0_initial = TH1F('h_z0_initial'+self.halftag,';z_{0} [mm]',50,-2.0e0,2.0e0)
            self.h_d0_gbl = TH1F('h_d0_gbl'+self.halftag,';d_{0} [mm]',50,-3.0e0,3.0e0)
            self.h_z0_gbl = TH1F('h_z0_gbl'+self.halftag,';z_{0} [mm]',50,-2.0e0,2.0e0)
            self.h_clPar_initial_xT = TH1F('h_clPar_initial_xT'+self.halftag,';x_{T} (mm)',50,-4.,4.)
            self.h_clPar_initial_yT = TH1F('h_clPar_initial_yT'+self.halftag,';y_{T} (mm)',50,-2.,2.)
            self.h_clPar_initial_lambda = TH1F('h_clPar_initial_lambda'+self.halftag,';#lambda (rad)',50,-0.1,0.1)
        self.h_clPar_initial_qOverP = TH1F('h_clPar_initial_qOverP'+self.halftag,';q/p (1/GeV)',50,-5.,5.)
        self.h_clPar_initial_phi = TH1F('h_clPar_initial_phi'+self.halftag,';#phi (rad)',50,-0.2,0.2)

        self.h_perPar_res_initial_d0 = TH1F('h_perPar_res_initial_d0'+self.halftag,';d0-d0_{truth}',50,-3.0e0,3.0e0)
        self.h_perPar_res_initial_phi0 = TH1F('h_perPar_res_initial_phi0'+self.halftag,';phi0-phi0_{truth}',50,-3.0e-2,3.0e-2)
        self.h_perPar_res_initial_kappa = TH1F('h_perPar_res_initial_kappa'+self.halftag,';kappa-kappa_{truth}',50,-3.0e-5,3.0e-5)
        self.h_perPar_res_initial_z0 = TH1F('h_perPar_res_initial_z0'+self.halftag,';z0-z0_{truth}',50,-3.0e0,3.0e0)
        self.h_perPar_res_initial_slope = TH1F('h_perPar_res_initial_slope'+self.halftag,';slope-slope_{truth}',50,-3.0e-2,3.0e-2)
                
        self.h_measMsCov = TH2F('h_measMsCov'+self.halftag,';Layer;Uncorrelated MS error in meas. dir. (mm)',12,1.,13.,100,0.,7.)
        self.h_xT_corr = TH2F('h_xT_corr'+self.halftag,';Point;x_{T} correction (mm)',24,1.,25.,100,-30.,30.)
        self.h_yT_corr = TH2F('h_yT_corr'+self.halftag,';Point;y_{T} correction (mm)',24,1.,25.,100,-0.3,0.3)
        self.h_clParGBL_res_xT = TH1F('h_clParGBL_res_xT'+self.halftag,';x_{T}-x_{T,truth}',50,-3.0e0,3.0e0)
        self.h_clParGBL_res_yT = TH1F('h_clParGBL_res_yT'+self.halftag,';y_{T}-y_{T,truth}',50,-3.0e0,3.0e0)
        self.h_clParGBL_res_qOverP = TH1F('h_clParGBL_res_qOverP'+self.halftag,';q/p-q/p_{truth}',50,-1.0e0,1.0e0)
        self.h_clParGBL_res_lambda = TH1F('h_clParGBL_res_lambda'+self.halftag,';#lambda-#lambda_{truth}',50,-1.0e-2,1.0e-2)
        self.h_clParGBL_res_phi = TH1F('h_clParGBL_res_phi'+self.halftag,';#phi-#phi_{truth}',50,-1.0e-2,1.0e-2)
        
        self.h_clParGBL_pull_xT = TH1F('h_clParGBL_pull_xT'+self.halftag,';x_{T}-x_{T,truth}',50,-5.0,5.0)
        self.h_clParGBL_pull_yT = TH1F('h_clParGBL_pull_yT'+self.halftag,';y_{T}-y_{T,truth}',50,-5.0,5.0)
        self.h_clParGBL_pull_qOverP = TH1F('h_clParGBL_pull_qOverP'+self.halftag,';q/p-q/p_{truth}',50,-5.0,5.0)
        self.h_clParGBL_pull_lambda = TH1F('h_clParGBL_pull_lambda'+self.halftag,';#lambda-#lambda_{truth}',50,-5.0,5.0)
        self.h_clParGBL_pull_phi = TH1F('h_clParGBL_pull_phi'+self.halftag,';#phi-#phi_{truth}',50,-5.0,5.0)
        self.h_map_res_layer = {}
        self.h_map_res_gbl_layer = {}
        self.h_map_corr_lambda_layer = {}
        self.h_map_corrdiff_lambda_layer = {}
        self.h_map_corr_phi_layer = {}
        self.h_map_corrdiff_phi_layer = {}
        self.h_map_res_truth_layer = {}
        self.h_map_pred_layer = {}
        self.h_map_iso_layer = {}

        self.gr_ures = TGraphErrors()
        self.gr_ures_truth = TGraph()
        self.gr_ures_simhit = TGraph()
        self.gr_corr_ures = TGraph()
        self.gr_ures_corr = TGraph()
    


    

        
    def fillSensorPlots(self,type,deName,val):
        h = None
        if type=="res":
            if deName in self.h_map_res_layer:
                h = self.h_map_res_layer[deName]
            else:
                l = getLayer(deName)
                xmax = 0.01 + (l-1)*0.6
                xmin = -1.*xmax
                h = TH1F('h_res_%s%s'%(deName,self.halftag),'%s;Residual in measurement direction (mm);Entries'%deName,50,xmin,xmax)
                self.h_map_res_layer[deName] = h
            h.Fill(val)
        elif type=="res_truth":
            if deName in self.h_map_res_truth_layer:
                h = self.h_map_res_truth_layer[deName]
            else:
                l = getLayer(deName)
                xmax = 0.01 + (l-1)*0.6
                xmin = -1.*xmax
                h = TH1F('h_res_truth_%s%s'%(deName,self.halftag),'%s;Residual in measurement direction (mm);Entries'%deName,50,xmin,xmax)
                self.h_map_res_truth_layer[deName] = h
            h.Fill(val)
        elif type=="res_gbl":
            if deName in self.h_map_res_gbl_layer:
                h = self.h_map_res_gbl_layer[deName]
            else:
                l = getLayer(deName)
                xmax = 0.01 
                xmin = -1.*xmax
                h = TH1F('h_res_gbl_%s%s'%(deName,self.halftag),'%s;Residual in measurement direction (mm);Entries'%deName,50,xmin,xmax)
                self.h_map_res_gbl_layer[deName] = h
            h.Fill(val)
        elif type=="corr_lambda":
            if deName in self.h_map_corr_lambda_layer:
                h = self.h_map_corr_lambda_layer[deName]
            else:
                l = getLayer(deName)
                xmax = 0.007
                xmin = -1.*xmax
                h = TH1F('h_corr_lambda_%s%s'%(deName,self.halftag),'%s;Lambda correction in CL frame;Entries'%deName,50,xmin,xmax)
                self.h_map_corr_lambda_layer[deName] = h
            h.Fill(val)
        elif type=="corrdiff_lambda":
            if deName in self.h_map_corrdiff_lambda_layer:
                h = self.h_map_corrdiff_lambda_layer[deName]
            else:
                l = getLayer(deName)
                xmax = 0.006
                xmin = -1.*xmax
                h = TH1F('h_corrdiff_lambda_%s%s'%(deName,self.halftag),'%s;#Delta#lambda correction in CL frame;Entries'%deName,50,xmin,xmax)
                self.h_map_corrdiff_lambda_layer[deName] = h
            h.Fill(val)
        elif type=="corr_phi":
            if deName in self.h_map_corr_phi_layer:
                h = self.h_map_corr_phi_layer[deName]
            else:
                l = getLayer(deName)
                xmax = 0.007
                xmin = -1.*xmax
                h = TH1F('h_corr_phi_%s%s'%(deName,self.halftag),'%s;Phi correction in CL frame;Entries'%deName,50,xmin,xmax)
                self.h_map_corr_phi_layer[deName] = h
            h.Fill(val)
        elif type=="corrdiff_phi":
            if deName in self.h_map_corrdiff_phi_layer:
                h = self.h_map_corrdiff_phi_layer[deName]
            else:
                l = getLayer(deName)
                xmax = 0.006
                xmin = -1.*xmax
                h = TH1F('h_corrdiff_phi_%s%s'%(deName,self.halftag),'%s;#Delta#phi correction in CL frame;Entries'%deName,50,xmin,xmax)
                self.h_map_corrdiff_phi_layer[deName] = h
            h.Fill(val)
        elif type=="pred_meas":
            if deName in self.h_map_pred_layer:
                h = self.h_map_pred_layer[deName]
            else:
                l = getLayer(deName)
                h = TH2F('h_pred_%s%s'%(deName,self.halftag),'%s;Hit pred. in v (mm);Hit pred. in u (mm)'%deName,20,-60,60,20,-25,25)
                self.h_map_pred_layer[deName] = h
            h.Fill(val[1],val[0])
        elif type=="iso":
            if deName in self.h_map_iso_layer:
                h = self.h_map_iso_layer[deName]
            else:
                l = getLayer(deName)
                h = TH1F('h_iso_%s%s'%(deName,self.halftag),'%s;Strip hit isolation in u (mm);Strip clusters'%deName,43,-2,41)
                self.h_map_iso_layer[deName] = h
            if val>999.9:
                h.Fill(-2.0)
            else:
                h.Fill(val)
        else:
            print "Thus type ius not defined ", type
    
                
    def show(self,save,nopause):
        #if not save:
        #    return
    

        gStyle.SetOptStat(111111)
        
        c_chi2_gbl = TCanvas('c_chi2_gbl'+self.halftag,'c_chi2_gbl'+self.halftag,10,10,690*2,500)
        c_chi2_gbl.Divide(3,1)
        c_chi2_gbl.cd(1)
        self.h_chi2.Draw()
        c_chi2_gbl.cd(2)
        self.h_chi2ndf.Draw()
        c_chi2_gbl.cd(3)
        self.h_chi2prob.Draw()

        
        c_chi2_initial = TCanvas('c_chi2_initial'+self.halftag,'c_chi2_initial'+self.halftag,10,10,690*2,500)
        c_chi2_initial.Divide(3,1)
        c_chi2_initial.cd(1)
        self.h_chi2_initial.Draw()
        c_chi2_initial.cd(2)
        self.h_chi2ndf_initial.Draw()
        c_chi2_initial.cd(3)
        self.h_chi2prob_initial.Draw()
        
        
        c_chi2_initial_truth = TCanvas('c_chi2_initial_truth'+self.halftag,'c_chi2_initial_truth'+self.halftag,10,10,690*2,500)
        c_chi2_initial_truth.Divide(3,1)
        c_chi2_initial_truth.cd(1)
        self.h_chi2_initial_truth.Draw()
        c_chi2_initial_truth.cd(2)
        self.h_chi2ndf_initial_truth.Draw()
        c_chi2_initial_truth.cd(3)
        self.h_chi2prob_initial_truth.Draw()
        
        
        c_chi2_gbl_truth = TCanvas('c_chi2_gbl_truth'+self.halftag,'c_chi2_gbl_truth'+self.halftag,10,10,690*2,500)
        c_chi2_gbl_truth.Divide(3,1)
        c_chi2_gbl_truth.cd(1)
        self.h_chi2_gbl_truth.Draw()
        c_chi2_gbl_truth.cd(2)
        self.h_chi2ndf_gbl_truth.Draw()
        c_chi2_gbl_truth.cd(3)
        self.h_chi2prob_gbl_truth.Draw()
        
        
        c_track_momentum = TCanvas('c_track_momentum'+self.halftag,'c_track_momentum'+self.halftag,10,10,690*2,500)
        c_track_momentum.Divide(3,1)
        c_track_momentum.cd(1)
        self.h_p.SetFillStyle(1001);
        self.h_p.SetFillColor(4);
        self.h_p.SetLineColor(1);
        self.h_p.Fit("gaus","R","",1.0,1.4)        
        self.h_p_gbl.SetFillStyle(3005);
        self.h_p_gbl.SetFillColor(2);
        self.h_p_gbl.SetLineColor(2);
        self.h_p_gbl.Draw("")
        self.h_p_gbl.Fit("gaus","R","",1.0,1.4)        
        self.h_p.Draw("same")
        self.h_p_gbl.Draw("same")
        self.h_p_truth.SetLineStyle(3)
        self.h_p_truth.SetLineColor(3)
        self.h_p_truth.Draw("same")
        f_gbl = self.h_p_gbl.GetFunction("gaus")
        if f_gbl != None:
            myText(0.7,0.63,"%.2f#sigma%.2f"% (f_gbl.GetParameter(1),f_gbl.GetParameter(2)),0.05,2)
        if self.h_p.GetFunction("gaus") != None:
            myText(0.7,0.55,"%.2f#sigma%.2f"% (self.h_p.GetFunction("gaus").GetParameter(1),self.h_p.GetFunction("gaus").GetParameter(2)),0.05,4)        
        c_track_momentum.cd(2)
        self.h_qOverP.SetFillStyle(1001);
        self.h_qOverP.SetFillColor(4);
        self.h_qOverP.SetLineColor(1);
        self.h_qOverP_gbl.SetFillStyle(3005);
        self.h_qOverP_gbl.SetFillColor(2);
        self.h_qOverP_gbl.SetLineColor(2);
        self.h_qOverP_gbl.Draw("")
        self.h_qOverP.Draw("same")
        self.h_qOverP_gbl.Draw("same")
        c_track_momentum.cd(3)
        self.h_qOverP_corr.Draw()
        
                        
        c_track_momentum_res = TCanvas('c_track_momentum_res'+self.halftag,'c_track_momentum_res'+self.halftag,10,10,690,490)
        c_track_momentum_res.Divide(2,2)
        c_track_momentum_res.cd(1)
        self.h_p_truth_res.Fit('gaus','Q')
        self.h_p_truth_res.Draw()
        if self.h_p_truth_res.GetFunction('gaus')!=None:
            myText(0.6,0.7,'#sigma=%.3fGeV'%self.h_p_truth_res.GetFunction('gaus').GetParameter(2),0.05,1)
        c_track_momentum_res.cd(2)
        self.h_p_truth_res_gbl.Fit('gaus','Q')
        self.h_p_truth_res_gbl.Draw()
        if self.h_p_truth_res_gbl.GetFunction('gaus') != None:
            myText(0.6,0.7,'#sigma=%.3fGeV'%self.h_p_truth_res_gbl.GetFunction('gaus').GetParameter(2),0.05,1)
        c_track_momentum_res.cd(3)
        self.h_qOverP_truth_res.Fit('gaus','Q')
        self.h_qOverP_truth_res.Draw()
        if self.h_qOverP_truth_res.GetFunction('gaus') != None:
            myText(0.6,0.7,'#sigma=%.3f'%self.h_qOverP_truth_res.GetFunction('gaus').GetParameter(2),0.05,1)
        c_track_momentum_res.cd(4)
        self.h_qOverP_truth_res_gbl.Fit('gaus','Q')
        self.h_qOverP_truth_res_gbl.Draw()
        if self.h_qOverP_truth_res_gbl.GetFunction('gaus') != None:
            myText(0.6,0.7,'#sigma=%.3f'%self.h_qOverP_truth_res_gbl.GetFunction('gaus').GetParameter(2),0.05,1)
        
        
        c_track_momentum_res_vs_p = TCanvas('c_track_momentum_res_vs_p'+self.halftag,'c_track_momentum_res_vs_p'+self.halftag,10,10,690,490)
        c_track_momentum_res_vs_p.Divide(2,2)
        c_track_momentum_res_vs_p.cd(1)
        self.h_p_truth_res_vs_p.Draw('colz')
        c_track_momentum_res_vs_p.cd(2)
        self.h_p_truth_res_gbl_vs_p.Draw('colz')    
        c_track_momentum_res_vs_p.cd(3)
        gr_vs_p = TGraphErrors()
        for b in range(1,self.h_p_truth_res_vs_p.GetNbinsX()+1):
            hprj = self.h_p_truth_res_vs_p.ProjectionY('%s_%d'%(self.h_p_truth_res_vs_p.GetName(),b),b,b,"")
            if hprj.Integral() > 20: 
                x = self.h_p_truth_res_vs_p.GetXaxis().GetBinCenter(b)
                hprj.Fit('gaus')
                y = hprj.GetFunction('gaus').GetParameter(2)
                y /= x
                dy = hprj.GetFunction('gaus').GetParError(2)
                dy /= x
                if (dy/y) < 1.0:
                    gr_vs_p.SetPoint(b-1,x,y)
                    gr_vs_p.SetPointError(b-1,0.,dy)
        gr_gbl_vs_p = TGraphErrors()
        for b in range(1,self.h_p_truth_res_gbl_vs_p.GetNbinsX()+1):
            hprj = self.h_p_truth_res_gbl_vs_p.ProjectionY('%s_%d'%(self.h_p_truth_res_gbl_vs_p.GetName(),b),b,b,"")
            if hprj.Integral() > 20: 
                x = self.h_p_truth_res_gbl_vs_p.GetXaxis().GetBinCenter(b)
                hprj.Fit('gaus')
                y = hprj.GetFunction('gaus').GetParameter(2)
                y /= x
                dy = hprj.GetFunction('gaus').GetParError(2)
                dy /= x
                if (dy/y) < 1.0:
                    gr_gbl_vs_p.SetPoint(b-1,x,y)
                    gr_gbl_vs_p.SetPointError(b-1,0.,dy)
        gr_vs_p.Draw('ALP')
        gr_vs_p.SetMarkerStyle(20)
        gr_vs_p.SetTitle(';Truth momentum (GeV);#sigma(p)/p')
        gr_gbl_vs_p.SetTitle(';Truth momentum (GeV);#sigma(p)/p')
        gr_gbl_vs_p.SetMarkerStyle(20)
        gr_gbl_vs_p.SetMarkerColor(2)
        gr_gbl_vs_p.SetLineColor(2)
        gr_gbl_vs_p.Draw('LP,same')
        leg_vs_p = TLegend(0.22,0.4,0.6,0.85)
        leg_vs_p.SetFillColor(0)
        leg_vs_p.SetBorderSize(0)
        leg_vs_p.AddEntry(gr_vs_p,'Initial fit', 'LP')
        leg_vs_p.AddEntry(gr_gbl_vs_p,'GBL refit', 'LP')
        leg_vs_p.Draw()
        
        c_track_momentum_res_vs_p.cd(4)
        gr_vs_p_ratio = TGraphErrors()
        for b in range(gr_vs_p.GetN()):
            x1 = ROOTDouble(0.)
            y1 = ROOTDouble(0.)
            gr_vs_p.GetPoint(b,x1,y1)
            dy1 = gr_vs_p.GetErrorY(b)
            x2 = ROOTDouble(0.)
            y2 = ROOTDouble(0.)
            gr_gbl_vs_p.GetPoint(b,x2,y2)
            dy2 = gr_gbl_vs_p.GetErrorY(b)
            r = y2/y1
            dr = (dy2/y1)*(dy2/y1) + (y2*dy1/(y1*y1))*(y2*dy1/(y1*y1))
            dr = sqrt(dr)
            gr_vs_p_ratio.SetPoint(b,x1,r)
            gr_vs_p_ratio.SetPointError(b,0.,dr)
        gr_vs_p_ratio.SetMarkerStyle(20)
        gr_vs_p_ratio.SetTitle(';Truth momentum (GeV);#sigma_{GBL}/#sigma_{initial}')
        gr_vs_p_ratio.Draw('ALP')
        
        
        
        
        
        c_vtx_corr = TCanvas('c_vtx_corr'+self.halftag,'c_vtx_corr'+self.halftag,10,10,690,490)
        c_vtx_corr.Divide(2,2)
        c_vtx_corr.cd(1)
        self.h_vtx_xT_corr.Draw()
        c_vtx_corr.cd(2)
        self.h_vtx_yT_corr.Draw()
        c_vtx_corr.cd(3)
        self.h_d0_corr.Draw()
        c_vtx_corr.cd(4)
        self.h_z0_corr.Draw()

        
        c_impactParameters_corr = TCanvas('c_impactParameters_corr'+self.halftag,'c_impactParameters_corr'+self.halftag,10,10,690,490)
        c_impactParameters_corr.Divide(2,2)
        c_impactParameters_corr.cd(1)
        self.h_d0_initial.Fit('gaus','Q')
        self.h_d0_initial.Draw()
        if self.h_d0_initial.GetFunction('gaus') != None:
            myText(0.6,0.76,'mean=%.3f'%self.h_d0_initial.GetFunction('gaus').GetParameter(1),0.05,1)
            myText(0.6,0.7,'#sigma=%.3f'%self.h_d0_initial.GetFunction('gaus').GetParameter(2),0.05,1)
        else:
            myText(0.6,0.7,'#sigma=?',0.05,1)            
        c_impactParameters_corr.cd(2)
        self.h_z0_initial.Fit('gaus','Q')
        self.h_z0_initial.Draw()
        if self.h_z0_initial.GetFunction('gaus') != None:
            myText(0.6,0.76,'mean=%.3f'%self.h_z0_initial.GetFunction('gaus').GetParameter(1),0.05,1)
            myText(0.6,0.7,'#sigma=%.3f'%self.h_z0_initial.GetFunction('gaus').GetParameter(2),0.05,1)
        else:
            myText(0.6,0.7,'#sigma=?',0.05,1)                        
        c_impactParameters_corr.cd(3)
        self.h_d0_gbl.Fit('gaus','Q')
        self.h_d0_gbl.Draw()
        if self.h_d0_gbl.GetFunction('gaus') != None:
            myText(0.6,0.76,'mean=%.3f'%self.h_d0_gbl.GetFunction('gaus').GetParameter(1),0.05,1)
            myText(0.6,0.7,'#sigma=%.3f'%self.h_d0_gbl.GetFunction('gaus').GetParameter(2),0.05,1)
        else:
            myText(0.6,0.7,'#sigma=?',0.05,1)                                    
        c_impactParameters_corr.cd(4)
        self.h_z0_gbl.Fit('gaus','Q')
        self.h_z0_gbl.Draw()
        if self.h_z0_gbl.GetFunction('gaus') != None:
            myText(0.6,0.76,'mean=%.3f'%self.h_z0_gbl.GetFunction('gaus').GetParameter(1),0.05,1)
            myText(0.6,0.7,'#sigma=%.3f'%self.h_z0_gbl.GetFunction('gaus').GetParameter(2),0.05,1)
        else:
            myText(0.6,0.7,'#sigma=?',0.05,1)                                                
        
        
        c_perPar_res = TCanvas('c_perPar_res_initial'+self.halftag,'c_perPar_res_initial'+self.halftag,10,10,690*2,390)
        c_perPar_res.Divide(5,1)
        c_perPar_res.cd(1)
        self.h_perPar_res_initial_d0.Draw()
        c_perPar_res.cd(2)
        self.h_perPar_res_initial_phi0.Draw()
        c_perPar_res.cd(3)
        self.h_perPar_res_initial_kappa.Draw()
        c_perPar_res.cd(4)
        self.h_perPar_res_initial_z0.Draw()
        c_perPar_res.cd(5)
        self.h_perPar_res_initial_slope.Draw()
        
        
        if self.doDetailed:
            c_measMsCov = TCanvas('c_measMsCov'+self.halftag,'c_measMsCov'+self.halftag,10,10,690*2,490)
            c_measMsCov.Divide(2,1)
            c_measMsCov.cd(1)
            self.h_measMsCov.SetStats(False)
            self.h_measMsCov.Draw('colz')
            c_measMsCov.cd(2)
            self.h_measMsCov_prf = self.h_measMsCov.ProfileX()
            self.h_measMsCov_prf.SetTitle(';Layer;Mean MS error in meas. dir (mm)')
            self.h_measMsCov_prf.SetStats(False)
            self.h_measMsCov_prf.Draw()
            gPad.SetLogy()
        
            if(save): c_measMsCov.SaveAs('measMsCov%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        
        if self.doDetailed:
            c_xT_corr_label = TCanvas('c_xT_corr_label'+self.halftag,'c_xT_corr_label'+self.halftag,10,10,690*2,390)
            c_xT_corr_label.Divide(2,1)
            c_xT_corr_label.cd(1)
            self.h_xT_corr_prf = self.h_xT_corr.ProfileX()
            for b in range(1,self.h_xT_corr.GetNbinsX()+1):
                label = int(self.h_xT_corr.GetBinLowEdge(b))
                if label % 2 == 0: 
                    lbl = -1*label/2
                else: 
                    lbl = (label+1)/2
                self.h_xT_corr.GetXaxis().SetBinLabel(b,'%d'%lbl)
                self.h_xT_corr_prf.GetXaxis().SetBinLabel(b,'%d'%lbl)
            self.h_xT_corr.SetStats(False)
            self.h_xT_corr.Draw('colz')
            c_xT_corr_label.cd(2)
            self.h_xT_corr_prf.SetStats(False)
            self.h_xT_corr_prf.SetTitle(';Point;Mean x_{T} correction (mm)')
            self.h_xT_corr_prf.Draw()
            
            if(save): c_xT_corr_label.SaveAs('xT_corr_label%s.%s'%(self.getTag(),self.picExt),self.picExt)
            
            c_yT_corr_label = TCanvas('c_yT_corr_label'+self.halftag,'c_yT_corr_label'+self.halftag,10,10,690*2,390)
            c_yT_corr_label.Divide(2,1)
            c_yT_corr_label.cd(1)
            self.h_yT_corr.SetStats(False)
            self.h_yT_corr_prf = self.h_yT_corr.ProfileX()
            for b in range(1,self.h_yT_corr.GetNbinsX()+1):
                label = int(self.h_yT_corr.GetBinLowEdge(b))
                if label % 2 == 0: 
                    lbl = -1*label/2
                else: 
                    lbl = (label+1)/2
                self.h_yT_corr.GetXaxis().SetBinLabel(b,'%d'%lbl)
                self.h_yT_corr_prf.GetXaxis().SetBinLabel(b,'%d'%lbl)
            self.h_yT_corr.Draw('colz')
            c_yT_corr_label.cd(2)
            self.h_yT_corr_prf.SetStats(False)
            self.h_yT_corr_prf.SetTitle(';Point;Mean y_{T} correction (mm)')
            self.h_yT_corr_prf.Draw()
            
            if(save): c_yT_corr_label.SaveAs('yT_corr_label%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        
        c_clPar_initial = TCanvas('c_clPar_initial'+self.halftag,'c_clPar_initial'+self.halftag,10,10,690*2,390)
        c_clPar_initial.Divide(5,1)
        c_clPar_initial.cd(1)
        self.h_clPar_initial_xT.Draw()
        c_clPar_initial.cd(2)
        self.h_clPar_initial_yT.Draw()
        c_clPar_initial.cd(3)
        self.h_clPar_initial_qOverP.Draw()
        c_clPar_initial.cd(4)
        self.h_clPar_initial_lambda.Draw()
        c_clPar_initial.cd(5)
        self.h_clPar_initial_phi.Draw()
        
        
        c_clParGBL_res = TCanvas('c_clParGBL_res'+self.halftag,'c_clParGBL_res'+self.halftag,10,10,690*2,390)
        c_clParGBL_res.Divide(5,1)
        c_clParGBL_res.cd(1)
        self.h_clParGBL_res_xT.Draw()
        c_clParGBL_res.cd(2)
        self.h_clParGBL_res_yT.Draw()
        c_clParGBL_res.cd(3)
        self.h_clParGBL_res_qOverP.Draw()
        c_clParGBL_res.cd(4)
        self.h_clParGBL_res_lambda.Draw()
        c_clParGBL_res.cd(5)
        self.h_clParGBL_res_phi.Draw()
        
        
        
        c_clParGBL_pull = TCanvas('c_clParGBL_pull'+self.halftag,'c_clParGBL_pull'+self.halftag,10,10,690*2,390)
        c_clParGBL_pull.Divide(5,1)
        c_clParGBL_pull.cd(1)
        self.h_clParGBL_pull_xT.Draw()
        c_clParGBL_pull.cd(2)
        self.h_clParGBL_pull_yT.Draw()
        c_clParGBL_pull.cd(3)
        self.h_clParGBL_pull_qOverP.Draw()
        c_clParGBL_pull.cd(4)
        self.h_clParGBL_pull_lambda.Draw()
        c_clParGBL_pull.cd(5)
        self.h_clParGBL_pull_phi.Draw()
        
        
        
        if self.doDetailed:
            c_res_example = TCanvas('c_res_example'+self.halftag,'c_res_example'+self.halftag,10,10,690*2,490)
            leg_res_example = TLegend(0.12,0.6,0.27,0.9)
            c_res_example.cd(1)
            self.gr_ures_truth.SetTitle('Example track fit;Path length (mm);Residual in measured direction (mm)')
            self.gr_ures_truth.SetLineColor(2)
            self.gr_ures_truth.SetLineStyle(2)
            self.gr_ures_truth.Draw('ALP')
            self.gr_ures.SetMarkerSize(0.8)
            self.gr_ures.SetMarkerStyle(20)
            self.gr_ures.Draw('LP,same')
            self.gr_ures_simhit.SetLineColor(2)
            self.gr_ures_simhit.SetMarkerStyle(5)
            self.gr_ures_simhit.Draw('P,same')
            leg_res_example.AddEntry(self.gr_ures_truth,'Truth (no MS)','L')
            leg_res_example.AddEntry(self.gr_ures_simhit,'Truth simhit','L')
            leg_res_example.AddEntry(self.gr_ures,'Initial','P')
            self.gr_corr_ures.SetTitle('Example track fit;Path length (mm);Correction in measured direction (mm)')
            self.gr_corr_ures.SetLineColor(3)
            self.gr_corr_ures.SetMarkerStyle(21)
            self.gr_corr_ures.SetMarkerSize(0.8)
            self.gr_corr_ures.SetMarkerColor(3)
            self.gr_ures_corr.SetLineColor(3)
            self.gr_ures_corr.SetMarkerColor(3)
            self.gr_ures_corr.SetMarkerStyle(20)
            self.gr_ures_corr.Draw('LP,same')
            leg_res_example.AddEntry(self.gr_ures_corr,'Corrected','LP')
            leg_res_example.SetFillStyle(0)
            leg_res_example.SetFillColor(0)
            leg_res_example.SetBorderSize(0)
            leg_res_example.Draw()
            
            if(save): c_res_example.SaveAs('res_example%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        
        c_res_initial_sensor = TCanvas('c_res_initial_sensor'+self.halftag,'c_res_initial_sensor'+self.halftag,10,10,690*2,390*2)
        c_res_initial_sensor.Divide(6,3)
        c_res_initial_sensor_mean = TCanvas('c_res_initial_sensor_mean'+self.halftag,'c_res_initial_sensor_mean'+self.halftag,10,10,690*2,390*2)
        c_res_initial_sensor_mean.SetBottomMargin(0.45)
        gr_res_initial_sensor_mean = TGraphErrors()
        gr_res_initial_sensor_mean.SetTitle(';Sensor;u residual (mm)')        
        c_res_initial_sensor_rms = TCanvas('c_res_initial_sensor_rms'+self.halftag,'c_res_initial_sensor_rms'+self.halftag,10,10,690*2,390*2)
        c_res_initial_sensor_rms.SetBottomMargin(0.45)
        gr_res_initial_sensor_rms = TGraphErrors()
        gr_res_initial_sensor_rms.SetTitle(';Sensor;u residual (mm)')        
        i = 1
        ms = sorted(self.h_map_res_layer,key=getLayer)
        idToSensor = {}
        for sensor in ms:
            c_res_initial_sensor.cd(i)
            h = self.h_map_res_layer[sensor] 
            if h.GetEntries() > 10.:
                h.Fit('gaus','Q')
            h.Draw()
            if h.GetEntries() > 0.:
                ii = gr_res_initial_sensor_mean.GetN()
                gr_res_initial_sensor_mean.SetPoint(ii,ii,h.GetMean())
                gr_res_initial_sensor_mean.SetPointError(ii,0.,h.GetMeanError())
                ii = gr_res_initial_sensor_rms.GetN()
                gr_res_initial_sensor_rms.SetPoint(ii,ii,h.GetRMS())
                gr_res_initial_sensor_rms.SetPointError(ii,0.,h.GetRMSError())
                idToSensor[ii] = sensor
            i=i+1
        
        setGraphXLabels(gr_res_initial_sensor_mean,idToSensor)
        setGraphXLabels(gr_res_initial_sensor_rms,idToSensor)


        c_res_initial_sensor2 = TCanvas('c_res_initial_sensor2'+self.halftag,'c_res_initial_sensor2'+self.halftag,10,10,690*2,390*2)
        c_res_initial_sensor2.Divide(6,3)
        i = 1
        ms = sorted(self.h_map_res_layer,key=getLayer)
        idToSensor = {}
        for sensor in ms:
            i = getCanvasIdxTwoCols(sensor)
            c_res_initial_sensor.cd(i)
            h = self.h_map_res_layer[sensor] 
            h.Draw()

        c_res_gbl_sensor = TCanvas('c_res_gbl_sensor'+self.halftag,'c_res_gbl_sensor'+self.halftag,10,10,690*2,390*2)
        c_res_gbl_sensor.Divide(6,3)
        c_res_gbl_sensor_mean = TCanvas('c_res_gbl_sensor_mean'+self.halftag,'c_res_gbl_sensor_mean'+self.halftag,10,10,690*2,390*2)
        c_res_gbl_sensor_mean.SetBottomMargin(0.45)
        gr_res_gbl_sensor_mean = TGraphErrors()
        gr_res_gbl_sensor_mean.SetTitle(';Sensor;u residual (mm)')        
        c_res_gbl_sensor_rms = TCanvas('c_res_gbl_sensor_rms'+self.halftag,'c_res_gbl_sensor_rms'+self.halftag,10,10,690*2,390*2)
        c_res_gbl_sensor_rms.SetBottomMargin(0.45)
        gr_res_gbl_sensor_rms = TGraphErrors()
        gr_res_gbl_sensor_rms.SetTitle(';Sensor;u residual (mm)')        
        i = 1
        ms = sorted(self.h_map_res_gbl_layer,key=getLayer)
        idToSensor = {}
        for sensor in ms:
            c_res_gbl_sensor.cd(i)
            h = self.h_map_res_gbl_layer[sensor] 
            if h.GetEntries() > 10.:
                h.Fit('gaus','Q')
            h.Draw()
            if h.GetEntries() > 0.:
                ii = gr_res_gbl_sensor_mean.GetN()
                gr_res_gbl_sensor_mean.SetPoint(ii,ii,h.GetMean())
                gr_res_gbl_sensor_mean.SetPointError(ii,0.,h.GetMeanError())
                ii = gr_res_gbl_sensor_rms.GetN()
                gr_res_gbl_sensor_rms.SetPoint(ii,ii,h.GetRMS())
                gr_res_gbl_sensor_rms.SetPointError(ii,0.,h.GetRMSError())
                idToSensor[ii] = sensor

            i=i+1

        setGraphXLabels(gr_res_gbl_sensor_mean,idToSensor)
        setGraphXLabels(gr_res_gbl_sensor_rms,idToSensor)



        c_res_gbl_sensor2 = TCanvas('c_res_gbl_sensor2'+self.halftag,'c_res_gbl_sensor2'+self.halftag,10,10,690,390*2)
        c_res_gbl_sensor2.Divide(2,12)
        i = 1
        ms = sorted(self.h_map_res_gbl_layer,key=getLayer)
        for sensor in ms:
            i = getCanvasIdxTwoCols(sensor)
            c_res_gbl_sensor2.cd(i)
            h = self.h_map_res_gbl_layer[sensor] 
            h.Draw()
        


        c_corr_lambda_sensor = TCanvas('c_corr_lambda_sensor'+self.halftag,'c_corr_lambda_sensor'+self.halftag,10,10,690,390*2)
        c_corr_lambda_sensor.Divide(2,12)
        c_corr_lambda_sensor_mean = TCanvas('c_corr_lambda_sensor_mean'+self.halftag,'c_corr_lambda_sensor_mean'+self.halftag,10,10,690*2,390*2)
        c_corr_lambda_sensor_mean.SetBottomMargin(0.45)
        gr_corr_lambda_sensor_mean = TGraphErrors()
        gr_corr_lambda_sensor_mean.SetTitle(';Sensor;#lambda correction mean')        
        c_corr_lambda_sensor_rms = TCanvas('c_corr_lambda_sensor_rms'+self.halftag,'c_corr_lambda_sensor_rms'+self.halftag,10,10,690*2,390*2)
        c_corr_lambda_sensor_rms.SetBottomMargin(0.45)
        gr_corr_lambda_sensor_rms = TGraphErrors()
        gr_corr_lambda_sensor_rms.SetTitle(';Sensor;#lambda correction width')        
        i = 1
        ms = sorted(self.h_map_corr_lambda_layer,key=getLayer)
        idToSensor = {}
        for sensor in ms:
            i = getCanvasIdxTwoCols(sensor)
            c_corr_lambda_sensor.cd(i)
            h = self.h_map_corr_lambda_layer[sensor]
            if h.GetEntries() > 10.0:
                h.Fit("gaus","Q") 
            h.Draw()
            if h.GetEntries() > 0.:
                f = h.GetFunction("gaus")
                if f != None:
                    ii = gr_corr_lambda_sensor_mean.GetN()
                    dm = f.GetParError(1)
                    if dm > 0.5:
                        print 'WARNING: fit for ', h.GetName(), ' has a large error ', dm, '. Set to zero.'
                        dm = 0.0
                    ds = f.GetParError(2)
                    if ds > 0.5:
                        print 'WARNING: fit for ', h.GetName(), ' has a large error ', ds, '. Set to zero.'
                        ds = 0.0
                    gr_corr_lambda_sensor_mean.SetPoint(ii,ii,f.GetParameter(1))
                    gr_corr_lambda_sensor_mean.SetPointError(ii,0.,dm)
                    gr_corr_lambda_sensor_rms.SetPoint(ii,ii,f.GetParameter(2))
                    gr_corr_lambda_sensor_rms.SetPointError(ii,0.,ds)
                    idToSensor[ii] = sensor

            i=i+1

        setGraphXLabels(gr_corr_lambda_sensor_mean,idToSensor)
        setGraphXLabels(gr_corr_lambda_sensor_rms,idToSensor)


        c_corrdiff_lambda_sensor = TCanvas('c_corrdiff_lambda_sensor'+self.halftag,'c_corrdiff_lambda_sensor'+self.halftag,10,10,690,390*2)
        c_corrdiff_lambda_sensor.Divide(2,12)
        c_corrdiff_lambda_sensor_mean = TCanvas('c_corrdiff_lambda_sensor_mean'+self.halftag,'c_corrdiff_lambda_sensor_mean'+self.halftag,10,10,690*2,390*2)
        c_corrdiff_lambda_sensor_mean.SetBottomMargin(0.45)
        gr_corrdiff_lambda_sensor_mean = TGraphErrors()
        gr_corrdiff_lambda_sensor_mean.SetTitle(';Sensor;#lambda correction mean')        
        c_corrdiff_lambda_sensor_rms = TCanvas('c_corrdiff_lambda_sensor_rms'+self.halftag,'c_corrdiff_lambda_sensor_rms'+self.halftag,10,10,690*2,390*2)
        c_corrdiff_lambda_sensor_rms.SetBottomMargin(0.45)
        gr_corrdiff_lambda_sensor_rms = TGraphErrors()
        gr_corrdiff_lambda_sensor_rms.SetTitle(';Sensor;#lambda correction width')        
        i = 1
        ms = sorted(self.h_map_corrdiff_lambda_layer,key=getLayer)
        idToSensor = {}
        for sensor in ms:
            i = getCanvasIdxTwoCols(sensor)
            c_corrdiff_lambda_sensor.cd(i)
            h = self.h_map_corrdiff_lambda_layer[sensor]
            if h.GetEntries() > 10.0:
                h.Fit("gaus","Q") 
            h.Draw()
            if h.GetEntries() > 0.:
                f = h.GetFunction("gaus")
                if f != None:                    
                    ii = gr_corrdiff_lambda_sensor_mean.GetN()
                    dm = f.GetParError(1)
                    if dm > 0.5:
                        print 'WARNING: fit for ', h.GetName(), ' has a large error ', dm, '. Set to zero.'
                        dm = 0.0
                    ds = f.GetParError(2)
                    if ds > 0.5:
                        print 'WARNING: fit for ', h.GetName(), ' has a large error ', ds, '. Set to zero.'
                        ds = 0.0

                    gr_corrdiff_lambda_sensor_mean.SetPoint(ii,ii,f.GetParameter(1))
                    gr_corrdiff_lambda_sensor_mean.SetPointError(ii,0.,dm)
                    gr_corrdiff_lambda_sensor_rms.SetPoint(ii,ii,f.GetParameter(2))
                    gr_corrdiff_lambda_sensor_rms.SetPointError(ii,0.,ds)
                    idToSensor[ii] = sensor
            i=i+1
        
        setGraphXLabels(gr_corrdiff_lambda_sensor_mean,idToSensor)
        setGraphXLabels(gr_corrdiff_lambda_sensor_rms,idToSensor)



        c_corr_phi_sensor = TCanvas('c_corr_phi_sensor'+self.halftag,'c_corr_phi_sensor'+self.halftag,10,10,690,390*2)
        c_corr_phi_sensor.Divide(2,12)
        c_corr_phi_sensor_mean = TCanvas('c_corr_phi_sensor_mean'+self.halftag,'c_corr_phi_sensor_mean'+self.halftag,10,10,690*2,390*2)
        c_corr_phi_sensor_mean.SetBottomMargin(0.45)
        gr_corr_phi_sensor_mean = TGraphErrors()
        gr_corr_phi_sensor_mean.SetTitle(';Sensor;#phi correction mean')        
        c_corr_phi_sensor_rms = TCanvas('c_corr_phi_sensor_rms'+self.halftag,'c_corr_phi_sensor_rms'+self.halftag,10,10,690*2,390*2)
        c_corr_phi_sensor_rms.SetBottomMargin(0.45)
        gr_corr_phi_sensor_rms = TGraphErrors()
        gr_corr_phi_sensor_rms.SetTitle(';Sensor;#phi correction width')        
        i = 1
        ms = sorted(self.h_map_corr_phi_layer,key=getLayer)
        idToSensor = {}
        for sensor in ms:
            i = getCanvasIdxTwoCols(sensor)
            c_corr_phi_sensor.cd(i)
            h = self.h_map_corr_phi_layer[sensor]
            if h.GetEntries() > 10.0:
                h.Fit("gaus","Q") 
            h.Draw()
            if h.GetEntries() > 0.:
                f = h.GetFunction("gaus")
                if f != None:
                    ii = gr_corr_phi_sensor_mean.GetN()
                    dm = f.GetParError(1)
                    if dm > 0.5:
                        print 'WARNING: fit for ', h.GetName(), ' has a large error ', dm, '. Set to zero.'
                        dm = 0.0
                    ds = f.GetParError(2)
                    if ds > 0.5:
                        print 'WARNING: fit for ', h.GetName(), ' has a large error ', ds, '. Set to zero.'
                        ds = 0.0
                    gr_corr_phi_sensor_mean.SetPoint(ii,ii,f.GetParameter(1))
                    gr_corr_phi_sensor_mean.SetPointError(ii,0.,dm)
                    gr_corr_phi_sensor_rms.SetPoint(ii,ii,f.GetParameter(2))
                    gr_corr_phi_sensor_rms.SetPointError(ii,0.,ds)
                    idToSensor[ii] = sensor

            i=i+1

        setGraphXLabels(gr_corr_phi_sensor_mean,idToSensor)
        setGraphXLabels(gr_corr_phi_sensor_rms,idToSensor)


        c_corrdiff_phi_sensor = TCanvas('c_corrdiff_phi_sensor'+self.halftag,'c_corrdiff_phi_sensor'+self.halftag,10,10,690,390*2)
        c_corrdiff_phi_sensor.Divide(2,12)
        c_corrdiff_phi_sensor_mean = TCanvas('c_corrdiff_phi_sensor_mean'+self.halftag,'c_corrdiff_phi_sensor_mean'+self.halftag,10,10,690*2,390*2)
        c_corrdiff_phi_sensor_mean.SetBottomMargin(0.45)
        gr_corrdiff_phi_sensor_mean = TGraphErrors()
        gr_corrdiff_phi_sensor_mean.SetTitle(';Sensor;#phi correction mean')        
        c_corrdiff_phi_sensor_rms = TCanvas('c_corrdiff_phi_sensor_rms'+self.halftag,'c_corrdiff_phi_sensor_rms'+self.halftag,10,10,690*2,390*2)
        c_corrdiff_phi_sensor_rms.SetBottomMargin(0.45)
        gr_corrdiff_phi_sensor_rms = TGraphErrors()
        gr_corrdiff_phi_sensor_rms.SetTitle(';Sensor;#phi correction width')        
        i = 1
        ms = sorted(self.h_map_corrdiff_phi_layer,key=getLayer)
        idToSensor = {}
        for sensor in ms:
            i = getCanvasIdxTwoCols(sensor)
            c_corrdiff_phi_sensor.cd(i)
            h = self.h_map_corrdiff_phi_layer[sensor]
            if h.GetEntries() > 10.0:
                h.Fit("gaus","Q") 
            h.Draw()
            if h.GetEntries() > 0.:
                f = h.GetFunction("gaus")
                if f != None:                    
                    ii = gr_corrdiff_phi_sensor_mean.GetN()
                    dm = f.GetParError(1)
                    if dm > 0.5:
                        print 'WARNING: fit for ', h.GetName(), ' has a large error ', dm, '. Set to zero.'
                        dm = 0.0
                    ds = f.GetParError(2)
                    if ds > 0.5:
                        print 'WARNING: fit for ', h.GetName(), ' has a large error ', ds, '. Set to zero.'
                        ds = 0.0
                    gr_corrdiff_phi_sensor_mean.SetPoint(ii,ii,f.GetParameter(1))
                    gr_corrdiff_phi_sensor_mean.SetPointError(ii,0.,dm)
                    gr_corrdiff_phi_sensor_rms.SetPoint(ii,ii,f.GetParameter(2))
                    gr_corrdiff_phi_sensor_rms.SetPointError(ii,0.,ds)
                    idToSensor[ii] = sensor
            i=i+1
        
        setGraphXLabels(gr_corrdiff_phi_sensor_mean,idToSensor)
        setGraphXLabels(gr_corrdiff_phi_sensor_rms,idToSensor)



        c_res_truth_sensor = TCanvas('c_res_truth_sensor'+self.halftag,'c_res_truth_sensor'+self.halftag,10,10,690*2,390*2)
        c_res_truth_sensor.Divide(6,3)
        c_res_truth_sensor_mean = TCanvas('c_res_truth_sensor_mean'+self.halftag,'c_res_truth_sensor_mean'+self.halftag,10,10,690*2,390*2)
        gr_res_truth_sensor_mean = TGraphErrors()
        gr_res_truth_sensor_mean.SetTitle(';Sensor;u residual (mm)')        
        i = 1
        ms = sorted(self.h_map_res_truth_layer,key=getLayer)
        for sensor in ms:
            c_res_truth_sensor.cd(i)
            h = self.h_map_res_truth_layer[sensor] 
            h.Draw()
            if h.GetEntries() > 0.:
                ii = gr_res_truth_sensor_mean.GetN()
                gr_res_truth_sensor_mean.SetPoint(ii,ii,h.GetMean())
                gr_res_truth_sensor_mean.SetPointError(ii,0.,h.GetMeanError())
            i=i+1



            
        c_res_initial_sensor_mean.cd()
        gr_res_initial_sensor_mean.SetMarkerStyle(20)
        gr_res_initial_sensor_mean.Draw('ALP')
        c_res_gbl_sensor_mean.cd()
        gr_res_gbl_sensor_mean.SetMarkerStyle(20)
        gr_res_gbl_sensor_mean.Draw('ALP')
        c_corr_lambda_sensor_mean.cd()
        gr_corr_lambda_sensor_mean.SetMarkerStyle(20)
        gr_corr_lambda_sensor_mean.Draw('ALP')
        c_corrdiff_lambda_sensor_mean.cd()
        gr_corrdiff_lambda_sensor_mean.SetMarkerStyle(20)
        gr_corrdiff_lambda_sensor_mean.Draw('ALP')
        c_corr_phi_sensor_mean.cd()
        gr_corr_phi_sensor_mean.SetMarkerStyle(20)
        gr_corr_phi_sensor_mean.Draw('ALP')
        c_corrdiff_phi_sensor_mean.cd()
        gr_corrdiff_phi_sensor_mean.SetMarkerStyle(20)
        gr_corrdiff_phi_sensor_mean.Draw('ALP')
        c_res_truth_sensor_mean.cd()
        gr_res_truth_sensor_mean.SetMarkerStyle(20)
        gr_res_truth_sensor_mean.Draw('ALP')

        c_res_initial_sensor_rms.cd()
        gr_res_initial_sensor_rms.SetMarkerStyle(20)
        gr_res_initial_sensor_rms.Draw('ALP')
        c_res_gbl_sensor_rms.cd()
        gr_res_gbl_sensor_rms.SetMarkerStyle(20)
        gr_res_gbl_sensor_rms.Draw('ALP')
        c_corr_lambda_sensor_rms.cd()
        gr_corr_lambda_sensor_rms.SetMarkerStyle(20)
        gr_corr_lambda_sensor_rms.Draw('ALP')
        c_corrdiff_lambda_sensor_rms.cd()
        gr_corrdiff_lambda_sensor_rms.SetMarkerStyle(20)
        gr_corrdiff_lambda_sensor_rms.Draw('ALP')
        c_corr_phi_sensor_rms.cd()
        gr_corr_phi_sensor_rms.SetMarkerStyle(20)
        gr_corr_phi_sensor_rms.Draw('ALP')
        c_corrdiff_phi_sensor_rms.cd()
        gr_corrdiff_phi_sensor_rms.SetMarkerStyle(20)
        gr_corrdiff_phi_sensor_rms.Draw('ALP')




        c_pred_sensor = TCanvas('c_pred_sensor'+self.halftag,'c_pred_sensor'+self.halftag,10,10,690*2,390*2)
        c_pred_sensor.Divide(6,3)
        i = 1
        ms = sorted(self.h_map_pred_layer,key=getLayer)
        for sensor in ms:
            c_pred_sensor.cd(i)
            h = self.h_map_pred_layer[sensor] 
            h.Draw("colz")
            i=i+1

        c_iso_sensor = TCanvas('c_iso_sensor'+self.halftag,'c_iso_sensor'+self.halftag,10,10,690*2,390*2)
        c_iso_sensor.Divide(6,3)
        i = 1
        ms = sorted(self.h_map_iso_layer,key=getLayer)
        for sensor in ms:
            c_iso_sensor.cd(i)
            h = self.h_map_iso_layer[sensor] 
            h.Draw("")
            i=i+1



        if(save):
            c_all = TCanvas('c_all'+self.halftag,'c_all'+self.halftag,10,10,690*2,390*2)
            c_all.Print('gbltst-hps-plots%s.ps['%self.getTag())
            c_res_initial_sensor2.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_res_initial_sensor.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_res_initial_sensor_mean.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_res_initial_sensor.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_res_initial_sensor_mean.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_res_initial_sensor_rms.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_res_truth_sensor.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_res_gbl_sensor2.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_res_gbl_sensor.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_res_gbl_sensor_mean.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_res_gbl_sensor_rms.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_corr_lambda_sensor.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_corr_lambda_sensor_mean.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_corr_lambda_sensor_rms.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_corrdiff_lambda_sensor.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_corrdiff_lambda_sensor_mean.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_corrdiff_lambda_sensor_rms.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_corr_phi_sensor.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_corr_phi_sensor_mean.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_corr_phi_sensor_rms.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_corrdiff_phi_sensor.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_corrdiff_phi_sensor_mean.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_corrdiff_phi_sensor_rms.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_res_truth_sensor.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_res_truth_sensor_mean.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_pred_sensor.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_iso_sensor.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_chi2_gbl.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_chi2_initial.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_chi2_initial_truth.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_chi2_gbl_truth.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_track_momentum.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_track_momentum_res.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_track_momentum_res_vs_p.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_vtx_corr.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_impactParameters_corr.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_perPar_res.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_clPar_initial.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_clParGBL_res.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_clParGBL_pull.Print('gbltst-hps-plots%s.ps'%self.getTag())
            c_all.Print('gbltst-hps-plots%s.ps]'%self.getTag())
            subprocess.call('ps2pdf gbltst-hps-plots%s.ps'%self.getTag(),shell=True)

            #c_res_initial_sensor.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_res_initial_sensor.GetName()))
            #c_res_initial_sensor2.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_res_initial_sensor2.GetName()))
            #c_res_initial_sensor_mean.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_res_initial_sensor_mean.GetName()))
            #c_res_initial_sensor_rms.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_res_initial_sensor_rms.GetName()))
            #c_res_gbl_sensor.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_res_gbl_sensor.GetName()))
            #c_res_gbl_sensor2.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_res_gbl_sensor2.GetName()))
            #c_res_gbl_sensor_mean.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_res_gbl_sensor_mean.GetName()))
            #c_res_gbl_sensor_rms.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_res_gbl_sensor_rms.GetName()))
            #c_corr_lambda_sensor.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_corr_lambda_sensor.GetName()))
            #c_corr_lambda_sensor_mean.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_corr_lambda_sensor_mean.GetName()))
            #c_corr_lambda_sensor_rms.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_corr_lambda_sensor_rms.GetName()))
            #c_corrdiff_lambda_sensor.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_corrdiff_lambda_sensor.GetName()))
            #c_corrdiff_lambda_sensor_mean.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_corrdiff_lambda_sensor_mean.GetName()))
            #c_corrdiff_lambda_sensor_rms.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_corrdiff_lambda_sensor_rms.GetName()))
            #c_corr_phi_sensor.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_corr_phi_sensor.GetName()))
            #c_corr_phi_sensor_mean.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_corr_phi_sensor_mean.GetName()))
            #c_corr_phi_sensor_rms.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_corr_phi_sensor_rms.GetName()))
            #c_corrdiff_phi_sensor.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_corrdiff_phi_sensor.GetName()))
            #c_corrdiff_phi_sensor_mean.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_corrdiff_phi_sensor_mean.GetName()))
            #c_corrdiff_phi_sensor_rms.SaveAs('gbltst-hps-plots%s-%s.png'%(self.getTag(),c_corrdiff_phi_sensor_rms.GetName()))


            saveHistosToFile(gDirectory,'gbltst-hps-plots%s.root'%self.getTag())
        
        if not nopause:
            ans = raw_input('press any key to continue')




def saveHistosToFile(direc,fileName):
    iter = TIter(direc.GetList())
    d = gDirectory
    outfile = TFile(fileName,'RECREATE')
    d.cd()
    print '%d items in collection' % iter.GetCollection().GetSize()
    n = 0
    while True:
        obj = iter.Next()
        #print obj
        if not obj:
            break
        if obj.InheritsFrom('TH1'):
            #print 'Writing %s' % obj.GetName()`
            outfile.cd()
            obj.Write()
            n=n+1
    print 'Saved %d histogram to %s' % (n,outfile.GetName())




