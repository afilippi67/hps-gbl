from ROOT import TH1F, TH2F, TGraph, TGraphErrors, TCanvas, TLegend, TLatex, gStyle, gDirectory, TIter, TFile, gPad
from ROOT import Double as ROOTDouble
from math import sqrt


class plotter:
    def __init__(self,tag,picExt,testRunFlag,isTop,isBot):
        self.picExt = picExt
        self.testRun = testRunFlag
        self.isTop = isTop
        self.isBot = isBot
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
        self.h_chi2_gbl_truth = TH1F('h_chi2_gbl_truth',';GBL truth #chi^{2};Entries',50,0.,100.)
        self.h_chi2ndf_gbl_truth = TH1F('h_chi2ndf_gbl_truth',';GBL truth #chi^{2}/ndf;Entries',50,0.,20.)
        self.h_chi2prob_gbl_truth = TH1F('h_chi2prob_gbl_truth',';GBL truth #chi^{2} prob.;Entries',50,0.,1.)
        self.h_chi2 = TH1F('h_chi2',';GBL #chi^{2};Entries',50,0.,20.)
        self.h_chi2ndf = TH1F('h_chi2ndf',';GBL #chi^{2}/ndf;Entries',50,0.,8.)
        self.h_chi2prob = TH1F('h_chi2prob',';GBL #chi^{2} prob;Entries',50,0.,1.)
        self.h_chi2_initial = TH1F('h_chi2_initial',';Track #chi^{2};Entries',50,0.,20.)
        self.h_chi2ndf_initial = TH1F('h_chi2ndf_initial',';Track #chi^{2}/ndf;Entries',50,0.,8.)
        self.h_chi2prob_initial = TH1F('h_chi2prob_initial',';Track #chi^{2} prob;Entries',50,0.,1.)
        self.h_chi2_initial_truth = TH1F('h_chi2_initial_truth',';Track #chi^{2};Entries',50,0.,20.)
        self.h_chi2ndf_initial_truth = TH1F('h_chi2ndf_initial_truth',';Track #chi^{2}/ndf;Entries',50,0.,8.)
        self.h_chi2prob_initial_truth = TH1F('h_chi2prob_initial_truth',';Track #chi^{2} prob.;Entries',50,0.,1.)
        self.h_p = TH1F('h_p',';Track momentum;Entries',50,0.,3.)
        self.h_p_gbl = TH1F('h_p_gbl',';Track momentum;Entries',50,0.,3.)
        self.h_p_truth = TH1F('h_p_truth',';Track momentum;Entries',50,0.,3.)
        self.h_p_truth_res = TH1F('h_p_truth_res',';Track - Truth momentum;Entries',50,-0.3,0.3)
        self.h_p_truth_res_vs_p = TH2F('h_p_truth_res_vs_p',';Truth momentum;Track - Truth momentum',6,0.4,1.6,50,-0.3,0.3)
        self.h_p_truth_res_gbl = TH1F('h_p_truth_res_gbl',';Track - Truth momentum;Entries',50,-0.3,0.3)
        self.h_p_truth_res_gbl_vs_p = TH2F('h_p_truth_res_gbl_vs_p',';Truth momentum;Track - Truth momentum',6,0.4,1.6,50,-0.3,0.3)
        self.h_qOverP_truth_res_gbl = TH1F('h_qOverP_truth_res_gbl',';q/p - Truth q/p;Entries',50,-0.3,0.3)
        self.h_qOverP_truth_res = TH1F('h_qOverP_truth_res',';q/p - Truth q/p;Entries',50,-0.3,0.3)
        self.h_qOverP_corr = TH1F('h_qOverP_corr',';q/p correction;Entries',50,-6.0e-1,6.0e-1)
        self.h_qOverP = TH1F('h_qOverP',';q/p;Entries',50,1.0e-3,1.0e-3)
        self.h_qOverP_gbl = TH1F('h_qOverP_gbl',';q/p;Entries',50,1.0e-3,1.0e-3)
        self.h_vtx_xT_corr = TH1F('h_vtx_xT_corr',';x_{T} vtx corr',50,-5.0e-1,5.0e-1)
        self.h_vtx_yT_corr = TH1F('h_vtx_yT_corr',';y_{T} vtx corr',50,-5.0e-1,5.0e-1)
        self.h_d0_corr = TH1F('h_d0_corr',';d_{0} corr [mm]',50,-5.0e-1,5.0e-1)
        self.h_z0_corr = TH1F('h_z0_corr',';z_{0} corr [mm]',50,-5.0e-1,5.0e-1)
        if self.testRun:
            self.h_d0_initial = TH1F('h_d0_initial',';d_{0} [mm]',50,-2.0e1,6.0e1)
            self.h_d0_gbl = TH1F('h_d0_gbl',';d_{0} [mm]',50,-2.0e1,6.0e1)
            self.h_clPar_initial_xT = TH1F('h_clPar_initial_xT',';x_{T} (mm)',50,-10.,40.)
            if self.isTop: 
                self.h_z0_initial = TH1F('h_z0_initial',';z_{0} [mm]',50,1.0e1,4.0e1)
                self.h_z0_gbl = TH1F('h_z0_gbl',';z_{0} [mm]',50,1.0e1,4.0e1)
                self.h_clPar_initial_yT = TH1F('h_clPar_initial_yT',';y_{T} (mm)',50,10.,40.)
                self.h_clPar_initial_lambda = TH1F('h_clPar_initial_lambda',';#lambda (rad)',50,0,0.08)
            elif self.isBot: 
                self.h_z0_initial = TH1F('h_z0_initial',';z_{0} [mm]',50,-4.0e1,-1.0e1)
                self.h_z0_gbl = TH1F('h_z0_gbl',';z_{0} [mm]',50,-4.0e1,-1.0e1)
                self.h_clPar_initial_yT = TH1F('h_clPar_initial_yT',';y_{T} (mm)',50,-40.,-10.)
                self.h_clPar_initial_lambda = TH1F('h_clPar_initial_lambda',';#lambda (rad)',50,-0.08,0.)
            else: 
                self.h_z0_initial = TH1F('h_z0_initial',';z_{0} [mm]',50,-4.0e1,4.0e1)
                self.h_z0_gbl = TH1F('h_z0_gbl',';z_{0} [mm]',50,-4.0e1,4.0e1)
                self.h_clPar_initial_yT = TH1F('h_clPar_initial_yT',';y_{T} (mm)',50,-40.,40.)
                self.h_clPar_initial_lambda = TH1F('h_clPar_initial_lambda',';#lambda (rad)',50,-0.1,0.1)
        else:
            self.h_d0_initial = TH1F('h_d0_initial',';d_{0} [mm]',50,-2.0e0,2.0e0)
            self.h_z0_initial = TH1F('h_z0_initial',';z_{0} [mm]',50,-2.0e0,2.0e0)
            self.h_d0_gbl = TH1F('h_d0_gbl',';d_{0} [mm]',50,-2.0e0,2.0e0)
            self.h_z0_gbl = TH1F('h_z0_gbl',';z_{0} [mm]',50,-2.0e0,2.0e0)
            self.h_clPar_initial_xT = TH1F('h_clPar_initial_xT',';x_{T} (mm)',50,-2.,2.)
            self.h_clPar_initial_yT = TH1F('h_clPar_initial_yT',';y_{T} (mm)',50,-2.,2.)
            self.h_clPar_initial_lambda = TH1F('h_clPar_initial_lambda',';#lambda (rad)',50,-0.1,0.1)
        self.h_clPar_initial_qOverP = TH1F('h_clPar_initial_qOverP',';q/p (1/GeV)',50,-3.,3.)
        self.h_clPar_initial_phi = TH1F('h_clPar_initial_phi',';#phi (rad)',50,-0.2,0.2)

        self.h_perPar_res_initial_d0 = TH1F('h_perPar_res_initial_d0',';d0-d0_{truth}',50,-3.0e0,3.0e0)
        self.h_perPar_res_initial_phi0 = TH1F('h_perPar_res_initial_phi0',';phi0-phi0_{truth}',50,-3.0e-2,3.0e-2)
        self.h_perPar_res_initial_kappa = TH1F('h_perPar_res_initial_kappa',';kappa-kappa_{truth}',50,-3.0e-5,3.0e-5)
        self.h_perPar_res_initial_z0 = TH1F('h_perPar_res_initial_z0',';z0-z0_{truth}',50,-3.0e0,3.0e0)
        self.h_perPar_res_initial_slope = TH1F('h_perPar_res_initial_slope',';slope-slope_{truth}',50,-3.0e-2,3.0e-2)
                
        self.h_measMsCov = TH2F('h_measMsCov',';Layer;Uncorrelated MS error in meas. dir. (mm)',12,1.,13.,100,0.,7.)
        self.h_xT_corr = TH2F('h_xT_corr',';Point;x_{T} correction (mm)',24,1.,25.,100,-30.,30.)
        self.h_yT_corr = TH2F('h_yT_corr',';Point;y_{T} correction (mm)',24,1.,25.,100,-0.3,0.3)
        self.h_clParGBL_res_xT = TH1F('h_clParGBL_res_xT',';x_{T}-x_{T,truth}',50,-3.0e0,3.0e0)
        self.h_clParGBL_res_yT = TH1F('h_clParGBL_res_yT',';y_{T}-y_{T,truth}',50,-3.0e0,3.0e0)
        self.h_clParGBL_res_qOverP = TH1F('h_clParGBL_res_qOverP',';q/p-q/p_{truth}',50,-1.0e0,1.0e0)
        self.h_clParGBL_res_lambda = TH1F('h_clParGBL_res_lambda',';#lambda-#lambda_{truth}',50,-1.0e-2,1.0e-2)
        self.h_clParGBL_res_phi = TH1F('h_clParGBL_res_phi',';#phi-#phi_{truth}',50,-1.0e-2,1.0e-2)
        
        self.h_clParGBL_pull_xT = TH1F('h_clParGBL_pull_xT',';x_{T}-x_{T,truth}',50,-5.0,5.0)
        self.h_clParGBL_pull_yT = TH1F('h_clParGBL_pull_yT',';y_{T}-y_{T,truth}',50,-5.0,5.0)
        self.h_clParGBL_pull_qOverP = TH1F('h_clParGBL_pull_qOverP',';q/p-q/p_{truth}',50,-5.0,5.0)
        self.h_clParGBL_pull_lambda = TH1F('h_clParGBL_pull_lambda',';#lambda-#lambda_{truth}',50,-5.0,5.0)
        self.h_clParGBL_pull_phi = TH1F('h_clParGBL_pull_phi',';#phi-#phi_{truth}',50,-5.0,5.0)
        self.h_res_layer = []
        self.h_res_truth_layer = []
        self.h_res_gbl_layer = []
        for l in range(1,13):
            l_tmp = l-1
            if l_tmp % 2 == 1:
                l_tmp = l_tmp - 1
            xmax = 0.02 + (l_tmp)*0.2
            xmin = -1.*xmax
            h = TH1F('h_res_layer%d'%l,'Initial trajectory;Residual in measurement direction layer %d (mm)'%l,50,xmin,xmax)
            self.h_res_layer.append(h)
            h = TH1F('h_res_gbl_layer%d'%l,'GBL trajectory;Residual in measurement direction layer %d (mm)'%l,50,-0.02,0.02)
            self.h_res_gbl_layer.append(h)
            h = TH1F('h_res_truth_layer%d'%l,'Truth trajectory;Residual in measurement direction layer %d (mm)'%l,50,xmin,xmax)
            self.h_res_truth_layer.append(h)
        
        self.gr_ures = TGraphErrors()
        self.gr_ures_truth = TGraph()
        self.gr_ures_simhit = TGraph()
        self.gr_corr_ures = TGraph()
        self.gr_ures_corr = TGraph()
    


    def show(self,save):
        #if not save:
        #    return
    

        gStyle.SetOptStat(111111)
        
        c_chi2_gbl = TCanvas('c_chi2_gbl','c_chi2_gbl',10,10,690*2,500)
        c_chi2_gbl.Divide(3,1)
        c_chi2_gbl.cd(1)
        self.h_chi2.Draw()
        c_chi2_gbl.cd(2)
        self.h_chi2ndf.Draw()
        c_chi2_gbl.cd(3)
        self.h_chi2prob.Draw()

        if(save): c_chi2_gbl.SaveAs('chi2_gbl%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        c_chi2_initial = TCanvas('c_chi2_initial','c_chi2_initial',10,10,690*2,500)
        c_chi2_initial.Divide(3,1)
        c_chi2_initial.cd(1)
        self.h_chi2_initial.Draw()
        c_chi2_initial.cd(2)
        self.h_chi2ndf_initial.Draw()
        c_chi2_initial.cd(3)
        self.h_chi2prob_initial.Draw()
        
        if(save): c_chi2_initial.SaveAs('chi2_initial%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        c_chi2_initial_truth = TCanvas('c_chi2_initial_truth','c_chi2_initial_truth',10,10,690*2,500)
        c_chi2_initial_truth.Divide(3,1)
        c_chi2_initial_truth.cd(1)
        self.h_chi2_initial_truth.Draw()
        c_chi2_initial_truth.cd(2)
        self.h_chi2ndf_initial_truth.Draw()
        c_chi2_initial_truth.cd(3)
        self.h_chi2prob_initial_truth.Draw()
        
        if(save): c_chi2_initial_truth.SaveAs('chi2_initial_truth%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        c_chi2_gbl_truth = TCanvas('c_chi2_gbl_truth','c_chi2_gbl_truth',10,10,690*2,500)
        c_chi2_gbl_truth.Divide(3,1)
        c_chi2_gbl_truth.cd(1)
        self.h_chi2_gbl_truth.Draw()
        c_chi2_gbl_truth.cd(2)
        self.h_chi2ndf_gbl_truth.Draw()
        c_chi2_gbl_truth.cd(3)
        self.h_chi2prob_gbl_truth.Draw()
        
        if(save): c_chi2_gbl_truth.SaveAs('chi2_gbl_truth%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        c_track_momentum = TCanvas('c_track_momentum','c_track_momentum',10,10,690*2,500)
        c_track_momentum.Divide(3,1)
        c_track_momentum.cd(1)
        self.h_p.Draw()
        self.h_p_gbl.SetLineStyle(2)
        self.h_p_gbl.Draw("same")
        self.h_p_truth.SetLineStyle(3)
        self.h_p_truth.SetLineColor(2)
        self.h_p_truth.Draw("same")
        c_track_momentum.cd(2)
        self.h_qOverP.Draw()
        self.h_qOverP_gbl.SetLineStyle(2)
        self.h_qOverP_gbl.Draw("same")
        c_track_momentum.cd(3)
        self.h_qOverP_corr.Draw()
        
        if(save): c_track_momentum.SaveAs('track_momentum%s.%s'%(self.getTag(),self.picExt),self.picExt)
                        
        c_track_momentum_res = TCanvas('c_track_momentum_res','c_track_momentum_res',10,10,690,490)
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
        
        if(save): c_track_momentum_res.SaveAs('track_momentum_resolution%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        c_track_momentum_res_vs_p = TCanvas('c_track_momentum_res_vs_p','c_track_momentum_res_vs_p',10,10,690,490)
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
        
        
        if(save): c_track_momentum_res_vs_p.SaveAs('track_momentum_resolution_vs_p%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        
        
        c_vtx_corr = TCanvas('c_vtx_corr','c_vtx_corr',10,10,690,490)
        c_vtx_corr.Divide(2,2)
        c_vtx_corr.cd(1)
        self.h_vtx_xT_corr.Draw()
        c_vtx_corr.cd(2)
        self.h_vtx_yT_corr.Draw()
        c_vtx_corr.cd(3)
        self.h_d0_corr.Draw()
        c_vtx_corr.cd(4)
        self.h_z0_corr.Draw()

        if(save): c_vtx_corr.SaveAs('vtx_correction%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        c_impactParameters_corr = TCanvas('c_impactParameters_corr','c_impactParameters_corr',10,10,690,490)
        c_impactParameters_corr.Divide(2,2)
        c_impactParameters_corr.cd(1)
        self.h_d0_initial.Fit('gaus','Q')
        self.h_d0_initial.Draw()
        if self.h_d0_initial.GetFunction('gaus') != None:
            myText(0.6,0.7,'#sigma=%.3f'%self.h_d0_initial.GetFunction('gaus').GetParameter(2),0.05,1)
        else:
            myText(0.6,0.7,'#sigma=?',0.05,1)            
        c_impactParameters_corr.cd(2)
        self.h_z0_initial.Fit('gaus','Q')
        self.h_z0_initial.Draw()
        if self.h_z0_initial.GetFunction('gaus') != None:
            myText(0.6,0.7,'#sigma=%.3f'%self.h_z0_initial.GetFunction('gaus').GetParameter(2),0.05,1)
        else:
            myText(0.6,0.7,'#sigma=?',0.05,1)                        
        c_impactParameters_corr.cd(3)
        self.h_d0_gbl.Fit('gaus','Q')
        self.h_d0_gbl.Draw()
        if self.h_d0_gbl.GetFunction('gaus') != None:
            myText(0.6,0.7,'#sigma=%.3f'%self.h_d0_gbl.GetFunction('gaus').GetParameter(2),0.05,1)
        else:
            myText(0.6,0.7,'#sigma=?',0.05,1)                                    
        c_impactParameters_corr.cd(4)
        self.h_z0_gbl.Fit('gaus','Q')
        self.h_z0_gbl.Draw()
        if self.h_z0_gbl.GetFunction('gaus') != None:
            myText(0.6,0.7,'#sigma=%.3f'%self.h_z0_gbl.GetFunction('gaus').GetParameter(2),0.05,1)
        else:
            myText(0.6,0.7,'#sigma=?',0.05,1)                                                
        
        if(save): c_impactParameters_corr.SaveAs('impactParameters%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        c_perPar_res = TCanvas('c_perPar_res_initial','c_perPar_res_initial',10,10,690*2,390)
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
        
        if(save): c_perPar_res.SaveAs('perPar_res%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        if self.doDetailed:
            c_measMsCov = TCanvas('c_measMsCov','c_measMsCov',10,10,690*2,490)
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
            c_xT_corr_label = TCanvas('c_xT_corr_label','c_xT_corr_label',10,10,690*2,390)
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
            
            c_yT_corr_label = TCanvas('c_yT_corr_label','c_yT_corr_label',10,10,690*2,390)
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
        
        
        c_clPar_initial = TCanvas('c_clPar_initial','c_clPar_initial',10,10,690*2,390)
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
        
        if(save): c_clPar_initial.SaveAs('clPar_initial%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        
        c_clParGBL_res = TCanvas('c_clParGBL_res','c_clParGBL_res',10,10,690*2,390)
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
        
        if(save): c_clParGBL_res.SaveAs('clParGBL_res%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        
        c_clParGBL_pull = TCanvas('c_clParGBL_pull','c_clParGBL_pull',10,10,690*2,390)
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
        
        if(save): c_clParGBL_pull.SaveAs('clParGBL_pull%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        
        if self.doDetailed:
            c_res_example = TCanvas('c_res_example','c_res_example',10,10,690*2,490)
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
        
        
        
        c_res_initial_layer = TCanvas('c_res_initial_layer','c_res_initial_layer',10,10,690*2,390*2)
        c_res_initial_layer.Divide(4,3)
        c_res_initial_layer_mean = TCanvas('c_res_initial_layer_mean','c_res_initial_layer_mean',10,10,690*2,390*2)
        gr_res_initial_layer_mean = TGraphErrors()
        gr_res_initial_layer_mean.SetTitle(';Layer;u residual (mm)')
        c_res_truth_layer = TCanvas('c_res_truth_layer','c_res_truth_layer',10,10,690*2,390*2)
        c_res_truth_layer.Divide(4,3)
        c_res_gbl_layer = TCanvas('c_res_gbl_layer','c_res_gbl_layer',10,10,690*2,390*2)
        c_res_gbl_layer.Divide(4,3)
        c_res_gbl_layer_mean = TCanvas('c_res_gbl_layer_mean','c_res_gbl_layer_mean',10,10,690*2,390*2)
        gr_res_gbl_layer_mean = TGraphErrors()
        gr_res_gbl_layer_mean.SetTitle(';Layer;u residual (mm)')
        for i in range(1,13):
            c_res_truth_layer.cd(i)
            self.h_res_truth_layer[i-1].Draw()
            c_res_initial_layer.cd(i)
            self.h_res_layer[i-1].Draw()
            if self.h_res_layer[i-1].GetEntries() > 0.:
                gr_res_initial_layer_mean.SetPoint(i-1,i,self.h_res_layer[i-1].GetMean())
                gr_res_initial_layer_mean.SetPointError(i-1,0.,self.h_res_layer[i-1].GetMeanError())
            c_res_gbl_layer.cd(i)
            self.h_res_gbl_layer[i-1].Draw()
            if self.h_res_gbl_layer[i-1].GetEntries() > 0.:
                gr_res_gbl_layer_mean.SetPoint(i-1,i,self.h_res_gbl_layer[i-1].GetMean())
                gr_res_gbl_layer_mean.SetPointError(i-1,0.,self.h_res_gbl_layer[i-1].GetMeanError())

            
        c_res_initial_layer_mean.cd()
        gr_res_initial_layer_mean.SetMarkerStyle(20)
        gr_res_initial_layer_mean.Draw('ALP')
        c_res_gbl_layer_mean.cd()
        gr_res_gbl_layer_mean.SetMarkerStyle(20)
        gr_res_gbl_layer_mean.Draw('ALP')
            
        if(save): c_res_initial_layer.SaveAs('res_initial_individual_layer%s.%s'%(self.getTag(),self.picExt),self.picExt)
        if(save): c_res_initial_layer_mean.SaveAs('res_initial_individual_layer_mean%s.%s'%(self.getTag(),self.picExt),self.picExt)
        if(save): c_res_truth_layer.SaveAs('res_truth_individual_layer%s.%s'%(self.getTag(),self.picExt),self.picExt)
        if(save): c_res_gbl_layer.SaveAs('res_gbl_individual_layer%s.%s'%(self.getTag(),self.picExt),self.picExt)
        if(save): c_res_gbl_layer_mean.SaveAs('res_gbl_individual_layer_mean%s.%s'%(self.getTag(),self.picExt),self.picExt)
        
        
        if(save): saveHistosToFile(gDirectory,'gbltst-hps-plots%s.root'%self.getTag())
        
        ans = raw_input('press any key to continue')



def myText(x,y,text, tsize,color):
    l = TLatex()
    l.SetTextSize(tsize); 
    l.SetNDC();
    l.SetTextColor(color);
    l.DrawLatex(x,y,text);

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
            #print 'Writing %s' % obj.GetName()
            outfile.cd()
            obj.Write()
            n=n+1
    print 'Saved %d histogram to %s' % (n,outfile.GetName())
