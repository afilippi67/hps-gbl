'''
Simple Test Program for General Broken Lines for HPS.

Created on Jul 27, 2011
Edited on Jun 20, 2013

@author: kleinwrt, phansson
'''
import numpy as np
import math
import time
import sys
import utils
#from gblpy.gblfit import GblPoint, GblTrajectory
sys.path.append('../GeneralBrokenLines/python')
from gblfit import GblPoint, GblTrajectory
#from toolspy.simpleHelix import SimpleHelix
from simpleHelix import SimpleHelix
#import hps_plots as plots
from ROOT import TH1F,TCanvas,TMath
#

debug = False
#debug = True
useUncorrMS = False # inflate MS errors instead of using scatterers
nEventsMax = 1000
inputfile = 'gblinput-proposal.txt'

 
def exampleHpsTest(inputfile):
  '''
  Read initial fit and points from  test file
  Create trajectory from points,
  fit and write trajectory to MP-II binary file,
  get track parameter corrections and covariance matrix at points.
  
  Detector arrangement from text file
  '''  

  Chi2Sum = 0.
  NdfSum = 0
  LostSum = 0.

  Bz = -0.5 # full detector  -0.491 for test run detector
  bfac = 0.0002998 * Bz # for Bz in Tesla, momentum in GeV and Radius in mm

  np.random.seed(47117)

  binaryFile = open("milleBinaryISN.dat", "wb")

  inputFile = open(inputfile, 'r')
  events = utils.readHPSEvents(inputFile, nEventsMax)
  
  print 'Read %d events from file' % len(events)

  
  #print " GblHpsTest $Rev: 234 $ ", nTry, nLayer
  nTry = 0
  start = time.clock()

  h_chi2prob_gbl_truth = TH1F('h_chi2prob_gbl_truth','h_chi2prob_gbl_truth',50,0,1)
  h_chi2prob_initial_truth = TH1F('h_chi2prob_initial_truth','h_chi2prob_initial_truth',50,0,1)


  for event in events:
    
    if debug: print '\nEvent %d has %d tracks ' % (event.id, len(event.tracks))
    

    for track in event.tracks:

      if debug: print '\nProcessing track %d \n p = %.3f p(truth)=%.3f' % (track.id, track.p(bfac), track.p_truth(bfac))

      # if there's no truth info -> skip it
      if track.p_truth(bfac) == 0.:
        print 'No truth info, skip track %d in event %d ' % (track.id, event.id)
        continue
      
      traj = GblTrajectory(True)

      #print " perPar ", track.perParTruth
#      hlxPar = [ cmp(bfac, 0.) * track.curvature(), track.phi0(), track.d0(), track.slope(), track.z0()]
      hlxPar = [ cmp(bfac, 0.) * track.curvature(), track.phi0(), -track.d0(), track.slope(), track.z0()]
      #hlxPar = [ cmp(bfac, 0.) * track.curvature_truth(), track.phi0_truth(), track.d0_truth(), track.slope_truth(), track.z0_truth()]
      #print " hlxPar ", hlxPar
      hlx = SimpleHelix(hlxPar)
      cosLambda = hlx.getZSDirection()[0]
      # stupid error for now
      #clCov = np.eye(5)
      #for i in range(5):
      #  clCov[i, i] = clErr[i] ** 2
      #clCov = track.clCov

      stripLabelMap = {}
      
      if debug: print 'Track has %d strip clusters' % len(track.strips)
      # arc length
      s = 0.
      # point-to-point jacobian (from previous point)    
      jacPointToPoint = np.eye(5)
      #start trajectory at reference point (defining s=0)
      point = GblPoint(jacPointToPoint)
      refLabel = traj.addPoint(point)

      # multiple scattering covariance matrix (for curvilinear track parameters)        
      msCov = np.zeros((5, 5))

      # store projection for later use
      proM2l_list = {} 
      proL2m_list = {} 

      for strip in track.strips:
        
        if debug: print '\nProcessing strip %d at layer %d ' % (strip.id, strip.layer)
 
        # direction in detector plane in XY: (xDir, yDir, 0.) 
        xDir = -strip.w[1]
        yDir = strip.w[0]
        # direction in detector plane in Z: (0., 0., 1.) 
        # position of detector (center)
        xDet = strip.origin[0]
        yDet = strip.origin[1]
        zDet = strip.origin[2]
        # get prediction along the 2 directions (in the detector plane)
        pred = hlx.getExpectedPlanePos(xDet, yDet, xDir, yDir, zDet)
        if pred is None:
          continue # no intersection of track and detector
        # prediction position (in global system)
        xPred = xDet + pred[0] * xDir
        yPred = yDet + pred[0] * yDir
        zPred = zDet + pred[1]

        # project onto u-direction
        uPred = pred[0] * (xDir * strip.u[0] + yDir * strip.u[1]) + pred[1] * strip.u[2]

        # stupid iterative intercept
        predIter = utils.getXPlanePositionIterative(track.perPar,strip.origin,strip.w,1.0e-8)
        diffIter = predIter - strip.origin
        uPredIter = np.dot(strip.u , diffIter.T )
        
        #xPred = predIter[0]
        #yPred = predIter[1]
        #zPred = predIter[2]
        #pred0 = (predIter[0] - xDet) / xDir
        #pred1 = (predIter[1] - zDet)
        #uPredIter = pred0 * (xDir * strip.u[0] + yDir * strip.u[1]) + pred1 * strip.u[2]
        

        # u residuum
        uRes = strip.meas - uPred
        uResIter = strip.meas - uPredIter
        #uRes = uResIter
        # (3D) arc-length
        sArc = pred[3] / cosLambda
        phi = pred[4]
        #print " pred ", sArc, xPred, yPred, zPred
        if nTry == 0:
          print " uRes ",strip.id, ' uRes ',  uRes, ' pred ', xPred, yPred, zPred, ' s(3D) ', sArc
          #print " uRes ", strip.id, 'uRes(Cl) ', uRes, ' uRes ', strip.ures, ' uResIter ', uResIter, ' pred ', xPred, yPred, zPred, ' predIter ', predIter, ' s(3D) ', sArc
         #print -1*track.d0(),track.z0(),track.phi0(),track.slope(),track.curvature()
          #print predIter, diffIter, strip.u, uPredIter, strip.meas
          #print " uRes ", strip.id, uRes, uResIter, strip.ures, strip.ures_err, ' pred ', xPred, yPred, zPred, ' s ', pred[3], ' s3D ', sArc, ' on plane ', np.dot(np.array([[xPred, yPred, zPred]]) - strip.origin , np.array([strip.w]).T )
          #print ' predIter ', predIter
          #print ' predIter diff ', (np.array([xPred,yPred,zPred]) - predIter) 

          #print " java ", strip.id , strip.ures, ' pred ', strip.tPos, ' on plane ', np.dot(np.array([strip.tPos]) - strip.origin , np.array([strip.w]).T ) 
        
        step = sArc - s

        if debug: print 'Step %f (s %f pathLen %f)' % (step, s, sArc)

        # measurement direction(s): (m[0]=u, m[1]=v) 
        if debug: print 'Strip udir', strip.u
        if debug: print 'Strip vdir', strip.v
        
        mDir = np.array([strip.u, strip.v])
        
        if debug: 
          print 'mDir:\n', mDir
        
        # track direction: in x directon
        sinLambda = strip.sinLambda
        cosLambda = math.sqrt(1.0 - sinLambda ** 2)
        sinPhi = strip.sinPhi
        cosPhi = math.sqrt(1.0 - sinPhi ** 2)
        
        if debug: print 'Track direction sinLambda=%f sinPhi=%f' % (sinLambda, sinPhi)
        
        #  tDir = np.array([cosLambda * cosPhi, cosLambda * sinPhi, sinLambda])
        # U = Z x T / |Z x T|, V = T x U
        uvDir = np.array([[-sinPhi, cosPhi, 0.], \
                            [-sinLambda * cosPhi, -sinLambda * sinPhi, cosLambda]])
        
        
        # projection measurement to local (curvilinear uv) directions (duv/dm)
        proM2l = np.dot(uvDir, mDir.T)

        proM2l_list[strip.id] = proM2l

        if debug: print 'proM2l:\n', proM2l

        # projection local (uv) to measurement directions (dm/duv)
        proL2m = np.linalg.inv(proM2l)

        proL2m_list[strip.id] = proL2m

        if debug: print 'proL2m:\n', proL2m

        # measurement/residual in the measurement system
        #meas = np.array([strip.ures, 0.])
        meas = np.array([uRes, 0.])
        #meas[0] += deltaU[iLayer] # misalignment
        measErr = np.array([strip.ures_err, strip.ures_err])
        measPrec = 1.0 / measErr ** 2
        measPrec[1] = 0. # 1D measurement perpendicular to strip direction
        
        if debug: print 'meas ', meas, ' measErr ', measErr, ' measPrec ', measPrec

        #propagate to this strip
        #jacPointToPoint = utils.gblSimpleJacobianLambdaPhi(step, cosLambda, bfac)
        jacPointToPoint = hlx.getPropagatorSimple(step, abs(bfac))
        #print jacPointToPoint
        
        point = GblPoint(jacPointToPoint)
        
        if debug: 
          print 'jacPointToPoint to extrapolate to this point:'
          print point.getP2pJacobian()
        
        #propagate MS covariance matrix        
        msCov = np.dot(jacPointToPoint, np.dot(msCov, jacPointToPoint.T))
        # MS covariance for measurements
        measMsCov = np.dot(proL2m, np.dot(msCov[3:, 3:], proL2m.T))
 
        if debug:
          print " uPred ", strip.id, pred[3], uPred, strip.meas, strip.ures, strip.ures_err, measMsCov[0, 0]
       
        #plots.h_measMsCov.Fill(float(strip.layer),measMsCov[0,0])
        
        if debug:
          print 'msCov propagated to this point:'
          print msCov
          print 'measMsCov at this point to be used in measPrec:'
          print measMsCov
        
        if useUncorrMS:
          # blow up measurement errors according to multiple scattering
          measPrec[0] = 1.0 / (measErr[0] ** 2 + measMsCov[0, 0])
        
        point.addMeasurement([proL2m, meas, measPrec])

        if debug:
          print 'measMsCov ', measMsCov[0, 0]
        
        scat = np.array([0., 0.])
        scatErr = np.array([ strip.scatAngle, strip.scatAngle / cosLambda]) 
        scatPrec = 1.0 / scatErr ** 2
        
        if not useUncorrMS:
          point.addScatterer([scat, scatPrec])
        
        #update MS covariance matrix
        msCov[1, 1] += scatErr[0] ** 2; msCov[2, 2] += scatErr[1] ** 2

        if debug:
          print 'adding scatError to the msCov from this point:'
          print scatErr
                
        addDer = np.array([[1.0], [0.0]])
        #top or bottom half        
        if math.copysign(1, sinLambda) > 0:
          offset = 11101
        else:
          offset = 21101
        labGlobal = np.array([[offset + strip.layer], [0]])
        point.addGlobals(labGlobal, addDer)
        
        # add point to trajectory
        iLabel = traj.addPoint(point)
        s += step
        stripLabelMap[strip] = iLabel
      
      
      if debug: print 'Do the fit'
      Chi2, Ndf, Lost = traj.fit()

      # write to millepede
      traj.milleOut(binaryFile)

      # sum up    
      Chi2Sum += Chi2
      NdfSum += Ndf
      LostSum += Lost
      # get corrections and covariance matrix at points 
      result = utils.GBLResults(track)
      #traj.dump()
        
      if nTry == 0:
        print 'fit result: Chi2=%f Ndf=%d Lost=%d' % (Chi2, Ndf, Lost)
        print 'get corrections and covariance matrix for %d points:' % 1 #traj.getNumPoints()
      
      for i in range(1, traj.getNumPoints() + 1):      
        # label start at 1
        locPar, locCov = traj.getResults(-i)
        if nTry < 0:
          print " >Point ", i
          print " locPar ", locPar
          #print " locCov ", locCov      
        result.addPoint(-i, locPar, locCov)
        locPar, locCov = traj.getResults(i)
        if nTry < 0:
          print " Point> ", i
          print " locPar ", locPar
          #print " locCov ", locCov
        result.addPoint(i, locPar, locCov)
      


      # calculate the truth chi2 from initial fit
      # get the truth and fitted params with indexes same as cov matrix of initial fit (dca,phi0,curv,z0,slope)
      perParVec = np.array([track.d0(), track.phi0(), track.curvature(), track.z0(), track.slope()])
      perParVecTruth = np.array([track.d0_truth(), track.phi0_truth(), track.curvature_truth(), track.z0_truth(), track.slope_truth()])
      perParVecRes = perParVec - perParVecTruth
      chi2_initial_truth = np.dot(perParVecRes, np.dot(np.linalg.inv(track.perCov) , perParVecRes))

      # calculate the truth chi2 from gbl fit at vertex
      clParVtx = np.array(track.clPar) + np.array(result.locPar[1])
      clParTruth = np.array(track.clParTruth)
      clParRes = clParVtx - clParTruth
      chi2_gbl_truth = np.dot(clParRes, np.dot(np.linalg.inv(result.locCov[1]), clParRes))

      #print " truth ", track.clParTruth
      #print " res ", refLabel, result.locPar[refLabel], result.locCov[refLabel]
      # calculate chi2 for seeding by truth 
      label = 1#refLabel
      chi2_res = np.dot(result.locPar[label], np.dot(np.linalg.inv(result.locCov[label]), result.locPar[label]))
      #chi2_res4 = np.dot(result.locPar[label][:4], np.dot(np.linalg.inv(result.locCov[label][:4, :4]), result.locPar[label][:4]))
      print " Chi2: ",
      #for i in range(5):
      #  print track.clParTruth[i], result.locPar[label][i] / math.sqrt(result.locCov[label][i][i]),
      print event.id, chi2_res, chi2_gbl_truth, chi2_initial_truth
      #print clParRes
      #print track.clPar
      #print result.locPar[1]
      #print clParVtx
      #print clParTruth
      #print result.locCov[label]
      h_chi2prob_gbl_truth.Fill(TMath.Prob(chi2_gbl_truth,5))
      h_chi2prob_initial_truth.Fill(TMath.Prob(chi2_initial_truth,5))

      '''
      print " clPar ", track.clPar
      print " clParTruth ", track.clParTruth
      print " clParVtx ", clParVtx
      print " clParRes ", clParRes
      print " res[1] ", np.array(result.locPar[1])
      print " cov[1] ", result.locCov[1]
      '''
      '''
      # plots

      plots.h_clPar_xT.Fill(track.clPar[3])
      plots.h_clPar_yT.Fill(track.clPar[4])
      plots.h_clPar_qOverP.Fill(track.clPar[0])
      plots.h_clPar_lambda.Fill(track.clPar[1])
      # transform phi to plot nicer
      if track.clPar[2]<math.pi:
        plots.h_clPar_phi.Fill(track.clPar[2])
      else:
        plots.h_clPar_phi.Fill(track.clPar[2]-math.pi*2)
      plots.h_clPar_res_qOverP.Fill(clParRes[0,0])
      plots.h_clPar_res_lambda.Fill(clParRes[0,1])
      plots.h_clPar_res_phi.Fill(clParRes[0,2])
      plots.h_clPar_res_xT.Fill(clParRes[0,3])
      plots.h_clPar_res_yT.Fill(clParRes[0,4])

      plots.h_clPar_pull_qOverP.Fill(clParRes[0,0]/math.sqrt(result.locCov[1][0,0]))
      plots.h_clPar_pull_lambda.Fill(clParRes[0,1]/math.sqrt(result.locCov[1][1,1]))
      plots.h_clPar_pull_phi.Fill(clParRes[0,2]/math.sqrt(result.locCov[1][2,2]))
      plots.h_clPar_pull_xT.Fill(clParRes[0,3]/math.sqrt(result.locCov[1][3,3]))
      plots.h_clPar_pull_yT.Fill(clParRes[0,4]/math.sqrt(result.locCov[1][4,4]))

      plots.h_perPar_res_d0.Fill(perParVecRes[0,0])
      plots.h_perPar_res_phi0.Fill(perParVecRes[0,1])
      plots.h_perPar_res_kappa.Fill(perParVecRes[0,2])
      plots.h_perPar_res_z0.Fill(perParVecRes[0,3])
      plots.h_perPar_res_slope.Fill(perParVecRes[0,4])
      plots.h_chi2_initial.Fill(track.chi2Initial)
      plots.h_chi2ndf_initial.Fill(track.chi2Initial/track.ndfInitial)
      plots.h_chi2_initial_truth.Fill(chi2_initial_truth)
      plots.h_chi2ndf_initial_truth.Fill(chi2_initial_truth/5.0)
      plots.h_chi2prob_initial_truth.Fill(utils.chi2Prob(chi2_initial_truth,5))
      plots.h_chi2_gbl_truth.Fill(chi2_gbl_truth)
      plots.h_chi2ndf_gbl_truth.Fill(chi2_gbl_truth/5.0)
      plots.h_chi2prob_gbl_truth.Fill(utils.chi2Prob(chi2_gbl_truth,5))
      plots.h_chi2.Fill(Chi2)
      plots.h_chi2ndf.Fill(Chi2/Ndf)
      plots.h_p.Fill(track.p(bfac))
      plots.h_qOverP.Fill(track.qOverP(bfac))
      plots.h_qOverP_truth_res.Fill(track.qOverP(bfac) - track.q()/track.p_truth(bfac))
      plots.h_p_truth.Fill(track.p_truth(bfac))
      plots.h_p_truth_res.Fill(track.p(bfac)-track.p_truth(bfac))
      plots.h_qOverP_corr.Fill(result.qOverPCorr())
      plots.h_qOverP_gbl.Fill(result.qOverP_gbl(bfac))
      plots.h_qOverP_truth_res_gbl.Fill(result.qOverP_gbl(bfac) - result.track.q()/result.track.p_truth(bfac))
      plots.h_p_corr.Fill(result.pCorr(bfac))
      plots.h_p_gbl.Fill(result.p_gbl(bfac))
      plots.h_p_truth_res_gbl.Fill(result.p_gbl(bfac) - result.track.p_truth(bfac))
      
      vtx_idx = 1 # first point is at s=0 (the "vtx" is at -670mm in test run)
      plots.h_vtx_xT_corr.Fill(result.xTCorr(vtx_idx))
      plots.h_vtx_yT_corr.Fill(result.yTCorr(vtx_idx))
      plots.h_d0_corr.Fill(result.d0Corr(vtx_idx))
      plots.h_z0_corr.Fill(result.z0Corr(vtx_idx))
      plots.h_d0.Fill(track.d0())
      plots.h_z0.Fill(track.z0())
      plots.h_d0_gbl.Fill(result.d0_gbl(vtx_idx))
      plots.h_z0_gbl.Fill(result.z0_gbl(vtx_idx))

      for label,corr in result.locPar.iteritems():
        if label>0:
          lbl = 2*(label-1) + 1
        else:
          lbl = -1*2*label
        plots.h_xT_corr.Fill(lbl, corr[result.idx_xT])
        plots.h_yT_corr.Fill(lbl, corr[result.idx_yT])
      
      

      for istrip in range(len(track.strips)):
        strip = track.strips[istrip]
          # find the label, if not found it's the vertex
        if strip in stripLabelMap:
          iLabel = stripLabelMap[strip]
        else:
          iLabel = 1

        #residuals 
        plots.h_res_layer.Fill(strip.layer,strip.ures)
        # correction to xT,yT from GBL fit
        corr = np.matrix( [result.locPar[iLabel][3], result.locPar[iLabel][4] ] )
        # project to measurement direction
        corr_meas = np.matrix( proL2m_list[strip.id] ) * np.transpose( np.matrix( corr ) )
        ures_gbl = strip.ures - corr_meas[0,0] # note minus sign due to definition of residual
        plots.h_res_gbl_layer.Fill(strip.layer,ures_gbl)
        
        # make plots for a given track only
        if nTry==0:
          plots.gr_ures.SetPoint(istrip,strip.pathLen,strip.ures)
          plots.gr_ures.SetPointError(istrip,0.,strip.ures_err)
          plots.gr_ures_truth.SetPoint(istrip,strip.pathLen,strip.uresTruth) 
          plots.gr_ures_simhit.SetPoint(istrip,strip.pathLen,strip.uresSimHit) 
          meas = np.array([strip.ures, 0.])
          #locRes = np.matrix(proM2l_list[strip.id]) *  np.transpose(np.matrix(meas))
          #xT_res = locRes[0,0]
          #yT_res = locRes[1,0]
          # find corrections to xT and yT
          plots.gr_corr_ures.SetPoint(istrip, strip.pathLen, corr_meas[0,0]) #u-direction
          ures_corr =  meas - corr_meas.T
          plots.gr_ures_corr.SetPoint(istrip, strip.pathLen, ures_corr[0,0]) #u-direction
      '''     
      nTry += 1

      
  #
  end = time.clock()
  print " Processed %d tracks " % nTry
  print " Time [s] ", end - start
  print " Chi2Sum/NdfSum ", Chi2Sum / NdfSum
  print " LostSum/nTry ", LostSum / nTry

  c = TCanvas('c','c',10,10,700,500)
  c.Divide(1,2)
  c.cd(1)
  h_chi2prob_initial_truth.Draw()
  c.cd(2)
  h_chi2prob_gbl_truth.Draw()
  ans = raw_input('kill...')
# 
  ''' 
  print " Make plots "
  savePlots = True
  plots.show(True)
  '''
 
def usage():
  print '%s: inputfile [-n nEventsMax] [-d|-dd debug level]' % sys.argv[0]
  sys.exit(0)


if __name__ == '__main__':
  '''
  if len(sys.argv)<2:
    usage()
  
  inputfile = sys.argv[1]
  for iparam in range(len(sys.argv)):
    s = sys.argv[iparam]
    if s=='-dd':
      debug=True
      utils.debug=True
    if s=='-d':
      debug=True
    if s=='-n':
      nEventsMax=int(sys.argv[iparam+1])
  '''
  inputfile = sys.argv[1]
  if len(sys.argv) > 2: 
    nEventsMax = int( sys.argv[2] )
  exampleHpsTest(inputfile)
