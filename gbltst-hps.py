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
import argparse
sys.path.append('../GeneralBrokenLines/python')
from gblfit import GblPoint, GblTrajectory
from hps_plots import plotter

#

debug = False
useUncorrMS = False # inflate MS errors instead of using scatterers
nEventsMax = -1
isTop=False
isBot=False
isTestRun = True
isMC = False
nametag = ''
savePlots = True;
 
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


  np.random.seed(47117)

  binaryFileName = "milleBinaryISN%s" % ("_%s" % nametag)
  if isTop: 
    binaryFileName += "_top"
  if isBot:
    binaryFileName += "_bot"
  binaryFile = open("%s.dat" % binaryFileName, "wb")

  inputFile = open(inputfile, 'r')
  events = utils.readHPSEvents(inputFile, nEventsMax)
  
  print 'Read %d events from file' % len(events)

  if len(events) > 0:
    Bz = events[0].Bz
  #Bz = -0.5 # full detector  -0.491 for test run detector
  bfac = 0.0002998 * Bz # for Bz in Tesla, momentum in GeV and Radius in mm
  print Bz, bfac
  
  plots = plotter(nametag,'pdf',isTestRun,isTop,isBot)
  
  #print " GblHpsTest $Rev: 234 $ ", nTry, nLayer
  nTry = 0
  start = time.clock()

  # loop over events
  for event in events:
    

    # loop over tracks in event
    for track in event.tracks:

      if ((nTry % 500) == 0 and nTry > 0) or debug: 
        print '\nProcessed %d events: now on event id %d and track id %d' % (nTry, event.id, track.id)

      # if there's no truth info -> skip the track
      if isMC and track.p_truth(bfac) == 0.:
        print 'No truth info, skip track %d in event %d ' % (track.id, event.id)
        continue

      # check if it is top or bottom
      if debug:
        print 'track with strip0 origin ', track.strips[0].origin
      if not isTop and not isBot:
        if debug:
          print 'track is not top and not bottom'
        continue
      else:
        if isBot and track.isTop():
          if debug:
            print 'not ok'
          continue
        if isTop and not track.isTop():
          if debug:
            print 'not ok'
          continue
      if debug:
        print 'ok'
      


      # create the trajectory
      traj = GblTrajectory(True)

      # save mapping between label and strip object 
      stripLabelMap = {}
      
      if debug: 
        print 'Track has %d strip clusters' % len(track.strips)

      # arc length
      s = 0.

      # point-to-point jacobian (from previous point)    
      jacPointToPoint = np.eye(5)

      #start trajectory at reference point (defining s=0)
      point = GblPoint(jacPointToPoint)
      refLabel = traj.addPoint(point)

      # multiple scattering covariance matrix (for curvilinear track parameters)
      msCov = np.zeros((5, 5))
      
      # store projections for later use
      proL2m_list = {} 

      # loop over strip clusters on the track

      for strip in track.strips:
        
        if debug: print '\nProcessing strip id %d at layer %d' % (strip.id, strip.layer)
        
        step = strip.pathLen3D - s

        if debug: print 'Path length step %f from %f to %f ' % (step, s, strip.pathLen3D)
        
        # measurement direction (in YZ plane: perpendicular/parallel to strip direction)
        mDir = np.array( [strip.u, strip.v] )
        
        if debug: print 'mDir:\n', mDir
        
        
        # track direction: in x directon
        sinLambda = strip.sinLambda
        cosLambda = math.sqrt(1.0 - sinLambda ** 2)
        sinPhi = strip.sinPhi
        cosPhi = math.sqrt(1.0 - sinPhi ** 2)
        
        if debug: print 'Track direction sinLambda=%f sinPhi=%f' % (sinLambda, sinPhi)

        # Track direction in the curvilinear frame (U,V,T)
        # U = Z x T / |Z x T|, V = T x U
        uvDir = np.array([[-sinPhi, cosPhi, 0.], \
                            [-sinLambda * cosPhi, -sinLambda * sinPhi, cosLambda]])
        
        
        # projection from  measurement to local (curvilinear uv) directions (duv/dm)
        proM2l = np.dot(uvDir, mDir.T)

        # projection from local (uv) to measurement directions (dm/duv)
        proL2m = np.linalg.inv(proM2l)
        proL2m_list[strip.id] = proL2m

        if debug: 
          print 'proM2l:\n', proM2l
          print 'proL2m:\n', proL2m

        # measurement/residual in the measurement system
        meas = np.array([strip.ures, 0.]) # only measurement in u-direction
        #meas[0] += deltaU[iLayer] # misalignment
        measErr = np.array([strip.ures_err, strip.ures_err])
        measPrec = 1.0 / measErr ** 2
        measPrec[1] = 0. # 1D measurement perpendicular to strip direction
        
        if debug: 
          print 'meas ', meas, ' measErr ', measErr, ' measPrec ', measPrec
      
        # cross-check track position and residuals
        if nTry < 10:
          predIter = utils.getXPlanePositionIterative(track.perPar,strip.origin,strip.w,1.0e-8)
          diffTrk = predIter - strip.origin
          uPredIter = np.dot(strip.u , diffTrk.T)
          uResIter = strip.meas - uPredIter
          if abs(uResIter - strip.ures) > 1.0e-6:
            print 'WARNING diff %.10f uResIter %.10f compared to %.10f' % (abs(uResIter - strip.ures),uResIter,strip.ures)
            print predIter,strip.origin,diffTrk,strip.u,diffTrk.T
            sys.exit(1)
        

        # Find the Jacobian to be able to propagate the covariance matrix to this strip position
        jacPointToPoint = utils.gblSimpleJacobianLambdaPhi(step, cosLambda, abs(bfac))
        
        if debug: 
          print 'jacPointToPoint to extrapolate to this point:'
          print jacPointToPoint
        
        # propagate MS covariance matrix (in the curvilinear frame) to this strip position
        msCov = np.dot(jacPointToPoint, np.dot(msCov, jacPointToPoint.T))

        # Get the MS covariance for measurements in the measurement frame
        measMsCov = np.dot(proL2m, np.dot(msCov[3:, 3:], proL2m.T))
        
        # Plot the MS variance in the u-direction
        plots.h_measMsCov.Fill(float(strip.layer),measMsCov[0,0])
        
        if debug:
          print 'msCov at this point:'
          print msCov
          print 'measMsCov at this point:'
          print measMsCov
        
        # Option to blow up measurement error according to multiple scattering
        if useUncorrMS:
          measPrec[0] = 1.0 / (measErr[0] ** 2 + measMsCov[0, 0])
          if debug:
            print 'Adding measMsCov ', measMsCov[0,0]
        
        # Create a GBL point
        point = GblPoint(jacPointToPoint)

        # Add measurement to the point
        point.addMeasurement([proL2m, meas, measPrec])

        # Add scatterer in curvilinear frame to the point
        # no direction in this frame
        scat = np.array([0., 0.])
        # Scattering angle in the curvilinear frame
        # Note the cosLambda to correct for the projection in the phi direction
        scatErr = np.array([ strip.scatAngle, strip.scatAngle / cosLambda]) 
        scatPrec = 1.0 / scatErr ** 2
        
        # add scatterer if not using the uncorrelated MS covariances
        if not useUncorrMS:
          point.addScatterer([scat, scatPrec])
          if debug:
            print 'adding scatError to this point:'
            print scatErr
        
        # Update MS covariance matrix 
        msCov[1, 1] += scatErr[0] ** 2; msCov[2, 2] += scatErr[1] ** 2
        
                
        ##### 
        ## Calculate global derivatives for this point
        # track direction in tracking/global frame
        tDirGlobal = np.array( [ [cosPhi * cosLambda, sinPhi * cosLambda, sinLambda] ] )        
        # Cross-check that the input is consistent
        if( np.linalg.norm( tDirGlobal - strip.tDir) > 0.00001):
          print 'ERROR: tDirs are not consistent!'
          sys.exit(1)
        # rotate track direction to measurement frame          
        tDirMeas = np.dot( tDirGlobal, np.array([strip.u, strip.v, strip.w]).T )
        #tDirMeas = utils.rotateGlToMeas(strip,tDirGlobal)
        normalMeas = np.dot( np.array([strip.w]) , np.array([strip.u, strip.v, strip.w]).T ) 
        #normalMeas = utils.rotateGlToMeas(strip,strip.w) 
        # non-measured directions 
        vmeas = 0.
        wmeas = 0.
        # calculate and add derivatives to point
        glDers = utils.globalDers(strip.layer,strip.meas,vmeas,wmeas,tDirMeas,normalMeas)
        ders = glDers.getDers(track.isTop())
        labGlobal = ders['labels']
        addDer = ders['ders']
        if debug:
          print 'global derivatives:'
          print labGlobal
          print addDer
        point.addGlobals(labGlobal, addDer)
        ##### 
        
        # add point to trajectory
        iLabel = traj.addPoint(point)

        if nTry==0:
          print 'uRes ', strip.id, ' uRes ', strip.ures, ' pred ', strip.tPos, ' s(3D) ', strip.pathLen3D
        
        # go to next point
        s += step

        # save strip and label map
        stripLabelMap[strip] = iLabel
      
      
      if debug: print 'Do the fit'
      Chi2, Ndf, Lost = traj.fit()

      # write to millepede
      traj.milleOut(binaryFile)

      # sum up    
      Chi2Sum += Chi2
      NdfSum += Ndf
      LostSum += Lost

      if nTry == 0 or debug:
        print 'fit result: Chi2=%f Ndf=%d Lost=%d' % (Chi2, Ndf, Lost)
      

      # get corrections and covariance matrix at points; collect the result in one object
      result = utils.GBLResults(track)

      #traj.dump()
      
      for i in range(1, traj.getNumPoints() + 1):      
        # label start at 1
        locPar, locCov = traj.getResults(-i)
        if nTry == 0:
          print " >Point ", i
          print " locPar ", locPar
          print " locCov ", locCov      
        result.addPoint(-i,locPar,locCov)
        locPar, locCov = traj.getResults(i)
        if nTry == 0:
          print " Point> ", i
          print " locPar ", locPar
          print " locCov ", locCov
        result.addPoint(i, locPar, locCov)
      


      # calculate the truth chi2 from initial fit
      # get the truth and fitted params with indexes same as cov matrix of initial fit (dca,phi0,curv,z0,slope)
      perParInitialVec = np.matrix([track.d0(), track.phi0(), track.curvature(), track.z0(), track.slope()])
      perParVecTruth = np.matrix([track.d0_truth(), track.phi0_truth(), track.curvature_truth(), track.z0_truth(), track.slope_truth()])
      perParInitialVecRes = perParInitialVec-perParVecTruth
      chi2_initial_truth = perParInitialVecRes * np.linalg.inv(track.perCov) * np.transpose(perParInitialVecRes)
      
      # calculate the truth chi2 from gbl fit at vertex
      clParGBLVtx = np.array(track.clPar) + np.array(result.locPar[1])
      clParTruth = np.array(track.clParTruth)
      clParGBLRes = clParGBLVtx - clParTruth
      chi2_gbl_truth = np.dot(clParGBLRes, np.dot(np.linalg.inv(result.locCov[1]), clParGBLRes))

      label = 1 #refLabel
      chi2_res = np.dot(result.locPar[label], np.dot(np.linalg.inv(result.locCov[label]), result.locPar[label]))
      if nTry == 0: 
        print " Chi2: ", event.id, chi2_res, chi2_gbl_truth, chi2_initial_truth
      




      # plots

      plots.h_clPar_initial_xT.Fill(track.clPar[3])
      plots.h_clPar_initial_yT.Fill(track.clPar[4])
      plots.h_clPar_initial_qOverP.Fill(track.clPar[0])
      plots.h_clPar_initial_lambda.Fill(track.clPar[1])
      # transform phi to plot nicer
      if track.clPar[2]<math.pi:
        plots.h_clPar_initial_phi.Fill(track.clPar[2])
      else:
        plots.h_clPar_initial_phi.Fill(track.clPar[2]-math.pi*2)
      plots.h_clParGBL_res_qOverP.Fill(clParGBLRes[0])
      plots.h_clParGBL_res_lambda.Fill(clParGBLRes[1])
      plots.h_clParGBL_res_phi.Fill(clParGBLRes[2])
      plots.h_clParGBL_res_xT.Fill(clParGBLRes[3])
      plots.h_clParGBL_res_yT.Fill(clParGBLRes[4])

      plots.h_clParGBL_pull_qOverP.Fill(clParGBLRes[0]/math.sqrt(math.fabs(result.locCov[1][0,0])))
      plots.h_clParGBL_pull_lambda.Fill(clParGBLRes[1]/math.sqrt(math.fabs(result.locCov[1][1,1])))
      plots.h_clParGBL_pull_phi.Fill(clParGBLRes[2]/math.sqrt(math.fabs(result.locCov[1][2,2])))
      plots.h_clParGBL_pull_xT.Fill(clParGBLRes[3]/math.sqrt(math.fabs(result.locCov[1][3,3])))
      plots.h_clParGBL_pull_yT.Fill(clParGBLRes[4]/math.sqrt(math.fabs(result.locCov[1][4,4])))
      
      plots.h_perPar_res_initial_d0.Fill(perParInitialVecRes[0,0])
      plots.h_perPar_res_initial_phi0.Fill(perParInitialVecRes[0,1])
      plots.h_perPar_res_initial_kappa.Fill(perParInitialVecRes[0,2])
      plots.h_perPar_res_initial_z0.Fill(perParInitialVecRes[0,3])
      plots.h_perPar_res_initial_slope.Fill(perParInitialVecRes[0,4])
      plots.h_chi2_initial.Fill(track.chi2Initial)
      if track.ndfInitial != 0: 
        plots.h_chi2ndf_initial.Fill(track.chi2Initial/track.ndfInitial)
        plots.h_chi2prob_initial.Fill(utils.chi2Prob(track.chi2Initial,track.ndfInitial))
      else:
        plots.h_chi2ndf_initial.Fill(0.)
        plots.h_chi2prob_initial.Fill(-1.)
      plots.h_chi2_initial_truth.Fill(chi2_initial_truth)
      plots.h_chi2ndf_initial_truth.Fill(chi2_initial_truth/5.0)
      plots.h_chi2prob_initial_truth.Fill(utils.chi2Prob(chi2_initial_truth,5))
      plots.h_chi2_gbl_truth.Fill(chi2_gbl_truth)
      plots.h_chi2ndf_gbl_truth.Fill(chi2_gbl_truth/5.0)
      plots.h_chi2prob_gbl_truth.Fill(utils.chi2Prob(chi2_gbl_truth,5))
      plots.h_chi2.Fill(Chi2)
      plots.h_chi2ndf.Fill(Chi2/Ndf)
      plots.h_chi2prob.Fill(utils.chi2Prob(Chi2,Ndf))
      plots.h_p.Fill(track.p(bfac))
      plots.h_qOverP.Fill(track.qOverP(bfac))
      if isMC:
        plots.h_qOverP_truth_res.Fill(track.qOverP(bfac) - track.q()/track.p_truth(bfac))   
        plots.h_p_truth.Fill(track.p_truth(bfac))
        plots.h_p_truth_res.Fill(track.p(bfac)-track.p_truth(bfac))
        plots.h_p_truth_res_vs_p.Fill(track.p_truth(bfac),track.p(bfac)-track.p_truth(bfac))
      
      
      plots.h_qOverP_corr.Fill(result.curvCorr())
      plots.h_qOverP_gbl.Fill(result.qOverP_gbl(bfac))
      plots.h_p_gbl.Fill(result.p_gbl(bfac))
      if isMC:
        plots.h_qOverP_truth_res_gbl.Fill(result.qOverP_gbl(bfac) - result.track.qOverP_truth(bfac))
        plots.h_p_truth_res_gbl.Fill(result.p_gbl(bfac) - result.track.p_truth(bfac))
        plots.h_p_truth_res_gbl_vs_p.Fill(result.track.p_truth(bfac), result.p_gbl(bfac) - result.track.p_truth(bfac))
      

      vtx_idx = 1 # first point is at s=0 (the "vtx" is at -670mm in test run)
      plots.h_vtx_xT_corr.Fill(result.xTCorr(vtx_idx))
      plots.h_vtx_yT_corr.Fill(result.yTCorr(vtx_idx))
      plots.h_d0_corr.Fill(result.d0Corr(vtx_idx))
      plots.h_z0_corr.Fill(result.z0Corr(vtx_idx))
      plots.h_d0_initial.Fill(track.d0())
      #print track.z0(), track.d0()
      plots.h_z0_initial.Fill(track.z0())
      plots.h_d0_gbl.Fill(result.d0_gbl(vtx_idx))
      plots.h_z0_gbl.Fill(result.z0_gbl(vtx_idx))

      if debug:
        print 'curvCorr ', result.curvCorr(), ' xT_corr ', result.xTCorr(vtx_idx), ' yT_corr ', result.yTCorr(vtx_idx)
        print 'd0_corr ', result.d0Corr(vtx_idx), ' z0_corr ', result.z0Corr(vtx_idx)
        print 'd0_gbl ', result.d0_gbl(vtx_idx), ' (', result.track.d0(), ') z0_gbl ' , result.z0_gbl(vtx_idx), ' (', result.track.z0(), ')' 
      
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
        plots.h_res_layer[strip.layer-1].Fill(strip.ures)
        plots.h_res_truth_layer[strip.layer-1].Fill(strip.uresTruth)
        # correction to xT,yT from GBL fit
        corr = np.matrix( [result.locPar[iLabel][3], result.locPar[iLabel][4] ] )
        # project to measurement direction
        corr_meas = np.matrix( proL2m_list[strip.id] ) * np.transpose( np.matrix( corr ) )
        ures_gbl = strip.ures - corr_meas[0,0] # note minus sign due to definition of residual
        plots.h_res_gbl_layer[strip.layer-1].Fill(ures_gbl)
        
        # make plots for a given track only
        if nTry==0:
          plots.gr_ures.SetPoint(istrip,strip.pathLen3D,strip.ures)
          plots.gr_ures.SetPointError(istrip,0.,strip.ures_err)
          plots.gr_ures_truth.SetPoint(istrip,strip.pathLen3D,strip.uresTruth) 
          plots.gr_ures_simhit.SetPoint(istrip,strip.pathLen3D,strip.uresSimHit) 
          meass = np.array([strip.ures, 0.])
          #locRes = np.matrix(proM2l_list[strip.id]) *  np.transpose(np.matrix(meas))
          #xT_res = locRes[0,0]
          #yT_res = locRes[1,0]
          # find corrections to xT and yT
          plots.gr_corr_ures.SetPoint(istrip, strip.pathLen3D, corr_meas[0,0]) #u-direction
          ures_corr =  meass - corr_meas.T
          plots.gr_ures_corr.SetPoint(istrip, strip.pathLen3D, ures_corr[0,0]) #u-direction
      
      nTry += 1

      
  #
  end = time.clock()
  print " Processed %d tracks " % nTry
  print " Time [s] ", end - start
  if nTry > 0:
    print " Chi2Sum/NdfSum ", Chi2Sum / NdfSum
    print " LostSum/nTry ", LostSum / nTry
  print " Make plots "
  if nTry > 0:
    plots.show(savePlots)


def getArgs():
  parser = argparse.ArgumentParser(description='Run HPS GBL code')
  parser.add_argument('file',help='Input file.')
  parser.add_argument('--debug','-d',type=int,default=0,help='Debug level.')
  parser.add_argument('--events','-n',type=int,default=-1,help='Max events to process.')
  parser.add_argument('--half',required=True,help='Top or bottom half tracks')
  parser.add_argument('--name',default='tmp',help='Name to add to results')
  parser.add_argument('--mc','-m',action='store_true',help='Simulation input')
  parser.add_argument('--save','-s',action='store_true',help='Save output')
  args = parser.parse_args();
  print args
  return args


if __name__ == '__main__':

  args = getArgs()

  nEventsMax = args.events
  nametag = args.name
  if args.debug > 0:
    debug=True
    if args.debug > 1:
      utils.debug = True
  if args.half == 'top' or args.half == 't':
    isTop=True
  if args.half == 'bottom' or args.half == 'b' or args.half == 'bot':
    isBot=True
  isMC = args.mc
  savePlots = args.save

  print 'top', isTop, ' bot ', isBot

  exampleHpsTest(args.file)

