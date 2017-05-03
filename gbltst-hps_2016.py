import numpy as np
import math
import time
import sys
import os
import utils
import hpsevent
import argparse
gblpythonpath = os.getenv('GBL','../GeneralBrokenLines/python')
sys.path.append(gblpythonpath)
from gblfit import GblPoint, GblTrajectory
import hps_plots
import simpleHelix
#import hps_utils
from ROOT import gROOT, gDirectory, TString

'''
Simple Test Program for General Broken Lines for HPS.

Created on Jul 27, 2011
Edited on Jun 20, 2013

@author: kleinwrt, phansson
'''



nametag = ''
 

def getArgs():
  parser = argparse.ArgumentParser(description='Run HPS GBL code')
  parser.add_argument('file',help='Input file.')
  parser.add_argument('--debug','-d',action='store_true',help='Debug output flag.')
  parser.add_argument('--nevents',type=int,default=-1,help='Max events to process.')
  parser.add_argument('--ntracks','-n',type=int,default=-1,help='Max tracks to process.')
  parser.add_argument('--notop',action='store_true',help='Reject top tracks.')
  parser.add_argument('--nobottom',action='store_true',help='Reject bottom tracks.')
  parser.add_argument('--name',help='Name to add to results')
  parser.add_argument('--mc','-m',action='store_true', help='Simulation input')
  parser.add_argument('--save','-s',action='store_true',help='Save output')
  parser.add_argument('--nopause',action='store_true',help='Require manual input to continue program.')
  parser.add_argument('--noshow',action='store_true',help='Do not show plots.')
  parser.add_argument('--useuncorrms',action='store_true',help='inflate MS errors instead of using scatterers')
  parser.add_argument('--testrun',action='store_true',help='Test Run input')
  parser.add_argument('--minStrips',type=int,default=0,help='Minimum number of strip clusters per track')
  parser.add_argument('--beamspot',action='store_true',help='Beamspot included as hit')
  parser.add_argument('--minP',type=float,help='Minimum track momentum in GeV/c')
  parser.add_argument('--batch','-b',action='store_true',help='Run ROOT in batch mode.')
  parser.add_argument('--elastic','-e',action='store_true',help='Select elastic events')
  parser.add_argument('--moller','-mol',action='store_true',help='Select Moller events')
  
  args = parser.parse_args();
  print args
  return args



def main(args):
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

  if args.save:
    binaryFileName = "milleBinaryISN" + "_" + nametag
    binaryFile = open("%s.dat" % binaryFileName, "wb")

  inputFile = open(args.file, 'r')
  events = hpsevent.readHPSEvents(inputFile, args.nevents, args.ntracks)
  
  print 'Read %d events from file' % len(events)

  if len(events) > 0:
    Bz = events[0].Bz
  bfac = 0.0002998 * Bz # for Bz in Tesla, momentum in GeV and Radius in mm
  print Bz, bfac
  
  plots = hps_plots.plotter(nametag,'pdf',args.testrun,True,True, args.beamspot)
  plotsTop = hps_plots.plotter(nametag,'pdf',args.testrun,True,False, args.beamspot)
  plotsBot = hps_plots.plotter(nametag,'pdf',args.testrun,False,True,args.beamspot)
  
  #print " GblHpsTest $Rev: 234 $ ", nTry, nLayer
  nTry = 0
  start = time.clock()


  # loop over events
  for event in events:
    
    # loop over tracks in event
    for track in event.tracks:

      if ((nTry % 500) == 0 and nTry > 0) or args.debug: 
        print '\nProcessed %d events: now on event id %d and track id %d' % (nTry, event.id, track.id)

      # if there's no truth info -> skip the track
      # use the track parameters and check curvature
      if args.mc and track.curvature_truth() == 0.:
        print 'Track curvature is zero for this MC track, skip track ', track.id, ' in event ', event.id, ' (perigee pars truth: ', track.perParTruth,')'
        continue

      # check if it is top or bottom
      if args.debug:
        print 'track with strip0 origin ', track.strips[0].origin

      if args.notop and track.isTop():
          if args.debug:
            print 'not ok'
          continue
      if args.nobottom and not track.isTop():
          if args.debug:
            print 'not ok'
          continue
      if args.debug:
        print 'ok'
      
      # check if it has enough hits
      if len(track.strips) < args.minStrips:
        if args.debug:
          print 'not enough strip clusters'
        continue

      if args.minP != None and track.p(bfac) < args.minP:
          continue

      # create the trajectory
      traj = GblTrajectory(True)

      # save mapping between label and strip object 
      stripLabelMap = {}
      
      if args.debug: 
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

      nhit = 0
      for strip in track.strips:
      
        if args.debug:
          print '\nProcessing strip id %d, millepedeId %d on sensor %s origin (%f,%f,%f)' % (strip.id, strip.millepedeId,strip.deName,strip.origin[0],strip.origin[1],strip.origin[2])
        
        step = strip.pathLen3D - s

        if args.debug: print 'Path length step %f from %f to %f ' % (step, s, strip.pathLen3D)
        
        # measurement direction (in YZ plane: perpendicular/parallel to strip direction)
        mDir = np.array( [strip.u, strip.v] )
        
        if args.debug: print 'mDir:\n', mDir
        
        
        # track direction: in x directon
        sinLambda = strip.sinLambda
        cosLambda = math.sqrt(1.0 - sinLambda ** 2)
        sinPhi = strip.sinPhi
        cosPhi = math.sqrt(1.0 - sinPhi ** 2)
        
        if args.debug: print 'Track direction sinLambda=%f sinPhi=%f' % (sinLambda, sinPhi)

        # Track direction in the curvilinear frame (U,V,T)
        # U = Z x T / |Z x T|, V = T x U
        uvDir = np.array([[-sinPhi, cosPhi, 0.], \
                            [-sinLambda * cosPhi, -sinLambda * sinPhi, cosLambda]])
        

        if args.debug: print 'Track direction in curvilinear frame\n',uvDir
        
        # projection from  measurement to local (curvilinear uv) directions (duv/dm)
        proM2l = np.dot(uvDir, mDir.T)

        # projection from local (uv) to measurement directions (dm/duv)
        proL2m = np.linalg.inv(proM2l)
        proL2m_list[strip.id] = proL2m

        if args.debug: 
          print 'proM2l:\n', proM2l
          print 'proL2m:\n', proL2m

        # measurement/residual in the measurement system
        meas = np.array([strip.ures, 0.]) # only measurement in u-direction
        #meas[0] += deltaU[iLayer] # misalignment
        measErr = np.array([strip.ures_err, strip.ures_err])
        measPrec = 1.0 / measErr ** 2
        measPrec[1] = 0. # 1D measurement perpendicular to strip direction
        
        if args.debug: 
          print 'meas ', meas, ' measErr ', measErr, ' measPrec ', measPrec
      
        # cross-check track position and residuals
        if nTry < 10:
          
          uResIter = utils.getMeasurementResidualIterative(track.perPar,strip.origin,strip.u,strip.w,strip.meas,1.0e-8)
          predIter = utils.getXPlanePositionIterative(track.perPar,strip.origin,strip.w,1.0e-8)
          diffTrk = predIter - strip.origin
          uPredIter = np.dot(strip.u , diffTrk.T)
          uResIter = strip.meas - uPredIter
          if abs(uResIter - strip.ures) > 1.0e-6:
            print 'WARNING diff %.10f uResIter %.10f compared to %.10f' % (uResIter - strip.ures,uResIter,strip.ures)
            #print 'predIter ', predIter, ' origin ', strip.origin, ' diffTrk ',diffTrk,' u ', strip.u, ' diffTrk ',diffTrk.T
            sys.exit(1)
        

        # Find the Jacobian to be able to propagate the covariance matrix to this strip position
        jacPointToPoint = utils.gblSimpleJacobianLambdaPhi(step, cosLambda, abs(bfac))
        
        if args.debug: 
          print 'jacPointToPoint to extrapolate to this point:'
          print jacPointToPoint
        
        # propagate MS covariance matrix (in the curvilinear frame) to this strip position
        msCov = np.dot(jacPointToPoint, np.dot(msCov, jacPointToPoint.T))

        # Get the MS covariance for measurements in the measurement frame
        measMsCov = np.dot(proL2m, np.dot(msCov[3:, 3:], proL2m.T))
        
        # Plot the MS variance in the u-direction
        plots.h_measMsCov.Fill(float(strip.id+1),measMsCov[0,0])
        if track.isTop():
          plotsTop.h_measMsCov.Fill(float(strip.id+1),measMsCov[0,0])
        else:
          plotsBot.h_measMsCov.Fill(float(strip.id+1),measMsCov[0,0])
        
        if args.debug:
          print 'msCov at this point:'
          print msCov
          print 'measMsCov at this point:'
          print measMsCov
        
        # Option to blow up measurement error according to multiple scattering
        if args.useuncorrms:
          measPrec[0] = 1.0 / (measErr[0] ** 2 + measMsCov[0, 0])
          if args.debug:
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
        if not args.useuncorrms:
          point.addScatterer([scat, scatPrec])
          if args.debug:
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


        # Projection matrix from tracking frame to measurement frame
        # t_u = dot(u,i)*t_i + dot(u,j)*t_j + dot(u,k)*t_k
        # where t_u is component in new u direction and t_i is in the i direction
        prjTrkToMeas = np.array([strip.u,strip.v,strip.w])
        # rotate to measurement frame          
        tDirMeas = np.dot( prjTrkToMeas, tDirGlobal.T) 
        normalMeas = np.dot( prjTrkToMeas, np.array([strip.w]).T )
        # vector coplanar with measurement plane from origin to prediction
        tDiff = np.array( [strip.tPos]) - np.array( [strip.origin] )
        # rotate to measurement frame          
        tPosMeas = np.dot( prjTrkToMeas, tDiff.T) 

        if args.debug: 
          print 'tDirGlobal ', tDirGlobal
          print 'rotation matrix to meas frame\n', prjTrkToMeas
          print 'tPosGlobal ', np.array( [strip.tPos]) , ' (origin ', np.array( [strip.origin] ),')'
          print 'tDiff ', tDiff
          print 'tPosMeas ', tPosMeas
        plots.fillSensorPlots("pred_meas", strip.deName, tPosMeas)
        if track.isTop():
          plotsTop.fillSensorPlots("pred_meas", strip.deName, tPosMeas)
        else:
          plotsBot.fillSensorPlots("pred_meas", strip.deName, tPosMeas)


# y vs x hit position in the sensor rf
        values = np.array([tPosMeas[0,0], tPosMeas[1,0]])
        plots.fillSensorPlots("xy", strip.deName, tPosMeas)
        if(track.qOverP(bfac)<0):
#            print strip.deName, event.id, track.id
            plots.fillSensorPlots("xy_neg", strip.deName, tPosMeas)
        else:
            plots.fillSensorPlots("xy_pos", strip.deName, tPosMeas)
        plots.fillSensorPlots("tpos", strip.deName, strip.tPos)
        plots.fillSensorPlots("hit3D", strip.deName, strip.hitPos3D)
          
        # rotate track direction to measurement frame          
        # non-measured directions         
        vmeas = 0.
        wmeas = 0.

        # calculate and add derivatives to point
        glDers = utils.globalDers(strip.millepedeId,strip.meas,vmeas,wmeas,tDirMeas,tPosMeas,normalMeas)
        if args.debug:
          glDers.dump()
        ders = glDers.getDers(track.isTop())
        labGlobal = ders['labels']
        addDer = ders['ders']
        if args.debug:
          print 'global derivatives for strip ', strip.id, ' which has millepede id ', strip.millepedeId
          print labGlobal.shape
          for ider in range(labGlobal.shape[1]):
            print labGlobal[0][ider], '\t', addDer[0][ider]
        point.addGlobals(labGlobal, addDer)
        ##### 
        
        # add point to trajectory
        iLabel = traj.addPoint(point)

        #if nTry==0:
        #print 'uRes ', strip.id, ' uRes ', strip.ures, ' pred ', tPosMeas, ' s(3D) ', strip.pathLen3D
        #print 'uRes ', strip.id, ' uRes ', strip.ures, ' pred ', strip.tPos, ' s(3D) ', strip.pathLen3D
        
        # go to next point
        s += step

        # save strip and label map
        stripLabelMap[strip] = iLabel
        nhit += 1

      if args.debug: print 'Do the fit'
      Chi2, Ndf, Lost = traj.fit()

      # write to millepede
      if args.save:
        traj.milleOut(binaryFile)

      # sum up    
      Chi2Sum += Chi2
      NdfSum += Ndf
      LostSum += Lost

      if nTry == 0 or args.debug:
        print 'fit result: Chi2=%f Ndf=%d Lost=%d' % (Chi2, Ndf, Lost)
      

      # get corrections and covariance matrix at points; collect the result in one object
      result = hpsevent.GBLResults(track)

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
        if nTry == -1:
          print " Point> ", i
          print " locPar ", locPar
          print " locCov ", locCov
        result.addPoint(i, locPar, locCov)
      
      if nTry == 0 or args.debug:
        result.printVertexCorr()
        result.printCorrection()
      
      # elastic cut selection
      if args.elastic:
        if(result.p_gbl(bfac) <2.0 or result.p_gbl(bfac)>2.6):
          continue
         
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
      for iplot in range(3):
        plot = None
        if iplot==0:
          plot = plots
        elif iplot==1 and track.isTop():
          plot = plotsTop
        elif iplot==2 and not track.isTop():
          plot = plotsBot

        if plot is None:
          continue

        # reject some tracks
        #if abs(track.d0())<2.0:
        #  continue
        #if(track.q<0):
        #  continue


        for i in range(0, traj.getNumPoints()-1):
          if(track.isTop()):
            plot.h_top_occupancy.Fill(track.strips[i].millepedeId)
          else:
            plot.h_bot_occupancy.Fill(track.strips[i].millepedeId)
          
        plot.h_nhit.Fill(nhit)
        plot.h_clPar_initial_xT.Fill(track.clPar[3])
        plot.h_clPar_initial_yT.Fill(track.clPar[4])
        plot.h_xT_vs_slope.Fill(track.slope(),track.clPar[3]) 
        plot.h_yT_vs_slope.Fill(track.slope(),track.clPar[4])
        if(track.slope()):
          plot.h_zTar.Fill(track.clPar[4]/track.slope())
        if track.isTop():
          plot.h_xT_vs_slope_top.Fill(track.slope(),track.clPar[3]) 
          plot.h_yT_vs_slope_top.Fill(track.slope(),track.clPar[4]) 
          plot.h_yT_vs_xT_top.Fill(track.clPar[3],track.clPar[4]) 
          if(track.slope()):
            plot.h_zTar_top.Fill(track.clPar[4]/track.slope())
        else:
          plot.h_xT_vs_slope_bot.Fill(track.slope(),track.clPar[3]) 
          plot.h_yT_vs_slope_bot.Fill(track.slope(),track.clPar[4]) 
          plot.h_yT_vs_xT_bot.Fill(track.clPar[3],track.clPar[4]) 
          if(track.slope()):
            plot.h_zTar_bot.Fill(track.clPar[4]/track.slope())
        
        plot.h_clPar_initial_qOverP.Fill(track.clPar[0])
        plot.h_clPar_initial_lambda.Fill(track.clPar[1])
        # transform phi to plot nicer
        if track.clPar[2]<math.pi:
          plot.h_clPar_initial_phi.Fill(track.clPar[2])
        else:
          plot.h_clPar_initial_phi.Fill(track.clPar[2]-math.pi*2)
        plot.h_clParGBL_res_qOverP.Fill(clParGBLRes[0])
        plot.h_clParGBL_res_lambda.Fill(clParGBLRes[1])
        plot.h_clParGBL_res_phi.Fill(clParGBLRes[2])
        plot.h_clParGBL_res_xT.Fill(clParGBLRes[3])
        plot.h_clParGBL_res_yT.Fill(clParGBLRes[4])
        if track.isTop():
          plot.h_clParGBL_res_xT_top.Fill(clParGBLRes[3])
          plot.h_clParGBL_res_yT_top.Fill(clParGBLRes[4])
        else: 
          plot.h_clParGBL_res_xT_bot.Fill(clParGBLRes[3])
          plot.h_clParGBL_res_yT_bot.Fill(clParGBLRes[4])        

        plot.h_clParGBL_pull_qOverP.Fill(clParGBLRes[0]/math.sqrt(math.fabs(result.locCov[1][0,0])))
        plot.h_clParGBL_pull_lambda.Fill(clParGBLRes[1]/math.sqrt(math.fabs(result.locCov[1][1,1])))
        plot.h_clParGBL_pull_phi.Fill(clParGBLRes[2]/math.sqrt(math.fabs(result.locCov[1][2,2])))
        plot.h_clParGBL_pull_xT.Fill(clParGBLRes[3]/math.sqrt(math.fabs(result.locCov[1][3,3])))
        plot.h_clParGBL_pull_yT.Fill(clParGBLRes[4]/math.sqrt(math.fabs(result.locCov[1][4,4])))

        plot.h_perPar_res_initial_d0.Fill(perParInitialVecRes[0,0])
        plot.h_perPar_res_initial_phi0.Fill(perParInitialVecRes[0,1])
        plot.h_perPar_res_initial_kappa.Fill(perParInitialVecRes[0,2])
        plot.h_perPar_res_initial_z0.Fill(perParInitialVecRes[0,3])
        plot.h_perPar_res_initial_slope.Fill(perParInitialVecRes[0,4])
        plot.h_chi2_initial.Fill(track.chi2Initial)
        if track.ndfInitial != 0: 
          plot.h_chi2ndf_initial.Fill(track.chi2Initial/track.ndfInitial)
          plot.h_chi2prob_initial.Fill(utils.chi2Prob(track.chi2Initial,track.ndfInitial))
        else:
          plot.h_chi2ndf_initial.Fill(0.)
          plot.h_chi2prob_initial.Fill(-1.)
        plot.h_chi2_initial_truth.Fill(chi2_initial_truth)
        plot.h_chi2ndf_initial_truth.Fill(chi2_initial_truth/5.0)
        plot.h_chi2prob_initial_truth.Fill(utils.chi2Prob(chi2_initial_truth,5))
        plot.h_chi2_gbl_truth.Fill(chi2_gbl_truth)
        plot.h_chi2ndf_gbl_truth.Fill(chi2_gbl_truth/5.0)
        plot.h_chi2prob_gbl_truth.Fill(utils.chi2Prob(chi2_gbl_truth,5))
        plot.h_chi2.Fill(Chi2)
        plot.h_chi2ndf.Fill(Chi2/Ndf)
        plot.h_chi2prob.Fill(utils.chi2Prob(Chi2,Ndf))
        plot.h_p.Fill(track.p(bfac))
#        print 'py', math.cos(track.theta())
        costhetax = math.sin(track.theta())*math.sin(track.phi0())
        py = track.p(bfac)*math.cos(track.theta())
#          if(abs(costhetax)<0.01 and py<0):
        if track.isTop():
          plot.h_p_top.Fill(track.p(bfac))
          plot.h_p_vs_costhetax_top.Fill(math.sin(track.theta())*math.sin(track.phi0()), track.p(bfac))
          plot.h_costhetay_vs_z0_top.Fill(track.z0(), math.cos(track.theta()))
        else: 
          plot.h_p_bot.Fill(track.p(bfac))
          plot.h_p_vs_costhetax_bot.Fill(math.sin(track.theta())*math.sin(track.phi0()), track.p(bfac))
          plot.h_costhetay_vs_z0_bot.Fill(track.z0(), math.cos(track.theta()))

        plot.h_qOverP.Fill(track.qOverP(bfac))
        if args.mc:
          plot.h_qOverP_truth_res.Fill(track.qOverP(bfac) - track.q()/track.p_truth(bfac))   
          plot.h_p_truth.Fill(track.p_truth(bfac))
          plot.h_p_truth_res.Fill(track.p(bfac)-track.p_truth(bfac))
          plot.h_p_truth_res_vs_p.Fill(track.p_truth(bfac),track.p(bfac)-track.p_truth(bfac))


        plot.h_qOverP_corr.Fill(result.curvCorr())
        plot.h_qOverP_gbl.Fill(result.qOverP_gbl(bfac))
        plot.h_p_gbl.Fill(result.p_gbl(bfac))
        if (track.isTop()):
          plot.h_p_gbl_top.Fill(result.p_gbl(bfac))
        else: 
          plot.h_p_gbl_bot.Fill(result.p_gbl(bfac))

        if args.mc:
          plot.h_qOverP_truth_res_gbl.Fill(result.qOverP_gbl(bfac) - result.track.qOverP_truth(bfac))
          plot.h_p_truth_res_gbl.Fill(result.p_gbl(bfac) - result.track.p_truth(bfac))
          plot.h_p_truth_res_gbl_vs_p.Fill(result.track.p_truth(bfac), result.p_gbl(bfac) - result.track.p_truth(bfac))


        vtx_idx = 1 # first point is at s=0 
        plot.h_vtx_xT_corr.Fill(result.xTCorr(vtx_idx))
        plot.h_vtx_yT_corr.Fill(result.yTCorr(vtx_idx))
        plot.h_d0_corr.Fill(result.d0Corr(vtx_idx))
        plot.h_z0_corr.Fill(result.z0Corr(vtx_idx))
        plot.h_d0_initial.Fill(track.d0())
        plot.h_z0_initial.Fill(track.z0())
        plot.h_d0_gbl.Fill(result.d0_gbl(vtx_idx))
        plot.h_z0_gbl.Fill(result.z0_gbl(vtx_idx))


        if args.debug:
          print 'curvCorr ', result.curvCorr(), ' xT_corr ', result.xTCorr(vtx_idx), ' yT_corr ', result.yTCorr(vtx_idx)
          print 'd0_corr ', result.d0Corr(vtx_idx), ' z0_corr ', result.z0Corr(vtx_idx)
          print 'd0_gbl ', result.d0_gbl(vtx_idx), ' (', result.track.d0(), ') z0_gbl ' , result.z0_gbl(vtx_idx), ' (', result.track.z0(), ')' 
          print 'locPar ', result.locPar[1]

        for label,corr in result.locPar.iteritems():
          if label>0:
            lbl = 2*(label-1) + 1
          else:
            lbl = -1*2*label
          plot.h_xT_corr.Fill(lbl, corr[result.idx_xT])
          plot.h_yT_corr.Fill(lbl, corr[result.idx_yT])

        for istrip in range(len(track.strips)):
          strip = track.strips[istrip]
            # find the label, if not found it's the vertex
          if strip in stripLabelMap:
            iLabel = stripLabelMap[strip]
          else:
            iLabel = 1

          prjTrkToMeas = np.array([strip.u,strip.v,strip.w])
          # rotate to measurement frame          
          tDirMeas = np.dot( prjTrkToMeas, tDirGlobal.T) 
          # vector coplanar with measurement plane from origin to prediction
          tDiff = np.array( [strip.tPos]) - np.array( [strip.origin] )
          # rotate to measurement frame          
          tPosMeas = np.dot( prjTrkToMeas, tDiff.T) 

          #residuals
          plot.fillSensorPlots("res", strip.deName, strip.ures)
          plot.fillSensorPlots("res_truth", strip.deName, strip.uresTruth)
          #track direction corrections
          point = istrip + 2
          plot.fillSensorPlots("corr_lambda",strip.deName, result.locPar[point][result.idx_lambda])
          plot.fillSensorPlots("corrdiff_lambda",strip.deName, result.locPar[point][result.idx_lambda]-result.locPar[point-1][result.idx_lambda])
          plot.fillSensorPlots("corr_phi",strip.deName, result.locPar[point][result.idx_phi])
          plot.fillSensorPlots("corrdiff_phi",strip.deName, result.locPar[point][result.idx_phi]-result.locPar[point-1][result.idx_phi])
          # correction to xT,yT from GBL fit
          corr = np.matrix( [result.locPar[iLabel][3], result.locPar[iLabel][4] ] )
          # project to measurement direction
          corr_meas = np.matrix( proL2m_list[strip.id] ) * np.transpose( np.matrix( corr ) )
          ures_gbl = strip.ures - corr_meas[0,0] # note minus sign due to definition of residual
          plot.fillSensorPlots("res_gbl", strip.deName, ures_gbl)
          plot.fillSensorPlots("res_gbl_vs_vpred", strip.deName, [ures_gbl,tPosMeas[1]])
          if abs(strip.meas) > 20.:
            print 'really, this shouldnt happen? ', strip.meas
            sys.exit(1)
          #if abs(tPosMeas[1]) > 50.:
          #  print 'really2? ', tPosMeas
          #  #sys.exit(1)
          plot.fillSensorPlots("res_gbl_vs_u", strip.deName, [ures_gbl, strip.meas] )
          plot.fillSensorPlots("iso", strip.deName, strip.iso)

          values = np.array([tPosMeas[0,0], ures_gbl])
          plot.fillSensorPlots("resVsPos", strip.deName, values)
          values = np.array([tPosMeas[1,0], ures_gbl])
          plot.fillSensorPlots("resUVsPosV", strip.deName, values)
# attach & exit angle
#          plot.fillSensorPlots("attachAngle", strip.deName, math.pi/2 - clParGBLRes[1])
          plot.fillSensorPlots("exitAngle", strip.deName, strip.scatAngle)

          # plot residuals of the seed vs the corrected seed
          if nTry < 999999:

            if args.debug: print '========= START DEBUG TRACK RESIDUALS ======== '

            # get the residual from the seed track perigee track parameters
            uResSeed = utils.getMeasurementResidualIterative(track.perPar,strip.origin,strip.u,strip.w,strip.meas,1.0e-8)

            # get the GBL corrections to the perigee track parameters at this point
            # NOTE: this is wrong!
            perParCorr = result.getPerParCorr(iLabel,bfac)
            if args.debug: print 'perPar     ', track.perPar
            if args.debug: print 'perParCorr ', perParCorr

            # get the residual from the *corrected* (see note above) seed track perigee track parameters
            uResSeedCorrWrong = utils.getMeasurementResidualIterative(perParCorr,strip.origin,strip.u,strip.w,strip.meas,1.0e-8)
            uResSeedCorrCmpWrong = abs(uResSeedCorrWrong) - abs(uResSeed)

            # plot the difference in residuals b/w the corrected and uncorrected track
            plot.fillSensorPlots("res_diff_wrong_gbl_seed", strip.deName, abs(uResSeedCorrWrong) - abs(uResSeed) )

            if args.debug: print 'WRONG diff ', uResSeedCorrCmpWrong, ' uResSeedCorrWrong', uResSeedCorrWrong, ' uResSeed ', uResSeed

            # This is the correct way of getting the corrected track parameters in perigee frame

            # Create a SimpleHelix object from the original seed track parameters
            # note that it uses slope instead of theta and difference ordering
            # [C,phi0,dca,slope,z0]
            helixSeed = simpleHelix.SimpleHelix([ track.perPar[0], track.perPar[2], track.perPar[3], math.tan(math.pi/2.0 - track.perPar[1]), track.perPar[4] ])

            if args.debug:
              print 'helixSeed '
              helixSeed.dump()

            # define reference points
            # global origin
            refPointAtOrg = [ 0., 0. ]
            # intersection of seed track with plane
            refPointAtPlane = [ strip.tPos[0],strip.tPos[1] ]

            if args.debug: print 'move helix to interception of seed track and plane which is at x,y ', refPointAtPlane, ' ( tPos ',strip.tPos,')'
            helixSeedParsAtPoint = helixSeed.moveToL3( refPointAtPlane )

            # create the helix at the new ref point
            helixSeedAtPoint = simpleHelix.SimpleHelix( helixSeedParsAtPoint, refPointAtPlane )

            if args.debug:
              print 'helixSeedAtPoint'
              helixSeedAtPoint.dump()
            
            # compare to the other propagation function
            if args.debug:
              helixSeedParsAtPointOther = helixSeed.moveTo( refPointAtPlane )
              helixSeedAtPointOther = simpleHelix.SimpleHelix( helixSeedParsAtPointOther, refPointAtPlane )
              print 'helixSeedAtPointOther'
              helixSeedAtPointOther.dump()
              print 'diff b/w propagations: ', np.array(helixSeedParsAtPoint) - np.array(helixSeedParsAtPointOther)
            
            # find the GBL corrections in perigee frame
            if args.debug: print 'get corrections in perFrame'
            perCorrections = result.getPerCorrections(iLabel, bfac)
            if args.debug: print 'perCorrection [C,theta,phi,d0,z0] ', perCorrections

            # get the corrections in the simpleHelix representation (slope instead of theta and different ordering)
            perCorrections = np.array( [ perCorrections[0], perCorrections[2], perCorrections[3], result.getSlopeCorrection(iLabel), perCorrections[4] ] )
            if args.debug: print 'perCorrection [C,phi,slope,d0,z0] ', perCorrections

            # cross-check the corrections in perigee frame with different formula
            perCorrectionsDiff = result.getPerCorrectionsOther(iLabel, bfac) - perCorrections
            #if args.debug: print 'perCorrectionOther ', perCorrectionsOther
            if args.debug: print 'Cross check perCorrections w/ different formula: ', perCorrectionsDiff
            if all( abs(i)>1e-5 for i in perCorrectionsDiff[:]):
                raise HpsGblException('Correction ', i, ' in perCorrections is different ', perCorrectionsDiff )
            
            # correct the perigee helix parameters at this point
            helixParsAtPoint = np.array( helixSeedAtPoint.getParameters() )            
            helixParsAtPointCorrected = helixParsAtPoint + perCorrections

            if args.debug:
              print 'helixParsAtPoint', helixParsAtPoint
              print 'perCorrections', perCorrections
              print 'helixParsAtPointCorrected', helixParsAtPointCorrected
            
            # create a SimpleHelix object from the corrected parameters
            helixCorrAtPoint = simpleHelix.SimpleHelix( helixParsAtPointCorrected, refPointAtPlane )

            if args.debug:
              print 'helixCorrAtPoint'
              helixCorrAtPoint.dump()

            # change reference point of the corrected helix to the original one
            #delta_refPointAtOrg = np.array( refPointAtOrg ) - np.array( refPointAtPlane )
            helixParsAtOrg = helixCorrAtPoint.moveToL3( refPointAtOrg )

            # create a SimpleHelix object from the corrected parameters at org
            helixCorr = simpleHelix.SimpleHelix( helixParsAtOrg, refPointAtOrg )

            if args.debug:
              print 'helixCorr'
              helixCorr.dump()

            # get the residual from the corrected track parameters at this point
            perParCorrSH = helixCorr.getParameters()
            if args.debug: print 'get residuals for parameters perParCorrSH ', perParCorrSH
            perParCorr = [ perParCorrSH[0], math.pi / 2.0 - math.atan( perParCorrSH[3] ), perParCorrSH[1], perParCorrSH[2], perParCorrSH[4] ]
            if args.debug: print ' and in [C,theta,phi,d0,z0] ', perParCorr
            uResSeedCorr = utils.getMeasurementResidualIterative(perParCorr,strip.origin,strip.u,strip.w,strip.meas,1.0e-8)
            uResSeedCorrCmp = abs(uResSeedCorr) - abs(uResSeed)

            # plot the difference in residuals b/w the corrected and uncorrected track
            plot.fillSensorPlots("res_diff_gbl_seed", strip.deName, abs(uResSeedCorr) - abs(uResSeed) )

            if args.debug: print 'CORR diff ', abs(uResSeedCorr) - abs(uResSeed), ' uResSeedCorr', uResSeedCorr, ' uResSeed ', uResSeed

            uResSeedCorrCmpVal = abs(uResSeedCorrCmp) - abs(uResSeedCorrCmpWrong)

            if args.debug: print 'CORR diff cmp ', uResSeedCorrCmpVal, ' uResSeedCorrCmp ', uResSeedCorrCmp, ' uResSeedCorrCmpWrong ', uResSeedCorrCmpWrong


            if args.debug: print '========= END DEBUG TRACK RESIDUALS ======== '


          # make plots for a given track only
          if nTry==0:
            plot.gr_ures.SetPoint(istrip,strip.pathLen3D,strip.ures)
            plot.gr_ures.SetPointError(istrip,0.,strip.ures_err)
            plot.gr_ures_truth.SetPoint(istrip,strip.pathLen3D,strip.uresTruth) 
            plot.gr_ures_simhit.SetPoint(istrip,strip.pathLen3D,strip.uresSimHit) 
            meass = np.array([strip.ures, 0.])
            #locRes = np.matrix(proM2l_list[strip.id]) *  np.transpose(np.matrix(meas))
            #xT_res = locRes[0,0]
            #yT_res = locRes[1,0]
            # find corrections to xT and yT
            plot.gr_corr_ures.SetPoint(istrip, strip.pathLen3D, corr_meas[0,0]) #u-direction
            ures_corr =  meass - corr_meas.T
            plot.gr_ures_corr.SetPoint(istrip, strip.pathLen3D, ures_corr[0,0]) #u-direction
      
      nTry += 1

      
  #
  end = time.clock()
  print " Processed %d tracks " % nTry
  print " Time [s] ", end - start
  if nTry > 0:
    print " Chi2Sum/NdfSum ", Chi2Sum / NdfSum
    print " LostSum/nTry ", LostSum / nTry
  print " Make plots "

  # Moller cut selection
  if args.moller:
    for event in events:
      if(len(event.tracks)!=2):
        continue
      track0 = event.tracks[0]
      track1 = event.tracks[1]
      sign0 = track0.qOverP(bfac)
      sign1 = track1.qOverP(bfac)
      if(sign0>0):
         continue
      if(sign1>0):
         continue
      p0 = track0.p(bfac)  
      p1 = track1.p(bfac)
      beamEnergy = 1.1
      if(abs(p0+p1- beamEnergy) < 0.1):
        print 'qui ho un moller!'
        p.h_thetay_vs_thetax_moller(math.acos(math.sin(track0.theta())*math.sin(track0.phi0())), track0.theta())
        p.h_thetay_vs_thetax_moller(math.acos(math.sin(track1.theta())*math.sin(track1.phi0())), track1.theta())

  if nTry > 0 and not args.noshow:
    plots.show(args.save,args.nopause)
    plotsTop.show(args.save,args.nopause)
    plotsBot.show(args.save,args.nopause)
    if args.save:
      hps_plots.saveHistosToFile(gDirectory,'gbltst-hps-plots-%s.root' % nametag)






if __name__ == '__main__':

  args = getArgs()

  gROOT.SetBatch(args.batch)


  if args.debug:
      hpsevent.debug = True
      utils.debug = True


      
  nametag = os.path.basename(args.file)
  if args.name:
    nametag = nametag + "-" + args.name

  if args.notop and args.nobottom:
    print 'Cannot reject both top and bottom'
    sys.exit(1)
  elif args.nobottom:
      nametag = nametag + "-top"
  elif args.notop:
      nametag = nametag + "-bot"
  

  main(args)

