import numpy as np
import math
import time
import sys
import os
import utils
import hpseventst
import argparse
gblpythonpath = os.getenv('GBL','../GeneralBrokenLines/python')
sys.path.append(gblpythonpath)
from gblfit import GblPoint, GblTrajectory
import hps_plots
from ROOT import gDirectory, gROOT
from hps_utils import getAxialStereo

'''
General Broken Lines script for HPS with no B-field

Created Dec 2, 2015

@author: Per Hansson Adrian <phansson@slac.stanford.edu>
'''

def getArgs():
  parser = argparse.ArgumentParser(description='Run HPS GBL code')
  parser.add_argument('file',help='Input file.')
  parser.add_argument('--debug','-d',action='store_true',help='Debug output flag.')
  parser.add_argument('--ntracks','-n',type=int,default=-1,help='Max tracks to process.')
  parser.add_argument('--nevents',type=int,default=-1,help='Max events to process.')
  parser.add_argument('--name',help='Name to add to results')
  parser.add_argument('--save','-s',action='store_true',help='Save output')
  parser.add_argument('--nopause',action='store_true',help='Require manual input to continue program.')
  parser.add_argument('--batch','-b',action='store_true',help='Run ROOT in batch mode.')
  args = parser.parse_args();
  #print args
  return args



def main():
    nametag = os.path.splitext( os.path.basename(args.file) )[0]
    if args.name:
        nametag = nametag + "-" + args.name

    Chi2Sum = 0.
    NdfSum  = 0.
    LostSum  = 0.
    nTracks = 0

    start = time.clock()

    np.random.seed(47117)

    # Open binary file
    binaryFile = None
    if args.save:
        binaryFileName = "milleBinaryISN" + "_" + nametag
        binaryFile = open("%s.dat" % binaryFileName, "wb")
    
    # Open input file
    events = hpseventst.readHPSEvents(args.file, args.nevents, args.ntracks)
    
    print 'Read %d events from file' % len(events)

    plotsTopBot = hps_plots.plotter(nametag,'pdf',False,True,True, False)
    plotsTop = hps_plots.plotter(nametag,'pdf',False,True,False, False)
    plotsBot = hps_plots.plotter(nametag,'pdf',False,False,True, False)


    # loop over all events
    for event in events:

        if args.ntracks > 0 and nTracks > args.ntracks:
            break

        if args.debug or nTracks % 1000 == 0:
            print 'Processed %d tracks, now at event id %d with %d tracks ' % (nTracks, event.id, len(event.tracks))


        # loop over all tracks in the event
        for track in event.tracks:

            if args.debug:
                print 'track %d in event %d has %d strips' % (track.id, event.id,len(track.strips))

            # create the GBL trajectory
            traj = GblTrajectory(False)
            
            # point-to-point Jacobian
            jacPointToPoint = np.eye(5)

            # store projections for later use
            proL2m_list = {} 

            # start trajectory at reference point
            # this defines zero path length (s=0)
            point = GblPoint(jacPointToPoint)
            iLabelRef = traj.addPoint(point)

            # save mapping between label and strip object 
            stripLabelMap = {}

            # track direction in global frame
            tDirGlobal = track.direction()

            sinLambda = math.sin( track.clPar[1] )
            sinPhi = math.sin( track.clPar[2] )
            cosLambda = math.sqrt(1 - sinLambda ** 2)
            cosPhi = math.sqrt(1 - sinPhi ** 2)

            if args.debug:
                print 'tDir', tDirGlobal
                print 'lambda ', track.clPar[1], ' phi ', track.clPar[2]
                print 'sinLambda ', sinLambda, ' sinPhi ', sinPhi
            

            # path length
            s = 0.


            # loop over all strip clusters on the track
            for strip in track.strips:
                
                if args.debug:
                    print 'Processing strip id %d, millepedeId %d on sensor %s ' % (strip.id, strip.millepedeId,strip.deName)

                # calculate step from previous point
                step = strip.pathLen - s

                if args.debug:
                    print 'step ', step , '     (path length ',strip.pathLen, ' s ', s, ')'
                

                # Find the projection from tracking to measurement frame
                mDir = np.array( [strip.u, strip.v] )

                if args.debug:
                    print 'u: ',strip.u
                    print 'v: ',strip.v
                    print 'w: ',strip.w
                    print 'mDir:\n', mDir

                # Find the projection from curvilinear to measurement frame

                # Track direction in the curvilinear frame (U,V,T)
                # U = Z x T / |Z x T|, V = T x U
                uvDir = np.array([[-sinPhi, cosPhi, 0.], \
                                    [-sinLambda * cosPhi, -sinLambda * sinPhi, cosLambda]])


                if args.debug:
                    print 'Track direction in curvilinear frame\n',uvDir

                # projection from  measurement to local (curvilinear uv) directions (duv/dm)
                proM2l = np.dot(uvDir, mDir.T)

                # projection from local (uv) to measurement directions (dm/duv)
                proL2m = np.linalg.inv(proM2l)
                proL2m_list[strip.id] = proL2m

                if args.debug: 
                  print 'proM2l:\n', proM2l
                  print 'proL2m:\n', proL2m


                    
                # first get the projection from curvilinear to XYZ (or "global") frame
                #prjGlobalToCl = track.perToClPrj
                #prjClToGlobal = np.linalg.inv(prjGlobalToCl)

                #if args.debug:
                #    print 'prjGlobalToCl\n',prjGlobalToCl
                #    print 'prjClToGlobal\n',prjClToGlobal

                # now get the projection from global frame to measurement (or "sensor") frame
                # in HPS we call the measurement frame the "local sensor frame" and denote by "uvw"
                #prjGlobalToMeas = np.array( [strip.u, strip.v, strip.w] )
            
                #if args.debug:
                #    print 'prjGlobalToMeas\n', prjGlobalToMeas
                #    print 'prjGlobalToMeas.T\n', prjGlobalToMeas.T

                # now calculate the projection of the curvilinear to measurement frame
                # this is the so-called local to measurement transformation
                #prjClToMeas = np.dot( prjClToGlobal, prjGlobalToMeas.T) 
                #prjMeasToCl = np.linalg.inv( prjClToMeas )
                # adjust dimension of the projection, track direction coordinate is always 0
                #proL2m = prjClToMeas[:2,:2]        
                
                #if args.debug:
                #    print 'prjClToMeas\n', prjClToMeas
                #    #print 'prjMeasToCl\n', prjMeasToCl

                # residual and errors in measurement frame
                meas = np.array( [strip.ures, 0.] )
                measErr = np.array( [strip.ures_err, strip.ures_err] )
                measPrec = 1.0 / measErr ** 2
                measPrec[1] = 0. # no weight for measurement along the strip direction

                # Find the Jacobian to be able to propagate the covariance matrix to this strip position
                # "cosLambda" is the projection from the difference views
                # note that for this Jacobian the cosLambda only enters if there is a B-field
                jacPointToPoint = utils.gblSimpleJacobian(step, cosLambda, 0.)

                # Create a GBL point
                point = GblPoint(jacPointToPoint)


                if args.debug: 
                    print 'meas ', meas, ' measErr ', measErr, ' measPrec ', measPrec
                    print 'jacPointToPoint\n', jacPointToPoint
                    print 'proL2m \n', proL2m

                # Add a measurement to the point                
                point.addMeasurement([proL2m, meas, measPrec])

                # Add scatterer in curvilinear frame to the point
                # no direction in this frame
                scat = np.array([0., 0.])
                # Scattering angle in the curvilinear frame
                scatErr = np.array([ strip.scatAngle, strip.scatAngle / cosLambda]) 
                scatPrec = 1.0 / scatErr ** 2

                if args.debug: 
                    print 'scatPrec ', scatPrec, ' scatErr ', scatErr, ' cosLambda ', cosLambda
                
                point.addScatterer( [scat, scatPrec] )

                # Calculate global derivatives for this point
                # needs a few vectors in measurement frame

                # Projection matrix from tracking frame to measurement frame
                # t_u = dot(u,i)*t_i + dot(u,j)*t_j + dot(u,k)*t_k
                # where t_u is component in new u direction and t_i is in the i direction
                prjTrkToMeas = np.array([strip.u,strip.v,strip.w])
                # rotate to track direction to measurement frame          
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
                    print 'normalMeas ', normalMeas

                # non-measured coordinates
                vmeas = 0.
                wmeas = 0.

                # actually calculate the derivatives
                glDers = utils.globalDers(strip.millepedeId, strip.meas, vmeas, wmeas, tDirMeas, tPosMeas, normalMeas)
                
                if args.debug:
                    glDers.dump()

                # restructure to two arrays to fit interface
                ders = glDers.getDers( track.isTop() )
                labGlobal = ders['labels']
                addDer = ders['ders']
                if args.debug or (1==1 and \
                                  not track.isTop() and \
                                  getAxialStereo(strip.deName) == 'stereo' and \
                                  math.tan( track.clPar[1] ) > -0.011 and \
                                  math.tan( track.clPar[1] ) < -0.01 and \
                                  track.clPar[2] >0.01 and \
                                  track.clPar[2] <0.011):
                    
                    print '=== Global derivatives ==='
                    tanLambda = math.tan ( track.clPar[1] )
                    phi0 = track.clPar[2]
                    print 'tanLambda ', tanLambda, ' phi0 ', phi0
                    print 'track dir tracking frame     ', tDirGlobal
                    print 'track dir measurement frame  ', tDirMeas
                    print 'track pred tracking frame    ', strip.tPos
                    print 'track pred measurement frame ', tPosMeas
                    print 'deName ', strip.deName
                    print 'normalMeas ', normalMeas
                    print 'strip.u ', strip.u
                    print 'strip.ures ', strip.ures
                    print 'global derivatives for ', strip.deName, ' with id ', strip.id, ' and millepede id ', strip.millepedeId
                    print labGlobal.shape
                    for ider in range(labGlobal.shape[1]):
                        print '%7d %10.3e  %s' % (labGlobal[0][ider], addDer[0][ider], strip.deName)
                    print '====================== ==='

                # actually add the global derivatives to the point
                point.addGlobals(labGlobal, addDer)

                # add point to trajectory
                iLabel = traj.addPoint(point)
                
                # save strip and label map
                stripLabelMap[strip] = iLabel
            
                # go to next point
                s += step

                if args.debug:
                    print 'Done processing strip %d for track %d in event %d' % (strip.id, track.id, event.id)

            if args.debug:
                print 'Done adding points to trjacetory for track %d in event %d' % (track.id, event.id)


            if args.debug:
                print 'Do the fit'

            Chi2, Ndf, Lost = traj.fit()


            #if utils.chi2Prob(Chi2,Ndf) < 0.1:
            #    continue

            if args.save:
                traj.milleOut( binaryFile )

            # sum up fits
            Chi2Sum += Chi2
            NdfSum += Ndf
            LostSum += Lost            

            # get corrections and covariance matrix at points; collect the result in one object
            result = hpseventst.GBLTrajectory(track,traj)

            if nTracks < 2 or args.debug:
                print 'fit result: Chi2=%f Ndf=%d Lost=%d' % (Chi2, Ndf, Lost)
                result.dump()
            



            # loop over the two halves and the combined to get all plots
            
            for iplot in range(3):
                if iplot == 0:
                    plots = plotsTopBot
                elif iplot == 1 and track.isTop():
                    plots = plotsTop
                elif iplot == 2 and not track.isTop():
                    plots = plotsBot
                else:
                    continue

                plots.h_chi2.Fill(Chi2)
                plots.h_chi2ndf.Fill(Chi2/Ndf)
                plots.h_chi2prob.Fill(utils.chi2Prob(Chi2,Ndf))

                # loop over all the strips
                for strip in track.strips:
                    if strip not in stripLabelMap:
                        raise HpsGblException('this strip is not in the label map?!')
                    iLabel = stripLabelMap[strip]
                    point = iLabel

                    #residual for initial fit
                    plots.fillSensorPlots("res", strip.deName, strip.ures)

                    locPar, locCov = result.traj.getResults(point)
                    kinkLambda = result.kink(point,result.idx_lambda)
                    kinkPhi = result.kink(point,result.idx_phi)

                    plots.fillSensorPlots("corr_lambda",strip.deName, locPar[result.idx_lambda])
                    plots.fillSensorPlots("corrdiff_lambda",strip.deName, kinkLambda)

                    plots.fillSensorPlots("corr_phi",strip.deName, locPar[result.idx_phi])
                    plots.fillSensorPlots("corrdiff_phi",strip.deName, kinkPhi)

                    # correction to xT,yT from GBL fit
                    plots.fillSensorPlots("xTcorr",strip.deName, locPar[result.idx_xT])
                    plots.fillSensorPlots("yTcorr",strip.deName, locPar[result.idx_yT])

                    corr = np.matrix( [ locPar[result.idx_xT], locPar[result.idx_yT] ] )

                    # project to measurement direction
                    corr_meas = np.matrix( proL2m_list[strip.id] ) * np.transpose( np.matrix( corr ) )
                    ures_gbl = strip.ures - corr_meas[0,0] # note minus sign due to definition of residual
                    plots.fillSensorPlots("res_gbl", strip.deName, ures_gbl)
                    plots.fillSensorPlots("res_gbl_vs_vpred", strip.deName, [ures_gbl,tPosMeas[1]])
                    if abs(strip.meas) > 20.:
                        raise HpsGblException('really, this shouldnt happen? meas= ' + str(strip.meas))
                    plots.fillSensorPlots("res_gbl_vs_u", strip.deName, [ures_gbl, strip.meas] )


            nTracks += 1  

            if args.debug:
                print 'Done processing track %d in event %d' % (track.id, event.id)

    
    if binaryFile != None:
        if not binaryFile.closed:
            binaryFile.close()

    end = time.clock()
    print " Processed %d tracks " % nTracks
    print " Time [s] ", end - start
    if nTracks > 0:
        print " Chi2Sum/NdfSum ", Chi2Sum / NdfSum
        print " LostSum/nTracks ", LostSum / nTracks
        plotsTopBot.show(args.save,args.nopause)
        plotsTop.show(args.save,args.nopause)
        plotsBot.show(args.save,args.nopause)
        if args.save:
            hps_plots.saveHistosToFile(gDirectory,'gbltst-hps-plots-%s.root' % nametag)
        
    else:
        print 'No tracks processed'
    
    
    


if __name__ == '__main__':
    args = getArgs()

    gROOT.SetBatch(args.batch)

    if args.debug:
        utils.debug = True
        
    main()
