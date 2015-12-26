import numpy as np
from hps_utils import getLayer
from math import copysign,cos,atan,tan,pi,sin,asin

debug = False

class Strip:
    def __init__(self,i,millepedeId,deName):
        self.id = i
        self.millepedeId = millepedeId
        self.deName = deName
        self.pathLen = 0.
        self.pathLen3D = 0.
        self.origin = []
        self.hitPos3D = []
        self.meas = 0.
        self.ures = 0.
        self.ures_err = 0.
        self.uresTruth = 0.
        self.uresTruth_err = 0.
        self.uresSimHit = 0.
        self.uresSimHit_err = 0.
        self.u = []
        self.v = []
        self.w = []
        self.scatAngle = 0.
        self.sinPhi = 0.
        self.sinLambda = 0.
        self.tDir = []
        self.tPos = []
        self.tPosMeas = []
        self.iso = 999999.9

class Track:
    def __init__(self,tracknr):
        self.id = tracknr
        self.strips = []
        self.perParDef = []
        self.perPar = []
        self.perParTruth = []
        self.clPar = []
        self.clParTruth = []
        self.perToClPrj = []
        self.clCov = []
        self.perCov = []
        self.chi2Initial = 0.
        self.ndfInitial = 0

    def isTop(self): 
        if copysign(1, self.strips[len(self.strips)-1].sinLambda) > 0:
            return True
        return False

    def pt(self,bfac):
        if self.curvature() == 0.:
            print 'curvature is zero? ', self.curvature(), ' ', self.perPar
            sys.exit(1)
            #print 'pt: bfac ', bfac, ' C ', self.curvature(), ' -> ', (bfac * 1.0/self.curvature()), ' abs(pt) ', TMath.Abs(bfac * 1.0/self.curvature()), ' ', abs(bfac * 1.0/self.curvature())
        return abs(bfac * 1.0/self.curvature())

    def p(self,bfac):
        p1 = self.pt(bfac) / sin(self.theta())
        if p1 == 0.0:
            print 'Warning p1 = ', p1, ' is zero, q/p params is ', self.clPar[0], ' pt ', self.pt(bfac),  ' sin(th)=', sin(self.theta())
        #cross-check
        p2 = abs(1.0/self.clPar[0])
        if (p1-p2)/p1 > 0.001:
            print 'Warning p1=%f and p2=%f!!' % ( p1, p2)
            #exit(1)
        return p1

    def pt_truth(self,bfac):
        if self.curvature_truth() == 0.:
            print 'curvature_truth is zero? ', self.curvature_truth(), ' ', self.perParTruth
            sys.exit(1)
        return abs(bfac * 1.0/self.curvature_truth())

    def p_truth(self,bfac):
        p1 = self.pt_truth(bfac) / sin(self.theta_truth())
        #cross-check
        p2 = abs(1.0/self.clParTruth[0])
        if (p1-p2)/p1 > 0.001:
            print 'Warning truth p1=%f and p2=%f!!' % (p1,p2)
            #exit(1)
        return p1

    def d0(self):
        return self.perPar[3]

    def d0_truth(self):
        return self.perParTruth[3]

    def z0(self):
        return self.perPar[4]

    def z0_truth(self):
        return self.perParTruth[4]

    def phi0(self):
        return self.perPar[2]

    def phi0_truth(self):
        return self.perParTruth[2]

    def q(self):
        if copysign(1,self.perPar[0])>0:
            return 1
        else:
            return -1
    def q_truth(self):
        if copysign(1,self.perParTruth[0])>0:
            return 1
        else:
            return -1
    def theta(self):
        return self.perPar[1]

    def theta_truth(self):
        return self.perParTruth[1]

    def slope(self):
        return tan(pi/2.0 - self.theta())

    def slope_truth(self):
        return tan(pi/2.0 - self.theta_truth())

    def qOverP(self,bfac):
        return float(self.q())/self.p(bfac)

    def qOverP_truth(self,bfac):
        return float(self.q_truth())/self.p_truth(bfac)

    def curvature(self):
        return self.perPar[0]

    def curvature_truth(self):
        return self.perParTruth[0]

    def prjClToPer(self,xT,yT):
        #print 'project UVT coordinates to IJK coordinates'
        vec_in = np.array([xT,yT,0.])
        prj = np.linalg.inv(self.perToClPrj)
        vec = np.dot(prj,vec_in)
        #print 'UVT coordinates ', vec_in, '\n perToClPrj \n', self.perToClPrj,'\n clToPerPrj \n', prj, ' \n to IJK (0,d0,z0)', vec
        return vec


class HPSEvent:
    """Simple event model."""
    def __init__(self,evtnr,Bz):
        self.id = evtnr
        self.Bz = Bz
        self.tracks = []

    def addTrack(self,track):
        """Add a track to list of tracks."""
        track.check()
        self.tracks.append(track)

class GBLResults:
    '''Class to encapsulate the GBL tracjectory and the seed track.'''

    def __init__(self,track):
        self.locPar = {}
        self.locCov = {}
        self.track = track
        self.idx_qoverp = 0
        self.idx_lambda = 1
        self.idx_phi = 2
        self.idx_xT = 3
        self.idx_yT = 4

    def addPoint(self,label,locPar,locCov):
        self.locPar[label] = locPar
        self.locCov[label] = locCov

    def curvCorr(self):
        label=1 # same correction for all points
        corr = self.locPar[label][self.idx_qoverp]
        return corr

    def p_gbl(self,bfac):
        return self.track.q()/self.qOverP_gbl(bfac)

    def qOverP_gbl(self,bfac):
        return self.track.qOverP(bfac) + self.curvCorr()

    def xTCorr(self,label):
            return self.locPar[label][self.idx_xT]

    def yTCorr(self,label):
            return self.locPar[label][self.idx_yT]

    def d0Corr(self,label):
        #correct for different sign convention of d0 in curvilinear frame
        corr = -1.0 * self.track.prjClToPer(self.xTCorr(label),self.yTCorr(label))[1]
        return corr

    def z0Corr(self,label):
        return self.track.prjClToPer(self.xTCorr(label),self.yTCorr(label))[2]

    def d0_gbl(self,label):
        return self.track.d0()+self.d0Corr(label) 

    def z0_gbl(self,label):
        return self.track.z0()+self.z0Corr(label)

    def printVertexCorr(self):
        print "qOverPCorr " + str(self.locPar[1][self.idx_qoverp]) + " xTCorr " + str(self.locPar[1][self.idx_xT]) + " yTCorr " +  str(self.locPar[1][self.idx_yT])  + " xtPrimeCorr " +  str(self.locPar[1][1])  + " yTPrimeCorr " +  str(self.locPar[1][2])

    def getCLParStr(self,point):
        s = ""
        for i in range(len(self.locPar[point])):
            s += "%5s %.10f " % (self.getIdxStr(i), self.locPar[point][i])
        return s

    def getIdxStr(self,i):
        s = ""
        if i==self.idx_qoverp:
            s = "q/p"
        elif i==self.idx_xT:
            s = "xT"
        elif i==self.idx_yT:
            s = "yT"
        elif i==self.idx_lambda:
            s = "lambda"
        elif i==self.idx_phi:
            s = "phi"
        else:
            print "ERROR, invalid idx", i, " quit"
            sys.exit(1)
        return s
    
    def printCorrection(self):
        print "\n==========="
        print "\nPrint corrections to curvilinear parameters on each side of scatterers"
        for i in range(1, len(self.locPar)/2 + 1):                  
            #for point,correction in self.locPar.iteritems():
            point = -i
            print "%3d %s" % (point, self.getCLParStr(point))
            point = i
            print "%3d %s" % (point, self.getCLParStr(point))
        print "\n==========="

    def getPerParCorr(self,label,bfac):
        '''Calculate new perigee frame track parameters - THIS IS WRONG'''

        # corrections to in perigee frame
        corrPer = self.track.prjClToPer(self.xTCorr(label),self.yTCorr(label))

        # corrections to phi
        xTPrimeCorr = self.locPar[label][self.idx_phi]

        # corrections to lambda
        yTPrimeCorr = self.locPar[label][self.idx_lambda]

        #calculate new slope
        lambda_gbl = atan( self.track.slope() ) + yTPrimeCorr
        slope_gbl = tan( lambda_gbl )
        theta_gbl = pi / 2.0 - atan( slope_gbl )

        # correction to curvature
        C_gbl = abs(bfac) * self.qOverP_gbl(bfac) / cos(lambda_gbl)

        # calculate new phi0
        phi0_gbl = self.track.phi0() + xTPrimeCorr - corrPer[0] * C_gbl

        return [ C_gbl, theta_gbl, phi0_gbl, self.d0_gbl(label), self.z0_gbl(label) ]

    def getSlopeCorrection(self,label):
        tmp_lambda_gbl = atan( self.track.slope() ) + self.locPar[label][self.idx_lambda]
        return tan( tmp_lambda_gbl )  - self.track.slope()
    
    def getPerCorrections(self,label,bfac):
        '''Calculate corrections in perigee frame'''

        # corrections to in perigee frame
        corrPer = self.track.prjClToPer(self.xTCorr(label),self.yTCorr(label))

        # corrections to phi
        xTPrimeCorr = self.locPar[label][self.idx_phi]

        # corrections to lambda
        yTPrimeCorr = self.locPar[label][self.idx_lambda]

        #calculate new theta
        #lambda_gbl = atan( self.track.slope() ) + yTPrimeCorr
        #slope_gbl = tan( lambda_gbl )
        slope_gbl = self.getSlopeCorrection(label) + self.track.slope()
        theta_gbl = pi / 2.0 - atan( slope_gbl )

        # correction to curvature
        C_gbl = abs(bfac) * self.qOverP_gbl(bfac) / cos( atan(slope_gbl) )

        # calculate new phi0
        phi0_gbl = self.track.phi0() + xTPrimeCorr - corrPer[0] * C_gbl

        # new parameters
        pars = [ C_gbl, theta_gbl, phi0_gbl, self.d0_gbl(label), self.z0_gbl(label) ]

        #print 'getPerCorrections: new perpars  ', np.array( pars )
        #print 'getPerCorrections: old perpars  ', np.array( self.track.perPar )

        # find the corrections from the difference compared to original parameters
        parsdiff = np.array( pars ) - np.array( self.track.perPar ) 

        #print 'getPerCorrections: diff perpars ', parsdiff

        return parsdiff


    
    def getPerCorrectionsOther(self,label,bfac):
        '''Calculate corrections in perigee frame'''
        locPars = self.locPar[label]
        la = atan( self.track.slope() )
        cosLambda = cos( la )
        sinLambda = sin( la )
        rInv = self.track.curvature()
        #print 'locPars ', locPars
        #print 'cosLambda ', cosLambda, ' sinLambda ', sinLambda, ' rInv ', rInv
        # transformation to perigee (-like) system
        # Pelle: note flipped sign for dca
        curv2per = np.zeros((5, 5))  
        curv2per[0][0] = -bfac / cosLambda
        curv2per[0][1] = sinLambda / cosLambda * rInv
        curv2per[3][1] = 1. / (cosLambda * cosLambda)  
        curv2per[1][2] = 1.  
        curv2per[1][4] = sinLambda * rInv  
        curv2per[2][3] = -1.  
        curv2per[4][4] = 1. / cosLambda
        aSolution = np.dot(curv2per, locPars)
        #print 'curv2per ', curv2per
        #print 'aSolution ', aSolution
        aSolution[2] = -1. * aSolution[2]
        #print 'aSolution dca fix ', np.array( [ aSolution ] )
        return aSolution




def readHPSEvents(infile,nEventsMax,nTracksMax):
    print 'Read max %d events or %d tracks from file' % (nEventsMax, nTracksMax)
    events = []
    event = None
    track = None
    strip = None
    ntracks = 0
    for line in infile:
        if debug: print 'Processing line \"%s\"' % line
        if 'New Event' in line:
            if event!=None:
                if strip!=None:
                    track.strips.append(strip)
                    if debug:  
                        print 'Added strip %d millepedeLayer %d to list of strips for track %d in event %d' % (strip.id,strip.millepedeId,track.id,event.id)
                if track!=None:
                    event.tracks.append(track)
                    if debug:  
                        print 'Added track %d to list of tracks in event %d' % (track.id,event.id)
                if len(event.tracks) > 0:
                        events.append(event)
                        ntracks += len( event.tracks )
                        if debug or (len(events) % 10000 == 0 and len(events) > 0):  
                            print 'Added event %d to list' % event.id
                        if debug or (ntracks % 1000 == 0 and ntracks > 0):  
                            print 'Added ', ntracks, ' so far'
                if nEventsMax!=-1 and len(events)>=nEventsMax:
                    return events
                if nTracksMax!=-1 and ntracks >= nTracksMax:
                    return events
                event = None
                track = None
                strip = None
            words = line.split('New Event')[1].split()
            if len(words) != 2: 
                print 'New event line is wrong: ', words
            event = HPSEvent(int(words[0]),float(words[1]))
        elif 'New Track' in line:
            if strip!=None:
                #print 'Adding strip id %d ' % strip.id
                track.strips.append(strip)
                if debug:  
                    print 'Added strip %d layer %d to list of strips for track %d in event %d' % (strip.id,strip.millepedeId,track.id,event.id)
            if track!=None:
                event.tracks.append(track)
                if debug:  
                    print 'Added track %d to list of tracks in event %d' % (track.id,event.id)
            track = None
            strip = None
            nr = int(line.split('New Track')[1])
            track = Track(nr)
        elif 'New Strip id layer' in line:
            if strip!=None:
                #print 'Adding strip id %d ' % strip.id
                track.strips.append(strip)
                if debug:  
                    print 'Added strip %d layer %d to list of strips for track %d in event %d' % (strip.id,strip.millepedeId,track.id,event.id)
            strip = None
            nr = int(line.split('New Strip id layer')[1].split()[0])
            layer = int(line.split('New Strip id layer')[1].split()[1])
            deName = line.split('New Strip id layer')[1].split()[2]
            #print 'check ', deName
            if getLayer(deName) == 0 and getHalf(deName)=='b' and track.slope()>0:
                print 'Fix deName for L0 for id ', nr, ' MP id ', layer, ' name ', deName
                deName = deName.replace('L0b_','L0t_')
                print 'new deName ' , deName
            strip = Strip(nr,layer,deName)
        elif 'Track perPar (R theta phi d0 z0)' in line:
            params = line.split('Track perPar (R theta phi d0 z0)')[1].split()
            if len(params)!=5:
                print 'Error: %d per params from \"%s\"' % (len(params),line)
                sys.exit(1)
            for p in range(len(params)):
                if p==0:
                    track.perPar.append(1.0/float(params[p]))
                else:
                    track.perPar.append(float(params[p]))
        elif 'Truth perPar (kappa theta phi d0 z0)' in line:
            params = line.split('Truth perPar (kappa theta phi d0 z0)')[1].split()
            if len(params)!=5:
                print 'Error: %d truth per params from \"%s\"' % (len(params),line)
                sys.exit(1)
            for p in params:
                track.perParTruth.append(float(p))
        elif 'Track perPar (R phi0 slope d0 z0)' in line:
            params = line.split('Track perPar (R phi0 slope d0 z0)')[1].split()
            if len(params)!=5:
                print 'Error: %d per params from \"%s\"' % (len(params),line)
                sys.exit(1)
            for p in params:
                track.perParDef.append(float(p))
        elif 'Track clPar (q/p lambda phi xT yT)' in line:
            params = line.split('Track clPar (q/p lambda phi xT yT)')[1].split()
            if len(params)!=5:
                print 'Error: %d cl params from \"%s\"' % (len(params),line)
                sys.exit(1)
            for p in params:
                track.clPar.append(float(p))
            # fix yT
            #if len(track.clPar)==5:
            #    track.clPar[4] = track.perPar[4] * cos( track.clPar[1] )
            # fix phi to be within -pi,pi
            #if track.clPar[2] > pi:
            #    track.clPar[2] -= 2 * pi
        elif 'Truth clPar (q/p lambda phi xT yT)' in line:
            params = line.split('Truth clPar (q/p lambda phi xT yT)')[1].split()
            if len(params)!=5:
                print 'Error: %d truth cl params from \"%s\"' % (len(params),line)
                sys.exit(1)
            for p in params:
                track.clParTruth.append(float(p))
        elif 'Track perToClPrj' in line:
            params = line.split('Track perToClPrj')[1].split()
            if len(params)!=9:
                print 'Error: %d perToClPrj from \"%s\"' % (len(params),line)
                sys.exit(1)
            track.perToClPrj = np.eye(3)
            i=0
            for irow in range(0,3):
                for icol in range(0,3):
                    track.perToClPrj[irow][icol] = float(params[i])
                    i=i+1
        elif 'Track clCov' in line:
            params = line.split('Track clCov')[1].split()
            if len(params)!=25:
                print 'Error: %d clCov from \"%s\"' % (len(params),line)
                sys.exit(1)
            track.clCov = np.eye(5)
            i=0
            for irow in range(0,5):
                for icol in range(0,5):
                    track.clCov[irow][icol] = float(params[i])
                    i=i+1
        elif 'Track perCov (idx: dca,phi0,curv,z0,slope)' in line:
            params = line.split('Track perCov (idx: dca,phi0,curv,z0,slope)')[1].split()
            if len(params)!=25:
                print 'Error: %d perCov from \"%s\"' % (len(params),line)
                sys.exit(1)
            track.perCov = np.eye(5)
            i=0
            for irow in range(0,5):
                for icol in range(0,5):
                    track.perCov[irow][icol] = float(params[i])
                    i=i+1
        elif 'Track chi2/ndf (circle,zfit)' in line:
            params = line.split('Track chi2/ndf (circle,zfit)')[1].split()
            if len(params)!=4:
                print 'Error: %d chi2 from \"%s\"' % (len(params),line)
                sys.exit(1)
            track.chi2Initial = float(params[0])+float(params[2])
            track.ndfInitial = int(params[1])+int(params[3])
        elif 'Strip pathLen3D ' in line:
            strip.pathLen3D = float(line.split('Strip pathLen3D ')[1])
        elif 'Strip pathLen ' in line:
            strip.pathLen = float(line.split('Strip pathLen ')[1])
        elif 'Strip meas dir' in line:
            params = line.split('Strip meas dir')[1].split()
            if len(params)!=3:
                print 'Error: %d meas dir from \"%s\"' % (len(params),line)
                sys.exit(1)
            strip.u = []
            for p in params:
                strip.u.append(float(p))
        elif 'Strip non-meas dir ' in line:
            params = line.split('Strip non-meas dir')[1].split()
            if len(params)!=3:
                print 'Error: %d non-meas dir from \"%s\"' % (len(params),line)
                sys.exit(1)
            strip.v = []
            for p in params:
                strip.v.append(float(p))
        elif 'Strip normal dir' in line:
            params = line.split('Strip normal dir')[1].split()
            if len(params)!=3:
                print 'Error: %d normal dir from \"%s\"' % (len(params),line)
                sys.exit(1)
            strip.w = []
            for p in params:
                strip.w.append(float(p))
        elif 'Strip ures' in line:
            params = line.split('Strip ures')[1].split()
            if len(params)!=2:
                print 'Error: %d ures from \"%s\"' % (len(params),line)
                sys.exit(1)
            strip.ures = float(params[0])
            strip.ures_err = float(params[1])
        elif 'Strip u ' in line:
            strip.meas = float( line.split('Strip u ')[1] )
        elif 'Strip origin pos' in line:
            params = line.split('Strip origin pos')[1].split()
            if len(params)!=3:
                print 'Error: %d strip origin from \"%s\"' % (len(params),line)
                sys.exit(1)
            for p in params:
                strip.origin.append(float(p))
        elif 'Strip truth ures' in line:
            params = line.split('Strip truth ures')[1].split()
            if len(params)!=2:
                print 'Error: %d truth ures from \"%s\"' % (len(params),line)
                sys.exit(1)
            strip.uresTruth = float(params[0])
            strip.uresTruth_err = float(params[1])
        elif 'Strip 3D hit pos' in line:
            params = line.split('Strip 3D hit pos')[1].split()
            if len(params)!=3:
                print 'Error: %d strip 3D hit pos from \"%s\"' % (len(params),line)
                sys.exit(1)
            for p in params:
                strip.hitPos3D.append(float(p))
        elif 'Strip simhit ures' in line:
            params = line.split('Strip simhit ures')[1].split()
            if len(params)!=2:
                print 'Error: %d simhit ures from \"%s\"' % (len(params),line)
                sys.exit(1)
            strip.uresSimHit = float(params[0])
            strip.uresSimHit_err = float(params[1])
        elif 'Strip scatangle' in line:
            scat = float(line.split('Strip scatangle')[1])
            strip.scatAngle = scat
        elif 'Strip sinPhi sinLambda' in line:
            uvDir = line.split('Strip sinPhi sinLambda')[1].split()
            if len(uvDir)!=2:
                print 'Error: %d tdir from \"%s\"' % (len(uvDir),line)
                sys.exit(1)
            strip.sinPhi = float(uvDir[0])
            strip.sinLambda = float(uvDir[1])
        elif 'Strip track dir' in line:
            params = line.split('Strip track dir')[1].split()
            if len(params)!=3:
                print 'Error: %d strip track dir from \"%s\"' % (len(params),line)
                sys.exit(1)
            for p in params:
                strip.tDir.append(float(p))
        elif 'Strip track pos' in line:
            params = line.split('Strip track pos')[1].split()
            if len(params)!=3:
                print 'Error: %d strip track pos from \"%s\"' % (len(params),line)
                sys.exit(1)
            for p in params:
                strip.tPos.append(float(p))
            #for p in params[3:]:
            #    strip.tPosMeas.append(float(p))
        elif 'Strip iso ' in line:
            strip.iso = float(line.split('Strip iso ')[1])
        else:
            print 'ERROR I should never get here line:\n%s' % line
            sys.exit()
    
    return events

