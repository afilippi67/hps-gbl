import sys
import numpy as np
import math
from hpsevent import HPSEvent
from hps_utils import getHalf
from utils import HpsGblException

debug = False





class Strip:
    """Simple strip cluster model."""
    

    def __init__(self,i,millepedeId,deName):
        self.id = i
        self.millepedeId = millepedeId
        self.deName = deName
        self.origin = []
        self.pathLen = 0.
        self.hitPos3D = []
        self.meas = 0.
        self.ures = 0.
        self.ures_err = 0.
        self.u = []
        self.v = []
        self.w = []
        self.scatAngle = 0.
        self.tDir = []
        self.tPos = []



class Track:
    """Simple track model."""

    # class variables
    idx_y0   = 0
    idx_z0   = 1
    idx_dydx = 2
    idx_dzdx = 3
    
    def __init__(self,tracknr):
        self.id = tracknr
        self.strips = []
        self.perPar = []
        self.clPar = []
        self.perToClPrj = None
    
    def direction(self):
        t = np.array( [ 1., self.perPar[Track.idx_dydx], self.perPar[Track.idx_dzdx]] )
        norm = np.linalg.norm(t)
        if norm == 0.:
            raise HpsGblException('can not normalize vector with no length.')
        #print 't ', t, ' norm ', norm, ' tNorm ', t/norm
        return t/norm

    def check(self):
        """Check consistency of this track. """

        # check that all hits are in the same half
        halfcount = {}
        for strip in self.strips:
            h = getHalf(strip.deName)
            if h in halfcount:
                halfcount[h] += 1
            else:
                halfcount[h] = 0
        if len(halfcount) != 1:
            raise HpsGblException('this track has strips from both sides?\n' + halfcount)

        # check that the correct nr of parameters exist
        if len(self.clPar) != 5 or len(self.perPar) != 4:
            raise HpsGblException('this track has wrong number of parameters?\n' + clPar + '\n' + perPar)

        # check dimenstion of transformation
        if self.perToClPrj == None:
            raise HpsGblException('this track has no CL projection matrix?')
        elif np.size( self.perToClPrj ) != 9:
            raise HpsGblException('this projection matrix has wrong dimenstions ' + np.size( self.perToClPrj ) )

        # check that direction of the track makes sense
        #t = self.direction()
        #theta = math.atan( 1. / math.tan( self.clPar[GBLTrajectory.idx_lambda] ) )
        #if theta < 0:
        #    theta = math.pi + theta
        #sinTheta = math.sin( theta )
        #cosTheta = math.sqrt ( 1 - sinTheta ** 2 )
        #sinPhi = math.sin(  self.clPar[GBLTrajectory.idx_phi] )
        #cosPhi = math.sqrt( 1 - sinPhi ** 2)
        #tTest = np.array( [ cosPhi*sinTheta, sinPhi*sinTheta, cosTheta] )
        #tDiff = t - tTest
        #print 'tDiff ', tDiff, ' t ', t, ' tTest ', tTest, ' lambda ',  self.clPar[GBLTrajectory.idx_lambda], ' theta ', theta, ' sinTheta ', sinTheta, ' cosTheta ', cosTheta, ' phi ',  self.clPar[GBLTrajectory.idx_phi], ' sinPhi ', sinPhi, ' cosPhi ' , cosPhi
        #if np.linalg.norm(tDiff) > 0.0001:
        #    raise HpsGblException('Track directions not consistent? ' + ' tDiff ' + np.array_str(tDiff) + ' t ' + np.array_str(t) + ' tTest ' + np.array_str(tTest))
    

        
    def isTop(self):
        """Check if this is a track in the top half of the detector."""

        if len(self.strips) == 0:
            raise HpsGblException('There are no strips on this track so this question makes no sense.')

        # the check method should have been called already so no point doing that. The below should be safe then.

        # get string representation
        h = getHalf(self.strips[0].deName)

        # check that it makes sense
        if h != 't' and h != 'b':
            raise HpsGblException('This half \"' + h + '\" is weird?')
        
        return h == 't'


    def prjClToPer(self,xT,yT):
        """Transform dimensions in curvilinear frame to perigee frame."""

        vec_in = np.array([xT,yT,0.])
        prj = np.linalg.inv(self.perToClPrj)
        vec = np.dot(prj,vec_in)
        return vec




class GBLTrajectory:
    """Encapsulates the trajectory and could be inherited I guess."""

    # class variables
    idx_qoverp = 0
    idx_lambda = 1
    idx_phi = 2
    idx_xT = 3
    idx_yT = 4
    
    def __init__(self,track, traj):
        self.track = track
        self.traj = traj

    def kink(self, point, idx):
        """Calculates the kink angle for a point."""

        kink = 0.
        if point > 1:
            locPar, locCov = self.traj.getResults(point)
            locParPrev, locCovPrev = self.traj.getResults(point-1)            
            kink = locPar[idx] - locParPrev[idx]
            if debug:
                print 'point ', point, ' kink ', kink, ' idx ', idx, ' locPar ', locPar[idx], ' locParPrev ', locParPrev[idx]
        return kink
    
    
    def dump(self):
        """Dump this trajectory."""
        
        print 'GBL trajectory:'
        for i in range(1, self.traj.getNumPoints() + 1):      
            # label start at 1
            locPar, locCov = self.traj.getResults(-i)
            print ' >Point ', i
            print ' locPar ', locPar
            print ' locCov ', locCov      
            locPar, locCov = self.traj.getResults(i)
            print ' Point> ', i
            print ' locPar ', locPar
            print ' locCov ', locCov
        print 'GBL kinks:'
        print '%8s %10s %10s' % ('Point','Lambda', 'Phi')
        for i in range(1, self.traj.getNumPoints() + 1):
            print '%8d %10f %10f' % (i, self.kink(i,GBLTrajectory.idx_lambda), self.kink(i,GBLTrajectory.idx_phi) )




def readHPSEvents(filename,nEventsMax):
    """Read and create the simple event model from a text file data format."""

    
    print 'Read max %d events from file' % nEventsMax
    events = []
    event = None
    track = None
    strip = None


    with open(filename, 'r') as f:

        for line in f:

            if debug:
                print 'Processing line \"%s\"' % line

            if 'New Event' in line:
                if event!=None:
                    if strip!=None:
                        track.strips.append(strip)
                        if debug:  
                            print 'Added strip %d millepedeLayer %d to list of strips for track %d in event %d' % (strip.id,strip.millepedeId,track.id,event.id)
                    if track!=None:
                        event.addTrack(track)
                        if debug:  
                            print 'Added track %d to list of tracks in event %d' % (track.id,event.id)
                    events.append(event)
                    if debug:  
                        print 'Added event %d to list' % event.id
                    if nEventsMax!=-1 and len(events)>=nEventsMax:
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
                    track.strips.append(strip)
                    if debug:  
                        print 'Added strip %d layer %d to list of strips for track %d in event %d' % (strip.id,strip.millepedeId,track.id,event.id)
                if track!=None:
                    event.addTrack(track)
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
                strip = Strip(nr,layer,deName)

            elif 'Track perPar (y0 z0 dydx dzdx)' in line:
                params = line.split('Track perPar (y0 z0 dydx dzdx)')[1].split()
                if len(params)!=4:
                    print 'Error: %d per params from \"%s\"' % (len(params),line)
                    sys.exit(1)
                for p in range(len(params)):
                    track.perPar.append(float(params[p]))

            elif 'Track clPar (q/p lambda phi xT yT)' in line:
                params = line.split('Track clPar (q/p lambda phi xT yT)')[1].split()
                if len(params)!=5:
                    print 'Error: %d cl params from \"%s\"' % (len(params),line)
                    sys.exit(1)
                for p in params:
                    track.clPar.append(float(p))

            elif 'Track clPrj' in line:
                params = line.split('Track clPrj')[1].split()
                if len(params)!=9:
                    print 'Error: %d clPrj from \"%s\"' % (len(params),line)
                    sys.exit(1)
                track.perToClPrj = np.eye(3)
                i=0
                for irow in range(0,3):
                    for icol in range(0,3):
                        track.perToClPrj[irow][icol] = float(params[i])
                        i=i+1

            #elif 'Strip jacobian' in line:
            #    params = line.split('Strip jacobian')[1].split()
            #    if len(params)!=16:
            #        print 'Error: %d jacobian from \"%s\"' % (len(params),line)
            #        sys.exit(1)
            #    strip.jacobian = np.eye(4)
            #    i=0
            #    for irow in range(0,4):
            #        for icol in range(0,4):
            #            strip.jacobian[irow][icol] = float(params[i])
            #            i=i+1
            
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

            elif 'Strip 3D hit pos' in line:
                params = line.split('Strip 3D hit pos')[1].split()
                if len(params)!=3:
                    print 'Error: %d strip 3D hit pos from \"%s\"' % (len(params),line)
                    sys.exit(1)
                for p in params:
                    strip.hitPos3D.append(float(p))

            elif 'Strip scatangle' in line:
                scat = float(line.split('Strip scatangle')[1])
                strip.scatAngle = scat

            elif 'Strip track pos' in line:
                params = line.split('Strip track pos')[1].split()
                if len(params)!=3:
                    print 'Error: %d strip track pos from \"%s\"' % (len(params),line)
                    sys.exit(1)
                for p in params:
                    strip.tPos.append(float(p))
            elif 'Strip origin' in line:
                params = line.split('Strip origin')[1].split()
                if len(params)!=3:
                    print 'Error: %d strip origin from \"%s\"' % (len(params),line)
                    sys.exit(1)
                for p in params:
                    strip.origin.append(float(p))

            else:
                print 'ERROR I should never get here line:\n%s' % line
                sys.exit()


    return events
