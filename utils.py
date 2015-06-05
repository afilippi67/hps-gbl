import sys
import numpy as np
from math import copysign,cos,atan,tan,pi,sin,asin,copysign
from ROOT import TMath


debug = False
global_params = {1:'u',2:'v',3:'w',4:'alpha',5:'beta',6:'gamma'}

def chi2Prob(chi2,ndf):
    return TMath.Prob(chi2,ndf)

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
        if copysign(1, self.strips[0].sinLambda) > 0:
            return True
        return False
    def pt(self,bfac):
        return TMath.Abs(bfac * 1.0/self.curvature())
    def p(self,bfac):
        p1 = self.pt(bfac) / sin(self.theta())
        #cross-check
        p2 = TMath.Abs(1.0/self.clPar[0])
        if (p1-p2)/p1 > 0.001:
            print 'Warning p1=%f and p2=%f!!' % ( p1, p2)
            #exit(1)
        return p1
    def pt_truth(self,bfac):
        return TMath.Abs(bfac * 1.0/self.curvature_truth())
    def p_truth(self,bfac):
        p1 = self.pt_truth(bfac) / sin(self.theta_truth())
        #cross-check
        p2 = TMath.Abs(1.0/self.clParTruth[0])
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
    def __init__(self,evtnr,Bz):
        self.id = evtnr
        self.Bz = Bz
        self.tracks = []

class GBLResults:
    def __init__(self,track):
        self.locPar = {}
        self.locCov = {}
        self.track = track
        self.idx_qoverp = 0
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
    

#def rotateGlToMeas(strip,vector):
#    rotMat = np.array([strip.u, strip.v, strip.w])
#    rotVector = np.dot(vector,rotMat)
#    return rotVector

    

def getBrackets(str):
    if not '<' in str:
        print 'Error no < in %s' % str
        sys.exit(1)
    str.split('<')[1]

def readHPSEvents(infile,nEventsMax):
    print 'Read max %d events from file' % nEventsMax
    events = []
    event = None
    track = None
    strip = None
    for line in infile:
        if debug: print 'Processing line \"%s\"' % line
        if 'New Event' in line:
            if event!=None:
                if strip!=None:
                #print 'Adding strip id %d ' % strip.id
                    track.strips.append(strip)
                    if debug:  
                        print 'Added strip %d millepedeLayer %d to list of strips for track %d in event %d' % (strip.id,strip.millepedeId,track.id,event.id)
                if track!=None:
                    event.tracks.append(track)
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
        else:
            print 'ERROR I should never get here line:\n%s' % line
            sys.exit()
    
    return events

        
            
            
        
def gblSimpleJacobianLambdaPhi(ds, cosl, bfac):
    '''
    Simple jacobian: quadratic in arc length difference.
    using lambda phi as directions
    
    @param ds: arc length difference
    @type ds: float
    @param cosl: cos(lambda)
    @type cosl: float
    @param bfac: Bz*c
    @type bfac: float
    @return: jacobian to move by 'ds' on trajectory
    @rtype: matrix(float)
       ajac(1,1)= 1.0D0
       ajac(2,2)= 1.0D0
       ajac(3,1)=-DBLE(bfac*ds)
       ajac(3,3)= 1.0D0
       ajac(4,1)=-DBLE(0.5*bfac*ds*ds*cosl)
       ajac(4,3)= DBLE(ds*cosl)
       ajac(4,4)= 1.0D0
       ajac(5,2)= DBLE(ds)
       ajac(5,5)= 1.0D0
    '''
    jac = np.eye(5)
    jac[2, 0] = -bfac * ds
    jac[3, 0] = -0.5 * bfac * ds * ds * cosl
    jac[3, 2] = ds * cosl
    jac[4, 1] = ds  
    return jac


def gblSimpleJacobian(ds, cosl, bfac):
    '''
    Simple jacobian: quadratic in arc length difference.
    
    @param ds: arc length difference
    @type ds: float
    @param cosl: cos(lambda)
    @type cosl: float
    @param bfac: Bz*c
    @type bfac: float
    @return: jacobian to move by 'ds' on trajectory
    @rtype: matrix(float)
    '''
    jac = np.eye(5)
    jac[1, 0] = -bfac * ds * cosl
    jac[3, 0] = -0.5 * bfac * ds * ds * cosl
    jac[3, 1] = ds
    jac[4, 2] = ds  
    return jac

class globalDers:
    def __init__(self,id,umeas,vmeas,wmeas, tDir, tPred, normal):
        self.millepedeId = id 
        self.umeas = umeas # measurement direction
        self.vmeas = vmeas # unmeasured direction
        self.wmeas = wmeas # normal to plane
        self.t = tDir # track direction
        self.p = tPred # track prediction
        self.n = normal # normal to plane
        # Global derivaties of the local measurements
        self.dm_dg = self.getMeasDers()
        # Derivatives of residuals w.r.t. measurement 
        self.dr_dm = self.getResDers()
        # Derivatives of residuals w.r.t. global parameters
        self.dr_dg = np.dot(self.dr_dm, self.dm_dg)
        #print 'dr_dm'
        #print self.dr_dm
        #print 'dm_dg'
        #print self.dm_dg
        #print 'dr_dg'
        #print self.dr_dg
    
    def dump(self):
        print 'globalDers:'
        print 'layer ', self.millepedeId
        print 'umeas ', self.umeas, ' vmeas ', self.vmeas, ' wmeas ', self.wmeas
        print 't ', self.t, ' p ', self.p, ' n ', self.n
        print 'dm_dg\n',self.dm_dg, '\ndr_dm\n',self.dr_dm,'\ndr_dg\n',self.dr_dg

    def getDers(self,isTop):
        half_offset = 10000 
        translation_offset = 1000
        direction_offset = 100
        topBot = 1
        transRot = 1
        direction = 1
        if not isTop:
            topBot = 2
        res = {}
        labels = []
        ders = []
        for ip, name in global_params.iteritems():
            if ip > 3:
                transRot = 2
                direction = ((ip-1) % 3) + 1
            else:
                direction = ip
            label = (int)(topBot * half_offset + transRot * translation_offset + direction * direction_offset + self.millepedeId)
            labels.append(label)
            ders.append(self.dr_dg[0,ip-1])
        
        return {'labels':np.array([labels]) , 'ders':np.array([ders])}
    
    

    
    def getResDers(self):
        # Derivatives of the local perturbed residual w.r.t. the measurements m (u,v,w)'
        tdotn = np.dot(self.t.T,self.n)
        drdg = np.eye(3)
        #print 't ', self.t, ' n ', self.n, ' dot(t,n) ', tdotn
        for i in range(3):
            for j in range(3):
                delta = 0.
                if i==j:
                    delta = 1.       
                drdg[i][j] = delta - self.t[i]*self.n[j]/tdotn[0]
        return drdg
    
    
    
    
    def getMeasDers(self):
        # Derivative of mt, the perturbed measured coordinate vector m 
        # w.r.t. to global parameters: u,v,w,alpha,beta,gamma

        # Derivative of the local measurement for a translation in u
        dmu_du = 1.
        dmv_du = 0.
        dmw_du = 0.
        # Derivative of the local measurement for a translation in v
        dmu_dv = 0.
        dmv_dv = 1.
        dmw_dv = 0.
        # Derivative of the local measurement for a translation in w
        dmu_dw = 0.
        dmv_dw = 0.
        dmw_dw = 1.
        # Derivative of the local measurement for a rotation around u-axis (alpha)
        dmu_dalpha = 0.
        dmv_dalpha = self.p[2] # self.wmeas
        dmw_dalpha = -1.0 * self.p[1] # -1.0 * self.vmeas
        # Derivative of the local measurement for a rotation around v-axis (beta)
        dmu_dbeta = -1.0 * self.p[2] #-1.0 * self.wmeas
        dmv_dbeta = 0.
        dmw_dbeta = self.p[0] #self.umeas
        # Derivative of the local measurement for a rotation around w-axis (gamma)
        dmu_dgamma = self.p[1] # self.vmeas
        dmv_dgamma = -1.0 * self.p[0]  # -1.0 * self.umeas 
        dmw_dgamma = 0.
        # put into matrix
        dmdg = np.array([[dmu_du, dmu_dv, dmu_dw, dmu_dalpha, dmu_dbeta, dmu_dgamma],[dmv_du, dmv_dv, dmv_dw, dmv_dalpha, dmv_dbeta, dmv_dgamma],[dmw_du, dmw_dv, dmw_dw, dmw_dalpha, dmw_dbeta, dmw_dgamma]])
        #print dmw_dbeta
        #dmdg = np.array([[dmu_du, dmu_dv],[dmu_dw, dmu_dalpha], [dmw_dbeta, dmw_dgamma]])
        return dmdg
    

def getHelixPathToX(parameters,x):
    dca = parameters[3]
    z0 = parameters[4]    
    phi0 = parameters[2]
    theta = parameters[1]    
    C = parameters[0]
    R = 1.0/C;
    
    xc = (R - dca) * sin(phi0)
    sinPhi = (xc - x)/R
    phi_at_x = asin(sinPhi)
    dphi_at_x = phi_at_x - phi0
    if dphi_at_x > pi: dphi_at_x -= 2.0*pi
    if dphi_at_x < -pi: dphi_at_x += 2.0*pi
    s_at_x = -1.0 * dphi_at_x * R    
    return s_at_x

def getHelixPosAtX(parameters,x):
    dca = parameters[3]
    z0 = parameters[4]    
    phi0 = parameters[2]
    theta = parameters[1]    
    C = parameters[0]
    R = 1.0/C;
    lam = pi/2.0 - theta
    slope = tan(lam)    
    #print '%.10f' % x
    #print '%.10f %.10f %.10f %.10f %.10f %.10f' % (dca, z0, phi0, slope, R, C)

    xc = (R - dca) * sin(phi0)
    sinPhi = (xc - x)/R
    phi_at_x = asin(sinPhi)
    dphi_at_x = phi_at_x - phi0
    if dphi_at_x > pi: 
        dphi_at_x -= 2.0*pi
    if dphi_at_x < -pi: 
        dphi_at_x += 2.0*pi
    s_at_x = -1.0 * dphi_at_x * R    
    y = dca * cos(phi0) - R * cos(phi0) + R * cos(phi_at_x)
    z = z0 + s_at_x * slope
    #print x,y,z,s_at_x,dphi_at_x,phi_at_x,sinPhi,xc
    return np.array([x,y,z])

def getXPlanePositionIterative(parameters,origin,normal,eps=0.0001):
    d0 = parameters[3]
    z0 = parameters[4]    
    phi0 = parameters[2]
    theta = parameters[1]    
    C = parameters[0]
    R = 1.0/C
    lam = pi/2.0 - theta
    slope = tan(lam)   
    #print 'Target origin ', origin, ' normal ', normal
    #print '%.10f %.10f %.10f %.10f %.10f %.10f' % (d0, z0, phi0, slope, R, C)
    #eps = 0.0001
    #print 'eps ', eps
    d = 9999.9
    x = origin[0]
    dx = 0.0
    nIter = 0
    pos = []
    while abs(d) > eps and nIter < 50:        
        # Calculate position on helix at x
        pos = getHelixPosAtX(parameters, x + dx)
        # Check if we are on the plane
        d = np.dot(pos-origin,np.array([normal]).T)
        # the direction to move depends on the direction of normal w.r.t. track
        # assume that the track moves in +x direction
        if np.sign(normal[0])==-1: 
            dx += 1.0*d[0]/2.0
        elif np.sign(normal[0])==1:
            dx += -1.0*d[0]/2.0
        else:
            print 'ERROR: normal is invalid, should move along x-axis'
            sys.exit(1)
        #print  nIter,' d ', d, ' pos ', pos, ' dx ', dx
        nIter +=  1
    return pos
