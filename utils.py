import sys, re, os
import numpy as np
from math import cos,atan,tan,pi,sin,asin,copysign
from ROOT import TMath
pyutilspath = os.getenv('PYTHONUTILS','pythonutils')
sys.path.append(pyutilspath)
import hps_utils


debug = False


class HpsGblException(Exception):
    def __init__(self, msg):
        Exception.__init__(self, msg)


def chi2Prob(chi2,ndf):
    return TMath.Prob(chi2,ndf)

    

def getBrackets(str):
    if not '<' in str:
        print 'Error no < in %s' % str
        sys.exit(1)
    str.split('<')[1]
        
            
            
        
def gblSimpleJacobianLambdaPhi(ds, cosl, bfac):
    '''
    Simple jacobian: quadratic in arc length difference.
    using lambda phi as directions

     curvilinear track parameter (q/p,lambda,phi,x_t,y_t)
    
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

    curvilinear local system (q/p,v',w',v,w), (v,w)=(x_t,y_t)
    
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
        global_params = {1:'u',2:'v',3:'w',4:'alpha',5:'beta',6:'gamma'}
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




def getMeasurementResidualIterative(perPar, origin, u ,w, meas, eps):
    '''Calculate the residual in the measurement direction for a set of track parameters.'''
    predIter = getXPlanePositionIterative(perPar,origin, w, eps)
    diffTrk = predIter - origin
    uPredIter = np.dot(u , diffTrk.T)
    uResIter = meas - uPredIter
    return uResIter





    
