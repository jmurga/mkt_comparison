from scipy.optimize import fsolve
from mpmath import gammainc
from scipy.special import polygamma
import mpmath
import numpy as np
class Analytical:

	def __init__(self,gam_neg=-40,N=250,theta_f=1e-5,theta_mid_neutral=1e-3,
                 alLow=0.2,alTot=0.2,gL=10,gH=200,al=0.184,be=0.000402,B=1.,
                 pposL=0.001,pposH=0.0,n=10,Lf=10**6,rho=0.001):
         self.gam_neg = gam_neg
         self.N = N
         self.NN = 2*N
         self.theta_f = theta_f
         self.theta_mid_neutral = theta_mid_neutral
         self.alLow = alLow
         self.alTot = alTot
         self.gL = gL
         self.gH = gH
         self.al = al
         self.be = be
         self.B = B
         self.pposL = pposL
         self.pposH = pposH
         self.n = n
         self.nn = 2*n
         self.Lf = Lf 


         self.rho = rho



	def fixNeut(self):
		return 0.25*(1./(self.B*self.NN))

	def fixNegB(self,ppos):
		return 0.75*(1-ppos)*(2**(-self.al))*(self.B**(-self.al))*(self.be**self.al)*(-mpmath.zeta(self.al,1.+self.be/(2.*self.B))+mpmath.zeta(self.al,0.5*(2-1./(self.N*self.B)+self.be/self.B)))

	def pFix(self,gamma):
		s = gamma/(self.NN+0.)
		pfix = (1.-mpmath.exp(-2.*s))/(1.-mpmath.exp(-2.*gamma))
		if s >= 0.1:
			pfix = mpmath.exp(-(1.+s))
			lim = 0
			while(lim < 200):
				pfix = mpmath.exp((1.+s)*(pfix-1.))
				lim +=1
			pfix = 1-pfix
		return pfix
		
	def fixPosSim(self,gamma,ppos):
		S = abs(self.gam_neg/(1.*self.NN))
		r = self.rho/(2.*self.NN)
		u = self.theta_f/(2.*self.NN)
		s = gamma/(self.NN*1.)
		p0 = polygamma(1,(s+S)/r)
		p1 = polygamma(1,(r+self.Lf*r+s+S)/r)
		#CC = 2*s/self.pFix(gamma)
		CC = 1.
		#print mpmath.exp(-2.*S*u*(p0-p1)/r**2)
		return 0.75*ppos*mpmath.exp(-2.*S*u*(p0-p1)*CC**2/r**2)*self.pFix(gamma)


	def alphaExpSimLow(self,pposL,pposH):
		return self.fixPosSim(self.gL,0.5*pposL)/(self.fixPosSim(self.gL,0.5*pposL)+self.fixPosSim(self.gH,0.5*pposH)+self.fixNegB(0.5*pposL+0.5*pposH))

	def alphaExpSimTot(self,pposL,pposH):
		return (self.fixPosSim(self.gL,0.5*pposL)+self.fixPosSim(self.gH,0.5*pposH))/(self.fixPosSim(self.gL,0.5*pposL)+self.fixPosSim(self.gH,0.5*pposH)+self.fixNegB(0.5*pposL+0.5*pposH))

	def solvEqns(self,params):
		pposL,pposH = params
		return (self.alphaExpSimTot(pposL,pposH)-self.alTot,self.alphaExpSimLow(pposL,pposH)-self.alLow)

	def setPpos(self):
		pposL,pposH =  fsolve(self.solvEqns,(0.001,0.001))
		#print self.alphaExpSimLow(pposL,pposH)
		#print self.alphaExpSimTot(pposL,pposH)
		if pposL < 0.:
			pposL = 0.
		if pposH < 0.:
			pposH = 0.
		self.pposH = pposH
		self.pposL = pposL



	def set_theta_f(self):
		theta_f  = fsolve(lambda theta: self.Br(self.Lf,theta)-self.B,0.00001)
		self.theta_f = theta_f[0]

	def Br(self,Lmax,theta):
		t = -1.*self.gam_neg/(self.NN+0.)
		u = theta/(2.*self.NN)
		r = self.rho/(2.*self.NN)
		#return np.exp(-2.*quad(lambda L: self.calcB(L,theta), 0., Lmax)[0])
		#print np.exp(-2.*t*u*(polygamma(1,(r+t)/r)-polygamma(1,1+Lmax+t/r))/r**2), np.exp(-u*2.*Lmax/(Lmax*r+2*t))
		return np.exp(-4*u*Lmax/(2*Lmax*r+t))