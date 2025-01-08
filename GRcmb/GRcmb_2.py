"calculate C_l"

import numpy as np
import timeit
from scipy.interpolate import interp1d as sp_interp1d
from scipy.interpolate import RectBivariateSpline
import matplotlib.pyplot as plt
import matplotlib
from scipy.special import spherical_jn
from scipy.integrate import trapezoid
from scipy.integrate import simpson




llp1c200 = 6.14403e-10  #for normalization


inputdata = np.genfromtxt('GR_input_data.txt')
inputns = inputdata[4]
#print(inputns)
 
 
t0 = timeit.time.time()
fontsize1 = 20  
 
   

class clcalc(object):
   """
   read and interpolate data 
   """

   def __init__(self, label=1, sol=0, lnacoanpt=4000, kcoanpt=200, nsnum=inputns ):  
      self.label, self.sol, self.kcoanpt, self.lnacoanpt, self.nsnum = label, sol, kcoanpt, lnacoanpt, nsnum

      bgdata = np.genfromtxt('GR_background_data.txt')       
      et, lna, calh, bt = bgdata[:, 0], bgdata[:, 1], bgdata[:, 2], bgdata[:, 3]    
      self.et2lna = sp_interp1d(et, lna, kind=3)
      self.lna2et = sp_interp1d(lna, et, kind=3)
      self.lna2calh = sp_interp1d(lna, calh, kind=3)
      self.lna2bt = sp_interp1d(lna, bt, kind=3)
      self.lnasta, self.lnaend, self.etsta, self.etend = lna[0], lna[-1], et[0], et[-1]
            
      lnacoaset = np.genfromtxt('GR_lna_data.txt')   
      knumcoaset = np.genfromtxt('GR_k_data.txt')
      source1data = np.genfromtxt('GR_sourcedis1set_data.txt')
      source2data = np.genfromtxt('GR_sourcedis2set_data.txt')             
      
      self.source1fit = RectBivariateSpline(lnacoaset, knumcoaset, source1data, kx=3, ky=3, s=0)
      self.source2fit = RectBivariateSpline(lnacoaset, knumcoaset, source2data, kx=3, ky=3, s=0)      
      #self.lnacoaset = lnacoaset
      #self.knumcoaset = knumcoaset
      self.lnamin, self.lnamax, self.kmin, self.kmax = lnacoaset[0], lnacoaset[-1], knumcoaset[0], knumcoaset[-1]
      
   
   
   def intglna(self, knum, lnum ):      
      lnafinnpt=8000
      #print(self.lnamin, self.lnamax, self.etsta, self.etend)
      lnafinset=np.linspace(self.lnamin, self.lnamax, lnafinnpt, endpoint=True )
      #self.lnafinset=lnafinset
      source1set = self.source1fit(lnafinset, knum)[:,0]
      source2set = self.source2fit(lnafinset, knum)[:,0]
      
      etfinset = self.lna2et(lnafinset) 
          
      jn_argset = knum*( self.etend - self.lna2et(lnafinset) )
      jnset = spherical_jn(lnum, jn_argset, derivative=0)
      jnpset = spherical_jn(lnum, jn_argset, derivative=1)
      
      #print(lnafinset)
      #print( source1set, source2set )
      #print( jnset, jnpset )
      #print( source1set*jnset, source2set*jnpset )
      testyset = source1set*jnset + source2set*jnpset
            
      return trapezoid(testyset, etfinset)
      #return simpson(testyset, etfinset) 
      #calhset = self.lna2calh(lnafinset)
      #testy1set = (source1set*jnset + source2set*jnpset)/calhset  
      #print(testy1set)         
      #return trapezoid( testy1set, lnafinset )


   def intgk(self, lnum ):      
      kfinnpt=2000
      #if lnum>300:
         #kfinnpt = 1500
      knumfinset=np.linspace(self.kmin, self.kmax, kfinnpt )
      transf = np.zeros(kfinnpt)
      testyset = np.zeros(kfinnpt)
      
      for i in range(kfinnpt):
         transf[i] = self.intglna(knumfinset[i], lnum)
         testyset[i] = ( self.intglna(knumfinset[i], lnum) )**2 * (knumfinset[i])**(self.nsnum-2)

      return trapezoid(testyset, knumfinset) 



      
      
      
     

t0 = timeit.time.time()
      
switch = 1   
if __name__ == '__main__':
   if switch:
      #lset = np.arange(10, 525, 25)
      lset = np.array([10, 30, 50, 70, 100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000 ])
      #lset = np.array([10, 30, 50, 70, 100, 150, 200, 250, 300, 350, 400])
      clset = np.zeros(len(lset))

      xx = clcalc()     
      for ii in range(len(lset)):
         clset[ii] = xx.intgk(lset[ii])
         print( lset[ii], clset[ii])
      
      clintp = sp_interp1d(lset, lset*(lset+1)*clset, kind=3)      
      lset1 = np.arange(min(lset), max(lset)+1, 5)
      clset1 = clintp(lset1)
      
      ind_norm = 6  
      normnum = llp1c200/(lset[ind_norm]*(lset[ind_norm]+1)*clset[ind_norm])       
      f = open('GR_cl_data.txt', 'w+')    
      for ii in range(len(lset1)):       
         f.write(('%d %.9e' % (lset1[ii], normnum*clset1[ii] ) ) + '\n')
         #f.write(('%d %.9e' % (lset1[ii], clintpset[ii] ) ) + '\n')         
      f.close()
      plt.plot(lset1, normnum*clset1)
      plt.ylabel(r'$l(l+1)C_l$', fontsize=fontsize1)
      plt.xlabel(r'$l$', fontsize=fontsize1)
      plt.savefig('GR_cl_fig.pdf', format='pdf', bbox_inches='tight')
      
         
   else:
      xx = clcalc()
      testkset = np.linspace(0.1, 1500, 2) 
      for ii in range(len(testkset)):
         print( testkset[ii], xx.intglna(testkset[ii], 200) )
      #print( xx.intgk(200) )
   
   
   
   
print( '\n *** use %.2f seconds\n' % (timeit.time.time() - t0))


#plt.show()     
