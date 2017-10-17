# -*- coding: utf-8 -*-
"""
Created on Wed Apr 12 10:23:39 2017
Updated on Tue Oct 17 10:54 2017

@author: ajb7
"""

import os
import csv
import glob
import re
import cmath
import ntpath
import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt, pi, exp, linspace
from numpy.linalg import inv
from scipy import optimize
from scipy.optimize import curve_fit

#turn interactive plotting off
plt.ioff()

#set IPython console to print arrays in their entirety (useful for debug)
np.set_printoptions(threshold=np.inf)

#------------------------------------------------------------------------------
#function definitions
def symLor(x, A, x0, d, C0, C1):
    return A*d**2/((x-x0)**2 + d**2) + C0 + C1*x

def MagFit(field, magS21, phaseS21):
    #generate initial guesses
    A_Lor = (max(magS21) - min(magS21))
    pkPos = np.argmax(magS21)
    x0_Lor = field[pkPos] #location of maximum
    
    #find full-width, half-max
    locs = np.where(magS21 < (max(magS21) - 1/2*A_Lor)) #returns tuple of indices where condition is true
    #convert loc to 1D array
    l = np.array(locs)
    l = l.T
    l = np.reshape(l,len(l))
    FWHMmin = field[int(min(np.where(l > pkPos, l, len(l)+1)))]
    FWHMmax = field[max(np.where(l < pkPos, l, -1))]    
    d_Lor = np.absolute(FWHMmax - FWHMmin)
    
    C0_Lor = 1/2*(magS21[0] + magS21[len(magS21)-1]) #average of first and last values of magS21
    C1_Lor = 0
    Lor_init = [A_Lor, x0_Lor, d_Lor, C0_Lor, C1_Lor]

    Lor_fit, Lor_covar = curve_fit(symLor, field, magS21, p0=Lor_init, bounds=([-np.inf,0,0,-np.inf,-np.inf],np.inf))
    #bounds prevent negative x0, d
    return Lor_fit

#------------------------------------------------------------------------------
def fitAll(f_min, f_max, f_min_lin, f_max_lin, field_col, reS21_col, imS21_col, window, windowSize, fname, resultsDir):
    
#define file directory and initialize arrays    
    
    filepath = ntpath.dirname(fname) #Returns the directory component of a pathname    
        
    freq = [] #frequency list
    Amp = [] #amplitude list
    Meff = [] #Meff list
    dH = [] #linewidth list
    phi = [] #phi list
    rC0 = [] #Re(C0) (background offset)
    iC0 = [] #Im(C0)
    rC1 = [] #Re(C1) (background slope)
    iC1 = [] #Im(C1)
    Hres = [] #resonant field list
    
    Amp_err = [] #amplitude list
    Meff_err = [] #Meff list
    dH_err = [] #linewidth list
    phi_err = [] #phi list
    rC0_err = [] #Re(C0) (background offset)
    iC0_err = [] #Im(C0)
    rC1_err = [] #Re(C1) (background slope)
    iC1_err = [] #Im(C1)
    Hres_err = [] #resonant field list
    
    p0a = [] #array of initial values
    
    #figure counter
    j = 0
    
    #set frequency range for susceptiblity fitting (in GHz)
    #f_min = 15
    #f_max = 35
    
    #field windowing
    #window = False
    #window = True
    #windowSize = 1.5 #no. of linewidths to fit
    
    #Make fit result directories
    fitDir = filepath+'\\'+resultsDir
    
    #Lorentzian pre-fits
    fitLor = fitDir+'\\'+'MagS21'
    
    #S21 Fits
    fitS21 = fitDir+'\\'+'S21'
    
    if not os.path.exists(fitDir):
        os.makedirs(fitDir)
    if not os.path.exists(fitLor):
        os.makedirs(fitLor)
    if not os.path.exists(fitS21):
        os.makedirs(fitS21)
    
    #------------------------------------------------------------------------------
    #iterate through VNA spectra, fitting complex susceptibility
    for file in glob.glob( os.path.join(filepath,'*Final*')):
        #set filename
        print(file)
        
        m = re.search(r"(\d+\.\d+) GHz", ntpath.basename(file))
        if m:
            freq.append(float(m.group(1)))
            n = len(freq)
            print ('Freq. is', freq[n-1])
        
        #only import data if frequency is in the desired range (set above)    
        if freq[n-1] >= f_min and freq[n-1] <= f_max:
               
            #load data file using numpy.loadtxt
            #yellow magnet data uses columns 6,7,8 (index starts at 0)
            
            data = np.loadtxt(file, delimiter='\t', dtype='float64',usecols=[field_col, reS21_col, imS21_col])
            
            #create individual arrays for field values, and real and imaginary S21 measurements
            field = data[0:,0]
            reS21 = data[0:,1]
            imS21 = data[0:,2]
            
            #for some reason, the averaged 'data' has a row of 0s in its last element
            #here, these 0s are removed
            reS21 = reS21[field>0]
            imS21 = imS21[field>0]
            field = field[field>0]
            
            #define constants
            g_i = 2.1
            mu_B = 9.274009994e-24 #Joules/Tesla
            hbar = 6.6260693e-34/(2*pi) #Joule*seconds
            gamma = g_i*mu_B/hbar
            f = freq[n-1]*1e9 #Hz
            mu_0 = 4*pi*1e-7 #H/m
            Heff = (2*pi*f)/(gamma*mu_0)*(4*pi/1e3) #(4*pi/1e3) needed to convert Oe to A/m
            
            #correct S21 for complex offset
            rC0_pre = np.mean(reS21)
            iC0_pre = np.mean(imS21)
            rC1_pre = (reS21[len(reS21)-1] - reS21[0])/(field[len(field)-1] - field[0]) #slope based on final and initial values of reS21, field
            iC1_pre = (imS21[len(reS21)-1] - imS21[0])/(field[len(field)-1] - field[0]) #slope based on final and initial values of imS21, field
                    
            reS21_corr = reS21 - ((rC0_pre - rC1_pre*(max(field) + min(field))/2) + rC1_pre*field)
            imS21_corr = imS21 - ((iC0_pre - iC1_pre*(max(field) + min(field))/2) + iC1_pre*field)
            
            #calculate magnitude and phase from raw Re and Im data to do Lorentzian pre-fit
            S21_corr = reS21_corr + 1j*imS21_corr
        
            phaseS21_corr = np.unwrap(np.angle(S21_corr), pi)
            magS21_corr = sqrt(reS21_corr**2 + imS21_corr**2)
            
            Lor_fit = MagFit(field,magS21_corr,phaseS21_corr)
                    
            plt.figure(j)
            j += 1
            plt.subplot(211)
            plt.plot(field,magS21_corr,'o') 
            #plot initial guess fit
            #plt.plot(field,symLor(field, Lor_init[0], Lor_init[1], Lor_init[2], Lor_init[3], Lor_init[4]))
            #plot final fit results
            plt.plot(field,symLor(field, Lor_fit[0], Lor_fit[1], Lor_fit[2], Lor_fit[3], Lor_fit[4]))
            plt.title('|S21| @ '+str(freq[-1])+'_GHz')
            plt.xlabel('Field (Oe)')
            plt.ylabel('|S21|')
            plt.subplot(212)
            plt.plot(field,phaseS21_corr,'go')
            plt.title('arg(S21) @'+str(freq[-1])+'_GHz')
            plt.xlabel('Field (Oe)')
            plt.ylabel('phase (rad)')
            plt.tight_layout()
            figName = 'MagFit_'+str(freq[-1])+'_GHz.png'
            figFile = os.path.join(fitLor, figName)
            plt.savefig(figFile)
            plt.close()
            
            #Use results from Lorentzian pre-fit to determine initial guesses for full susceptibility fit
            A_i = 2*Lor_fit[0]
            Meff_i = Lor_fit[1]-(2*pi*f)/(gamma*mu_0)*(4*pi/1e3) #Meff = Hres -  Heff. Conversion factor (4*pi/1e3) to convert A/m --> Oe
            dH_i = 2*Lor_fit[2]
            phi_i = (phaseS21_corr[np.argmax(magS21_corr)] + pi/2)
                       
            #rescale A_i, since chi_yy = Meff/dH on resonance
            A_i = A_i*dH_i/Meff_i            
            
            #set phi_i to a value between 0 and 2*pi
            while phi_i < 0:
                phi_i += 2*pi
            while phi_i > 2*pi:
                phi_i -= 2*pi
            
            p0 = [A_i, Meff_i, dH_i, phi_i, rC0_pre, iC0_pre, 0, 0]
            
            #array of initial values
            p0a.append(p0)            
        
            #lower bounds
            A_lb = -np.inf
            Meff_lb = -np.inf
            dH_lb = 0
            phi_lb = 0
            rC0_lb = -np.inf
            iC0_lb = -np.inf
            rC1_lb = -np.inf
            iC1_lb = -np.inf
            
            #upper bounds
            A_ub = np.inf
            Meff_ub = np.inf
            dH_ub = np.inf
            phi_ub = 2*pi
            rC0_ub = np.inf
            iC0_ub = np.inf
            rC1_ub = np.inf
            iC1_ub = np.inf
            bnds = ([A_lb, Meff_lb, dH_lb, phi_lb, rC0_lb, iC0_lb, rC1_lb, iC1_lb], [A_ub, Meff_ub, dH_ub, phi_ub, rC0_ub, iC0_ub, rC1_ub, iC1_ub])
            
            if A_i < A_lb or A_i > A_ub:
                print('A_i out of bounds')
            if Meff_i < Meff_lb or Meff_i > Meff_ub:
                print ('Meff_i out of bounds')
            if dH_i < dH_lb or dH_i > dH_ub:
                print ('dH_i out of bounds')
            if phi_i < phi_lb or phi_i > phi_ub:
                print ('phi_i out of bounds')
                
            chiyy = lambda Meff, dH, x: Meff*(x - Meff)/((x-Meff)**2 - Heff**2 + 1j*dH*Heff)  
            
            SusFit = lambda p, x: (p[4] + 1j*p[5]) + (p[6] + 1j*p[7])*x + p[0]*cmath.exp(1j*p[3])*chiyy(p[1], p[2], x)
            #where
            #A = p[0]
            #Meff = p[1]
            #dH = p[2]
            #phi = p[3]
            #rC0 = p[4]
            #iC0 = p[5]
            #rC1 = p[6]
            #iC1 = p[7]
            
            #define error function for leastsq fit (minimizes difference between data and target function)
            errSusFit = lambda p, x, y1, y2: np.r_[np.real(SusFit(p,x)) - y1 , np.imag(SusFit(p,x)) - y2]
            
            #add fitting constraints (A > 0)
            #perform the fit using leastsq
            #fit_result, cov_result, infodict, mesg, ier = optimize.leastsq(errSusFit, p0, args=(field[:], reS21[:], imS21[:]), full_output = True, maxfev = 500, xtol = 0.0000000001, ftol = 0.000000000001)
            #without bounds
            #fit_result = optimize.least_squares(errSusFit, p0, args=(field[:], reS21[:], imS21[:]), max_nfev = 500, xtol = 0.0000000001, ftol = 0.000000000001)
            #with bounds
            fit_result = optimize.least_squares(errSusFit, p0, bounds = bnds, args=(field[:], reS21[:], imS21[:]), max_nfev = 1000, xtol = 0.0000000001, ftol = 0.000000000001)
            
            if window:
                #in v10 of the code, I moved the windowing to use SusFit results instead of Lorentzian pre-fit results
                H_r = 2*pi*f/(gamma*mu_0)*(4*pi/1e3) + fit_result.x[1]
                deltaH = fit_result.x[2]
                print("dH (Oe) = "+str(deltaH))
                #if fitting only a smaller field range, crop field, reS21, and imS21 arrays
                keepIndx = (field < (H_r + windowSize*deltaH)) & (field > (H_r - windowSize*deltaH))
                field  = field[keepIndx]
                print("dH (# of linewidths) = "+str((field[0] - field[len(field)-1])/deltaH))
                reS21 = reS21[keepIndx]
                imS21 = imS21[keepIndx]
                
                #redo fit with windowed data
                fit_result = optimize.least_squares(errSusFit, p0, bounds = bnds, args=(field[:], reS21[:], imS21[:]), max_nfev = 1000, xtol = 0.0000000001, ftol = 0.000000000001)
            
            #fill result arrays with fit_result
            #calculate Hres = 2*pi*f/(gamma*mu0) + Meff
            H_res = 2*pi*f/(gamma*mu_0)*(4*pi/1e3) + fit_result.x[1]
            
            Amp.append(fit_result.x[0])
            Meff.append(fit_result.x[1])
            dH.append(fit_result.x[2])
            phi.append(fit_result.x[3])
            rC0.append(fit_result.x[4])
            iC0.append(fit_result.x[5])
            rC1.append(fit_result.x[6])
            iC1.append(fit_result.x[7])
            Hres.append(H_res)
            
            #calculate fit errors
            #hessian = J^T*J
            hessian = np.dot(np.transpose(fit_result.jac),fit_result.jac)
            h_inv = inv(hessian)
            #variance of residual
            sigma_res = np.var(fit_result.fun)
            #covariance matrix
            cov = h_inv*sigma_res
            #standard deviation is square-root of variance
            std_err = np.sqrt(np.absolute(cov))
            
            Amp_err.append(np.diag(std_err)[0])
            Meff_err.append(np.diag(std_err)[1])
            dH_err.append(np.diag(std_err)[2])
            phi_err.append(np.diag(std_err)[3])
            rC0_err.append(np.diag(std_err)[4])
            iC0_err.append(np.diag(std_err)[5])
            rC1_err.append(np.diag(std_err)[6])
            iC1_err.append(np.diag(std_err)[7])
            
            #error in Hres is same as error in Meff
            Hres_err.append(np.diag(std_err)[1])
            
            #final result fit
            plt.figure(j)
            j += 1
            plt.subplot(211)        
            plt.plot(field,reS21,'bo',fillstyle='none', label='data')
            plt.plot(field,np.real(SusFit(p0,field)),'k--', label='initial guess')
            plt.plot(field,np.real(SusFit(fit_result.x,field)),'m', label='final fit')
            plt.title('Re(S_21) @ f=' + str(freq[n-1]) + ' GHz')
            plt.legend()
            plt.xlabel('Field (Oe)')
            plt.ylabel('Re(S_21)')
            plt.subplot(212)
            plt.plot(field,imS21,'ro',fillstyle='none', label='data')
            plt.plot(field,np.imag(SusFit(p0,field)),'k--', label='initial guess')
            plt.plot(field,np.imag(SusFit(fit_result.x,field)),'c', label='final fit')
            plt.title('Im(S_21) @ f=' + str(freq[n-1]) + ' GHz')
            plt.legend()
            plt.xlabel('Field (Oe)')
            plt.ylabel('Im(S_21)')
            plt.tight_layout()
            figName = 'S21Fit_'+str(freq[-1])+'_GHz.png'
            figFile = os.path.join(fitS21, figName)
            plt.savefig(figFile)
            plt.close()
            
        else:
            del freq[n-1] #don't include frequencies outside the specified range in the final list of frequencies
    
    #------------------------------------------------------------------------------    
    #after completion of FOR loop and VNA spectra fits, plot data as a function of frequency,
    #and perform various fits (FMR dispersion, damping, inductance, ...)
    
    #define a finely spaced frequency array for plotting purposes
    freqscale = linspace(0,max(freq))
    
    #sort result arrays in order of increasing frequency
    #               0    1     2     3       4     5   6     7     8     9     10    11    12    13  14      15   16      17    18     19
    out = list(zip(freq,Amp,Amp_err,Meff,Meff_err,dH,dH_err,phi,phi_err,rC0,rC0_err,iC0,iC0_err,rC1,rC1_err,iC1,iC1_err,Hres,Hres_err,p0a))
    out.sort()
    
    #reassign individual result arrays from sorted "out" list of tuples
    freq = [i[0] for i in out]
    Amp = [i[1] for i in out]
    Amp_err = [i[2] for i in out]
    Meff = [i[3] for i in out]
    Meff_err = [i[4] for i in out]
    dH = [i[5] for i in out]
    dH_err = [i[6] for i in out]
    phi = [i[7] for i in out]
    phi_err = [i[8] for i in out]
    rC0 = [i[9] for i in out]
    rC0_err = [i[10] for i in out]
    iC0 = [i[11] for i in out]
    iC0_err = [i[12] for i in out]
    rC1 = [i[13] for i in out]
    rC1_err = [i[14] for i in out]
    iC1 = [i[15] for i in out]
    iC1_err = [i[16] for i in out]
    Hres = [i[17] for i in out]
    Hres_err = [i[18] for i in out]
    p0a = [i[19] for i in out]
    
    #convert lists to arrays
    freq = np.array(freq)
    Amp = np.array(Amp)
    Amp_err = np.array(Amp_err)
    Meff = np.array(Meff)
    Meff_err = np.array(Meff_err)
    dH = np.array(dH)
    dH_err = np.array(dH_err)
    phi = np.array(phi)
    phi_err = np.array(phi_err)
    rC0 = np.array(rC0)
    rC0_err = np.array(rC0_err)
    iC0 = np.array(iC0)
    iC0_err = np.array(iC0_err)
    rC1 = np.array(rC1)
    rC1_err = np.array(rC1_err)
    iC1 = np.array(iC1)
    iC1_err = np.array(iC1_err)
    Hres = np.array(Hres)
    Hres_err = np.array(Hres_err)
    p0a = np.array(p0a)
    
    #calculate insertion loss = 20*log_10(sqrt(Re(Z_BG)^2 + Im(Z_BG)^2))
    IL = 20*np.log10(sqrt((rC0 + rC1*Hres)**2 + (iC0 + iC1*Hres)**2))
    IL_err = 20/np.log(10)*sqrt(((iC0*iC1 + rC0*rC1 + Hres*(iC1**2 + rC1**2))**2 * Hres_err**2 + 
                                (iC0 + Hres*iC1)**2 * iC0_err**2 +
                                Hres**2 * (iC0 + Hres*iC1)**2 * iC1_err**2 + 
                                (rC0 + Hres*rC1)**2 * rC0_err**2 +
                                Hres**2 * (rC0 + Hres*rC1)**2 * rC1_err**2) / ((iC0 + Hres*iC1)**2 + (rC0 + Hres*rC1)**2)**2)
    #save fit results
    txtName = 'Susceptibility Fit Results.csv'
    txtFile = os.path.join(fitDir, txtName)
    
    with open(txtFile,'w',newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["Freq (GHz)","H_res (Oe)", "err H_res (Oe)", "dH (Oe)", "err dH (Oe)", "A", "err A", "phi (rad)", "err phi (rad)",
                         "Re(C_0)","err Re(C_O)","Im(C_0)","err Im(C_0)","Re(C_1)","err Re(C_1)","Im(C_1)", "err Im(C_1)","IL (dB)","err IL (dB)"])
        writer.writerows(zip(freq,Hres,Hres_err,dH,dH_err,Amp,Amp_err,phi,phi_err,rC0,rC0_err,iC0,iC0_err,rC1,rC1_err,iC1,iC1_err,IL,IL_err))
    
    #--------------------------------------------------------------------------
    #linear range for Kittel, damping, and inductance fits
    linRegime_min = f_min_lin #GHz
    linRegime_max = f_max_lin
    mask = np.logical_and(freq >= linRegime_min, freq <= linRegime_max)
    #--------------------------------------------------------------------------
    #Kittel fit    
      
    def linearFit(x, m, b):
        return m*x + b
    
    def KittelFit(freq, Hres, Hres_err):
        #initial guesses
        g_i = 2.1
        Meff_i = 10000 #1 Tesla
        
        Kittel_init = [2*pi/((g_i*mu_B/hbar)*mu_0)*(4*pi/1e3), Meff_i]    
    
        result = curve_fit(linearFit, freq, Hres, p0=Kittel_init, sigma=Hres_err)
        #Kittel_fit[0] = (1/linear_result[0])*2*pi*hbar/(mu_0*mu_B)*(4*pi/1e3)
        #Kittel_fit[1] = linear_result[1]
        return result
    
    freq_lin = freq[mask]
    Hres_lin = Hres[mask]
    Hres_err_lin = Hres_err[mask]
    
    Kittel_fit = KittelFit(freq_lin*1e9, Hres_lin, Hres_err_lin)
    g_final = 2*pi*hbar/(mu_0*mu_B*Kittel_fit[0][0])*(4*pi/1e3)
    Meff_final = Kittel_fit[0][1]
    Kittel_err = np.sqrt(np.diag(Kittel_fit[1]))
    g_final_err = 2*pi*hbar/(mu_0*mu_B*Kittel_fit[0][0]**2)*(4*pi/1e3)*Kittel_err[0]
    Meff_final_err = Kittel_err[1]
    print("Kittel fit results:")
    print("g =",g_final,"+/-",g_final_err)
    print("M_eff =",Meff_final,"+/-",Meff_final_err,"Oe")
    
    plt.figure(j)
    j += 1
    plt.errorbar(freq,Hres,yerr=Hres_err,fmt='bo',fillstyle='none',capsize=4)
    plt.plot(freqscale,linearFit(freqscale,Kittel_fit[0][0]*1e9,Kittel_fit[0][1]),'k-')
    plt.title("Kittel fit")
    plt.xlabel("Freq (GHz)")
    plt.ylabel("H_res (Oe)")
    figName = "Kittel Fit.png"
    figFile = os.path.join(fitDir, figName)
    plt.savefig(figFile)
    plt.close()
        
    #------------------------------------------------------------------------------
    #Linewidth fit
    def dampingFit(freq, dH, dH_err):
        #initial guesses
        alpha_i = .01
        dH_0_i = 0
            
        damping_init = [4*pi*alpha_i/(gamma*mu_0)*(4*pi/1e3), dH_0_i]
        
        result = curve_fit(linearFit, freq, dH, p0=damping_init, sigma=dH_err)
        return result
    
    dH_lin = dH[mask]
    dH_err_lin = dH_err[mask]
    
    damping_fit = dampingFit(freq_lin*1e9, dH_lin, dH_err_lin)
    alpha_final = damping_fit[0][0]*gamma*mu_0/(4*pi)*(1e3/(4*pi))
    dH_0_final = damping_fit[0][1]
    damping_err = np.sqrt(np.diag(damping_fit[1]))
    alpha_final_err = damping_err[0]*gamma*mu_0/(4*pi)*(1e3/(4*pi))
    dH_0_final_err = damping_err[1]
    print("Damping fit results:")
    print("alpha =",alpha_final,"+/-",alpha_final_err)
    print("dH_0 =",dH_0_final,"+/-",dH_0_final_err,"Oe")
    
    plt.figure(j)
    j += 1
    plt.errorbar(freq,dH,yerr=dH_err,fmt='go',fillstyle='none',capsize=4)
    plt.plot(freqscale,linearFit(freqscale,damping_fit[0][0]*1e9,damping_fit[0][1]),'k-')
    plt.title("Damping fit")
    plt.xlabel("Freq (GHz)")
    plt.ylabel("dH (Oe)")
    figName = "Damping Fit.png"
    figFile = os.path.join(fitDir, figName)
    plt.savefig(figFile)
    plt.close()
    
    #------------------------------------------------------------------------------
    #save spectroscopy data
    txtName = 'Spectroscopy Fit Results.csv'
    txtFile = os.path.join(fitDir, txtName)
    
    with open(txtFile,'w',newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["M_eff (Oe)", "err M_eff (Oe)", "g", "err g", "dH0 (Oe)","err dH0 (Oe)","alpha", "err alpha"])
        writer.writerow([Meff_final,Meff_final_err,g_final,g_final_err,dH_0_final,dH_0_final_err,alpha_final,alpha_final_err])
    
    #------------------------------------------------------------------------------
    #Complex Inductance Analysis
    def ReLvsf_Fit(p, x): 
        ReL0 = p[0]
        dReL_df = p[1]
        dImL_df = p[2]
        phi_anom = p[3]
        return (ReL0 + dReL_df*x)*np.cos(phi_anom) - dImL_df*x*np.sin(phi_anom)
    
    def ImLvsf_Fit(p, x):
        ReL0 = p[0]
        dReL_df = p[1]
        dImL_df = p[2]
        phi_anom = p[3]
        return dImL_df*x*np.cos(phi_anom) + (ReL0 + dReL_df*x)*np.sin(phi_anom)
    
    def errLvsf(p, f, Lreal, Limag, err_Lreal, err_Limag):
        #concatenate real and imaginary errors into 1D array
        return np.r_[(ReLvsf_Fit(p,f) - Lreal)/err_Lreal, (ImLvsf_Fit(p,f) - Limag)/err_Limag]
    
    #######
    #drop ZBG, Amp, and phi points for frequencies in the non-linear regime
    
    Amp = Amp[mask]
    Amp_err = Amp_err[mask]
    phi = phi[mask]
    phi_err = phi_err[mask]
    
    rC0 = rC0[mask]
    rC0_err = rC0_err[mask]
    iC0 = iC0[mask]
    iC0_err = iC0_err[mask]
    rC1 = rC1[mask]
    rC1_err = rC1_err[mask]
    iC1 = iC1[mask]
    iC1_err = iC1_err[mask]
    
    Hres = Hres[mask]
    Hres_err = Hres_err[mask]
    
    freq = freq[mask]
    #######
    
    ZBG = (rC0 + rC1*Hres) + 1j*(iC0 + iC1*Hres)
    rZBG = np.real(ZBG)
    iZBG = np.imag(ZBG)
    ABG = np.absolute(ZBG)
    phi_BG = np.angle(ZBG)
     
    ABG_err = np.sqrt((rZBG)**2*rC0_err**2/(ABG**2) + \
                      (iZBG)**2*iC0_err**2/(ABG**2) + \
                      (Hres**2*(rZBG)**2)*rC1_err**2/(ABG**2) + \
                      (Hres**2*(iZBG)**2)*iC1_err**2/(ABG**2) + \
                      (2*rC1*(rZBG) + 2*iC1*(iZBG))**2*Hres_err**2/(4*ABG**2))
    
    phi_BG_err = np.sqrt( ((iC1*rC0 - iC0*rC1)**2*Hres_err**2 + \
                           (rZBG**2*iC0_err**2 + iC0**2*rC0_err**2) + \
                           Hres*(Hres*(rZBG)**2*iC1_err**2 + iC1*(2*iC0 + iC1*Hres)*rC0_err**2 + Hres*(iZBG)**2*rC1_err**2)) / \
                            (iZBG**2 + rC0**2 + 2*Hres*rC0*rC1 + Hres**2*rC1**2)**2 )
    
    
    Z = Amp*np.exp(1j*phi)
    
    #need to add pi/2 phase to convert Z to L
    L = 2*50/(2*pi*freq*1e9)*Z/ZBG*exp(1j*pi/2)
    rL = np.real(L)
    iL = np.imag(L)
    
    #calculate uncertainties in Re(L) and Im(L) vs. f
    #see "Inductance Uncertainties.nb"
    rL_err = 50/(2*pi*freq*1e9)*np.sqrt( (np.cos(phi-phi_BG)/ABG)**2*Amp_err**2 + \
                      (Amp*np.cos(phi-phi_BG)/ABG**2)**2*ABG_err**2 + \
                      (Amp*np.sin(phi-phi_BG)/ABG)**2*phi_err**2 + \
                      (Amp*np.sin(phi-phi_BG)/ABG)**2*phi_BG_err**2 )
    iL_err = 50/(2*pi*freq*1e9)*np.sqrt( (np.sin(phi-phi_BG)/ABG)**2*Amp_err**2 + \
                      (Amp*np.sin(phi-phi_BG)/ABG**2)**2*ABG_err**2 + \
                      (Amp*np.cos(phi-phi_BG)/ABG)**2*phi_err**2 + \
                      (Amp*np.cos(phi-phi_BG)/ABG)**2*phi_BG_err**2 )
    #L_err = np.sqrt(rL_err**2 + iL_err**2)
    
    #do quick linear fit to uncorrected inductances in order to generate initial guesses for L_corr fit
    ReL_fit, ReL_covar = curve_fit(linearFit,freq,np.real(L),sigma=rL_err)
    ImL_fit, ImL_covar = curve_fit(linearFit,freq,np.imag(L),sigma=iL_err) 
    
    dReL_df_i = ReL_fit[0]
    dImL_df_i = ImL_fit[0]
    ReL0_i = ReL_fit[1] #intercept of Re(L) vs. f
    ImL0_i = ImL_fit[1] #intercept of Im(L) vs. f
    
    phi_anom_i = np.arctan2(ImL0_i,ReL0_i)    
  
    #------------------------------------------------------------------------------
    #simultaneous least_squares fit of Re(L) and Im(L) in order to correct for anomalous phase
    Lcorr_init = [ReL0_i, dReL_df_i, dImL_df_i, phi_anom_i]
    
    Lcorr_result = optimize.least_squares(errLvsf, Lcorr_init, args=(freq[:], rL[:], iL[:], rL_err[:], iL_err[:]), method='lm')
    
    ReL0_final = Lcorr_result.x[0]
    dReL_df_final = Lcorr_result.x[1]
    dImL_df_final = Lcorr_result.x[2]
    phi_anom_final = Lcorr_result.x[3]
    
    #calculate fit errors
    #hessian = J^T*J
    Hess_L = np.dot(np.transpose(Lcorr_result.jac),Lcorr_result.jac)
    Hess_inv_L = inv(Hess_L)
    #variance of residual
    sigma_res_L = np.var(Lcorr_result.fun)
    #covariance matrix
    cov_L = Hess_inv_L*sigma_res_L
    #standard deviation is square-root of variance
    stderr_L = np.sqrt(np.absolute(cov_L))
    
    ReL0_final_err = np.diag(stderr_L)[0]
    dReL_df_final_err = np.diag(stderr_L)[1]
    dImL_df_final_err = np.diag(stderr_L)[2]
    phi_anom_final_err = np.diag(stderr_L)[3]
    
    print("Simultaneous inductance fit results:")
    print("L_0 =",ReL0_final,"+/-",ReL0_final_err," H")
    print("d(Re(L))/df =",dReL_df_final,"+/-",dReL_df_final_err," H/GHz")
    print("d(Im(L))/df =",dImL_df_final,"+/-",dImL_df_final_err," H/GHz")
    print("phi_anom =",phi_anom_final,"+/-",phi_anom_final_err," rad")
    
    #rotate inductance by anomalous phase
    Lcorr = L*np.exp(-1j*phi_anom_final)
    rLcorr = np.real(Lcorr)
    iLcorr = np.imag(Lcorr)
    
    rLcorr_err = np.sqrt( np.absolute((np.cos(phi_anom_final))**2*rL_err**2 + (np.sin(phi_anom_final))**2*iL_err**2 + \
                           (iL*np.cos(phi_anom_final) - rL*np.sin(phi_anom_final))**2*phi_anom_final_err**2) )
    iLcorr_err = np.sqrt( np.absolute((np.cos(phi_anom_final))**2*iL_err**2 - (np.sin(phi_anom_final))**2*rL_err**2 + \
                           (-iL*np.sin(phi_anom_final) - rL*np.cos(phi_anom_final))**2*phi_anom_final_err**2) )
    
    #------------------------------------------------------------------------------
    #save Inductance results (from Simultaneous fit)
    txtName = 'Inductance Results.csv'
    txtFile = os.path.join(fitDir, txtName)
    
    with open(txtFile,'w',newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["Freq (GHz)","Re(L)","err Re(L)","Im(L)","err Im(L)","Re(L_corr)","err Re(L_corr)","Im(L_corr)","err Im(L_corr)"])
        writer.writerows(zip(freq,rL,rL_err,iL,iL_err,rLcorr,rLcorr_err,iLcorr,iLcorr_err))
    
    #plot corrected Re(L) and Im(L) vs. f
    plt.figure(j)
    j += 1
    #uncorrected Re(L)
    plt.errorbar(freq,np.real(L),yerr=rL_err,fmt='bs',fillstyle='none',capsize=4,label='Re(L)')
    plt.plot(freqscale,linearFit(freqscale,ReL_fit[0],ReL_fit[1]),'b--')
    #corrected Re(L), from seperate fits
    #plt.errorbar(freq,rLcorr_i ,yerr=rLcorr_i_err,fmt='kv',fillstyle='none',capsize=4,label='Re(L_corr)_sep')
    #plt.plot(freqscale,linearFit(freqscale,ReLcorr_fit[0],ReLcorr_fit[1]),'k--')
    #corrected Re(L), from simultaneous fit
    plt.errorbar(freq,rLcorr,yerr=rLcorr_err,fmt='bo',capsize=4,label='Re(L_corr)')
    plt.plot(freqscale,linearFit(freqscale,dReL_df_final,ReL0_final),'b')
    plt.legend()
    plt.title('Re(L)')
    plt.ylabel('Re(L) (H)')
    plt.xlabel('Freq (GHz)')
    figName = "Re(L).png"
    figFile = os.path.join(fitDir, figName)
    plt.savefig(figFile)
    plt.close()
    
    plt.figure(j)
    j += 1
    #uncorrected Im(L)
    plt.errorbar(freq,np.imag(L),yerr=iL_err,fmt='rs',fillstyle='none',capsize=4,label='Im(L)')
    plt.plot(freqscale,linearFit(freqscale,ImL_fit[0],ImL_fit[1]),'r--')
    #corrected Im(L), from seperate fits
    #plt.errorbar(freq,iLcorr_i ,yerr=iLcorr_i_err,fmt='kv',fillstyle='none',capsize=4,label='Im(L_corr)_sep')
    #plt.plot(freqscale,linearFit(freqscale,ImLcorr_fit[0],ImLcorr_fit[1]),'k--')
    #corrected Im(L), from simultaneous fit
    plt.errorbar(freq,iLcorr,yerr=iLcorr_err,fmt='ro',capsize=4,label='Im(L_corr)')
    plt.plot(freqscale,linearFit(freqscale,dImL_df_final,0),'r')
    plt.legend()
    plt.title('Im(L)')
    plt.ylabel('Im(L) (H)')
    plt.xlabel('Freq (GHz)')
    figName = "Im(L).png"
    figFile = os.path.join(fitDir, figName)
    plt.savefig(figFile)
    plt.close()
    
    #save Inductance vs. f linear fit results
    txtName = 'L_corr vs f Fit Results.csv'
    txtFile = os.path.join(fitDir, txtName)
    
    with open(txtFile,'w',newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["L0 (H)","err L0 (H)","dRe(L)/df (H/GHz)","err dRe(L)/df (H/GHz)","dIm(L)/df (H/GHz)","err dIm(L)/df (H/GHz)","phi_anom (rad)","err phi_anom"])
        writer.writerow([ReL0_final,ReL0_final_err,dReL_df_final,dReL_df_final_err,dImL_df_final,dImL_df_final_err,phi_anom_final,phi_anom_final_err])
        
    #COMPLETE results
    txtName = 'All Results.csv'
    txtFile = os.path.join(fitDir, txtName)
    
    with open(txtFile,'w',newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["M_eff (Oe)", "err M_eff (Oe)", "g", "err g", "dH0 (Oe)","err dH0 (Oe)","alpha", "err alpha","L0 (H)","err L0 (H)","dRe(L)/df (H/GHz)","err dRe(L)/df (H/GHz)","dIm(L)/df (H/GHz)","err dIm(L)/df (H/GHz)","phi_anom (rad)","err phi_anom"])
        writer.writerow([Meff_final,Meff_final_err,g_final,g_final_err,dH_0_final,dH_0_final_err,alpha_final,alpha_final_err,ReL0_final,ReL0_final_err,dReL_df_final,dReL_df_final_err,dImL_df_final,dImL_df_final_err,phi_anom_final,phi_anom_final_err])
    
