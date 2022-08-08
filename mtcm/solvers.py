from . import functions

import pandas as pd
import numpy as np

from math import pi as pi
from math import sqrt as sqrt

class mtcm():
    
    def __init__(self,
        Lc: float,
        As_tot: float,
        n_phi_s: float,
        phi_s: float,
        rho_s: float,
        Es: float,
        Ecm: float,
        fctm: float,
        fs_yield: float=500,
        fs_ult: float=640,
        eps_ult: float=100e-03,
        zeta: float=1.0,
        psi: float=0.70,
        u1: float=0.1,
        tau_max: float=5.0,
        alpha: float=0.35,
        beta_sm: float=1.0
    ) -> None:
        """Method for instantiating RC tie
        
        Args:    
            eps_sr: steel strain at the crack, e.g. sigma_sr/Es = 2.0*1e-3
            Lc: Member length [mm]
            As_tot: Total steel rebar area [mm2]
            n_phi_s: Number of reinforcing steel bars 
            phi_s: Steel rebar diameter [mm] 
            rho_s: Reinforcement ratio, e.g. As_tot/Ac_ef
            Es: Youngs Modulus for steel [MPa]
            Ecm: Youngs modulus for concrete [MPa]
            fctm: Tensile strength [MPa]
        
        Kwargs:
            zeta: MTCM parameter related to uniform bond stresses, default is 1.0 
            psi: MTCM parameter related to ratio between strains at rebar level 
                and the mean strains over the cover, default is 0.7
            u1: MTCM bond slip parameter [mm], default is 0.1
            tau_max: MTCM bond slip parameter [MPa], default is 5.0
            alpha: MTCM bond slip parameter, default is 0.35
        """        
        
        # RC tie attributes
        self.Lc = Lc
        self.As_tot = As_tot
        self.n_phi_s = n_phi_s
        self.phi_s = phi_s
        self.rho_s = rho_s
        self.Es = Es
        self.Ecm = Ecm
        self.fctm = fctm
        
        # Steel nonlinearity
        self.fs_yield = fs_yield
        self.fs_ult = fs_ult
        self.eps_ult = eps_ult
        
        # Bond slip parameters
        self.zeta = zeta
        self.psi = psi
        self.u1 = u1
        self.tau_max = tau_max
        self.alpha = alpha
        self.beta_sm = beta_sm
        
    def stress(self,
        eps_sr: float,
    ):
        """Method for calculating cracks based on steel stresses in cracks in a member
        
        Args:    
            eps_sr: Steel strain at crack, i.e. sigma_sr/Es
        """

        L_calc = self.Lc
        alpha_E = self.Es/self.Ecm                                                            # Modular ratio
        eps_ctm = self.fctm/self.Ecm                                                          # Cracking strain

        beta = 1+self.alpha                                                              # Eq. (36)
        delta = (1-self.alpha)/2                                                         # Eq. (41)
        xi = alpha_E*self.rho_s/self.psi                                                   # Eq. (27)    
        chi = (1+xi)*(self.zeta*self.n_phi_s*pi*self.phi_s/(self.As_tot*self.Es))                                  # Eq. (31)
        gamma = chi*self.tau_max/(beta*self.u1**self.alpha)                                        # Eq. (37)
        xcr0 = (1/delta)*((1/self.psi)*((1+xi)/xi)*eps_ctm*(1/(2*gamma))
            **(1/(2*delta)))**(2*delta/beta)                                     # Eq. (75) Crack spacing [mm]
        eps_sr_cr = (1/self.psi)*(1 + (1/xi))*eps_ctm                                    # Eq. (74) Cracking strain at the end of transfer length  

        eps_sr_S = (2*gamma)**(1/(2*delta))*((L_calc/2)*delta)**(beta/(2*delta))        # Eq. (72) Limit steel strain at symmetry section

        if eps_sr <= 0:
            condition = 'Compression or zero'
            concept = 'None'
            Lt = 0
            eps_sm = eps_sr
            eps_cm = eps_sr
            wcr = (eps_sm-eps_cm)*Lt
            
        else:
            if eps_sr_S > eps_sr_cr:
                condition = 'Condition 1'
                if eps_sr < eps_sr_cr:
                    concept = 'CLLM'
                    (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CLLM(eps_sr,delta,gamma,beta,xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                elif eps_sr >= eps_sr_cr:
                    condition = 'Condition 1 and Condition 2 for new cracked member'
                    L_calc = xcr0
                    eps_sr_S = (2*gamma)**(1/(2*delta))*((L_calc/2)*delta)**(beta/(2*delta)) 
                    if eps_sr < eps_sr_S:
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CLLM(eps_sr,delta,gamma,beta,xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                    elif eps_sr >= eps_sr_S:
                        concept = 'CHLM'
                        # (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CLLM(eps_sr,delta,gamma,beta,xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                        while eps_sr:
                            (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L_calc,delta,gamma,beta,xi,eps_sr_cr,self.psi,self.tau_max,self.u1,self.alpha,xcr0)
                            if eps_cm_cover_max >= eps_ctm:
                                # print('MEMBER CRACKED i.e. NEW MEMBER LENGTH L_calc = L_calc/2 CHOSEN')
                                L_calc = L_calc/2
                            elif eps_cm_cover_max < eps_ctm:
                                break
            elif eps_sr_S <= eps_sr_cr:
                condition = 'Condition 2'
                if eps_sr < eps_sr_S:
                    concept = 'CLLM'
                    (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CLLM(eps_sr,delta,gamma,beta,xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                elif eps_sr >= eps_sr_S:
                    concept = 'CHLM'
                    # (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L_calc,delta,gamma,beta,xi,eps_sr_cr,psi,tau_max,u1,alpha,xcr0)
                    while eps_sr:
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L_calc,delta,gamma,beta,xi,eps_sr_cr,self.psi,self.tau_max,self.u1,self.alpha,xcr0)
                        if eps_cm_cover_max >= eps_ctm:
                            # print('MEMBER CRACKED i.e. NEW MEMBER LENGTH L_calc = L_calc/2 CHOSEN')
                            L_calc = L_calc/2
                        elif eps_cm_cover_max < eps_ctm:
                            break
            
            # Output for plotting           
            self.plot_dict = {
                'xcoord':xcoord[0],
                'u':u[0],
                'tau':tau[0],
                'eps_s':eps_s[0]*1000,
                'eps_c':eps_c[0]*1000,
                'eps_sm':eps_sm_list[0]*1000,
                'eps_cm':eps_cm_list[0]*1000
            }
            self.df = pd.DataFrame.from_dict(self.plot_dict)
            
        # Output for accessing class attributes
        self.condition = condition
        self.concept = concept
        self.eps_sm = eps_sm
        self.eps_cm = eps_cm
        self.xcr0 = xcr0
        self.Lt = Lt
        self.wcr = wcr
        
    def strain(self,
            eps_m: float,
        ):
        """Method for calculating cracks based on imposed mean strains in a member
        
        Args:
            eps_m: mean imposed strain for the member, e.g. from shrinkage or temperature loads
        """
        
        L_calc = self.Lc
        alpha_E = self.Es/self.Ecm                                                            # Modular ratio
        eps_ctm = self.fctm/self.Ecm                                                          # Cracking strain

        beta = 1+self.alpha                                                              # Eq. (36)
        delta = (1-self.alpha)/2                                                         # Eq. (41)
        xi = alpha_E*self.rho_s/self.psi                                                      # Eq. (27)    
        chi = (1+xi)*(self.zeta*self.n_phi_s*pi*self.phi_s/(self.As_tot*self.Es))                                  # Eq. (31)
        gamma = chi*self.tau_max/(beta*self.u1**self.alpha)                                        # Eq. (37)
        xcr0 = (1/delta)*((1/self.psi)*((1+xi)/xi)*eps_ctm*(1/(2*gamma))
            **(1/(2*delta)))**(2*delta/beta)                                     # Eq. (75) Crack spacing [mm]
        eps_sr_cr = (1/self.psi)*(1 + (1/xi))*eps_ctm                                    # Eq. (74) Cracking strain at the end of transfer length  
        
        eps_sr_S = (2*gamma)**(1/(2*delta))*((L_calc/2)*delta)**(beta/(2*delta))        # Eq. (72) Limit steel strain at symmetry section
        
        beta_sm = 1.0
        
        self.stress(self.fs_yield/self.Es)
        L_yield = self.Lt
        x_yield = self.plot_dict['xcoord']
        tau_yield = self.plot_dict['tau']
        tau_m_yield = abs(np.trapz(tau_yield,x_yield)/max(x_yield))
        tau_b0 = tau_m_yield
        tau_b1 = tau_b0/2
        Esh = (self.fs_ult-self.fs_yield)/(self.eps_ult-self.fs_yield/self.Es)
        
        if eps_m <= 0:
            condition = 'Compression or zero'
            concept = 'None'
            (sigma_sr, sigma_sm, eps_sr) = functions.nakedsteel(eps_m,self.Es,self.fs_yield,self.fs_ult,self.eps_ult)
            sigma_cm = 0
            eps_sm = eps_m
            eps_cm = eps_m
            Lt = 0
            wcr = (eps_sm-eps_cm)*Lt
            
        elif (eps_m > self.fs_yield/self.Es-tau_b0*L_yield/(self.phi_s*self.Es) and 
            eps_m <= (self.fs_yield/self.Es+tau_b1*L_yield/(self.phi_s*Esh))):
            condition = 'Regime 2'
            concept = 'None'
            sigma_sr = self.fs_yield + 2*(((tau_b0*L_yield/self.phi_s) - sqrt((self.fs_yield-
                self.Es*eps_m)*tau_b1*L_yield/self.phi_s*(tau_b0/tau_b1-self.Es/Esh)+
                self.Es/Esh*tau_b0*tau_b1*L_yield**2/self.phi_s**2))/
                (tau_b0/tau_b1 - self.Es/Esh))
            sigma_sm = self.fs_yield-((sigma_sr-self.fs_yield)**2*self.phi_s/(4*tau_b1*
                                L_yield))*(tau_b0/tau_b1-1)+(sigma_sr-
                                self.fs_yield)*(tau_b0/tau_b1)-(tau_b0*L_yield/self.phi_s)
            sigma_cm = self.rho_s/(1-self.rho_s)*(sigma_sr-sigma_sm)
            eps_sm = eps_m
            eps_cm = sigma_cm/self.Ecm
            eps_sr = sigma_sr/self.Es
            Lt = max([L_yield,xcr0])
            wcr = (eps_sm-eps_cm)*Lt
            
        elif (eps_m > (self.fs_yield/self.Es + tau_b1*L_yield/(self.phi_s*Esh)) and 
            eps_m <= self.fs_yield/self.Es+(self.fs_ult-self.fs_yield)/Esh-tau_b1*L_yield/(self.phi_s*Esh)):
            condition = 'Regime 3'
            concept = 'None'
            sigma_sr = self.fs_yield+Esh*(eps_m-self.fs_yield/self.Es)+(tau_b1*L_yield/self.phi_s)
            sigma_sm = self.fs_yield + Esh*(eps_m-self.fs_yield/self.Es)
            sigma_cm = self.rho_s/(1-self.rho_s)*(sigma_sr-sigma_sm)
            eps_sm = eps_m
            eps_cm = sigma_cm/self.Ecm
            eps_sr = sigma_sr/self.Es
            Lt = max([L_yield,xcr0])
            wcr = (eps_sm-eps_cm)*Lt
            
        else:
            eps_sr = (1+xi)/(delta+xi)*eps_m       
            if eps_sr_S > eps_sr_cr:
                condition = 'Regime 1 - Condition 1'
                if eps_sr < eps_sr_cr:
                    if eps_sr >= self.fs_yield/self.Es:
                        concept = 'CLLM_yielding'
                        (eps_sr, wcr) = functions.CLLM_yield(eps_m,L_calc,self.phi_s,self.rho_s,self.Es,self.Ecm,Esh,
                                                alpha_E,delta,gamma,beta,
                                                self.fs_yield,self.tau_max,self.u1,self.alpha)
                    else: 
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CLLM(eps_sr,delta,gamma,beta,xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                
                elif eps_sr >= eps_sr_cr:
                    condition = 'Regime 1 - Condition 1 and Condition 2 for new cracked member'
                    L_calc = xcr0
                    eps_sr_S = (2*gamma)**(1/(2*delta))*((L_calc/2)*delta)**(beta/(2*delta)) # Eq. (72) Limit steel strain at symmetry section
                    if eps_sr < eps_sr_S:
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CLLM(eps_sr,delta,gamma,beta,xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                    elif eps_sr >= eps_sr_S:
                        concept = 'CHLM'
                        for i in range(0,50):
                            eps_sr = eps_m/beta_sm
                            # (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L_calc,delta,gamma,beta,xi,eps_sr_cr,self.psi,self.tau_max,self.u1,self.alpha,xcr0)
                            while eps_sr:
                                (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L_calc,delta,gamma,beta,xi,eps_sr_cr,self.psi,self.tau_max,self.u1,self.alpha,xcr0)
                                if eps_cm_cover_max >= eps_ctm:
                                    # print('MEMBER CRACKED i.e. NEW MEMBER LENGTH L_calc = L_calc/2 CHOSEN')
                                    L_calc = L_calc/2
                                elif eps_cm_cover_max < eps_ctm:
                                    break
                            if abs(eps_m - eps_sm) < 1e-10:
                                break
                            beta_sm = eps_sm/eps_sr
            
            elif eps_sr_S <= eps_sr_cr:
                condition = 'Regime 1 - Condition 2'
                if eps_sr < eps_sr_S:
                    if eps_sr >= self.fs_yield/self.Es:
                        concept = 'CLLM_yielding'
                        (eps_sr, wcr) = functions.CLLM_yield(eps_m,L_calc,self.phi_s,self.rho_s,self.Es,self.Ecm,Esh,
                                                alpha_E,delta,gamma,beta,
                                                self.fs_yield,self.tau_max,self.u1,self.alpha)
                    else:
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CLLM(eps_sr,delta,gamma,beta,xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                
                elif eps_sr >= eps_sr_S:
                    concept = 'CHLM'
                    for i in range(0,50):
                        eps_sr = eps_m/beta_sm
                        # (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L_calc,delta,gamma,beta,xi,eps_sr_cr,self.psi,self.tau_max,self.u1,self.alpha,xcr0)
                        while eps_sr:
                            (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L_calc,delta,gamma,beta,xi,eps_sr_cr,self.psi,self.tau_max,self.u1,self.alpha,xcr0)
                            if eps_cm_cover_max >= eps_ctm:
                                # print('MEMBER CRACKED i.e. NEW MEMBER LENGTH L_calc = L_calc/2 CHOSEN')
                                L_calc = L_calc/2
                            elif eps_cm_cover_max < eps_ctm:
                                break
                        if abs(eps_m - eps_sm) < 1e-10:
                            break
                        beta_sm = eps_sm/eps_sr

            # Output for plotting           
            plot_dict = {
                'xcoord':xcoord[0],
                'u':u[0],
                'tau':tau[0],
                'eps_s':eps_s[0]*1000,
                'eps_c':eps_c[0]*1000,
                'eps_sm':eps_sm_list[0]*1000,
                'eps_cm':eps_cm_list[0]*1000
            }
            self.df = pd.DataFrame.from_dict(plot_dict)

        # Output for accessing class attributes
        self.condition = condition
        self.concept = concept
        self.eps_sm = eps_sm
        self.eps_cm = eps_cm
        self.eps_sr = eps_sr
        self.sigma_sr = eps_sr*self.Es
        self.xcr0 = xcr0
        self.Lt = Lt
        self.wcr = wcr

