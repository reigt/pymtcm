from . import functions

import pandas as pd
import numpy as np

from math import pi as pi
from math import sqrt as sqrt

class mtcm():
    
    def __init__(self) -> None:
        self._contains_data = False
        
    def stress(self,
        eps_sr: float,
        Lc: float,
        As_tot: float,
        n_phi: float,
        phi_s: float,
        rho_s: float,
        Es: float,
        Ecm: float,
        fctm: float,
        zeta: float=1.0,
        psi: float=0.70,
        u1: float=0.1,
        tau_max: float=5.0,
        alpha: float=0.35
    ):
        """Method for calculating cracks based on steel stresses in cracks in a member
        
        Args:    
            eps_sr: steel strain at the crack, e.g. sigma_sr/Es = 2.0*1e-3
            Lc: Member length [mm]
            As_tot: Total steel rebar area [mm2]
            n_phi: Number of reinforcing steel bars 
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

        L_calc = Lc
        alpha_E = Es/Ecm                                                            # Modular ratio
        eps_ctm = fctm/Ecm                                                          # Cracking strain

        beta = 1+alpha                                                              # Eq. (36)
        delta = (1-alpha)/2                                                         # Eq. (41)
        xi = alpha_E*rho_s/psi                                                   # Eq. (27)    
        chi = (1+xi)*(zeta*n_phi*pi*phi_s/(As_tot*Es))                                  # Eq. (31)
        gamma = chi*tau_max/(beta*u1**alpha)                                        # Eq. (37)
        xcr0 = (1/delta)*((1/psi)*((1+xi)/xi)*eps_ctm*(1/(2*gamma))
            **(1/(2*delta)))**(2*delta/beta)                                     # Eq. (75) Crack spacing [mm]
        eps_sr_cr = (1/psi)*(1 + (1/xi))*eps_ctm                                    # Eq. (74) Cracking strain at the end of transfer length  

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
                    (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CLLM(eps_sr,delta,gamma,beta,xi,psi,Lc,tau_max,u1,alpha)
                elif eps_sr >= eps_sr_cr:
                    condition = 'Condition 1 and Condition 2 for new cracked member'
                    L_calc = xcr0
                    eps_sr_S = (2*gamma)**(1/(2*delta))*((L_calc/2)*delta)**(beta/(2*delta)) 
                    if eps_sr < eps_sr_S:
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CLLM(eps_sr,delta,gamma,beta,xi,psi,Lc,tau_max,u1,alpha)
                    elif eps_sr >= eps_sr_S:
                        concept = 'CHLM'
                        # (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L_calc,delta,gamma,beta,xi,eps_sr_cr,psi,tau_max,u1,alpha,xcr0)
                        while eps_sr:
                            (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L_calc,delta,gamma,beta,xi,eps_sr_cr,psi,tau_max,u1,alpha,xcr0)
                            if eps_cm_cover_max >= eps_ctm:
                                # print('MEMBER CRACKED i.e. NEW MEMBER LENGTH L_calc = L_calc/2 CHOSEN')
                                L_calc = L_calc/2
                            elif eps_cm_cover_max < eps_ctm:
                                break
            elif eps_sr_S <= eps_sr_cr:
                condition = 'Condition 2'
                if eps_sr < eps_sr_S:
                    concept = 'CLLM'
                    (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CLLM(eps_sr,delta,gamma,beta,xi,psi,Lc,tau_max,u1,alpha)
                elif eps_sr >= eps_sr_S:
                    concept = 'CHLM'
                    # (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L_calc,delta,gamma,beta,xi,eps_sr_cr,psi,tau_max,u1,alpha,xcr0)
                    while eps_sr:
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L_calc,delta,gamma,beta,xi,eps_sr_cr,psi,tau_max,u1,alpha,xcr0)
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
            Lc: float,
            As_tot: float,
            n_phi: float,
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
        ):
        """
        Function for calculating cracks based on imposed mean strains in a member
        
        Syntax: strain(eps_m,Lc,As_tot,n_phi,phi_s,rho_s,Es,Ecm,fctm,
                    zeta=1.0,psi=0.70,u1=0.1,tau_max=5.0,alpha=0.35)
        
        eps_m: mean imposed strain for the member, e.g. from shrinkage or temperature loads
        Lc: Member length [mm]
        As_tot: Total steel rebar area [mm2]
        n_phi: Number of reinforcing steel bars 
        phi_s: Steel rebar diameter [mm] 
        rho_s: Reinforcement ratio, e.g. As_tot/Ac_ef
        Es: Youngs Modulus for steel [MPa]
        Ecm: Youngs modulus for concrete [MPa]
        fctm: Tensile strength [MPa]
        fs_yield: Steel stress at yielding, default 500 [MPa]
        fs_ult: Steel stress at ultimate strength, default 640 [MPa]
        eps_ult: Steel strain at ultimate strength, default 100e-3
        zeta: MTCM parameter related to uniform bond stresses, default is 1.0 
        psi: MTCM parameter related to ratio between strains at rebar level 
            and the mean strains over the cover, default is 0.7
        u1: MTCM bond slip parameter [mm], default is 0.1
        tau_max: MTCM bond slip parameter [MPa], default is 5.0
        alpha: MTCM bond slip parameter, default is 0.35
        """

        L = Lc
        alpha_E = Es/Ecm                                                            # Modular ratio
        eps_ctm = fctm/Ecm                                                          # Cracking strain

        beta = 1+alpha                                                              # Eq. (36)
        delta = (1-alpha)/2                                                         # Eq. (41)
        xi = alpha_E*rho_s/psi                                                      # Eq. (27)    
        chi = (1+xi)*(zeta*n_phi*pi*phi_s/(As_tot*Es))                                  # Eq. (31)
        gamma = chi*tau_max/(beta*u1**alpha)                                        # Eq. (37)
        xcr0 = (1/delta)*((1/psi)*((1+xi)/xi)*eps_ctm*(1/(2*gamma))
            **(1/(2*delta)))**(2*delta/beta)                                     # Eq. (75) Crack spacing [mm]
        eps_sr_cr = (1/psi)*(1 + (1/xi))*eps_ctm                                    # Eq. (74) Cracking strain at the end of transfer length  
        
        eps_sr_S = (2*gamma)**(1/(2*delta))*((L/2)*delta)**(beta/(2*delta))        # Eq. (72) Limit steel strain at symmetry section
        
        self.stress(fs_yield/Es,Lc,As_tot,n_phi,phi_s,rho_s,Es,Ecm,fctm)
        L_yield = self.Lt
        x_yield = self.plot_dict['xcoord']
        tau_yield = self.plot_dict['tau']
        tau_m_yield = abs(np.trapz(tau_yield,x_yield)/max(x_yield))
        tau_b0 = tau_m_yield
        tau_b1 = tau_b0/2
        Esh = (fs_ult-fs_yield)/(eps_ult-fs_yield/Es)
        
        if eps_m <= 0:
            condition = 'Compression or zero'
            concept = 'None'
            (sigma_sr, sigma_sm, eps_sr) = functions.nakedsteel(eps_m,Es,fs_yield,fs_ult,eps_ult)
            sigma_cm = 0
            eps_sm = eps_m
            eps_cm = eps_m
            Lt = 0
            wcr = (eps_sm-eps_cm)*Lt
            
        elif (eps_m > fs_yield/Es-tau_b0*L_yield/(phi_s*Es) and 
            eps_m <= (fs_yield/Es+tau_b1*L_yield/(phi_s*Esh))):
            condition = 'Regime 2'
            concept = 'None'
            sigma_sr = fs_yield + 2*(((tau_b0*L_yield/phi_s) - sqrt((fs_yield-
                Es*eps_m)*tau_b1*L_yield/phi_s*(tau_b0/tau_b1-Es/Esh)+
                Es/Esh*tau_b0*tau_b1*L_yield**2/phi_s**2))/
                (tau_b0/tau_b1 - Es/Esh))
            sigma_sm = fs_yield-((sigma_sr-fs_yield)**2*phi_s/(4*tau_b1*
                                L_yield))*(tau_b0/tau_b1-1)+(sigma_sr-
                                fs_yield)*(tau_b0/tau_b1)-(tau_b0*L_yield/phi_s)
            sigma_cm = rho_s/(1-rho_s)*(sigma_sr-sigma_sm)
            eps_sm = eps_m
            eps_cm = sigma_cm/Ecm
            eps_sr = sigma_sr/Es
            Lt = max([L_yield,xcr0])
            wcr = (eps_sm-eps_cm)*Lt
            
        elif (eps_m > (fs_yield/Es + tau_b1*L_yield/(phi_s*Esh)) and 
            eps_m <= fs_yield/Es+(fs_ult-fs_yield)/Esh-tau_b1*L_yield/(phi_s*Esh)):
            condition = 'Regime 3'
            concept = 'None'
            sigma_sr = fs_yield+Esh*(eps_m-fs_yield/Es)+(tau_b1*L_yield/phi_s)
            sigma_sm = fs_yield + Esh*(eps_m-fs_yield/Es)
            sigma_cm = rho_s/(1-rho_s)*(sigma_sr-sigma_sm)
            eps_sm = eps_m
            eps_cm = sigma_cm/Ecm
            eps_sr = sigma_sr/Es
            Lt = max([L_yield,xcr0])
            wcr = (eps_sm-eps_cm)*Lt
            
        else:
            eps_sr = (1+xi)/(delta+xi)*eps_m       
            if eps_sr_S > eps_sr_cr:
                condition = 'Regime 1 - Condition 1'
                if eps_sr < eps_sr_cr:
                    if eps_sr >= fs_yield/Es:
                        concept = 'CLLM_yielding'
                        (eps_sr, wcr) = functions.CLLM_yield(eps_m,L,phi_s,rho_s,Es,Ecm,Esh,
                                                alpha_E,delta,gamma,beta,
                                                fs_yield,tau_max,u1,alpha)
                    else: 
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CLLM(eps_sr,delta,gamma,beta,xi,psi,Lc,tau_max,u1,alpha)
                
                elif eps_sr >= eps_sr_cr:
                    condition = 'Regime 1 - Condition 1 and Condition 2 for new cracked member'
                    L = xcr0
                    eps_sr_S = (2*gamma)**(1/(2*delta))*((L/2)*delta)**(beta/(2*delta)) # Eq. (72) Limit steel strain at symmetry section
                    if eps_sr < eps_sr_S:
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CLLM(eps_sr,delta,gamma,beta,xi,psi,Lc,tau_max,u1,alpha)
                    elif eps_sr >= eps_sr_S:
                        concept = 'CHLM'
                        for i in range(0,50):
                            eps_sr = eps_m/beta_sm
                            # (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L,delta,gamma,beta,xi,eps_sr_cr,psi,tau_max,u1,alpha,xcr0)
                            while eps_sr:
                                (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L,delta,gamma,beta,xi,eps_sr_cr,psi,tau_max,u1,alpha,xcr0)
                                if eps_cm_cover_max >= eps_ctm:
                                    # print('MEMBER CRACKED i.e. NEW MEMBER LENGTH L = L/2 CHOSEN')
                                    L = L/2
                                elif eps_cm_cover_max < eps_ctm:
                                    break
                            if abs(eps_m - eps_sm) < 1e-10:
                                break
                            beta_sm = eps_sm/eps_sr
            
            elif eps_sr_S <= eps_sr_cr:
                condition = 'Regime 1 - Condition 2'
                if eps_sr < eps_sr_S:
                    if eps_sr >= fs_yield/Es:
                        concept = 'CLLM_yielding'
                        (eps_sr, wcr) = functions.CLLM_yield(eps_m,L,phi_s,rho_s,Es,Ecm,Esh,
                                                alpha_E,delta,gamma,beta,
                                                fs_yield,tau_max,u1,alpha)
                    else:
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CLLM(eps_sr,delta,gamma,beta,xi,psi,Lc,tau_max,u1,alpha)
                
                elif eps_sr >= eps_sr_S:
                    concept = 'CHLM'
                    for i in range(0,50):
                        eps_sr = eps_m/beta_sm
                        # (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L,delta,gamma,beta,xi,eps_sr_cr,psi,tau_max,u1,alpha,xcr0)
                        while eps_sr:
                            (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list) = functions.CHLM(eps_sr,L,delta,gamma,beta,xi,eps_sr_cr,psi,tau_max,u1,alpha,xcr0)
                            if eps_cm_cover_max >= eps_ctm:
                                # print('MEMBER CRACKED i.e. NEW MEMBER LENGTH L = L/2 CHOSEN')
                                L = L/2
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
        self.sigma_sr = eps_sr*Es
        self.xcr0 = xcr0
        self.Lt = Lt
        self.wcr = wcr

