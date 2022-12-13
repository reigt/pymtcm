from . import functions

import pandas as pd
import numpy as np

from math import pi as pi
from math import sqrt as sqrt

class mtcm():
    
    def __init__(self,
        phi_s: float,
        n_phi_s: float,
        hc_ef: float,
        wc_ef: float,
        Es: float,
        Ecm: float,
        fctm: float,
        Lc: float=10000,
        fs_yield: float=500,
        fs_ult: float=640,
        eps_ult: float=100e-3,
        zeta: float=1.0,
        psi: float=0.70,
        u1: float=0.1,
        tau_max: float=5.0,
        alpha: float=0.35,
    ) -> None:
        """Method for instantiating RC tie
        
        Args:    
            phi_s: Steel rebar diameter [mm] 
            n_phi_s: Number of reinforcing steel bars 
            hc_ef: Effective height of RC tie [mm]
            wc_ef: Effective width of RC tie [mm]
            Es: Youngs Modulus for steel [MPa]
            Ecm: Youngs modulus for concrete [MPa]
            fctm: Tensile strength [MPa]
        
        Kwargs:
            Lc: Uncracked length of member [mm], default is 10000
            fs_yield: Yield strength of steel [MPa], defaults is 500
            fs_ult: Ultimate strength of steel [MPa], default is 640
            eps_ult: Ultimate strain of steel, default is 100e-3
            zeta: MTCM parameter related to uniform bond stresses, default is 1.0 
            psi: MTCM parameter related to ratio between strains at rebar level 
                and the mean strains over the cover, default is 0.7
            u1: MTCM bond slip parameter [mm], default is 0.1
            tau_max: MTCM bond slip parameter [MPa], default is 5.0
            alpha: MTCM bond slip parameter, default is 0.35
        """        
        
        # RC tie geometry
        self.phi_s = phi_s
        self.n_phi_s = n_phi_s
        self.hc_ef = hc_ef
        self.wc_ef = wc_ef
        self.Lc = Lc

        # RC tie materials
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
        
        # RC tie parameters
        self.As_tot = n_phi_s*pi*phi_s**2/4
        self.Ac_ef = hc_ef*wc_ef
        self.rho_s = self.As_tot/self.Ac_ef
        
        # MTCM parameters
        self.alpha_E = Es/Ecm
        self.eps_ctm = fctm/Ecm
        self.beta = 1+alpha
        self.delta = (1-alpha)/2
        self.xi = self.alpha_E*self.rho_s/psi
        self.chi = (1+self.xi)*(self.zeta*self.n_phi_s*pi*self.phi_s/(self.As_tot*self.Es))
        self.gamma = self.chi*self.tau_max/(self.beta*self.u1**self.alpha)
        self.xcr0 = (1/self.delta)*((1/self.psi)*((1+self.xi)/self.xi)*self.eps_ctm*(1/(2*self.gamma))**(1/(2*self.delta)))**(2*self.delta/self.beta)

        # Limit state
        self.eps_sr_cr = (1/self.psi)*(1 + (1/self.xi))*self.eps_ctm  

    def stress(self,
        sigma_sr: float,
    ):
        """Method for calculating cracks based on steel stresses in cracks in a member
        
        Args:    
            sigma_sr: Steel stress at crack [MPa]
        """
        
        # Steel strains
        eps_sr = sigma_sr/self.Es
        
        # Steel strains at symmetry section
        L_calc = self.Lc
        eps_sr_S = (2*self.gamma)**(1/(2*self.delta))*((L_calc/2)*self.delta)**(self.beta/(2*self.delta))        # Eq. (72) Limit steel strain at symmetry section

        if eps_sr <= 0:
            condition = 'Compression or zero'
            concept = 'None'
            Lt = 0
            eps_sm = eps_sr
            eps_cm = eps_sr
            wcr = (eps_sm-eps_cm)*Lt
            
        else:
            if eps_sr_S > self.eps_sr_cr:
                condition = 'Regime 1 - Condition 1'
                if eps_sr < self.eps_sr_cr:
                    concept = 'CLLM'
                    (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list, tau_m) = functions.CLLM(eps_sr,self.delta,self.gamma,self.beta,self.xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                elif eps_sr >= self.eps_sr_cr:
                    condition = 'Condition 1 and Condition 2 for new cracked member'
                    L_calc = self.xcr0
                    eps_sr_S = (2*self.gamma)**(1/(2*self.delta))*((L_calc/2)*self.delta)**(self.beta/(2*self.delta)) 
                    if eps_sr < eps_sr_S:
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list, tau_m) = functions.CLLM(eps_sr,self.delta,self.gamma,self.beta,self.xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                    elif eps_sr >= eps_sr_S:
                        concept = 'CHLM'
                        while eps_sr:
                            (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list, tau_m) = functions.CHLM(eps_sr,L_calc,self.delta,self.gamma,self.beta,self.xi,self.eps_sr_cr,eps_sr_S,self.psi,self.tau_max,self.u1,self.alpha,self.xcr0,condition)
                            if eps_cm_cover_max >= self.eps_ctm:
                                L_calc = L_calc/2
                            elif eps_cm_cover_max < self.eps_ctm:
                                break
            elif eps_sr_S <= self.eps_sr_cr:
                condition = 'Regime 1 - Condition 2'
                if eps_sr < eps_sr_S:
                    concept = 'CLLM'
                    (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list, tau_m) = functions.CLLM(eps_sr,self.delta,self.gamma,self.beta,self.xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                elif eps_sr >= eps_sr_S:
                    concept = 'CHLM'
                    while eps_sr:
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list, tau_m) = functions.CHLM(eps_sr,L_calc,self.delta,self.gamma,self.beta,self.xi,self.eps_sr_cr,eps_sr_S,self.psi,self.tau_max,self.u1,self.alpha,self.xcr0,condition)
                        if eps_cm_cover_max >= self.eps_ctm:
                            L_calc = L_calc/2
                        elif eps_cm_cover_max < self.eps_ctm:
                            break
            
            # Output for plotting           
            plot_dict = {
                'xcoord':xcoord,
                'u':u,
                'tau':tau,
                'eps_s':eps_s*1000,
                'eps_c':eps_c*1000,
                'eps_sm':eps_sm_list*1000,
                'eps_cm':eps_cm_list*1000
            }
            self.df = pd.DataFrame.from_dict(plot_dict)
            
            # Mean bond stress
            self.tau_m = tau_m
            
        # Output for accessing class attributes
        self.condition = condition
        self.concept = concept
        self.eps_sm = eps_sm
        self.eps_cm = eps_cm
        self.sigma_sr = eps_sr*self.Es
        self.Lt = Lt
        self.wcr = wcr
        
    def strain(self,
            eps_m: float,
        ):
        """Method for calculating cracks based on imposed mean strains in a member
        
        Args:
            eps_m: mean imposed strain for the member, e.g. from shrinkage or temperature loads
        """

        # Steel strains at symmetry section
        L_calc = self.Lc
        eps_sr_S = (2*self.gamma)**(1/(2*self.delta))*((L_calc/2)*self.delta)**(self.beta/(2*self.delta))        # Eq. (72) Limit steel strain at symmetry section
        
        # Initial values
        beta_sm = 1.0
        
        # Bond stresses at yielding
        self.stress(self.fs_yield)
        L_yield = self.Lt
        x_yield = self.df['xcoord']
        tau_yield = self.df['tau']
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
            Lt = max([L_yield,self.xcr0])
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
            Lt = max([L_yield,self.xcr0])
            wcr = (eps_sm-eps_cm)*Lt
            
        else:
            eps_sr = (1+self.xi)/(self.delta+self.xi)*eps_m       
            if eps_sr_S > self.eps_sr_cr:
                condition = 'Regime 1 - Condition 1'
                if eps_sr < self.eps_sr_cr:
                    if eps_sr >= self.fs_yield/self.Es:
                        concept = 'CLLM_yielding'
                        (eps_sr, eps_sm, eps_cm, Lt, wcr) = functions.CLLM_yield(eps_m,L_calc,self.phi_s,self.rho_s,self.Es,self.Ecm,Esh,
                                                self.alpha_E,self.delta,self.gamma,self.beta,
                                                self.fs_yield,self.tau_max,self.u1,self.alpha)
                    else: 
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list, tau_m) = functions.CLLM(eps_sr,self.delta,self.gamma,self.beta,self.xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                
                elif eps_sr >= self.eps_sr_cr:
                    condition = 'Regime 1 - Condition 1 and Condition 2 for new cracked member'
                    L_calc = self.xcr0
                    eps_sr_S = (2*self.gamma)**(1/(2*self.delta))*((L_calc/2)*self.delta)**(self.beta/(2*self.delta)) 
                    if eps_sr < eps_sr_S:
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list, tau_m) = functions.CLLM(eps_sr,self.delta,self.gamma,self.beta,self.xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                    elif eps_sr >= eps_sr_S:
                        concept = 'CHLM'
                        for i_ in range(0,50):
                            eps_sr = eps_m/beta_sm
                            while eps_sr:
                                (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list, tau_m) = functions.CHLM(eps_sr,L_calc,self.delta,self.gamma,self.beta,self.xi,self.eps_sr_cr,eps_sr_S,self.psi,self.tau_max,self.u1,self.alpha,self.xcr0,condition)
                                if eps_cm_cover_max >= self.eps_ctm:
                                    L_calc = L_calc/2
                                elif eps_cm_cover_max < self.eps_ctm:
                                    break
                            if abs(eps_m - eps_sm) < 1e-10:
                                break
                            beta_sm = eps_sm/eps_sr
            
            elif eps_sr_S <= self.eps_sr_cr:
                condition = 'Regime 1 - Condition 2'
                if eps_sr < eps_sr_S:
                    if eps_sr >= self.fs_yield/self.Es:
                        concept = 'CLLM_yielding'
                        (eps_sr, eps_sm, eps_cm, Lt, wcr) = functions.CLLM_yield(eps_m,L_calc,self.phi_s,self.rho_s,self.Es,self.Ecm,Esh,
                                                self.alpha_E,self.delta,self.gamma,self.beta,
                                                self.fs_yield,self.tau_max,self.u1,self.alpha)
                    else:
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list, tau_m) = functions.CLLM(eps_sr,self.delta,self.gamma,self.beta,self.xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                
                elif eps_sr >= eps_sr_S:
                    concept = 'CHLM'
                    for i_ in range(0,50):
                        eps_sr = eps_m/beta_sm
                        while eps_sr:
                            (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list, tau_m) = functions.CHLM(eps_sr,L_calc,self.delta,self.gamma,self.beta,self.xi,self.eps_sr_cr,eps_sr_S,self.psi,self.tau_max,self.u1,self.alpha,self.xcr0,condition)
                            if eps_cm_cover_max >= self.eps_ctm:
                                L_calc = L_calc/2
                            elif eps_cm_cover_max < self.eps_ctm:
                                break
                        if abs(eps_m - eps_sm) < 1e-10:
                            break
                        beta_sm = eps_sm/eps_sr

            # Output for plotting
            try:
                
                # Output for plotting           
                plot_dict = {
                    'xcoord':xcoord,
                    'u':u,
                    'tau':tau,
                    'eps_s':eps_s*1000,
                    'eps_c':eps_c*1000,
                    'eps_sm':eps_sm_list*1000,
                    'eps_cm':eps_cm_list*1000
                }
                df = pd.DataFrame.from_dict(plot_dict)
                
            except UnboundLocalError:
                
                df = None
                tau_m = None
                
            self.df = df
            self.tau_m = tau_m

        # Output for accessing class attributes
        self.condition = condition
        self.concept = concept
        self.eps_sm = eps_sm
        self.eps_cm = eps_cm
        self.eps_sr = eps_sr
        self.sigma_sr = eps_sr*self.Es
        self.Lt = Lt
        self.wcr = wcr

class smtcm(mtcm):

    def strain(self,
            eps_m: float,
        ):
        """Method for calculating cracks based on imposed mean strains in a member
        
        Args:
            eps_m: mean imposed strain for the member, e.g. from shrinkage or temperature loads
        """

        # Steel strains at symmetry section
        L_calc = self.Lc
        eps_sr_S = (2*self.gamma)**(1/(2*self.delta))*((L_calc/2)*self.delta)**(self.beta/(2*self.delta))

        # Strains at yielding        
        self.stress(self.fs_yield)
        L_yield = self.Lt
        x_yield = self.df['xcoord']
        tau_yield = self.df['tau']
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
            Lt = max([L_yield,self.xcr0])
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
            Lt = max([L_yield,self.xcr0])
            wcr = (eps_sm-eps_cm)*Lt
            
        else:
            eps_sr = (1+self.xi)/(self.delta+self.xi)*eps_m       
            if eps_sr_S > self.eps_sr_cr:
                condition = 'Regime 1 - Condition 1'
                if eps_sr < self.eps_sr_cr:
                    if eps_sr >= self.fs_yield/self.Es:
                        concept = 'CLLM_yielding'
                        (eps_sr, eps_sm, eps_cm, Lt, wcr) = functions.CLLM_yield(eps_m,L_calc,self.phi_s,self.rho_s,self.Es,self.Ecm,Esh,
                                                self.alpha_E,self.delta,self.gamma,self.beta,
                                                self.fs_yield,self.tau_max,self.u1,self.alpha)
                    else: 
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list, tau_m) = functions.CLLM(eps_sr,self.delta,self.gamma,self.beta,self.xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                
                elif eps_sr >= self.eps_sr_cr:
                    condition = 'Regime 1 - Condition 1 and Condition 2 for new cracked member'
                    L_calc = self.xcr0
                    eps_sr_S = (2*self.gamma)**(1/(2*self.delta))*((L_calc/2)*self.delta)**(self.beta/(2*self.delta))
                    if eps_sr < eps_sr_S:
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list, tau_m) = functions.CLLM(eps_sr,self.delta,self.gamma,self.beta,self.xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                    elif eps_sr >= eps_sr_S:
                        concept = 'SCHLM'
                        (eps_sr, eps_sm, eps_cm, Lt, wcr) = functions.SCHLM_strain(
                            eps_m,
                            L_calc,
                            self.delta,
                            eps_sr_S,
                            self.gamma,
                            self.beta,
                            self.xi,
                            self.eps_sr_cr,
                            self.psi,
                            self.fs_yield,
                            self.Es
                        )
            
            elif eps_sr_S <= self.eps_sr_cr:
                condition = 'Regime 1 - Condition 2'
                if eps_sr < eps_sr_S:
                    if eps_sr >= self.fs_yield/self.Es:
                        concept = 'CLLM_yielding'
                        (eps_sr, eps_sm, eps_cm, Lt, wcr) = functions.CLLM_yield(eps_m,L_calc,self.phi_s,self.rho_s,self.Es,self.Ecm,Esh,
                                                self.alpha_E,self.delta,self.gamma,self.beta,
                                                self.fs_yield,self.tau_max,self.u1,self.alpha)
                    else:
                        concept = 'CLLM'
                        (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list, tau_m) = functions.CLLM(eps_sr,self.delta,self.gamma,self.beta,self.xi,self.psi,self.Lc,self.tau_max,self.u1,self.alpha)
                
                elif eps_sr >= eps_sr_S:
                    concept = 'SCHLM'
                    (eps_sr, eps_sm, eps_cm, Lt, wcr) = functions.SCHLM_strain(
                        eps_m,
                        L_calc,
                        self.delta,
                        eps_sr_S,
                        self.gamma,
                        self.beta,
                        self.xi,
                        self.eps_sr_cr,
                        self.psi,
                        self.fs_yield,
                        self.Es
                    )

        # Output for accessing class attributes
        self.condition = condition
        self.concept = concept
        self.eps_sm = eps_sm
        self.eps_cm = eps_cm
        self.eps_sr = eps_sr
        self.sigma_sr = eps_sr*self.Es
        self.Lt = Lt
        self.wcr = wcr

class miso(mtcm):
        
    def ansys_input(self,
        ansfilename: str=None,
        mat_id: int=None,
        gamma_s: float=7.775,
        v_s: float=0.2,
        alpha_s: float=1.2e-5,
        eps_m_max: float=2.5e-3,
        fully_discretized: bool=False,
    ):
        """Method for generating Multilinear Isotropic model for MTCM in ANSYS. 

        Kwargs:
            ansfilename (str): ANSYS filename, default is None.
            mat_id (int): Material ID number, default is None. 
            gamma_s (float): Density steel, default is 7.775 [t/m3].
            v_s (float): Poisson's ratio, default is 0.2. 
            alpha_s (float): Temperature dilation coefficient, default is 1.2e-5.
            eps_m_max (float): Max mean strain in stress-strain curve, default is 2.5e-3. 
            fully_discretized (bool): Discretize MTCM at each discrete point, default is False
        """

        # Set defaults
        if ansfilename is None:
            ansfilename = 'miso_mtcm.ans'
        if mat_id is None:
            mat_id = 1000
            
        # Tension chord responses
        eps_m_vector = np.linspace(0,eps_m_max)
        sigma_sr_mtcm_strain = []
        sigma_sr_naked_steel = []
        scr_mtcm = []
        wcr_mtcm = []
        for eps_m in eps_m_vector:
            
            # MTCM
            self.strain(eps_m)
            sigma_sr_mtcm_strain.append(self.sigma_sr) 
            wcr_mtcm.append(self.wcr)
            scr_mtcm.append(self.Lt)
            
            # Naked steel
            (sigma_sr, sigma_sm, eps_sr_naked_steel) = functions.nakedsteel(
                eps_m,self.Es,fs_yield=self.fs_yield
            )
            sigma_sr_naked_steel.append(sigma_sr)

        # Discretize MTCM
        eps_m_list_cllm = []
        eps_m_list_chlm = []
        sigma_sr_list_cllm = []
        sigma_sr_list_chlm = []
        
        # Prior to yielding
        for eps_m,sigma_sr_mtcm in zip(eps_m_vector,sigma_sr_mtcm_strain):
            eps_sr = sigma_sr_mtcm/self.Es
            if eps_sr < self.eps_sr_cr:
                eps_m_list_cllm.append(eps_m)
                sigma_sr_list_cllm.append(sigma_sr_mtcm)
            elif eps_sr >= self.eps_sr_cr and eps_sr < self.fs_yield/self.Es:
                eps_m_list_chlm.append(eps_m)
                sigma_sr_list_chlm.append(sigma_sr_mtcm)
                
        # Post yielding
        eps_m_yield = self.fs_yield/self.Es
        self.strain(eps_m_yield)
        eps_sr_yield = self.sigma_sr
        eps_m_list_yield = [eps_m_yield,min([self.eps_ult,eps_m])]
        sigma_sr_list_yield = [eps_sr_yield,min([self.fs_ult,sigma_sr_mtcm_strain[-1]])]

        # Generate MISO for MTCM
        eps_m_miso = []
        sigma_sr_miso = []
        try:
            
            if fully_discretized:
                
                # Strains
                for this_eps in eps_m_list_cllm:
                    eps_m_miso.append(this_eps)
                
                for this_eps in eps_m_list_chlm:
                    eps_m_miso.append(this_eps)
                    
                for this_eps in eps_m_list_yield:
                    eps_m_miso.append(this_eps)
                
                # Steel stress
                for this_sigma in sigma_sr_list_cllm:
                    sigma_sr_miso.append(this_sigma)
                
                for this_sigma in sigma_sr_list_chlm:
                    sigma_sr_miso.append(this_sigma)
                    
                for this_sigma in sigma_sr_list_yield:
                    sigma_sr_miso.append(this_sigma)
                    
            else:
                
                # Strains
                eps_m_miso.append(eps_m_list_cllm[0])
                eps_m_miso.append(eps_m_list_cllm[round(len(eps_m_list_cllm)/2)])
                eps_m_miso.append(eps_m_list_cllm[-1])
                eps_m_miso.append(eps_m_list_chlm[0])
                eps_m_miso.append(eps_m_list_chlm[round(len(eps_m_list_chlm)/2)])
                eps_m_miso.append(eps_m_list_chlm[-1])
                eps_m_miso.append(eps_m_list_yield[0])
                eps_m_miso.append(eps_m_list_yield[-1])

                # Steel stress
                sigma_sr_miso.append(sigma_sr_list_cllm[0])
                sigma_sr_miso.append(sigma_sr_list_cllm[round(len(sigma_sr_list_cllm)/2)])
                sigma_sr_miso.append(sigma_sr_list_cllm[-1])
                sigma_sr_miso.append(sigma_sr_list_chlm[0])
                sigma_sr_miso.append(sigma_sr_list_chlm[round(len(sigma_sr_list_chlm)/2)])
                sigma_sr_miso.append(sigma_sr_list_chlm[-1])
                sigma_sr_miso.append(sigma_sr_list_yield[0])
                sigma_sr_miso.append(sigma_sr_list_yield[-1])
            
        except IndexError:
            
            print('IndexError: CLLM yielding, please increase reinforcement ratio!')
            pass

        # Initial values
        eps_pl_miso = eps_m_miso[1]
        sigma_sr_pl_miso = sigma_sr_miso[1]
        Es_miso = sigma_sr_pl_miso/eps_pl_miso
        sigma_sr_pl_miso_list = [sigma_sr_pl_miso]
        eps_pl_miso_list = [0]
        eps_el_miso_list = [eps_pl_miso]
        
        # Elastic and plastic strains
        for eps_tot,sigma_sr in zip(eps_m_miso[2:],sigma_sr_miso[2:]):
            eps_el_miso = sigma_sr/Es_miso
            eps_el_miso_list.append(eps_el_miso)
            eps_pl_miso_list.append(eps_tot-eps_el_miso)
            sigma_sr_pl_miso_list.append(sigma_sr)

        # Write out as ANSYS inputfile
        with open(ansfilename,'w') as file:
            file.write('! ---------------------------------------------------------------- \n')
            file.write('! INPUT FOR MULTLINEAR ISOTROPIC HARDENING CURVE FOR MTCM \n')
            file.write('! -----------------------------START------------------------------ \n \n')
            file.write('    ! Material ID \n')
            file.write(f'    mtcm_mat_id = {mat_id} \n \n')
            file.write('    ! Linear elastic parameters \n')
            file.write(f'    mp,ex  ,mtcm_mat_id,{int(Es_miso*1000)} \n')
            file.write(f'    mp,dens,mtcm_mat_id,{gamma_s} \n')
            file.write(f'    mp,nuxy,mtcm_mat_id,{v_s} \n')
            file.write(f'    mp,alpx,mtcm_mat_id,{alpha_s} \n \n')
            file.write('    ! Multilinear isotropic hardening parameters \n')
            file.write(f'    tb,plastic,mtcm_mat_id,1,{len(eps_pl_miso_list)},miso \n')
            file.write('    tbtemp,0.0 \n')
            for eps_pl,sigma_sr in zip(eps_pl_miso_list,sigma_sr_pl_miso_list):
                file.write(f'    tbpt,defi,{eps_pl:.7f},{sigma_sr*1000:.2f} \n')

            file.write('! -----------------------------END------------------------------ \n \n')
        
        # Output for accessing class attributes
        self.eps_m_vector = eps_m_vector
        self.eps_el_miso_list = eps_el_miso_list
        self.eps_pl_miso_list = eps_pl_miso_list
        self.eps_m_miso = eps_m_miso
        self.sigma_sr_mtcm_strain = sigma_sr_mtcm_strain
        self.sigma_sr_naked_steel = sigma_sr_naked_steel
        self.sigma_sr_miso = sigma_sr_miso
        self.sigma_sr_pl_miso_list = sigma_sr_pl_miso_list
        self.scr_mtcm = scr_mtcm
        self.wcr_mtcm = wcr_mtcm
        