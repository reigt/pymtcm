# Import python built-in modules
import numpy as np
   
# Import math functions
from math import pi as pi
from math import sqrt as sqrt
from math import factorial as factorial

def CLLM(eps_sr,delta,gamma,beta,xi,psi,Lc,tau_max,u1,alpha):
    """Function comparatively lightly loaded member behaviour
    
    Args:    
        eps_sr: steel strain at the crack, e.g. sigma_sr/Es = 2.0*1e-3
        delta: Eq. (41)
        gamma: Eq. (37)
        beta: Eq. (36)
        xi: Eq. (27)
        psi: MTCM parameter related to ratio between strains at rebar level 
            and the mean strains over the cover, default is 0.7
        Lc: Uncracked length of member [mm]
        tau_max: MTCM bond slip parameter [MPa], default is 5.0
        u1: MTCM bond slip parameter [mm], default is 0.1
        alpha: MTCM bond slip parameter, default is 0.35
    """
    
    # Parameters
    u0 = (eps_sr**2/(2*gamma))**(1/beta)                                        # Eq. (50) Slip at the loaded end [mm]
    xr = (1/delta)*(eps_sr*(1/(2*gamma))**(1/(2*delta)))**(2*delta/beta)        # Eq. (51) Transfer length [mm]  
    eps_sm = ((xi*eps_sr*xr+u0)/(1+xi))/xr                                      # Mean steel strains over bar length
    eps_cm = (psi*xi*(eps_sr*xr-u0)/(1+xi))/xr                                  # Mean steel strains over bar length
    eps_cm_cover_max = psi*xi/(1+xi)*eps_sr                                     # Eq. (73) Mean concrete strains over cover 
    Lt = 2*xr                                                                   # Transfer length [mm]
    wcr = Lt*(eps_sm-eps_cm)                                                    # Crack width [mm]
    
    # Plotting
    STEPS = 30
    xcoord = np.zeros((1,STEPS))
    xcoord[0,:] = np.linspace(0,Lc/2,STEPS)
    u = np.zeros((1,STEPS))
    tau = np.zeros((1,STEPS))
    eps_s = np.zeros((1,STEPS))
    eps_c = np.zeros((1,STEPS))
    eps_sm_list = np.zeros((1,STEPS))
    eps_cm_list = np.zeros((1,STEPS))
    
    
    for i in range(0,STEPS):
            if xcoord[0,i] < xr: 
                u[0,i] = (delta*sqrt(2*gamma)*(xr-xcoord[0,i]))**(1/delta) # Eq. (49) Slip [mm]
                tau[0,i] = tau_max*(u[0,i]/u1)**alpha                           # Eq. (33) Bond stresses [MPa]
                eps_s[0,i] = ((xi*eps_sr + (2*gamma)**(1/(2*delta))*(delta*
                    (xr-xcoord[0,i]))**(beta/(2*delta)))/(1+xi))                # Eq. (52) Steel strains                 
                eps_c[0,i] = (xi*(eps_sr - (2*gamma)**(1/(2*delta))*(delta*
                    (xr-xcoord[0,i]))**(beta/(2*delta)))/(1+xi))                # Eq. (53) Concrete strains at interface
                eps_sm_list[0,i] = eps_sm
                eps_cm_list[0,i] = eps_cm

            elif xcoord[0,i] >= xr:
                u[0,i] = 0
                eps_s[0,i] = xi/(1+xi)*eps_sr                                   # Eq. (52) Steel strains beyond xr
                eps_c[0,i] = xi/(1+xi)*eps_sr                                   # Eq. (53) Concrete strains at interface beyond xr
                eps_sm_list[0,i] = eps_sm
                eps_cm_list[0,i] = eps_cm

    return (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list)
        
def CLLM_yield(eps_m,L,phi_s,rho_s,Es,Ecm,Esh,alpha_E,delta,gamma,beta,
               fs_yield,tau_max,u1,alpha):
    """Function for yielding of comparatively lightly loaded member behaviour
    
    Args:
    """
    
    xr = (1/delta)*(fs_yield/Es*(1/(2*gamma))**(1/(2*delta)))**(2*delta/beta)
    tau_b0 = ((tau_max/u1**alpha)*(delta*sqrt(2*gamma))**
              (alpha/delta)*(delta/(alpha+delta))*xr**
              ((alpha+delta)/delta))/xr
    tau_b1 = tau_b0/2
    alpha_yield = 1+alpha_E*rho_s
    
    L_yield = ((phi_s*fs_yield*Esh)/(4*tau_b1*alpha_yield*Es)*
                (sqrt(1+4*alpha_yield*Es/Esh*(L*tau_b1/
                (phi_s*fs_yield)*(alpha_yield*Es*eps_m/fs_yield-
                alpha_E*rho_s)-tau_b1/(4*alpha_yield*tau_b0)))-1))
    if not np.isreal(L_yield):
        L_yield = 0
    L_yield = min(max(L_yield,0),0.5*L)
    
    sigma_sr = fs_yield+L_yield*4*tau_b1/phi_s
    sigma_sm = Es*eps_m
    sigma_cm = rho_s/(1-rho_s)*(sigma_sr-sigma_sm)
    eps_sm = eps_m
    eps_cm = sigma_cm/Ecm
    eps_sr = sigma_sr/Es
    Lt = L_yield
    wcr = (eps_sm-eps_cm)*Lt
    
    return (eps_sr, wcr)
        
def CHLM(eps_sr,L,delta,gamma,beta,xi,eps_sr_cr,psi,tau_max,u1,alpha,xcr0):
    """Function comparatively heavily loaded member behaviour
    
    Args:    
        eps_sr: steel strain at the crack, e.g. sigma_sr/Es = 2.0*1e-3
        L: Member length [mm]
        delta: Eq. (41)
        gamma: Eq. (37)
        beta: Eq. (36)
        eps_sr_cr: Eq. (74) Cracking strain at the end of transfer length
        psi: MTCM parameter related to ratio between strains at rebar level 
            and the mean strains over the cover, default is 0.7
        tau_max: MTCM bond slip parameter [MPa], default is 5.0
        u1: MTCM bond slip parameter [mm], default is 0.1
        alpha: MTCM bond slip parameter, default is 0.35
        xcr0: Crack spacing [mm]
    """
    
    udCASE = (eps_sr**2/(4*gamma))**(1/beta)                                    # Eq. (62) or Eq. (69)
    u0max = (eps_sr**2/(2*gamma))**(1/beta)                                     # Eq. (57)
            
    Delta_x = 0.1                                                               # Case 2 parameter dx1 + dx2 (See section 4.4.3 and Fig. 7(b)) 
    Delta_u = 5.8000e-05                                                        # Case 2 parameter du (See section 4.4.3 and Fig. 7(b))
    m = 10                                                                      # Number of chosen terms in Eq. (39) and Eq. (40) (MAX 170 TERMS BECAUSE factorial(171) = infty)
    R = -1/2                                                                    # Falling factorial (Pocchammer) coefficient
    
    ## Case check
    # We need to determine wether Case 1 or Case 2 occurs according to the 
    # discussions in Section 4.4.4. In general, each term in the series is 
    # calculated before summed in the end.
    
    u0_CHECK = (eps_sr**2/(4*gamma)-Delta_u)**(1/beta)                          # Eq. (62) Choose a value of u0 close to the limit value discriminating Case 1 and Case 2
    C_CHECK = eps_sr**2/2 - gamma*u0_CHECK**beta                                # Eq. (56) Integration constant for case check     

    r = list(range(0,m+1))                                                      # Vector for calculating the falling factorial 
    n = list(range(0,m+1))                                                      # Vector for calculating the factorial k!
    bink = list(range(0,m+1))                                                   # Vector for the binkomial coefficients
    FCASE_CHECK = list(range(0,m+1))                                            # Vector for function in Case 1 check
    fCASE_CHECK = list(range(0,m+1))                                            # Vector for function in Case 1 check
    
    for k in range(0,m):
        r[0] = 1
        r[k+1] = r[k]*(R-(k+1)+1)                                               # Falling factorial  
    
    for k in range(0,m+1):             
        n[k] = factorial(k)                                                # The factorial k!    
        bink[k] = r[k]/n[k]                                                     # The binkomial coefficients in Eq. (61)
        FCASE_CHECK[k] = (gamma/C_CHECK)**k*(u0_CHECK**(1+k*beta)/(1+k*beta))   # The function term in Eq. (61)        
        fCASE_CHECK[k] = bink[k]*FCASE_CHECK[k]                                 # Eq. (61) The binkomial coefficients multiplied with the function terms     

    fCASE1_CHECK = L/2 - 1/sqrt(2*C_CHECK)*sum(fCASE_CHECK)                # Eq. (61) 
    
    
    ## Determine u0
    # Now that we know which case that occurs, we can determine u0. 
    FCASE1 = list(range(0,m+1))                                                 # Eq. (61) Vector for function in Case 1     
    dFCASE1 = list(range(0,m+1))                                                # Eq. (78) Vector for the derivatives of function in Case 1    
    F1CASE2 = list(range(0,m+1))                                               # Eq. (66) Vector for function 1 in Case 2 
    F2CASE2 = list(range(0,m+1))                                               # Eq. (67) Vector for function 1 in Case 2     
    F3CASE2 = list(range(0,m+1))                                               # Eq. (68) Vector for function 1 in Case 2 
    
    dF1CASE2 = list(range(0,m+1))                                              # Eq. (80) Vector for the derivatives of function 1 in Case 2    
    dF2CASE2 = list(range(0,m+1))                                              # Eq. (81) Vector for the derivatives of function 2 in Case 2    
    dF3CASE2 = list(range(0,m+1))                                              # Eq. (82) Vector for the derivatives of function 3 in Case 2
    
    FB2 = list(range(0,m+1))                                                   # Eq. (60) Vector for function in integration constant B2
    
    if fCASE1_CHECK < 0:
        
        CASE = 'CASE 1'
        
        u0itr = udCASE - Delta_u                                                # Initial value for u0 in Case 1
        
        for j in range(0,20):
            
            u0 = min(u0itr,u0max-Delta_u)                                      # Iterated value for u0, however not exceeding u0max in Eq. (58) due to Eq. (56)            
            C = eps_sr**2/2 - gamma*u0**beta                                    # Eq. (56) 
            
            for k in range(0,m+1):
                FCASE1[k] = bink[k]*(gamma**k*(1/C)**(1/2+k)*(u0**(1+k*beta))/
                    (1+k*beta))                                               # Eq. (61)
                    
                dFCASE1[k] = bink[k]*(gamma**k*((gamma*beta*u0**(beta-1)*(1/2+k)*
                            C**(-3/2-k)*u0**(1+k*beta)/(1+k*beta))+
                            (C**(-(1/2+k))*u0**(k*beta))))                  # Eq. (78)
                
            fCASE1_u0 = L/2 - 1/sqrt(2)*sum(FCASE1)                       # Eq. (61)
            dfCASE1_u0 = -1/sqrt(2)*sum(dFCASE1)                          # Eq. (78)
            u0itr = u0 - fCASE1_u0/dfCASE1_u0                                  # Eq. (70)
            
            if abs(fCASE1_u0) < 1e-4:
                break                                   
        
    elif fCASE1_CHECK > 0:
        
        CASE = 'CASE 2'
        
        u0itr = udCASE - Delta_u                                                # Initial value for u0 in Case 1
        
        for j in range(0,20):
            
            u0 = min(u0itr,u0max-Delta_u)                                      # Iterated value for u0, however not exceeding u0max in Eq. (58) due to Eq. (56)
            C = eps_sr**2/2 - gamma*u0**beta                                    # Eq. (56) 
            
            for k in range(0,m+1): 
                F1CASE2[k] = (bink[k]*(C/gamma)**k*(u0**(delta-k*beta)/
                            (delta-k*beta)))                                  # Eq. (66)
                    
                F2CASE2[k] = (bink[k]*(((C/gamma)**(k/(delta-k*beta)+1/beta)+
                            Delta_u*(C/gamma)**(k/(delta-k*beta)))**
                            (delta-k*beta))/(delta-k*beta))                   # Eq. (67)
                
                F3CASE2[k] = (bink[k]*(gamma**k*(((1/gamma)**(1/beta)*C**
                            ((2-beta)/(2*beta*(1+k*beta))))-(Delta_u*C**
                            (-(1/2+k)/(1+k*beta))))**(1+k*beta))/(1+k*beta))  # Eq. (68)
                
                dF1CASE2[k] = bink[k]*((1/gamma)**k*(C**k*u0**(delta-k*beta-1)-
                            ((gamma*beta*k*C**(k-1))/(delta-k*beta))*u0**
                            (beta*(1-k)+delta-1)))                           # Eq. (80)
                
                dF2CASE2[k] = bink[k]*((((C/gamma)**(k/(delta-k*beta)+
                            (1/beta))+Delta_u*(C/gamma)**(k/(delta-
                            k*beta)))**(delta-k*beta-1))*
                            (-gamma*beta*u0**(beta-1))*(((1/gamma)**
                            (k/(delta-k*beta)+1/beta)*(k/(delta-k*beta)+
                            1/beta)*C**(k/(delta-k*beta)+1/beta-1))+
                            (Delta_u*(1/gamma)**(k/(delta-k*beta))*
                            (k/(delta-k*beta))*C**(k/(delta-k*beta)-1))))    # Eq. (81)
                
                dF3CASE2[k] = bink[k]*(gamma**k*((((1/gamma)**(1/beta)*C**
                            ((2-beta)/(2*beta*(1+k*beta))))-(Delta_u*C**
                            (-(1/2+k)/(1+k*beta))))**(k*beta))*
                            (-gamma*beta*u0**(beta-1))*(((1/gamma)**
                            (1/beta)*((2-beta)/(2*beta*(1+k*beta)))*C**
                            ((2-beta)/(2*beta*(1+k*beta))-1))+(Delta_u*
                            ((1/2+k)/(1+k*beta))*C**(-((1/2+k)/
                            (1+k*beta)+1)))))                                # Eq. (82)
                
                FB2[k] = bink[k]*((C/gamma)**k*(u0**(delta-k*beta)/
                        (delta-k*beta)))                                      # Eq. (60)
                
            fCASE2 = (L/2 - (1/sqrt(2*gamma))*(sum(F1CASE2)-sum(F2CASE2))-
                    (1/sqrt(2))*sum(F3CASE2) - Delta_x)                  # Eq. (65)       
            
            dfCASE2 = (-1/sqrt(2*gamma)*(sum(dF1CASE2)-sum(dF2CASE2))-
                    1/sqrt(2)*sum(dF3CASE2))                            # Eq. (79)
            
            u0itr = u0 - fCASE2/dfCASE2                                        # Eq. (70)
            
            if abs(fCASE2) < 1e-4:
                break
    
    B1 = L/2                                                                   # Eq. (59) 
    B2 = (1/sqrt(2*gamma))*sum(FB2)                                       # Eq. (60)
    ud = (eps_sr**2/(2*gamma) - u0**beta)**(1/beta)                            # Eq. (58)
    eps_c_max = xi*(eps_sr-sqrt(2*C))/(1+xi)                              # Eq. (71)
    eps_cm_cover_max = psi*eps_c_max                                           # Eq. (71)   
    eps_sm = ((xi*eps_sr*(L/2)+u0)/(1+xi))/(L/2)                               # Mean steel strains
    eps_cm = (psi*xi*(eps_sr*(L/2)-u0)/(1+xi))/(L/2)                           # Mean concrete strains    
    
    # Plotting
    STEPS = 30
    u = np.zeros((1,STEPS))
    u[0,:] = np.linspace(0,u0,STEPS)
    FX = np.zeros((1,STEPS))
    xcoord = np.zeros((1,STEPS))
    tau = np.zeros((1,STEPS))
    eps_s = np.zeros((1,STEPS))
    eps_c = np.zeros((1,STEPS))
    eps_sm_list = np.zeros((1,STEPS))
    eps_cm_list = np.zeros((1,STEPS))

    for i in range(0,STEPS):
        for k in range(0,len(bink)):
            if u[0,i] < ud:
                FX[0,i] = FX[0,i] + bink[k]*(gamma**k*(1/C)**(1/2+k)*(u[0,i]**(1+k*beta)/(1+k*beta))) 
            elif u[0,i] > ud:
                FX[0,i] = FX[0,i] + bink[k]*((C/gamma)**k*(u[0,i]**(delta-k*beta)/(delta-k*beta)))
        if u[0,i] < ud:
            xcoord[0,i] = B1 - 1/sqrt(2)*FX[0,i] 
        elif u[0,i] > ud:
            xcoord[0,i] = B2 - 1/sqrt(2*gamma)*FX[0,i] 
        
        tau[0,i] = tau_max*(u[0,i]/u1)**alpha                                  # Eq. (33) Bond stresses [MPa]
        eps_s[0,i] = ((xi*eps_sr + sqrt(2*(gamma*u[0,i]**beta+C)))/
                        (1+xi))                                               # Eq. (44) Steel strains
        eps_c[0,i] = (xi*(eps_sr - sqrt(2*(gamma*u[0,i]**beta+C)))/
                        (1+xi))                                               # Eq. (45) Steel strains       
        eps_sm_list[0,i] = eps_sm
        eps_cm_list[0,i] = eps_cm
    
    # Calculate crack widths
    u0_cllm = (eps_sr_cr**2/(2*gamma))**(1/beta)                                        # Eq. (50) Slip at the loaded end [mm]
    xr_cllm = (1/delta)*(eps_sr_cr*(1/(2*gamma))**(1/(2*delta)))**(2*delta/beta)        # Eq. (51) Transfer length [mm]  
    eps_sm_cllm = ((xi*eps_sr_cr*xr_cllm+u0_cllm)/(1+xi))/xr_cllm                                      # Mean steel strains over bar length
    eps_cm_cllm = (psi*xi*(eps_sr_cr*xr_cllm-u0_cllm)/(1+xi))/xr_cllm                                  # Mean steel strains over bar length
    Lt_cllm = 2*xr_cllm                                                                   # Transfer length [mm]
    wcr_cllm = Lt_cllm*(eps_sm_cllm-eps_cm_cllm)                                                    # Crack width [mm]

    Lt = xcr0                                                                      # Crack spacing [mm]    
    wcr = max(Lt*(eps_sm-eps_cm),wcr_cllm)                                                    # Crack width [mm]

    return (u0, u0, eps_sm, eps_cm, eps_cm_cover_max, Lt, wcr, xcoord, u, tau, eps_s, eps_c, eps_sm_list, eps_cm_list)

def nakedsteel(eps_m,Es,fs_yield=500,fs_ult=640,eps_ult=100e-03):
    """Bilinear material model for naked steel
    
    Args:
        eps_m: steel strain at the crack, e.g. sigma_sr/Es = 2.0*1e-3
    
    Kwargs:
        Es: Youngs Modulus for steel [MPa]
        fs_yield: Steel stress at yielding
        fs_ult: Steel stress at ultimate strength [MPa]
        eps_ult: Steel strain at ultimate strength
    """ 
    
    Esh = (fs_ult-fs_yield)/(eps_ult-fs_yield/Es)

    if abs(eps_m) <= fs_yield/Es:
        sigma_sr = Es*eps_m
        sigma_sm = sigma_sr
        eps_sr = eps_m
    elif abs(eps_m) > fs_yield/Es and eps_m > 0:
        sigma_sr = fs_yield+Esh*(eps_m-fs_yield/Es)
        sigma_sm = sigma_sr
        eps_sr = sigma_sr/Es
    elif abs(eps_m) > fs_yield/Es and eps_m < 0:
        sigma_sr = -fs_yield+Esh*(eps_m+fs_yield/Es)
        sigma_sm = sigma_sr
        eps_sr = sigma_sr/Es
            
    return (sigma_sr, sigma_sm, eps_sr)

def prestressingsteel(eps_m,eps_p0,Ep,fp_yield,fp_ult,eps_p_ult):
    """Bilinear material model for prestressing steel
    
    Args:    
        eps_m: steel strain at the crack, e.g. sigma_sr/Es = 2.0*1e-3
        Es: Youngs Modulus for steel [MPa]
        fp_yield: Steel stress at yielding
        fp_ult: Steel stress at ultimate strength [MPa]
        eps_p_ult: Steel strain at ultimate strength
    """ 
    
    eps_p = eps_m+eps_p0
    Eph = (fp_ult-fp_yield)/(eps_p_ult-fp_yield/Ep)

    if abs(eps_p) <= fp_yield/Ep:
        sigma_pm = Ep*eps_m
        sigma_pr = Ep*eps_p
    elif abs(eps_p) > fp_yield/Ep:
        sigma_pm = (fp_yield/Ep-eps_p0)*Ep+((eps_p0+eps_m)-fp_yield/Ep)*Eph      
        sigma_pr = fp_yield+Eph*(eps_m-fp_yield/Ep)
        
    return (sigma_pm, sigma_pr, eps_p)
