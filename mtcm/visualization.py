from . import functions
from .solvers import mtcm

import numpy as np
import plotly.express as px
import plotly.graph_objects as go

class PlotMTCM():
    
    def __init__(self, mtcm_obj) -> None:
        
        for attr in dir(mtcm_obj):
            value = getattr(mtcm_obj, attr)
            if (not callable(value)) and (not attr.startswith('_')):
                setattr(self, attr, value)
                
    def chord_distribution(self):
        """Method for plotting distributions of slip, bond stress and strains in chord
        """
        
        try: 
        
            # Plot slip
            fig = px.line(
                self.df,x="xcoord",y="u",
                title='Slip',
                template='plotly_dark',
                labels={
                    'xcoord': u'$x [mm]$',
                    'u': u'$u [mm]$'
                }
            )
            fig.show()

            # Plot bond stresses
            fig = px.line(
                self.df,x="xcoord",y="tau",
                title='Bond stress',
                template='plotly_dark',
                labels={
                    'xcoord': u'$x [mm]$',
                    'tau': u'$\u03C4 [MPa]$'
                }
            )
            fig.show()
            
            # Plot strains
            fig = go.Figure()

            fig.add_trace(go.Scatter(
                x=self.df['xcoord'],
                y=self.df['eps_s'],
                name=u"$\u03B5_{s}$"
            ))

            fig.add_trace(go.Scatter(
                x=self.df['xcoord'],
                y=self.df['eps_c'],
                name=u"$\u03B5_{c}$"
            ))

            fig.add_trace(go.Scatter(
                x=self.df['xcoord'],
                y=self.df['eps_sm'],
                name=u"$\u03B5_{sm}$"
            ))

            fig.add_trace(go.Scatter(
                x=self.df['xcoord'],
                y=self.df['eps_cm'],
                name=u"$\u03B5_{cm}$"
            ))

            fig.update_layout(
                title="Stress vs. Strain",
                xaxis_title=u"$x [m]$",
                yaxis_title=u"$\u03B5 [‰]$",
                template='plotly_dark',
            )

            fig.show()
            
        except AttributeError:
            
            print(f'ERROR: The RC tie is in "{self.condition} condition". Plot of chord distributions is only possible for Regime 1, i.e. for tensile stresses at crack below yielding. ')
        
    def stress_vs_strain(self,
        eps_m_max: float = 3.0*1e-3,
        eps_m_min: float = 0*1e-3,
        show_stress_from_calc: bool = False,
    ):
        """Method for plotting stress vs. strain 
        
        Kwargs:
            eps_m_max: Maximum mean strain to plot, deafult is 3.0*1e-3
            eps_m_min: Minimum mean strain to plot, default is 0
            show_stress_from_calc: Show corresponding stress to a mean strain from a calculation, default is False
        """
        
        # Instantiate object
        mtmc_bar = mtcm(
            self.Lc,self.As_tot,self.n_phi_s,self.phi_s,self.rho_s,
            self.Es,self.Ecm,self.fctm,self.fs_yield,self.fs_ult,self.eps_ult,
            self.zeta,self.psi,self.u1,self.tau_max,self.alpha,self.beta_sm
        )
        
        # Calculate stresses
        eps_m = np.linspace(eps_m_min,eps_m_max)
        sigma_sr_mtcm = []
        sigma_sr_nakedsteel = []
        for eps in eps_m:
            
            # MTCM
            mtmc_bar.strain(eps)
            sigma_sr_mtcm.append(mtmc_bar.sigma_sr)
            
            # Naked steel
            (sigma_sr, sigma_sm, eps_sr) = functions.nakedsteel(eps,self.Es,self.fs_yield,self.fs_ult,self.eps_ult)
            sigma_sr_nakedsteel.append(sigma_sr)
        
        # Plot stress vs. strain        
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=eps_m,
            y=sigma_sr_mtcm,
            name="MTCM"
        ))
        
        fig.add_trace(go.Scatter(
            x=eps_m,
            y=sigma_sr_nakedsteel,
            name="Naked steel"
        ))
        
        if show_stress_from_calc:
            fig.add_trace(go.Scatter(
                x=[self.eps_sm],
                y=[self.sigma_sr],
                name="Stress level"
            ))

        fig.update_layout(
            title="Stress vs. Strain",
            xaxis_title=u"$\u03B5_{m}$",
            yaxis_title=u"$\u03C3_{sr} [MPa]$",
            template='plotly_dark',
        )

        fig.show()

    