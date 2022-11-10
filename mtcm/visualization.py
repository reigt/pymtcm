from . import functions
from .solvers import mtcm
from .solvers import smtcm

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

class PlotBase():
    
    def __init__(self, mtcm_obj) -> None:
        
        # Obtain attributes from object argument
        for attr in dir(mtcm_obj):
            value = getattr(mtcm_obj, attr)
            if (not callable(value)) and (not attr.startswith('_')):
                setattr(self, attr, value)
                
        # Reset plotting dictionary
        self._reset_plot_dict()

    def _reset_plot_dict(self):
        
        # Create empty dict for plotting
        self.plot_dict = {
            'mtcm': [],
            'smtcm': [],
            'naked_steel': [],
        }

    def _gen_data(
        self,
        eps_m_min: float=0,
        eps_m_max: float=None,
    ):
        
        # Set defaults
        if eps_m_max is None:
            eps_m_max = 2*self.fs_yield/self.Es
        
        # Instantiate MTCM object
        mtmc_bar = mtcm(
            self.phi_s,self.n_phi_s,self.hc_ef,self.wc_ef,
            self.Es,self.Ecm,self.fctm,self.Lc,self.fs_yield,self.fs_ult,self.eps_ult,
            self.zeta,self.psi,self.u1,self.tau_max,self.alpha
        )
        
        # Instantiate SMTCM object
        smtmc_bar = smtcm(
            self.phi_s,self.n_phi_s,self.hc_ef,self.wc_ef,
            self.Es,self.Ecm,self.fctm,self.Lc,self.fs_yield,self.fs_ult,self.eps_ult,
            self.zeta,self.psi,self.u1,self.tau_max,self.alpha
        )

        # Calculate stresses
        eps_m = np.linspace(eps_m_min,eps_m_max)

        # Dictionary for MTCM
        mtcm_dict = {
            'eps_m': eps_m,
            'sigma_sr': [],
            'sigma_sr_nakedsteel': [],
            'Lt': [],
            'wcr': []
        }

        # Lists for SMTCM
        smtcm_dict = {
            'eps_m': eps_m,
            'sigma_sr': [],
            'sigma_sr_nakedsteel': [],
            'Lt': [],
            'wcr': []
        }

        for eps in eps_m:
            
            # Naked steel
            (sigma_sr, sigma_sm, eps_sr) = functions.nakedsteel(eps,self.Es,self.fs_yield,self.fs_ult,self.eps_ult)
            
            # MTCM
            mtmc_bar.strain(eps)
            mtcm_dict['sigma_sr'].append(mtmc_bar.sigma_sr)
            mtcm_dict['Lt'].append(mtmc_bar.Lt)
            mtcm_dict['wcr'].append(mtmc_bar.wcr)
            mtcm_dict['sigma_sr_nakedsteel'].append(sigma_sr)

            # SMTCM
            smtmc_bar.strain(eps)
            smtcm_dict['sigma_sr'].append(smtmc_bar.sigma_sr)
            smtcm_dict['Lt'].append(smtmc_bar.Lt)
            smtcm_dict['wcr'].append(smtmc_bar.wcr)
            smtcm_dict['sigma_sr_nakedsteel'].append(sigma_sr)
            
        # Generate dataframes
        self.df_mtcm = pd.DataFrame.from_dict(mtcm_dict)
        self.df_smtcm = pd.DataFrame.from_dict(smtcm_dict)       
        
class PlotMTCM(PlotBase):
    
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
            fig = go.Figure()

            fig.add_trace(go.Scatter(
                x=self.df['xcoord'],
                y=self.df['tau'],
                name=u"$\u03C4 [MPa]$"
            ))

            fig.add_trace(go.Scatter(
                x=[min(self.df['xcoord']),max(self.df['xcoord'])],
                y=[self.tau_m,self.tau_m],
                name=u"$\u03C4_{m} [MPa]$"
            ))
            
            fig.update_layout(
                title='Bond stress',
                xaxis_title=u'$x [mm]$',
                yaxis_title=u'$\u03C4 [MPa]$',
                template='plotly_dark',
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
                title="Strains",
                xaxis_title=u"$x [m]$",
                yaxis_title=u"$\u03B5 [‰]$",
                template='plotly_dark',
            )

            fig.show()
            
        except AttributeError:
            
            print(f'ERROR: The RC tie is in "{self.condition} condition". Plot of chord distributions is only possible for Regime 1, i.e. for tensile stresses at crack below yielding. ')
        
    def tie_response(self,
        eps_m_min: float = 0*1e-3,
        eps_m_max: float = None,
        add_stress_strain_plot: bool = True,
        add_crack_width_plot: bool = True,
        add_transfer_length_plot: bool = True,
        save_html: bool = False,
        add_smtcm: bool = False,
    ):
        """Method for plotting stress vs. strain 
        
        Kwargs:
            eps_m_max: Maximum mean strain to plot, deafult is 3.0*1e-3
            eps_m_min: Minimum mean strain to plot, default is 0
            show_stress_from_calc: Show corresponding stress to a mean strain from a calculation, default is False
        """
        
        # Set defaults
        if eps_m_max is None:
            eps_m_max = 2*self.fs_yield/self.Es
        
        # Generate plotting data
        self._gen_data(
            eps_m_min = eps_m_min,
            eps_m_max = eps_m_max
        )
        
        # Plot stress vs. strain
        if add_stress_strain_plot:
            fig = go.Figure()
            
            fig.add_trace(go.Scatter(
                x=self.df_mtcm['eps_m'],
                y=self.df_mtcm['sigma_sr'],
                name="MTCM"
            ))
            
            fig.add_trace(go.Scatter(
                x=self.df_mtcm['eps_m'],
                y=self.df_mtcm['sigma_sr_nakedsteel'],
                name="Naked steel"
            ))

            if add_smtcm:
                fig.add_trace(go.Scatter(
                    x=self.df_smtcm['eps_m'],
                    y=self.df_smtcm['sigma_sr'],
                    name="SMTCM"
                ))

            try:
                fig.add_trace(go.Scatter(
                    x=[self.eps_sm],
                    y=[self.sigma_sr],
                    name="Stress level"
                ))
            
            except AttributeError:
                None

            fig.update_layout(
                title="Stress vs. Strain",
                xaxis_title=u"$\u03B5_{m} [‰]$",
                yaxis_title=u"$\u03C3_{sr} [MPa]$",
                template='plotly_dark',
            )

            fig.show()
            
            # Save HTML
            if save_html:
                fig.write_html("stress_vs_strain.html")

        # Plot mean strains vs. crack widths
        if add_crack_width_plot:
            
            fig = go.Figure()

            fig.add_trace(go.Scatter(
                x=self.df_mtcm['eps_m'],
                y=self.df_mtcm['wcr'],
                name='MTCM',
                # showlegend=False
            ))

            if add_smtcm:
                fig.add_trace(go.Scatter(
                    x=self.df_smtcm['eps_m'],
                    y=self.df_smtcm['wcr'],
                    name='SMTCM',
                    # showlegend=False
                ))

            try:
                fig.add_trace(go.Scatter(
                    x=[self.eps_sm],
                    y=[self.wcr],
                    name="Crack width at stress level"
                ))
            
            except AttributeError:
                None

            fig.update_layout(
                title="Crack widths",
                xaxis_title=u"$\u03B5_m [‰]$",
                yaxis_title=u"$w [mm]$",
                template='plotly_dark',
                xaxis_range=[0,self.fs_yield/self.Es],
            )

            fig.show()
            
            # Save HTML
            if save_html:
                fig.write_html("crack_widths.html")

        # Plot mean strains vs. transfer length
        if add_transfer_length_plot:
            
            fig = go.Figure()

            fig.add_trace(go.Scatter(
                x=self.df_mtcm['eps_m'],
                y=self.df_mtcm['Lt'],
                name='MTCM'
                # showlegend=False
            ))

            if add_smtcm:
                fig.add_trace(go.Scatter(
                    x=self.df_smtcm['eps_m'],
                    y=self.df_smtcm['Lt'],
                    name='SMTCM'
                ))

            try:
                fig.add_trace(go.Scatter(
                    x=[self.eps_sm],
                    y=[self.Lt],
                    name="Transfer length at stress level"
                ))
            
            except AttributeError:
                None

            fig.update_layout(
                title="Transfer lengths",
                xaxis_title=u"$\u03B5_m [‰]$",
                yaxis_title=u"$2L_{t} [mm]$",
                template='plotly_dark',
                xaxis_range=[0,self.fs_yield/self.Es],
            )

            fig.show()
            
            # Save HTML
            if save_html:
                fig.write_html("transfer_length.html")

    def plot_miso(self):
        
        # Plot MTCM
        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=self.eps_m_vector,
            y=self.sigma_sr_mtcm_strain,
            name="MTCM"
        ))

        fig.add_trace(go.Scatter(
            x=self.eps_m_vector,
            y=self.sigma_sr_naked_steel,
            name="Naked steel"
        ))

        fig.add_trace(go.Scatter(
            x=self.eps_m_miso,
            y=self.sigma_sr_miso,
            name="MISO"
        ))

        fig.update_layout(
            title="Stress vs. Strain MTCM MISO",
            xaxis_title=u"$\u03B5 [‰]$",
            yaxis_title=u"$\u03C3_{sr} [MPa]$",
            template='plotly_dark',
            xaxis_range=[0,1.0*max(self.eps_m_vector)],
        )

        fig.show()

        # Plot discretized MTCM
        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=self.eps_el_miso_list,
            y=self.sigma_sr_pl_miso_list,
            name="Elastic strains"
        ))

        fig.add_trace(go.Scatter(
            x=self.eps_pl_miso_list,
            y=self.sigma_sr_pl_miso_list,
            name="Plastic strains"
        ))

        fig.add_trace(go.Scatter(
            x=self.eps_m_miso,
            y=self.sigma_sr_miso,
            name="Total strains"
        ))

        fig.add_trace(go.Scatter(
            x=self.eps_m_vector,
            y=self.sigma_sr_naked_steel,
            name="Naked steel"
        ))

        fig.update_layout(
            title="Discretized stress vs. strain for MTCM MISO",
            xaxis_title=u"$\u03B5 [‰]$",
            yaxis_title=u"$\u03C3_{sr} [MPa]$",
            template='plotly_dark',
            xaxis_range=[0,1.0*max(self.eps_m_vector)],
        )

        fig.show()
