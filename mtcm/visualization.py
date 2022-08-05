import numpy as np
import plotly.express as px

class PlotMTCM():
    
    def __init__(self, mtcm_obj) -> None:
        for attr in dir(mtcm_obj):
            value = getattr(mtcm_obj, attr)
            if (not callable(value)) and (not attr.startswith('_')):
                setattr(self, attr, value)
                
    def chord_distribution(self):
        """Method for plotting distributions of slip, bond stress and strains in chord
        """
        
        # Plot slip
        fig = px.line(
            self.df,x="xcoord",y="u",
            title='Slip',
            labels={
                'xcoord': 'x [mm]',
                'u': 'u [mm]'
            }
        )
        fig.show()

        # Plot bond stresses
        fig = px.line(
            self.df,x="xcoord",y="tau",
            title='Bond stress',
            labels={
                'xcoord': 'x [mm]',
                'tau': '\u03C4 [MPa]'
            }
        )
        fig.show()
        
        # Plot strains
        fig = px.line(
            self.df,x="xcoord",y=["eps_s","eps_c","eps_sm","eps_cm"],
            title='Strains',
            labels={
                'xcoord': 'x [mm]',
            }
        )
        fig.update_layout(yaxis_title='\u03B5 [1e-3]')
        fig.show()
        
    def stress_from_strain(self,
        # eps_m_max: float = 1.0
    ):
        pass
        # eps_m_list = np.linspace(0,eps_m_max)
        # for eps_m in np.linspace(0,eps_m_max):
            
    