# -*- coding: utf-8 -*-
"""
Created on Fri Jun  6 11:38:40 2025

@author: Michele Schiavo
"""


from python_anesthesia_simulator.pd_models import BIS_model

# BIS_model object to test the Bouillon propofol-remifentanil interaction model
bouillon_bis_model = BIS_model('Bouillon')

bouillon_bis_model.plot_surface()


# BIS_model object to test the Vanluchene propofol model
vanluchene_bis_model = BIS_model('Vanluchene')

