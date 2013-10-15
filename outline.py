import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib

from copy import deepcopy


class planet
    ''' Represents a snapshot of an evoloving planet, with methods
    for comparing different snapshots and relating thermal evolution to 
    time.'''

    def __init__(self, params=None, pcopy=None):
        ''' Create planet snapshot, either with new parameters or with 
        another planet object.'''
        self.layers = []
        self.boundaries = []
        self.mass = None
        self.radius = None
        # tracks whether structure (radii,T-profiles) have been calculated since modifications were made
        self.structured = False 
        if not pcopy is none:
            self.layers = deepcopy(pcopy.layers)
            self.boundaries = deepcopy(pcopy.boundaries)
            self.mass = deepcopy(pcopy.mass)
            self.radius = deepcopy(pcopy.layers)

    def setLayers(self,layers):
        ''' Define the layers to be included in the planet snapshot.
        Note: incormation not copied'''
        self.layers = layers

        # add layers consistant with the boundaries
        self.boundaries = []
        for layer in layers:
            lbounds = layer.boundaries
            self.boundaries += [ lb for lb in lbounds if lb not in self.boundaries ]

    def setBoundaries(self,boundaries):
        ''' Define the layers to be included in the planet snapshot.
        Note: Information not copied'''
        self.boundaries = boundaries

        # add layers consistant with the boundaries
        self.layers = []
        for bound in boundaries:
            blayers = bound.layers
            self.layers += [ bl for bl in blayers if bl not in self.layers ]

      # Following functions will probably be feature of the more specific model,
      # might want to add placeholder functions, however.

#     def stepT(dT,boundary):
#         ''' Simulated timestep in which a boundaries lower temperature is
#         changed by dT'''

#     def QtoTime(self,dQ):
#         ''' Placeholder, should be defined by specific model '''

#     def diff(comp_snapshot):
#         ''' Compare energy and entropy budgets of various layers. Relates the
#         difference to a heat flow for one snapshot to evolve into the previous one.'''

#     def structure(T,boundary):
#         ''' Integrate Strucuture, specifing the temperature at a boundary '''

    def budgets():
        ''' Report integrated energy and entropy budgets for layers in the planet'''
        pass

class xlznCorePlanet(planet):
    ''' Represent particular situation of a crystallizing metallic core with a 
    simple model for mantle evoloution'''
    
    def __init__(self,planet_file=None,pcopy=None):
        ''' Initizialize xlznCorePlanet object. Will eventually want to take input
        file names as arguments.'''
        
        if not pcopy is None:
            self.mantle = deepcopy(pcopy.mantle)
            self.core = deepcopy(pcopy.core)
            self.layers = deepcopy(pcopy.layers)
            self.boundaries = deepcopy(pcopy.boundaries)
            self.surface = deepcopy(pcopy.surface)
            self.cmb = deepcopy(pcopy.cmb)
            self.icb = deepcopy(pcopy.icb)
            self.center = deepcopy(pcopy.center)
            self.mass = pcopy.mass
            self.radius = pcopy.radius
            self.mantle_density = pcopy.mantle_density
            self.core_mass = pcopy.core_mass
            self.structured = pcopy.structured
            return None

        # Define materials from text files
        rock = material('simpleSilicateMantle.txt')
        liqFeS = material('binaryFeS.txt',mtype='liquid')
        solFe = material('simpleSolidIron.txt')
        
        # Define melting/solidifcation relationships under consideration
        liqFeS.interpLiquidus('FeS',solFe)

        self.mantle = layer('mantle',material=rock)
        self.core = xlznLayer('core', material=liqFeS, layernames=['outer_core','inner_core'],
                boundnames=['icb'])

        # set list of layers (self.core.layers should initially return a sinlge liquid layer
        self.layers = [self.mantle] + self.core.layers
    
        # set list of boundaries
        self.boundaries = GenerateBoundaries(self.layers,['surface','cmb','center'])
        self.surface = boundaries[0], self.cmb = boundaries[1], self.center = boundaries[-1]
        self.icb = None # indicates icb has yet to form

        # read in parameters from file and decide whether this is new planet or a continuation
        # at a different condition
        try:
            params = parseParamFile(open(planet_file,'r'))
        except:
            raise Exception('File not found: {}'.format(planet_file))

        try:
            cmb.r = params.Rcore
            self.mass = params.Mplanet
            self.radius = params.Rplanet
            core_radius = params.Rcore

            self.cmb.r = core_radius
            self.surface.r = self.radius
            self.mantle_density = None
            self.core_mass = None

            self.center.P = P0
            self.cmb.T = params.Tcmb
            smode = 'initial'
        except:
            try:
                self.mass = params.Mplanet
                self.mantle_density = params.rhomantle
                self.core_mass = params.Mcore
                self.core.M = core_mass

                self.center.P = P0
                self.cmb.T = params.Tcmb

                self.radius = None

                smode = 'cont'
            except:
                raise Exception('Invalid input file.')

        # Integrate structure (making entire stucture consistant with starting values
        success = self.structure(mode=smode)

        # should check to make sure the integration succeeded
        if success:

            self.radius = self.surface.r
            #self.mantle_density = 
            self.core_mass = pcopy.core.M
            self.structured = pcopy.structured

            Mplanet = self.mass
            Mcore = self.core_mass
            rhomantle = self.mantle_density
            PO = self.center.P
            params = [ Mplanet Mcore rhomantle P0]

            # write new parameter file to run a planet with consistent mass and mantle density
            writeParamFile(open('./xlzncoreplanet_cont.txt','w'),params)
            return None
        else
            raise Exception



            

        


                


class layer(object):
    ''' A layer is a portion of a planet with an adiabatic temperature profile,
    composed of a single material '''

    def __init__(self,name='',mass=None,material=None,comp=None):
        self.name = name
        self.boundaries = None
        self.mass = mass
        self.material = material
        self.comp = comp

    def specEnergy():
        pass




class liquidLayer(layer):
    ''' A liquid layer has a specific entropy which can be related to viscous and/or
    ohmic dissipation.'''

    def specEntropy():
        pass


class xlznLayer(object):
    ''' Defines a layer of liquid material that is free to crystalize upon cooling.
    Contains a list of solids with a corresponding liquidus. Upon intersection 
    of the liquidus, there are three possible occurances upon intersection with 
    a liquidus.
        1)  solid more dense and adiabat remains below liquidus to the bottom of the
            layer, forming a settled region at the bottom.
        2)  Identical case with less dense settling to the top
        3)  'Snow' regime, where sink/floating crystals would remelt before settling
    for 1) and 2) a separate solid layer is formed. For 3) the liquid adiabat
    is instead constrained to follow the liquidus'''
    self.name
    self.liquid
    self.solids
    self.comp # Mass ratio of different components
    self.mass
    self.liquidi # a liquidi corresponding to each solid phase
    self.adiabat # a modified adiabat, following the liquidus in a 'snow' region.


class boundary(object):

    self.T # upper and lower temperature
    self.d
    self.layers

    def calcEnergy():
        pass

    def calcEntropy():
        pass


class Material(object):
    '''Class for keeping track of various physical properties of a material and how
    these vary as a function of P, T and composition.'''

    self.liquidus
    self.components

    self.td_params # holds functions for returning thermodynamic parameters

    def interp_liquidus(data_file,solid):
        pass

    def set_td_params(self,param_file):
        pass

def shellIntegral(funcs,r0=0.,r1=1.,tols=[],limits=[]):
    ''' Integrate an arbitrary number of function
