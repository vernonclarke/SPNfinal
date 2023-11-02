'''
The MSN class defining the cell
'''

from neuron import h
import numpy as np
import json

# Distributions:
'''
T-type Ca: g = 1.0/( 1 +np.exp{(x-70)/-4.5} )
naf (den): (0.1 + 0.9/(1 + np.exp((x-60.0)/10.0)))

'''

def calculate_distribution(d3, dist, a4, a5,  a6,  a7, g8):
    '''
    Used for setting the maximal conductance of a segment.
    Scales the maximal conductance based on somatic distance and distribution type.

    Parameters:
    d3   = distribution type:
         0 linear, 
         1 sigmoidal, 
         2 exponential
         3 step function
    dist = somatic distance of segment
    a4-7 = distribution parameters 
    g8   = base conductance (similar to maximal conductance)

    '''

    if   d3 == 0: 
        value = a4 + a5*dist
    elif d3 == 1: 
        value = a4 + a5/(1 + np.exp((dist-a6)/a7) )
    elif d3 == 2: 
        value = a4 + a5*np.exp((dist-a6)/a7)
    elif d3 == 3:
        if (dist > a6) and (dist < a7):
            value = a4
        else:
            value = a5

    if value < 0:
        value = 0

    value = value*g8
    return value 




# ======================= the MSN class ==================================================

class MSN:
    def __init__(self,  params=None,                                        \
                        morphology=None,     \
                        variables=None,                                     \
                        section=None                                        ):
        Import = h.Import3d_SWC_read()
        Import.input(morphology)
        imprt = h.Import3d_GUI(Import, 0)
        imprt.instantiate(None)
        h.define_shape()
        # h.cao0_ca_ion = 2  # default in nrn
        h.celsius = 35
        self._create_sectionlists()
        self._set_nsegs(section=section)
        self.v_init = -80

        self.dendritic_channels =   [
                    "naf",      
                    "kaf",
                    "kas",
                    "kdr",
                    "kir",
                    "cal12",
                    "cal13",
                    "can",
                    "car",
                    "cav32",
                    "cav33",
                    "sk",
                    "bk",
                    "gaba1",
                    "gaba2"         ]

        self.somatic_channels = [
                    "naf",
                    "kaf",
                    "kas",
                    "kdr",
                    "kir",
                    "cal12",
                    "cal13",
                    "cav32",
                    "cav33",
                    "can",
                    "car",
                    "sk",
                    "nap",
                    "bk",
                    "gaba1",
                    "gaba2"     ]

        self.axonal_channels = [
                    "naf",
                    "kas" ,
                    "Im"       ]


        # insert active mechanisms (related to channels) -------------
        for sec in self.somalist:
            for mech in self.somatic_channels+["cadyn", "caldyn"]:
                sec.insert(mech)

        for sec in self.axonlist:
            for mech in self.axonal_channels:
                sec.insert(mech)

        for sec in self.dendlist:
            for mech in self.dendritic_channels+["cadyn", "caldyn"]:
                sec.insert(mech)

        with open(params) as file:
            par = json.load(file)

        # set passive parameters --------------------------------------------        
        for sec in self.allseclist:
            if not sec.name().find('spine') >= 0: # some mech don't exist at spine, when remaking cell causes error
                sec.Ra = 200 # increased membrane resistance 150 -> 200
                sec.cm = 1.0
                sec.insert('pas')
                #sec.g_pas = 1e-5 # set using json file
                sec.e_pas = -85 # changed to resting potential of spn
                sec.g_pas = float(par['g_pas_all']['Value'])
                sec.ena = 50
                sec.ek = -85 

        self.distribute_channels("soma", "gbar_naf",   0, 1, 0, 0, 0, float(par['gbar_naf_somatic']['Value']))
        self.distribute_channels("soma", "gbar_kaf",   0, 1, 0, 0, 0, float(par['gbar_kaf_somatic']['Value']))
        self.distribute_channels("soma", "gbar_kas",   0, 1, 0, 0, 0, float(par['gbar_kas_somatic']['Value']))
        self.distribute_channels("soma", "gbar_kdr",   0, 1, 0, 0, 0, float(par['gbar_kdr_somatic']['Value']))
        self.distribute_channels("soma", "gbar_bk",    0, 1, 0, 0, 0, float(par['gbar_bk_somatic' ]['Value']))
        self.distribute_channels("soma", "pbar_cal12", 0, 1, 0, 0, 0, 1.34e-5)
        self.distribute_channels("soma", "pbar_cal13", 0, 1, 0, 0, 0, 1.34e-6)
        self.distribute_channels("soma", "pbar_cav32", 0, 1, 0, 0, 0, 0) # off
        self.distribute_channels("soma", "pbar_cav33", 0, 1, 0, 0, 0, 0) # off
        self.distribute_channels("soma", "pbar_car",   0, 1, 0, 0, 0, 1.34e-4)
        self.distribute_channels("soma", "pbar_can",   0, 1, 0, 0, 0,    4e-5)
        self.distribute_channels("soma", "gbar_nap",   0, 1, 0, 0, 0, 0.0007)
        self.distribute_channels("dend", "gbar_kdr", 1,   0.25, 1,  50.0,  30.0, float(par['gbar_kdr_basal']['Value'])) 
        
        # use gaba1 for ohmic extrasynaptic GABA; use gaba2 for outward rectificying GABA
        # default both off; if employed, choose 1
        self.distribute_channels("soma", "gbar_gaba1", 0, 1, 0, 0, 0, float(par['gbar_gaba']['Value']))
        self.distribute_channels("dend", "gbar_gaba1", 0, 1, 0, 0, 0, float(par['gbar_gaba']['Value']))
        self.distribute_channels("soma", "gbar_gaba2", 0, 1, 0, 0, 0, float(par['gbar_gaba']['Value']))
        self.distribute_channels("dend", "gbar_gaba2", 0, 1, 0, 0, 0, float(par['gbar_gaba']['Value']))
        
        self.distribute_channels("dend", "gbar_bk",    0, 1, 0, 0, 0, float(par['gbar_bk_basal' ]['Value']))
        self.distribute_channels("dend", "pbar_cal12", 0, 1, 0, 0, 0, 1e-5)
        self.distribute_channels("dend", "pbar_cal13", 0, 1, 0, 0, 0, 1e-6)
        self.distribute_channels("dend", "pbar_car",   0, 1, 0, 0, 0, 5e-4) #1e-4
        self.distribute_channels("axon", "gbar_kas",   0, 1, 0, 0, 0,      float(par['gbar_kas_axonal']['Value']))
        self.distribute_channels("axon", "gbar_naf",   3, 1, 1.1, 30, 500, float(par['gbar_naf_axonal']['Value']))
        self.distribute_channels("axon", "gbar_Im",   0, 1, 0, 0, 0, 1.0e-3)

        if variables: # currently not using variables
            self.distribute_channels("dend", "gbar_naf", 1,   1.0-variables['naf'][1],  \
                                                              variables['naf'][1],      \
                                                              variables['naf'][2],      \
                                                              variables['naf'][3],      \
                                                              np.power(10,variables['naf'][0])*float(par['gbar_naf_basal']['Value']))
            self.distribute_channels("dend", "gbar_kaf", 1,   1.0,                      \
                                                              variables['kaf'][1],      \
                                                              variables['kaf'][2],      \
                                                              variables['kaf'][3],      \
                                                              np.power(10,variables['kaf'][0])*float(par['gbar_kaf_basal']['Value']))
            self.distribute_channels("dend", "gbar_kas", 1,   0.1,                      \
                                                              0.9,                      \
                                                              variables['kas'][1],      \
                                                              variables['kas'][2],      \
                                                              np.power(10,variables['kas'][0])*float(par['gbar_kas_basal']['Value']))

            self.distribute_channels("dend", "gbar_kir", 0,   np.power(10,variables['kir'][0]), 0, 0, 0,    float(par['gbar_kir_basal'  ]['Value']))
            self.distribute_channels("soma", "gbar_kir", 0,   np.power(10,variables['kir'][0]), 0, 0, 0,    float(par['gbar_kir_somatic']['Value']))
            self.distribute_channels("dend", "gbar_sk",  0,   np.power(10,variables['sk' ][0]), 0, 0, 0,    float(par['gbar_sk_basal'   ]['Value']))
            self.distribute_channels("soma", "gbar_sk",  0,   np.power(10,variables['sk' ][0]), 0, 0, 0,    float(par['gbar_sk_somatic' ]['Value']))

            self.distribute_channels("dend", "pbar_can",   1, 1.0-variables['can'][1],  \
                                                              variables['can'][1],      \
                                                              variables['can'][2],      \
                                                              variables['can'][3],      \
                                                              np.power(10,variables['can'][0]))
            self.distribute_channels("dend", "pbar_cav32", 1, 0,                        \
                                                              1,                        \
                                                              variables['c32'][1],      \
                                                              variables['c32'][2],      \
                                                              np.power(10,variables['c32'][0]))
            self.distribute_channels("dend", "pbar_cav33", 1, 0,                        \
                                                              1,                        \
                                                              variables['c33'][1],      \
                                                              variables['c33'][2],      \
                                                              np.power(10,variables['c33'][0]))
        else: # channel variables set here
            self.distribute_channels("dend", "gbar_naf", 1, 0, 1,   50.0,   10.0, float(par['gbar_naf_basal']['Value']))
            # self.distribute_channels("dend", "gbar_kaf", 1,   0.25, 1,  120.0,  30.0, float(par['gbar_kaf_basal']['Value'])) 

            self.distribute_channels("dend", "gbar_kaf", 1,   0.5, 0.25,  120.0,  30.0, float(par['gbar_kaf_basal']['Value'])) 
            self.distribute_channels("dend", "gbar_kas", 2,   0.25, 5,  0.0, -10.0, float(par['gbar_kas_basal']['Value']))
            self.distribute_channels("dend", "gbar_kir", 0, 1, 0, 0, 0, float(par['gbar_kir_basal']['Value']) * 2) # 2 scaling value for kir to accurately fit below-V-threshold experimental data 
            self.distribute_channels("soma", "gbar_kir", 0, 1, 0, 0, 0, float(par['gbar_kir_somatic']['Value']))
            self.distribute_channels("dend", "gbar_sk",  0, 1, 0, 0, 0, float(par['gbar_sk_basal']['Value']))
            self.distribute_channels("soma", "gbar_sk",  0, 1, 0, 0, 0, float(par['gbar_sk_basal']['Value']))
            self.distribute_channels("dend", "pbar_can", 0, 1, 0, 0, 0, 1e-7) 
            self.distribute_channels("dend", "pbar_cav32", 1, 0, 1.0, 100.0, -30.0, 1e-5) # 100, 1e-7
            self.distribute_channels("dend", "pbar_cav33", 1, 0, 1.0, 100.0, -30.0, 2.5e-6) # midpoint was 100, 1e-8, 2.5 scaling factor of cav33 used to more accurately produce up-states


    def _create_sectionlists(self):
        self.allsecnames = []
        self.allseclist  = h.SectionList()
        for sec in h.allsec():
            self.allsecnames.append(sec.name())
            self.allseclist.append(sec=sec)
        self.nsomasec = 0
        self.somalist = h.SectionList()
        for sec in h.allsec():
            if sec.name().find('soma') >= 0:
                self.somalist.append(sec=sec)
                if self.nsomasec == 0:
                    self.soma = sec
                self.nsomasec += 1
        self.axonlist = h.SectionList()
        for sec in h.allsec():
            if sec.name().find('axon') >= 0:
                self.axonlist.append(sec=sec)
        self.dendlist = h.SectionList()
        for sec in h.allsec():
            if sec.name().find('dend') >= 0:
                self.dendlist.append(sec=sec)


    def _set_nsegs(self, section=None, N=20):
        """ def seg/sec. if section: set seg ~= 1/um  """
        if section:
            dend_name = 'dend[' + str(int(section)) + ']'
            for sec in self.allseclist:
                if sec.name() == dend_name:
                    # TODO: this needs some thinking; how to best set number of segments
                    n = 2*int(sec.L/2.0)+1
                    if n > N:
                        sec.nseg = n
                    else:
                        sec.nseg = 2*(N/2) + 1 # odd number of seg
                else:
                    sec.nseg = 2*int(sec.L/40.0)+1
        else:
            for sec in self.allseclist:
                sec.nseg = 2*int(sec.L/40.0)+1
        for sec in self.axonlist:
            sec.nseg = 2  # two segments in axon initial segment

    def distribute_channels(self, as1, as2, d3, a4, a5, a6, a7, g8):
        h.distance(sec=self.soma)

        for sec in self.allseclist:

            # if right cellular compartment (axon, soma or dend)
            if sec.name().find(as1) >= 0:
                for seg in sec:
                    dist = h.distance(seg.x, sec=sec)
                    val = calculate_distribution(d3, dist, a4, a5, a6, a7, g8)
                    cmd = 'seg.%s = %g' % (as2, val)
                    exec(cmd)

class Spine():
    """
    Spine class. Create a spine with neck and head.
    inspired by Mattioni and Le Novere, (2013).
    https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=150284&file=/TimeScales-master/neuronControl/spine.py#tabs-2

    if a parent section is passed as argument, the spine will be connected to that section (using the default orientation)
        to connect a spine at a later stage use the connect_spine(parent) method, where parent is a section.

    Since default orientation is used, to place multiple spines spread over one section, first use function split_section(sec) in common_functions.
        split_section(sec) splits one section into multiple sections with retained total surface area (and close conductances).
    """
    #TODO: 
    # -test spine (one and multiple): rheobase etc retained?
    # -run in vivo with spines

    def __init__(self, id,              \
                       parent=None,     \
                       x=None,          \
                       neck_L=1.0,      \
                       neck_dia=0.1,    \
                       head_L=0.5,      \
                       head_dia=0.5,    \
                       Ra=150.0,        \
                       Cm=1.0          
                                ):
        """ Create a spine with geometry given by the arguments"""

        self.id         =   id
        self.neck       =   self.create_neck(neck_L=neck_L, neck_dia=neck_dia, Ra=Ra, Cm=Cm)
        self.head       =   self.create_head(head_L=head_L, head_dia=head_dia, Ra=Ra, Cm=Cm)
        self.parent     =   parent  # the parent section connected to the neck
        self.stim       =   None    # attribute for saving spike apparatus (netStim, synapse and
        self.x          =   x

        # connect spine parts
        self.connect_head2neck()

        # connect spine to parent
        self.connect_spine(parent, x)



    def create_neck(self, neck_L=1.0, neck_dia=0.1, Ra=150.0, Cm=1.0):
        """ Create a spine neck"""

        sec_name        =   'spine_%d_neck' % (self.id)
        neck            =   h.Section(name=sec_name)
        neck.nseg       =   1
        neck.L          =   neck_L 
        neck.diam       =   neck_dia
        neck.Ra         =   Ra 
        neck.cm         =   Cm

        for mech in [   'pas',      \
                        'cav32',    \
                        'cav33',    \
                        'caldyn'     ]:
            neck.insert(mech)

        neck(0.5).pbar_cav33 = 1e-7
        neck(0.5).pbar_cav33 = 1e-8


        neck.g_pas      =  1.25e-7
        neck.e_pas      =   -70 

        return neck



    def create_head(self, head_L=0.5, head_dia=0.5, Ra=150.0, Cm=1.0):
        """Create the head of the spine and populate it with channels"""

        sec_name        =   'spine_%d_head' % (self.id)
        head            =   h.Section(name=sec_name)

        head.nseg       =   1
        head.L          =   head_L
        head.diam       =   head_dia
        head.Ra         =   Ra
        head.cm         =   1.0

        for mech in [   'pas',      \
                        'kir',      \
                        'cav32',    \
                        'cav33',    \
                        'car',      \
                        'cal12',    \
                        'cal13',    \
                        'cadyn',    \
                        'sk',       \
                        'caldyn'    ]:
            head.insert(mech)

        head(0.5).pbar_cav32 = 1e-7
        head(0.5).pbar_cav33 = 1e-8
        head(0.5).pbar_car   = 1e-8
        head(0.5).pbar_cal12 = 1e-7
        head(0.5).pbar_cal13 = 1e-6 #1e-8
        head(0.5).gbar_kir   = 1e-7
        head(0.5).gbar_sk    = 0.00002

        head.g_pas      =   1.25e-7
        head.e_pas      =   -70 

        return head


    def connect_head2neck(self, parent=None):
        ''' connect spine head to neck and if parent is not None connect spine to parent 
        connection is hardcoded to use default connection orientation (0 end connected to 1 end on parent)
        To use other orientation, first create and then connect spine, 
            using the connect spine method with updated arguments
        '''
        self.head.connect(self.neck(1),0)
        if not parent is None:
            self.neck.connect(parent(1),0) 

    def connect_spine(self, parent, x=1, end=0):
        ''' connect spine to parent sec.
        '''
        self.neck.connect(parent(x),end)
        self.parent = parent

    def move_spine(self, new_parent):
        ''' move spine from one section to another (using default orientation)
        '''
        h.disconnect(sec=self.neck)
        self.neck.connect(new_parent(1),0)
        self.parent = new_parent

def add_spines(cell, spines_per_sec): 
    SPINES = {}
    ID = 0
    for sec in cell.dendlist:
        SPINES[sec.name()] = {}
        for i,seg in enumerate(sec):

            x = h.distance(seg) # distance to soma(0.5)

            if 80 > x > 30: # none 0-30, 30+ linear increase
                b = 0.2 # spines_per_sec / (80-30)
                n_spines = int(spines_per_sec + (x-80)*b)

                for j in range(n_spines):
                    SPINES[sec.name()][ID] = Spine(ID,sec,seg.x)
                    ID += 1

            elif x >= 80: # peak constant 
                n_spines = spines_per_sec

                for j in range(n_spines):
                    SPINES[sec.name()][ID] = Spine(ID,sec,seg.x)
                    ID += 1
    print("Number of spines added:{}.".format(ID))
    return SPINES
                   
