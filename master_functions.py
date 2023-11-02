'''
functions for model
'''
from   neuron           import h
import pandas as pd
import math as math
import numpy as np
import random 
import os
from tqdm import tqdm
import seaborn as sns
import matplotlib as mpl
import plotly.graph_objects as go
    
def cell_build(cell_type='dspn', specs=None, addSpines=False, branch=False,spine_per_length=1.711):
    model = specs[cell_type]['model']
    morphology = specs[cell_type]['morph']
    if model == 0:
        import MSN_builder0          as build
        if cell_type == 'dspn':
            params='params_dMSN0.json'
        elif cell_type == 'ispn':
            params='params_iMSN0.json'
        cell = build.MSN(params=params, morphology=morphology, variables=None)
        dend_tree = get_root_branches(cell)
        if branch: branch_groups = get_root_groups(cell)
        if addSpines: spines = build.add_spines(cell, spines_per_sec=30)
    elif model == 1:
        import MSN_builder1          as build
        if cell_type == 'dspn':
            params='params_dMSN1.json'
        elif cell_type == 'ispn':
            params='params_iMSN1.json'
        cell = build.MSN(params=params, morphology=morphology, variables=None)
        dend_tree = get_root_branches(cell)
        if branch: branch_groups = get_root_groups(cell)
        if addSpines: spines = build.add_spines(cell=cell, spine_per_length=spine_per_length)
    elif model == 2:
        import MSN_builder2          as build
        if cell_type == 'dspn':
            params='params_dMSN2.json'
        elif cell_type == 'ispn':
            params='params_iMSN2.json'
        cell = build.MSN(params=params, morphology=morphology, variables=None)
        dend_tree = get_root_branches(cell)
        if branch: branch_groups = get_root_groups(cell)
        if addSpines: spines = build.add_spines(params=params, cell=cell, spine_per_length=spine_per_length)
    if addSpines and branch: return(cell, spines, dend_tree, branch_groups)
    elif addSpines and not branch: return(cell, spines, dend_tree)
    elif not addSpines and branch: return(cell, dend_tree, branch_groups)
    elif not addSpines and not branch: return(cell, dend_tree)

def params_selector(cell_type, specs):
    model = specs[cell_type]['model']
    morphology = specs[cell_type]['morph']
    if cell_type == 'dspn':        
        if model == 0:
            params='params_dMSN0.json'
        elif model == 1:
            params='params_dMSN1.json'
        elif model == 2:
            params='params_dMSN2.json'
            
    if cell_type == 'ispn':        
        if model == 0:
            params='params_iMSN0.json'
        elif model == 1:
            params='params_iMSN1.json'
        elif model == 2:
            params='params_iMSN2.json'
    return(params)

def spines_per_dend(cell, spines):
    for sec in cell.dendlist:
        print(sec.name(), len(spines[sec.name()])) 
        
# function gives the distance of EVERY segment within a given dendrite from soma 
def seg_dist(cell, dend):
    dist = []
    for sec in cell.dendlist:
        if sec.name() == dend:
            for i,seg in enumerate(sec):
                dist.append(h.distance(seg))
    return(dist)

# get idxs for spines with UNIQUE locations on a particular dendrite 
# then use to find a unique spines in a given dendrite
def spine_idx(cell, spines, dend):
    for sec in cell.dendlist:
        if sec.name() == dend:
            Nseg = sec.nseg
    spine_locs = (2*np.linspace(1, Nseg, Nseg)-1)/Nseg/2
#     spine_loc = spine_locs[0]
#     spine_loc
    # Get possible spines from section
    candidate_spines = []
    sec_spines = list(spines[dend].items())

    for spine_i, spine_obj in sec_spines: 
        candidate_spines.append(spine_obj)

    # len(sec_spines)
    candidate_spines_locs = []
    for spine in candidate_spines:
        candidate_spines_locs.append(spine.x)
    # candidate_spines_locs
    # spine_idxs = []
    output = []
    for ii in range(Nseg):
        spine_loc = spine_locs[ii]
        a = abs(candidate_spines_locs - spine_loc)
        idx = np.argmin(a)
        # spine_idxs.append(idx)  
        output.append(candidate_spines[idx])
        # only return unique spines
        output = list(dict.fromkeys(output))
    return(output)

def calculate_dist(d3, dist, a4, a5,  a6,  a7, g8):
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

# function to alter a particular conductance g_name in all spines 
def spine_alter(cell, spines, g_name, d3, a4, a5, a6, a7, g8, cell_type='dspn', model=1):

    if model == 2:
        if cell_type == 'dspn':
            par ='params_dMSN2.json'
        elif cell_type == 'dspn':
            par ='params_iMSN2.json'
    if g_name == 'kir':
        for sec in cell.dendlist:
            sec_spines = list(spines[sec.name()].items())
            for spine_i, spine_obj in sec_spines: 
                spine_obj.head.gbar_kir = g8
                if model not in [0,1]:
                    spine_obj.neck.gbar_kir = g8

    if g_name == 'cav32':
        if model == 2:
            for sec in cell.dendlist:
                sec_spines = list(spines[sec.name()].items())
                for spine_i, spine_obj in sec_spines: 
                    dist = h.distance(sec(spine_obj.x))
                    spine_obj.head.pbar_cav32 = calculate_dist(d3=d3, dist=dist, a4=a4, a5=a5, a6=a6, a7=a7, g8=g8)
                    spine_obj.neck.pbar_cav32 = calculate_dist(d3=d3, dist=dist, a4=a4, a5=a5, a6=a6, a7=a7, g8=g8)
        else:
            for sec in cell.dendlist:
                sec_spines = list(spines[sec.name()].items())
                for spine_i, spine_obj in sec_spines: 
                    spine_obj.head.pbar_cav32 = g8
                    spine_obj.neck.pbar_cav32 = g8


    if g_name == 'cav33':
        if model == 2:
            for sec in cell.dendlist:
                sec_spines = list(spines[sec.name()].items())
                for spine_i, spine_obj in sec_spines: 
                    dist = h.distance(sec(spine_obj.x))
                    spine_obj.head.pbar_cav33 = calculate_dist(1, dist, 0, 1.0, 100.0, -30.0, g8)
                    spine_obj.neck.pbar_cav33 = calculate_dist(1, dist, 0, 1.0, 100.0, -30.0, g8)
        else:
            for sec in cell.dendlist:
                sec_spines = list(spines[sec.name()].items())
                for spine_i, spine_obj in sec_spines: 
                    spine_obj.head.pbar_cav33 = g8
                    spine_obj.neck.pbar_cav33 = g8

    if g_name == 'car':
        for sec in cell.dendlist:
            sec_spines = list(spines[sec.name()].items())
            for spine_i, spine_obj in sec_spines: 
                spine_obj.head.pbar_car = g8

    if g_name == 'cal12':
        for sec in cell.dendlist:
            sec_spines = list(spines[sec.name()].items())
            for spine_i, spine_obj in sec_spines: 
                spine_obj.head.pbar_cal12 = g8


    if g_name == 'cal13':
        for sec in cell.dendlist:
            sec_spines = list(spines[sec.name()].items())
            for spine_i, spine_obj in sec_spines: 
                spine_obj.head.pbar_cal13 = g8


    if g_name == 'sk':
        for sec in cell.dendlist:
            sec_spines = list(spines[sec.name()].items())
            for spine_i, spine_obj in sec_spines: 
                spine_obj.head.gbar_sk = g8

    if g_name == 'bk':
        if model == 2:
            for sec in cell.dendlist:
                sec_spines = list(spines[sec.name()].items())
                for spine_i, spine_obj in sec_spines: 
                    spine_obj.head.gbar_bk = g8
                    

# finds dendrites with at least 3 spines
def dend_spine_selector(cell, spines, branch_groups):
    dends_with_spines = []
    # Make list of dendrite sections with at least 2 spines 
    for dend in cell.dendlist:
        sec_spines = list(spines[dend.name()].items())
        for group in branch_groups: # for each nrn dendrite sec, one plot per branch
            if dend in group:
                if len(group) > 2:
                    if len(sec_spines) > 2:
                        dends_with_spines.append(dend)
    return dends_with_spines

# function to alter a given conductance in any location (other than spine) in cell                
def conductance_alter(cell, g_name, d3, a4, a5, a6, a7, g8, g8_somatic_scale, g8_axonal_scale):
    if g_name in ['naf', 'kaf', 'kas', 'kdr', 'kir', 'sk', 'bk']:
        gbar = 'gbar_{}'.format(g_name)
    else:
        gbar = 'pbar_{}'.format(g_name)        
    cell.distribute_channels('dend', gbar, d3, a4, a5, a6, a7, g8) 
    if g_name in ['naf', 'nap', 'kaf', 'kas', 'kdr', 'kir', 'cal12', 'cal13', 'can', 'car', 'sk', 'bk']:
        cell.distribute_channels('soma', gbar, 0, 1, 0, 0, 0, g8*g8_somatic_scale)    
    if g_name == 'naf':
        cell.distribute_channels('axon', gbar, 3, 1, 1.1, 30, 500, g8*g8_axonal_scale)
    if g_name == 'kas':
        cell.distribute_channels('axon', gbar, 0, 1, 0, 0, 0, g8*g8_axonal_scale)   

# this function will alter the relevant conductances in all sections and if present in spine heads and necks
def g_alter(cell, spines, g_name, g8, specs, cell_type='dspn'):
    g8_orig=g8
    model = specs[cell_type]['model'] 
    params =  params_selector(cell_type, specs)          
    # run once:
    [d3, a4, a5, a6, a7, g8, g8_somatic_scale, g8_axonal_scale] = scaling_factors(g_name=g_name, params = params) # maintains original ratios between dend, soma (and axon)
    [d3, a4, a5, a6, a7, g8, g8_somatic_scale, g8_axonal_scale]
    g8 = g8_orig
    conductance_alter(cell=cell, g_name=g_name, d3=d3, a4=a4, a5=a5, a6=a6, a7=a7, g8=g8, g8_somatic_scale=g8_somatic_scale, g8_axonal_scale=g8_axonal_scale)
    spine_alter(cell=cell, spines=spines, g_name=g_name, d3=d3, a4=a4, a5=a5, a6=a6, a7=a7, g8=g8, cell_type=cell_type, model=model)

# rectification is False then ohmic else rectification Pavlov 
def tonic_gaba(cell, gaba_reversal, gbar_gaba, d3=0, a4=1, a5=0, a6=0, a7=0, rectification=False):
    if rectification:
        for sec in cell.dendlist:
            sec.e_gaba2 = gaba_reversal
        for sec in cell.somalist:
            sec.e_gaba2 = gaba_reversal
        cell.distribute_channels('dend', 'gbar_gaba2', d3, a4, a5, a6, a7, gbar_gaba)
        
        g_name = 'gaba2'
        g = []
        gbar = 'gbar_{}'.format(g_name)  
        for sec in cell.dendlist:
            g.append((eval('sec.{}'.format(gbar))))        
        
        if g[0] < g[-1]:
            cell.distribute_channels('soma', 'gbar_gaba2', 0, 1, 0, 0, 0, g[0])
        else:
            cell.distribute_channels('soma', 'gbar_gaba2', 0, 1, 0, 0, 0, gbar_gaba)

    else:
        for sec in cell.dendlist:
            sec.e_gaba1 = gaba_reversal
        for sec in cell.somalist:
            sec.e_gaba1 = gaba_reversal
        cell.distribute_channels('dend', 'gbar_gaba1', d3, a4, a5, a6, a7, gbar_gaba)
        cell.distribute_channels('soma', 'gbar_gaba1', 0, 1, 0, 0, 0, gbar_gaba)

        g_name = 'gaba1'
        g = []
        gbar = 'gbar_{}'.format(g_name)  
        for sec in cell.dendlist:
            g.append((eval('sec.{}'.format(gbar))))        
        
        if g[0] < g[-1]:
            cell.distribute_channels('soma', 'gbar_gaba1', 0, 1, 0, 0, 0, g[0])
        else:
            cell.distribute_channels('soma', 'gbar_gaba1', 0, 1, 0, 0, 0, gbar_gaba)

    
# Get dendrite branches, list for each unique branch structure (TODO: there's probably a neuron func for this)
def get_children(dend, branch_list):
    branch_list.append(dend)
    branches = []

    for child in dend.children():
        branch_list_cpy = branch_list.copy()
        branches.append(get_children(child, branch_list_cpy))

    if len(branches) == 0:
        return branch_list
    else:
        return branches
    

# Parser helper func
def branch_parser_helper(tree):
    for branch in tree:
        if all(type(b) == list for b in branch):
            # need to keep parsing
            for b in branch:
                tree.append(b)
            tree.remove(branch)
        # done parsing branch
    branch_parser(tree)
        
# Parses children into list format
def branch_parser(tree):
    for branch in tree:
        if all(type(b) == list for b in branch):
            branch_parser_helper(tree)
    return

# Takes nrn cell and int for origin dendrite segment index that branches occur from
# Returns parsed list with each entry a list of each unique branch path from origin dendrite segment to termination
def get_dend_branches_from(cell, origin):
    i = 0 
    for dend in cell.dendlist:
        if i == origin: # origin dendrite number to get branches from
            dend_tree = []
            dend_tree = get_children(dend, dend_tree)
            branch_parser(dend_tree)
            return dend_tree
        i += 1
        
def get_root_branches(cell):
    sref_soma = h.SectionRef(sec=cell.soma)

    # Get sec roots (excluding axon)
    roots = []

    for child in sref_soma.child:
        roots.append(child)

    roots = roots[1:]

    # Get dend tree from all roots
    root_tree = []

    for root in roots:
        dend_branch = []
        branch = get_children(root, dend_branch)
        branch_parser(branch)
        root_tree.append(branch)
    
    return root_tree

# gets path from dend to soma
def path_finder2(cell, dend_tree, dend):
    dend_tree2 = [num for sublist in dend_tree for num in sublist]
    for XX in dend_tree2:
        if not isinstance(XX, list):
            XX = [XX]
        for XXX in XX:
            if XXX == dend:
                return XX

def include_upto(iterable, value):
    for it in iterable:
        yield it
        if it == value:
            return

def path_finder(cell, dend_tree, dend):               
    pathlist = []
    pathlist = path_finder2(cell=cell, dend_tree=dend_tree, dend=dend)
    pathlist =  [cell.soma] + pathlist
    return list(include_upto(pathlist, dend))      
            
# Takes cell
# Return the dendrites that are in a branch (with root dendrite first followed by children in that branch ordered by first instance in tree)
# Useful if you want to do something to all dendrites in a branch, without the ordered duplication of root branches 
def get_root_groups(cell):
    root_tree = get_root_branches(cell)
    branch_groups = []
    for branch in root_tree:
        dend_list = []

        for dend in branch:
            
            if isinstance(dend, list):
                for d in dend:
                    if d not in dend_list:
                        dend_list.append(d)
            else:
                if dend not in dend_list:
                    dend_list.append(dend)
                    
        branch_groups.append(dend_list)
    return branch_groups

def nsegs(cell):
    nsegs =[]
    for sec in cell.dendlist:
        nsegs.append(sec.nseg)
    N = sum(nsegs)
    return(N)

def extract(d):
    lists = sorted(d.items()) # sorted by key, return a list of tuples
    x, y = zip(*lists) # unpack a list of pairs into two tuples
    return y

def extract2(d):
    out = []
    for x in d:
        out.append(x)
    return out

def list2df(lst):
    df = pd.DataFrame() 
    df['time'] = extract2(lst[0])
    df['pas'] = extract2(lst[1])
    df['kdr'] = extract2(lst[2])
    df['naf'] = extract2(lst[3])
    df['kaf'] = extract2(lst[4])
    df['kas'] = extract2(lst[5])
    df['kir'] = extract2(lst[6])
    df['cal12'] = extract2(lst[7])
    df['cal13'] = extract2(lst[8])
    df['can'] = extract2(lst[9])
    df['car'] = extract2(lst[10])
    df['cav32'] = extract2(lst[11])
    df['cav33'] = extract2(lst[12])
    #     df['kcnq'] = extract2(lst[xx])
    df['sk'] = extract2(lst[13])
    df['bk'] = extract2(lst[14])  
    return df

def plot_mech(d, mech_name):
    lists = sorted(d.items()) # sorted by key, return a list of tuples
    x, y = zip(*lists) # unpack a list of pairs into two tuples
    plt.title(mech_name)
    plt.plot(x, y)
    plt.show()
    
# return all dendritic inserted mechanisms
def mechanisms(cell):
    d_ = {}
    df = pd.DataFrame()
    # mechs = ['kdr', 'naf', 'kaf', 'kas', 'kdr', 'kir', 'cal12', 'cal13', 'can', 'car', 'cav32', 'cav33', 'sk', 'bk']

    for sec in cell.dendlist:
         d_[h.distance(sec(0.5))] = sec.gbar_kdr
    lists = sorted(d_.items()) # sorted by key, return a list of tuples
    x, y = zip(*lists) # unpack a list of pairs into two tuples

    df['dist'] = x
    df['kdr'] = y

    for sec in cell.dendlist:
         d_[h.distance(sec(0.5))] = sec.gbar_naf
    df['naf'] = extract(d_)

    for sec in cell.dendlist:
         d_[h.distance(sec(0.5))] = sec.gbar_kaf
    df['kaf'] = extract(d_)

    for sec in cell.dendlist:
         d_[h.distance(sec(0.5))] = sec.gbar_kas
    df['kas'] = extract(d_)


    for sec in cell.dendlist:
         d_[h.distance(sec(0.5))] = sec.gbar_kdr
    df['kdr'] = extract(d_)


    for sec in cell.dendlist:
         d_[h.distance(sec(0.5))] = sec.gbar_kir
    df['kir'] = extract(d_)


    for sec in cell.dendlist:
         d_[h.distance(sec(0.5))] = sec.pbar_cal12
    df['cal12'] = extract(d_)


    for sec in cell.dendlist:
         d_[h.distance(sec(0.5))] = sec.pbar_cal13
    df['cal13'] = extract(d_)


    for sec in cell.dendlist:
         d_[h.distance(sec(0.5))] = sec.pbar_can
    df['can'] = extract(d_)


    for sec in cell.dendlist:
         d_[h.distance(sec(0.5))] = sec.pbar_car
    df['car'] = extract(d_)


    for sec in cell.dendlist:
         d_[h.distance(sec(0.5))] = sec.pbar_cav32
    df['cav32'] = extract(d_)


    for sec in cell.dendlist:
         d_[h.distance(sec(0.5))] = sec.pbar_cav33
    df['cav33'] = extract(d_)


    # for sec in cell.dendlist:
    #      d_[h.distance(sec(0.5))] = sec.gbar_kcnq
    # mechanisms.append(extract(d_))
    for sec in cell.dendlist:      
         d_[h.distance(sec(0.5))] = sec.gbar_sk
    df['sk'] = extract(d_)


    for sec in cell.dendlist:
         d_[h.distance(sec(0.5))] = sec.gbar_bk
    df['bk'] = extract(d_)

    return(df)


def scaling_factors(g_name, params='params_dMSN.json'):
    # Parameters:
    # d3   = distribution type:
    #      0 linear, 
    #      1 sigmoidal, 
    #      2 exponential
    #      3 step function
    # dist = somatic distance of segment
    # a4-7 = distribution parameters 
    # g8   = base conductance (similar to maximal conductance)

    import json
    with open(params) as file:
        par = json.load(file)

    # cell, SPINES, branch_groups, dend_tree = utils.make_cell()
    g8_axonal = 0
    g8_somatic = 0
    if g_name == 'naf':
        d3 = 1
        a4 = 0  # ALT
        a5 = 1  # SPREAD
        a6 = 50
        a7 = 10
        g8_basal = float(par['gbar_{}_basal'.format(g_name)]['Value'])
        g8_somatic = float(par['gbar_{}_somatic'.format(g_name)]['Value'])
        g8_axonal = float(par['gbar_{}_axonal'.format(g_name)]['Value'])

    if g_name == 'kaf':
        d3 = 1
        a4 = 0.5   # ALT
        a5 = 0.25  # SPREAD
        a6 = 120
        a7 = 30
        g8_basal = float(par['gbar_{}_basal'.format(g_name)]['Value'])
        g8_somatic = float(par['gbar_{}_somatic'.format(g_name)]['Value'])

    if g_name == 'kas':
        d3 = 2
        a4 = 0.25  # ALT
        a5 = 5     # SPREAD
        a6 = 0
        a7 = -10
        g8_basal = float(par['gbar_{}_basal'.format(g_name)]['Value'])
        g8_somatic = float(par['gbar_{}_somatic'.format(g_name)]['Value'])
        g8_axonal = float(par['gbar_{}_axonal'.format(g_name)]['Value'])

    if g_name == 'kdr':
        d3 = 1
        a4 = 0.25  # ALT
        a5 = 1     # SPREAD
        a6 = 50
        a7 = 30
        g8_basal = float(par['gbar_{}_basal'.format(g_name)]['Value'])
        g8_somatic = float(par['gbar_{}_somatic'.format(g_name)]['Value']) # build.calculate_distribution(d3, -cell.soma.diam/2, a4, a5,  a6,  a7, g8_basal)

    if g_name == 'kir':
        d3 = 0
        a4 = 1  # ALT
        a5 = 0  # SPREAD
        a6 = 0
        a7 = 0
        if params == 'params_dMSN.json':
            g8_basal = float(par['gbar_{}_basal'.format(g_name)]['Value'])*2
            g8_somatic = float(par['gbar_{}_somatic'.format(g_name)]['Value'])*2
        else:
            g8_basal = float(par['gbar_{}_basal'.format(g_name)]['Value'])
            g8_somatic = float(par['gbar_{}_somatic'.format(g_name)]['Value'])

    if g_name == 'sk':
        d3 = 0
        a4 = 1  # ALT
        a5 = 0  # SPREAD
        a6 = 0
        a7 = 0
        g8_basal = float(par['gbar_{}_basal'.format(g_name)]['Value'])
        g8_somatic = float(par['gbar_{}_somatic'.format(g_name)]['Value'])

    if g_name == 'bk':
        d3 = 0
        a4 = 1  # ALT
        a5 = 0  # SPREAD
        a6 = 0
        a7 = 0
        g8_basal = float(par['gbar_{}_basal'.format(g_name)]['Value'])
        g8_somatic = float(par['gbar_{}_somatic'.format(g_name)]['Value'])

    if g_name == 'car':
        d3 = 0
        a4 = 1  # ALT
        a5 = 0  # SPREAD
        a6 = 0
        a7 = 0
        g8_basal = float(par['pbar_{}_basal'.format(g_name)]['Value'])
        g8_somatic = float(par['pbar_{}_somatic'.format(g_name)]['Value'])

    if g_name == 'can':
        d3 = 0
        a4 = 1  # ALT
        a5 = 0  # SPREAD
        a6 = 0
        a7 = 0
        g8_basal = float(par['pbar_{}_basal'.format(g_name)]['Value'])
        g8_somatic = float(par['pbar_{}_somatic'.format(g_name)]['Value'])
        

    if g_name == 'cal12':
        d3 = 0
        a4 = 1  # ALT
        a5 = 0  # SPREAD
        a6 = 0
        a7 = 0
        g8_basal = float(par['pbar_{}_basal'.format(g_name)]['Value'])
        g8_somatic = float(par['pbar_{}_somatic'.format(g_name)]['Value'])
        
    if g_name == 'cal13':
        d3 = 0
        a4 = 1  # ALT
        a5 = 0  # SPREAD
        a6 = 0
        a7 = 0
        g8_basal = float(par['pbar_{}_basal'.format(g_name)]['Value'])
        g8_somatic = float(par['pbar_{}_somatic'.format(g_name)]['Value'])

    if g_name == 'cav32':
        d3 = 1
        a4 = 0  # ALT
        a5 = 1  # SPREAD
        a6 = 100
        a7 = -30
        g8_basal = float(par['pbar_{}_basal'.format(g_name)]['Value'])
 
    if g_name == 'cav33':
        d3 = 1
        a4 = 0  # ALT
        a5 = 1  # SPREAD
        a6 = 100
        a7 = -30
        g8_basal = float(par['pbar_{}_basal'.format(g_name)]['Value'])
 
    g8_somatic_scale = g8_somatic / g8_basal
    g8_axonal_scale = g8_axonal / g8_basal
    
    out = [d3, a4, a5, a6, a7, g8_basal, g8_somatic_scale, g8_axonal_scale]
    return(out)

# Set up branch assignment
def branch_selection(cell, cell_type='dpsn'):
    branch1_dends = [None] * 2 
    branch2_dends = [None] * 2 
    branch3_dends = [None] * 2 
    branch4_dends = [None] * 2 
    branch5_dends = [None] * 2 

    # Define dendrite sections for each predefined branch (chosen as sterotypical primary dendrites)
    if cell_type == 'dspn':    
        for dend in cell.dendlist:
            if dend.name() == 'dend[28]':
                branch1_dends[-1] = dend
            if dend.name() == 'dend[25]':
                branch1_dends[0] = dend

            if dend.name() == 'dend[15]':
                branch2_dends[-1] = dend
            if dend.name() == 'dend[13]':
                branch2_dends[0] = dend

            if dend.name() == 'dend[46]':
                branch3_dends[-1] = dend
            if dend.name() == 'dend[43]':
                branch3_dends[0] = dend

            if dend.name() == 'dend[36]':
                branch4_dends[-1] = dend
            if dend.name() == 'dend[32]':
                branch4_dends[0] = dend

            if dend.name() == 'dend[57]':
                branch5_dends[-1] = dend
            if dend.name() == 'dend[55]':
                branch5_dends[0] = dend

    elif cell_type == 'ispn':
        for dend in cell.dendlist:
            if dend.name() == 'dend[29]':
                branch1_dends[-1] = dend
            if dend.name() == 'dend[27]':
                branch1_dends[0] = dend

            if dend.name() == 'dend[15]':
                branch2_dends[-1] = dend
            if dend.name() == 'dend[13]':
                branch2_dends[0] = dend

            if dend.name() == 'dend[17]':
                branch3_dends[-1] = dend
            if dend.name() == 'dend[12]':
                branch3_dends[0] = dend

            if dend.name() == 'dend[45]':
                branch4_dends[-1] = dend
            if dend.name() == 'dend[41]':
                branch4_dends[0] = dend

            if dend.name() == 'dend[36]':
                branch5_dends[-1] = dend
            if dend.name() == 'dend[32]':
                branch5_dends[0] = dend
    
    # For sparse plotting
    return [branch1_dends] + [branch2_dends] + [branch3_dends] + [branch4_dends] + [branch5_dends]

# change all spine neck diameters
def spine_neck_diameter(cell, spines, diam):
    for sec in cell.dendlist:
        sec_spines = list(spines[sec.name()].items())
        for spine_i, spine_obj in sec_spines: 
            spine_obj.neck.diam = diam

def spine_neck_length(cell, spines, length):
    for sec in cell.dendlist:
        sec_spines = list(spines[sec.name()].items())
        for spine_i, spine_obj in sec_spines: 
            spine_obj.neck.L = length
            
def spine_head_diameter(cell, spines, diam, length):
    for sec in cell.dendlist:
        sec_spines = list(spines[sec.name()].items())
        for spine_i, spine_obj in sec_spines: 
            spine_obj.head.diam = diam
            spine_obj.head.L = length
            
# Set up branch assignment and add glutamate
def glut_add(cell=None,
               branch1_glut = False, 
               branch2_glut = True, 
               branch3_glut = False, 
               branch4_glut = False, 
               branch5_glut = False, 
               num_gluts = 15,
               glut_placement = 'distal',
               glut = True,
               cell_type='dspn'):
    [branch1_dends, branch2_dends, branch3_dends, branch4_dends, branch5_dends] = branch_selection(cell, cell_type=cell_type) 
    glut_secs = []
    glut_secs_orig = []
    # Define placement on dendritic branch (prox/dist)
    if 'proximal' in glut_placement:
        glut_site = 0
    else:
        glut_site = -1

    # Define branch for glutamate (multiple possible)
    if branch1_glut:
        glut_secs.append(branch1_dends[glut_site])
        glut_secs_orig.append(branch1_dends[glut_site])

    if branch2_glut:
        glut_secs.append(branch2_dends[glut_site])
        glut_secs_orig.append(branch2_dends[glut_site])

    if branch3_glut:
        glut_secs.append(branch3_dends[glut_site])
        glut_secs_orig.append(branch3_dends[glut_site])

    if branch4_glut:
        glut_secs.append(branch4_dends[glut_site])
        glut_secs_orig.append(branch4_dends[glut_site])

    if branch5_glut:
        glut_secs.append(branch5_dends[glut_site])
        glut_secs_orig.append(branch5_dends[glut_site])

    # Number of glutamatergic inputs per section is num_gluts
    glut_secs *= num_gluts     

    if not glut:
        glut_secs = []
    return glut_secs, glut_secs_orig

# def glut_place(spines,
#                method=0, 
#                physiological=True, 
#                AMPA=True, 
#                g_AMPA = 0.001,
#                NMDA=True,
#                ratio = 2,
#                glut_time = 200,
#                glut_secs = None,
#                glut_onsets=None,
#                num_gluts=15,
#                return_currents = True):
#     nmda_currents = [None]*len(glut_secs)
#     ampa_currents = [None]*len(glut_secs)
#     glut_synapses = [0]*len(glut_secs)
#     glut_stimulator = {}
#     glut_connection = {}
#     if len(glut_secs) > 0: 
#         glut_id = 0 # index used for glut_synapse list and printing
#         final_spine_locs = []
#         random.seed(42)
#         for dend_glut in glut_secs:
#             # Get possible spines from section
#             candidate_spines = []
#             sec_spines = list(spines[dend_glut.name()].items())

#             for spine_i, spine_obj in sec_spines: 
#                 candidate_spines.append(spine_obj)

#             if len(glut_secs) < len(sec_spines):
#                 if method==1:
#                     # reversed order so activate along dendrite towards soma
#                     spine_idx = 2*len(candidate_spines)//3-1 # arbitrary start point at 2/3 of spines
#                     spine = candidate_spines[spine_idx - glut_id] 
#                 else:
#                     spine_idx = 2*len(candidate_spines)//3 - num_gluts # arbitrary start point at 1/3 of spines
#                     if spine_idx < 0:
#                         if len(candidate_spines) >= num_gluts:
#                             spine_idx = len(candidate_spines) - num_gluts
#                         else:
#                             spine_idx = 0        
#                     spine = candidate_spines[spine_idx + glut_id] 
#             else:
#                 spine = random.choice(candidate_spines)

#             spine_loc = spine.x
#             spine_head = spine.head
#             final_spine_locs.append(spine_loc) 

#             # Define glutamate syn 
#             glut_synapses[glut_id] = h.glutsynapse(spine_head(0.5))
#             if physiological:
#                 if AMPA:
#                     glut_synapses[glut_id].gmax_AMPA = g_AMPA
#                 else:
#                     glut_synapses[glut_id].gmax_AMPA = 0
#                 if NMDA:
#                     glut_synapses[glut_id].gmax_NMDA = g_AMPA*ratio # 
#                 else:
#                     glut_synapses[glut_id].gmax_NMDA = 0 # NMDA:AMPA ratio is 0.5
#                 # values from Ding et al., 2008; AMPA decay value similar in Kreitzer & Malenka, 2007
#                 glut_synapses[glut_id].tau1_ampa = 0.86 # 10-90% rise 1.9; tau = 1.9/2.197
#                 glut_synapses[glut_id].tau2_ampa = 4.8                
#                 # physiological kinetics for NMDA from Chapman et al. 2003, 
#                 # NMDA decay is weighted average of fast and slow 231 +- 5 ms
#                 # rise time 10-90% is 12.13 ie tau = 12.13 / 2.197 
#                 # values from Kreitzer & Malenka, 2007 are 2.5 and 50 
#                 glut_synapses[glut_id].tau1_nmda = 5.52
#                 glut_synapses[glut_id].tau2_nmda = 231   
#                 # alpha and beta determine neg slope of Mg block for NMDA
#                 glut_synapses[glut_id].alpha = 0.096
#                 glut_synapses[glut_id].beta = 17.85  # ie 5*3.57  
#             else:
#                 glut_synapses[glut_id].gmax_AMPA = 0.001 
#                 glut_synapses[glut_id].gmax_NMDA = 0.007
#                 # physiological kinetics for NMDA from Chapman et al. 2003, 
#                 # NMDA decay is weighted average of fast and slow 231 +- 5 ms
#                 # rise time 10-90% is 12.13 ie tau = 12.13 / 2.197 
#                 glut_synapses[glut_id].tau1_nmda = 5.52
#                 glut_synapses[glut_id].tau2_nmda = 231            

#             # Stim to play back spike times as defined by onsets
#             glut_stimulator[glut_id] = h.VecStim()
#             glut_stimulator[glut_id].play(h.Vector(1, glut_onsets[glut_id]))

#             # Connect stim and syn
#             glut_connection[glut_id] = h.NetCon(glut_stimulator[glut_id], glut_synapses[glut_id])
#             glut_connection[glut_id].weight[0] = 0.35

#             if return_currents:
#                 # Record NMDA current for synapse
#                 nmda_currents[glut_id] = h.Vector().record(glut_synapses[glut_id]._ref_i_nmda)
#                 ampa_currents[glut_id] = h.Vector().record(glut_synapses[glut_id]._ref_i_ampa)

#             glut_id += 1 # Increment glutamate counter

#         print("# glutamate added:{}, on sections:{}, with final spine locs:{}".format(glut_id, glut_secs, final_spine_locs))
#     return glut_synapses, glut_stimulator, glut_connection, ampa_currents, nmda_currents




def glut_place(spines,
               method=0, 
               physiological=True, 
               AMPA=True, 
               g_AMPA = 0.001,
               NMDA=True,
               ratio = 2,
               glut_time = 200,
               glut_secs = None,
               glut_onsets=None,
               num_gluts=15,
               return_currents = True,
               model = 1):
    nmda_currents = [None]*len(glut_secs)
    ampa_currents = [None]*len(glut_secs)
    glut_synapses = [0]*len(glut_secs)
    glut_stimulator = {}
    glut_connection = {}
    if len(glut_secs) > 0: 
        glut_id = 0 # index used for glut_synapse list and printing
        final_spine_locs = []
        random.seed(42)
        for dend_glut in glut_secs:
            # Get possible spines from section
            candidate_spines = []
            sec_spines = list(spines[dend_glut.name()].items())

            if model in [1,2]:
            
                for spine_i, spine_obj in sec_spines: 
                    candidate_spines.append(spine_obj)

                if len(glut_secs) < len(sec_spines):
                    if method==1:
                        # reversed order so activate along dendrite towards soma
                        spine_idx = 2*len(candidate_spines)//3-1 # arbitrary start point at 2/3 of spines
                        spine = candidate_spines[spine_idx - glut_id] 
                    else:
                        spine_idx = 2*len(candidate_spines)//3 - num_gluts # arbitrary start point at 1/3 of spines
                        if spine_idx < 0:
                            if len(candidate_spines) >= num_gluts:
                                spine_idx = len(candidate_spines) - num_gluts
                            else:
                                spine_idx = 0        
                        spine = candidate_spines[spine_idx + glut_id] 
                else:
                    spine = random.choice(candidate_spines)
                    
            else:
            
                for spine_i, spine_obj in sec_spines: 
                    candidate_spines.append(spine_obj)
                if len(glut_secs) < len(sec_spines):
                    spine_idx = len(candidate_spines)//3-1 # arbitrary start point at 1/3 of spines
                    spine = candidate_spines[spine_idx + glut_id] 

                else:
                    spine = random.choice(candidate_spines)


            spine_loc = spine.x
            spine_head = spine.head
            final_spine_locs.append(spine_loc) 

            # Define glutamate syn 
            glut_synapses[glut_id] = h.glutsynapse(spine_head(0.5))
            if physiological:
                if AMPA:
                    glut_synapses[glut_id].gmax_AMPA = g_AMPA
                else:
                    glut_synapses[glut_id].gmax_AMPA = 0
                if NMDA:
                    glut_synapses[glut_id].gmax_NMDA = g_AMPA*ratio # 
                else:
                    glut_synapses[glut_id].gmax_NMDA = 0 # NMDA:AMPA ratio is 0.5
                # values from Ding et al., 2008; AMPA decay value similar in Kreitzer & Malenka, 2007
                glut_synapses[glut_id].tau1_ampa = 0.86 # 10-90% rise 1.9; tau = 1.9/2.197
                glut_synapses[glut_id].tau2_ampa = 4.8                
                # physiological kinetics for NMDA from Chapman et al. 2003, 
                # NMDA decay is weighted average of fast and slow 231 +- 5 ms
                # rise time 10-90% is 12.13 ie tau = 12.13 / 2.197 
                # values from Kreitzer & Malenka, 2007 are 2.5 and 50 
                glut_synapses[glut_id].tau1_nmda = 5.52
                glut_synapses[glut_id].tau2_nmda = 231   
                # alpha and beta determine neg slope of Mg block for NMDA
                glut_synapses[glut_id].alpha = 0.096
                glut_synapses[glut_id].beta = 17.85  # ie 5*3.57  
            else:
                glut_synapses[glut_id].gmax_AMPA = 0.001 
                glut_synapses[glut_id].gmax_NMDA = 0.007
                # physiological kinetics for NMDA from Chapman et al. 2003, 
                # NMDA decay is weighted average of fast and slow 231 +- 5 ms
                # rise time 10-90% is 12.13 ie tau = 12.13 / 2.197 
                glut_synapses[glut_id].tau1_nmda = 5.52
                glut_synapses[glut_id].tau2_nmda = 231            

            # Stim to play back spike times as defined by onsets
            glut_stimulator[glut_id] = h.VecStim()
            glut_stimulator[glut_id].play(h.Vector(1, glut_onsets[glut_id]))

            # Connect stim and syn
            glut_connection[glut_id] = h.NetCon(glut_stimulator[glut_id], glut_synapses[glut_id])
            glut_connection[glut_id].weight[0] = 0.35

            if return_currents:
                # Record NMDA current for synapse
                nmda_currents[glut_id] = h.Vector().record(glut_synapses[glut_id]._ref_i_nmda)
                ampa_currents[glut_id] = h.Vector().record(glut_synapses[glut_id]._ref_i_ampa)

            glut_id += 1 # Increment glutamate counter

        print("# glutamate added:{}, on sections:{}, with final spine locs:{} with timing onsets:{}".format(glut_id, glut_secs, final_spine_locs, glut_onsets))
    return glut_synapses, glut_stimulator, glut_connection, ampa_currents, nmda_currents

# finds distances of syanpses from soma
def synapse_dist(spines,
               method=0,
               glut_secs = None,
               num_gluts=15):
    if len(glut_secs) > 0: 
        glut_id = 0 # index used for glut_synapse list and printing
        final_spine_dists = []
        random.seed(42)
        for dend_glut in glut_secs:
            # Get possible spines from section
            candidate_spines = []
            sec_spines = list(spines[dend_glut.name()].items())

            for spine_i, spine_obj in sec_spines: 
                candidate_spines.append(spine_obj)

            if len(glut_secs) < len(sec_spines):
                if method==1:
                    # reversed order so activate along dendrite towards soma
                    spine_idx = 2*len(candidate_spines)//3-1 # arbitrary start point at 2/3 of spines
                    spine = candidate_spines[spine_idx - glut_id] 
                else:
                    spine_idx = 2*len(candidate_spines)//3 - num_gluts # arbitrary start point at 1/3 of spines
                    if spine_idx < 0:
                        if len(candidate_spines) >= num_gluts:
                            spine_idx = len(candidate_spines) - num_gluts
                        else:
                            spine_idx = 0        
                    spine = candidate_spines[spine_idx + glut_id] 
            else:
                spine = random.choice(candidate_spines)

            spine_loc = spine.x
            final_spine_dists.append(h.distance(dend_glut(spine_loc))) 
            glut_id += 1 # Increment glutamate counter
    return final_spine_dists

def gaba_onset(gaba_time, num_gabas, num_branch2, model=1):
    if model == 0:
        gaba_onsets = list(range(gaba_time, gaba_time + int(num_gabas/3)+1)) * 3 * num_branch2
        gaba_onsets = gaba_onsets[:num_gabas]
    else:
        if num_branch2 in [0,1]:
            if (num_gabas < 4):
                gaba_onsets = list(range(gaba_time, gaba_time + num_gabas)) 
            else:
                if num_gabas % 3 == 0:
                    gaba_onsets = list(range(gaba_time, gaba_time + int(num_gabas/3))) * 3 * num_branch2
                else:
                    gaba_onsets = list(range(gaba_time, gaba_time + int(num_gabas/3)+1)) * 3 * num_branch2
            gaba_onsets = gaba_onsets[:num_gabas]
        else:
            onsets = list(range(gaba_time, gaba_time + num_gabas)) 
            gaba_onsets = [x for x in onsets for _ in range(num_branch2)]
    return gaba_onsets


def gaba_add(cell=None,
               gaba=True, 
               branch1_gaba = False, 
               branch2_gaba = True, 
               branch3_gaba = False, 
               branch4_gaba = False, 
               branch5_gaba = False, 
               gaba_placement = 'distal',
               num_gabas=15,
               show=True,
               cell_type='dspn'):

    if gaba > 0: 
        [branch1_dends, branch2_dends, branch3_dends, branch4_dends, branch5_dends] = branch_selection(cell, cell_type) 

        gaba_secs = []

        # Define gaba spatial placement 
        if 'soma' in gaba_placement:
            gaba_secs.append(cell.soma)

            gaba_secs *= num_gabas # need to duplicate sections to place synapses 

        elif 'everywhere' in gaba_placement: # append to every dendrite section
            for dend in cell.dendlist:
                gaba_secs.append(dend)

            gaba_secs *= num_gabas # need to duplicate sections to place synapses 


        elif 'distributed_branch' in gaba_placement: # append to specific branches

            # Define placement on dendritic branch (prox/dist)
            if 'proximal' in gaba_placement:
                gaba_site = 0
            else:
                gaba_site = -1

            # Define branch for gaba (multiple possible)
            if branch1_gaba:
                gaba_secs.append(branch1_dends[gaba_site])

            if branch2_gaba:
                gaba_secs.append(branch2_dends[gaba_site])

            if branch3_gaba:
                gaba_secs.append(branch3_dends[gaba_site])

            if branch4_gaba:
                gaba_secs.append(branch4_dends[gaba_site])

            if branch5_gaba:
                gaba_secs.append(branch5_dends[gaba_site])

            gaba_secs *= num_gabas # need to duplicate sections to place synapses 

    else:
        # No gaba
        gaba_secs = []

    if show:
        print("gaba:{}".format(gaba_secs))
    return gaba_secs

def gaba_place(physiological=True,
               gaba_reversal = -60,
               gaba_weight = 0.001,
               gaba_time = 200,
               gaba_secs = None,
               gaba_onsets=None,
               gaba_locations = None,
               num_gabas=15,
               return_currents = True,
               show=True):
    
    gaba_conductances = [0] * len(gaba_secs)
    gaba_currents = [0] * len(gaba_secs)
    gaba_synapses = [0]*len(gaba_secs) # list of gaba synapses
    gaba_stimulator = {}
    gaba_connection = {}
    if gaba_locations is None:
        gaba_locations = [0.5] * len(gaba_secs)
        
    # Place gabaergic synapses
    if len(gaba_secs) > 0:

        gaba_id = 0 # index used for gaba_synapse list and printing
        gaba_locs = []

        for dend_gaba in gaba_secs:

            # For now, just assign to middle of section instead of uniform random
            gaba_loc = gaba_locations[gaba_id]

            # Choose random location along section
    #                 gaba_loc = round(random.uniform(0, 1), 2)

            gaba_locs.append(gaba_loc)

            # Define gaba synapse
            gaba_synapses[gaba_id] = h.gabasynapse(dend_gaba(gaba_loc)) 
            if physiological:
                gaba_synapses[gaba_id].tau1 = 0.9 
                gaba_synapses[gaba_id].tau2 = 18
            else:
                gaba_synapses[gaba_id].tau2 = 0.9 # TODO: Tune tau2 further for accurate response 
            gaba_synapses[gaba_id].erev = gaba_reversal

            # Stim to play back spike times
            gaba_stimulator[gaba_id] = h.VecStim()

            # Use with deterministic onset times
            gaba_stimulator[gaba_id].play(h.Vector(1, gaba_onsets[gaba_id]))

            # Connect stim and syn
            gaba_connection[gaba_id] = h.NetCon(gaba_stimulator[gaba_id], gaba_synapses[gaba_id])
            gaba_connection[gaba_id].weight[0] = gaba_weight # Depending on desired EPSP response at soma, tune this

            if return_currents:
                # Measure conductance and current
                gaba_currents[gaba_id] = h.Vector().record(gaba_synapses[gaba_id]._ref_i)
                gaba_conductances[gaba_id] = h.Vector().record(gaba_synapses[gaba_id]._ref_g)

            gaba_id += 1 # increment gaba counter

        if show:
            print("# gaba synapses added:{} on:{} with locs:{} with timing onsets:{}".format(gaba_id, gaba_secs, gaba_locs, gaba_onsets))
    return gaba_synapses, gaba_stimulator, gaba_connection, gaba_currents, gaba_conductances

# def glut_place2(cell,
#                spines,
#                method=0, 
#                physiological=True, 
#                AMPA=True, 
#                g_AMPA = 0.001,
#                NMDA=True,
#                ratio = 2,
#                glut_time = 200,
#                glut_secs = None,
#                glut_onsets=None,
#                glut_locs = None,
#                num_gluts=15,
#                return_currents = True):
#     nmda_currents = [None]*len(glut_secs)
#     ampa_currents = [None]*len(glut_secs)
#     glut_synapses = [0]*len(glut_secs)
#     glut_stimulator = {}
#     glut_connection = {}
#     final_spine_locs = []
#     final_spines = []
#     if len(glut_secs) > 0: 
#         glut_id = 0 # index used for glut_synapse list and printing

#         for ii in range(len(glut_secs)):
#             synapse_loc = glut_locs[ii]
#             # Get candidate spines from section
#             candidates = []
#             sec_spines = list(spines[glut_secs[ii].name()].items())

#             for spine_i, spine_obj in sec_spines: 
#                 candidates.append(spine_obj)

#             locs = []
#             for spine in candidates:
#                 locs.append(spine.x)
#             loc, idx = find_closest_value(locs, synapse_loc)
#             spine = candidates[idx] # choose last spine

#             spine_loc = spine.x
#             spine_head = spine.head
#             final_spine_locs.append(spine_loc) 
#             final_spines.append(spine)
#             # Define glutamate syn 
#             glut_synapses[glut_id] = h.glutsynapse(spine_head(0.5))
#             if physiological:
#                 if AMPA:
#                     glut_synapses[glut_id].gmax_AMPA = g_AMPA
#                 else:
#                     glut_synapses[glut_id].gmax_AMPA = 0
#                 if NMDA:
#                     glut_synapses[glut_id].gmax_NMDA = g_AMPA*ratio # 
#                 else:
#                     glut_synapses[glut_id].gmax_NMDA = 0 # NMDA:AMPA ratio is 0.5
#                 # values from Ding et al., 2008; AMPA decay value similar in Kreitzer & Malenka, 2007
#                 glut_synapses[glut_id].tau1_ampa = 0.86 # 10-90% rise 1.9; tau = 1.9/2.197
#                 glut_synapses[glut_id].tau2_ampa = 4.8                
#                 # physiological kinetics for NMDA from Chapman et al. 2003, 
#                 # NMDA decay is weighted average of fast and slow 231 +- 5 ms
#                 # rise time 10-90% is 12.13 ie tau = 12.13 / 2.197 
#                 # values from Kreitzer & Malenka, 2007 are 2.5 and 50 
#                 glut_synapses[glut_id].tau1_nmda = 5.52
#                 glut_synapses[glut_id].tau2_nmda = 231   
#                 # alpha and beta determine neg slope of Mg block for NMDA
#                 glut_synapses[glut_id].alpha = 0.096
#                 glut_synapses[glut_id].beta = 17.85  # ie 5*3.57  
#             else:
#                 glut_synapses[glut_id].gmax_AMPA = 0.001 
#                 glut_synapses[glut_id].gmax_NMDA = 0.007
#                 # physiological kinetics for NMDA from Chapman et al. 2003, 
#                 # NMDA decay is weighted average of fast and slow 231 +- 5 ms
#                 # rise time 10-90% is 12.13 ie tau = 12.13 / 2.197 
#                 glut_synapses[glut_id].tau1_nmda = 5.52
#                 glut_synapses[glut_id].tau2_nmda = 231            

#             # Stim to play back spike times as defined by onsets
#             glut_stimulator[glut_id] = h.VecStim()
#             glut_stimulator[glut_id].play(h.Vector(1, glut_onsets[glut_id]))

#             # Connect stim and syn
#             glut_connection[glut_id] = h.NetCon(glut_stimulator[glut_id], glut_synapses[glut_id])
#             glut_connection[glut_id].weight[0] = 0.35

#             if return_currents:
#                 # Record NMDA current for synapse
#                 nmda_currents[glut_id] = h.Vector().record(glut_synapses[glut_id]._ref_i_nmda)
#                 ampa_currents[glut_id] = h.Vector().record(glut_synapses[glut_id]._ref_i_ampa)

#             glut_id += 1 # Increment glutamate counter

# #         rounded_locs = [round(value, 4) for value in glut_locs]
#         rounded_locs = [round(value, 4) for value in final_spine_locs]
#         print("# glutamate added:{}, on sections:{}, with final spine locs:{} with timing onsets:{}".format(glut_id, glut_secs, rounded_locs, glut_onsets))
#     return glut_synapses, glut_stimulator, glut_connection, ampa_currents, nmda_currents, final_spines, final_spine_locs

def glut_place2(cell,
               spines,
               method=0, 
               physiological=True, 
               AMPA=True, 
               g_AMPA = 0.001,
               NMDA=True,
               ratio = 2,
               glut=True,
               glut_time = 200,
               glut_secs = None,
               glut_onsets=None,
               glut_locs = None,
               num_gluts=15,
               return_currents = True):
    nmda_currents = [None]*len(glut_secs)
    ampa_currents = [None]*len(glut_secs)
    glut_synapses = [0]*len(glut_secs)
    glut_stimulator = {}
    glut_connection = {}
    final_spine_locs = []
    final_spines = []
    if num_gluts > 0: 
        glut_id = 0 # index used for glut_synapse list and printing

        for ii in list(range(0,num_gluts)):
            synapse_loc = glut_locs[ii]
            # Get candidate spines from section
            candidates = []
            sec_spines = list(spines[glut_secs[ii].name()].items())

            for spine_i, spine_obj in sec_spines: 
                candidates.append(spine_obj)

            locs = []
            for spine in candidates:
                locs.append(spine.x)
            loc, idx = find_closest_value(locs, synapse_loc)
            spine = candidates[idx] # choose last spine

            spine_loc = spine.x
            spine_head = spine.head
            final_spine_locs.append(spine_loc) 
            final_spines.append(spine)
            if glut:
                # Define glutamate syn 
                glut_synapses[glut_id] = h.glutsynapse(spine_head(0.5))
                if physiological:
                    if AMPA:
                        glut_synapses[glut_id].gmax_AMPA = g_AMPA
                    else:
                        glut_synapses[glut_id].gmax_AMPA = 0
                    if NMDA:
                        glut_synapses[glut_id].gmax_NMDA = g_AMPA*ratio # 
                    else:
                        glut_synapses[glut_id].gmax_NMDA = 0 # NMDA:AMPA ratio is 0.5
                    # values from Ding et al., 2008; AMPA decay value similar in Kreitzer & Malenka, 2007
                    glut_synapses[glut_id].tau1_ampa = 0.86 # 10-90% rise 1.9; tau = 1.9/2.197
                    glut_synapses[glut_id].tau2_ampa = 4.8                
                    # physiological kinetics for NMDA from Chapman et al. 2003, 
                    # NMDA decay is weighted average of fast and slow 231 +- 5 ms
                    # rise time 10-90% is 12.13 ie tau = 12.13 / 2.197 
                    # values from Kreitzer & Malenka, 2007 are 2.5 and 50 
                    glut_synapses[glut_id].tau1_nmda = 5.52
                    glut_synapses[glut_id].tau2_nmda = 231   
                    # alpha and beta determine neg slope of Mg block for NMDA
                    glut_synapses[glut_id].alpha = 0.096
                    glut_synapses[glut_id].beta = 17.85  # ie 5*3.57  
                else:
                    glut_synapses[glut_id].gmax_AMPA = 0.001 
                    glut_synapses[glut_id].gmax_NMDA = 0.007
                    # physiological kinetics for NMDA from Chapman et al. 2003, 
                    # NMDA decay is weighted average of fast and slow 231 +- 5 ms
                    # rise time 10-90% is 12.13 ie tau = 12.13 / 2.197 
                    glut_synapses[glut_id].tau1_nmda = 5.52
                    glut_synapses[glut_id].tau2_nmda = 231            

                # Stim to play back spike times as defined by onsets
                glut_stimulator[glut_id] = h.VecStim()
                glut_stimulator[glut_id].play(h.Vector(1, glut_onsets[glut_id]))

                # Connect stim and syn
                glut_connection[glut_id] = h.NetCon(glut_stimulator[glut_id], glut_synapses[glut_id])
                glut_connection[glut_id].weight[0] = 0.35

                if return_currents:
                    # Record NMDA current for synapse
                    nmda_currents[glut_id] = h.Vector().record(glut_synapses[glut_id]._ref_i_nmda)
                    ampa_currents[glut_id] = h.Vector().record(glut_synapses[glut_id]._ref_i_ampa)

            glut_id += 1 # Increment glutamate counter

#         rounded_locs = [round(value, 4) for value in glut_locs]
        rounded_locs = [round(value, 4) for value in final_spine_locs]
        if glut:
            print("# glutamate added:{}, on sections:{}, with final spine locs:{} with timing onsets:{}".format(glut_id, glut_secs, rounded_locs, glut_onsets))
    return glut_synapses, glut_stimulator, glut_connection, ampa_currents, nmda_currents, final_spines, final_spine_locs


def gaba_place2(physiological=True,
               gaba_reversal = -60,
               gaba_weight = 0.001,
               gaba_time = 200,
               gaba_secs = None,
               gaba_onsets=None,
               gaba_locations = None,
               num_gabas=15,
               return_currents = True,
               show=True):
    
    gaba_conductances = [0] * len(gaba_secs)
    gaba_currents = [0] * len(gaba_secs)
    gaba_synapses = [0]*len(gaba_secs) # list of gaba synapses
    gaba_stimulator = {}
    gaba_connection = {}
    if gaba_locations is None:
        gaba_locations = [0.5] * len(gaba_secs)
    gaba_locs = []    
    # Place gabaergic synapses
    if len(gaba_secs) > 0:

        gaba_id = 0 # index used for gaba_synapse list and printing
        for dend_gaba in gaba_secs:

            # For now, just assign to middle of section instead of uniform random
            gaba_loc = gaba_locations[gaba_id]

            # Choose random location along section
    #                 gaba_loc = round(random.uniform(0, 1), 2)

            gaba_locs.append(gaba_loc)

            # Define gaba synapse
            gaba_synapses[gaba_id] = h.gabasynapse(dend_gaba(gaba_loc)) 
            if physiological:
                gaba_synapses[gaba_id].tau1 = 0.9 
                gaba_synapses[gaba_id].tau2 = 18
            else:
                gaba_synapses[gaba_id].tau2 = 0.9 # TODO: Tune tau2 further for accurate response 
            gaba_synapses[gaba_id].erev = gaba_reversal

            # Stim to play back spike times
            gaba_stimulator[gaba_id] = h.VecStim()

            # Use with deterministic onset times
            gaba_stimulator[gaba_id].play(h.Vector(1, gaba_onsets[gaba_id]))

            # Connect stim and syn
            gaba_connection[gaba_id] = h.NetCon(gaba_stimulator[gaba_id], gaba_synapses[gaba_id])
            gaba_connection[gaba_id].weight[0] = gaba_weight # Depending on desired EPSP response at soma, tune this

            if return_currents:
                # Measure conductance and current
                gaba_currents[gaba_id] = h.Vector().record(gaba_synapses[gaba_id]._ref_i)
                gaba_conductances[gaba_id] = h.Vector().record(gaba_synapses[gaba_id]._ref_g)

            gaba_id += 1 # increment gaba counter

        rounded_locs = [round(value, 4) for value in gaba_locs]
        print("# gaba synapses added:{} on:{} with locs:{} with timing onsets:{}".format(gaba_id, gaba_secs, gaba_locs, gaba_onsets))
    return gaba_synapses, gaba_stimulator, gaba_connection, gaba_currents, gaba_conductances, gaba_locs

def count_unique_dends(input_list):
    unique_names = set(input_list)
    count = len(unique_names)
    return count

# Records voltage across all sections
def record_all_v(cell, loc=0.4):
    all_v = {}
    for sec in cell.allseclist:
        all_v[sec.name()] = h.Vector()
        all_v[sec.name()].record(sec(loc)._ref_v) # given a sec with multiple seg, 
    return all_v

# Records voltage across selected sections
def record_v(cell, seclist, loc=0.4):
    all_v = {}
    for sec in seclist:
        all_v[sec.name()] = h.Vector()
        all_v[sec.name()].record(sec(loc)._ref_v) # given a sec with multiple seg, 
    return all_v


# function gets all unique locations on path to soma
def record_all_path_secs_v(cell, dend_tree, dendrite):
    all_v = {}
    dists = []
    for dend in cell.allseclist:
        if dend.name() == dendrite:
            dendrite = dend
    # get path to soma
    if dendrite.name() != 'soma[0]':
        pathlist = path_finder(cell=cell, dend_tree=dend_tree, dend=dendrite)
    else:
        pathlist = [dendrite]
    # for each dendrite in path find unique locations corresponding to each seg of that dendrite
    i=0
    for sec in pathlist:
        for seg in sec:
            dist = h.distance(seg)
            dists.append(dist)
            loc = seg.x
            all_v[i] = h.Vector()
            all_v[i].record(sec(loc)._ref_v) # given a sec with multiple seg
            i = i + 1
    return all_v, dists

# Records Cai across selected sections
def record_cai(cell, seclist, loc=0.4, return_Ca=True):
    all_cai = {}
    if return_Ca:
        for sec in seclist:
            all_cai[sec.name()] = h.Vector()
            all_cai[sec.name()].record(sec(loc)._ref_cai) # given a sec with multiple seg, 
    return all_cai

# Returns vectors for impedance recording
def record_impedance(dend, loc=0.4):
    imp = h.Impedance()
    # imp.loc(0.5, sec=cell.soma) 
    # define location either for current stim or voltage measuring electrode
    # this is needed for the transfer impedance calculation
    imp.loc(loc, sec=dend) # location of interest; nb voltages are measured at 0.4 ; not necessary if computing imp.input()  
    zvec1 = h.Vector()  
    zvec1.append(0)
    zvec2 = h.Vector()  
    zvec2.append(0)
    return imp, zvec1, zvec2

def record_i_mechs(cell, dend, loc=0.4, return_currents=True):
    i_mechs_dend = []
    if return_currents:
        t = h.Vector().record(h._ref_t) 
        # record for ONE glut input ONLY
#         if (len(dend) != 1): dend = dend[0]
        print("i_mechanisms recorded in {}".format(dend))
        ipas = h.Vector().record(dend(loc)._ref_i_pas)
        ikdr = h.Vector().record(dend(loc)._ref_ik_kdr)
        inaf = h.Vector().record(dend(loc)._ref_ina_naf)
        ikaf = h.Vector().record(dend(loc)._ref_ik_kaf)
        ikas = h.Vector().record(dend(loc)._ref_ik_kas)
        ikir = h.Vector().record(dend(loc)._ref_ik_kir)    
        ical12 = h.Vector().record(dend(loc)._ref_ical_cal12)    
        ical13 = h.Vector().record(dend(loc)._ref_ical_cal13)
        ican = h.Vector().record(dend(loc)._ref_ica_can)
        icar = h.Vector().record(dend(loc)._ref_ica_car)
        icav32 = h.Vector().record(dend(loc)._ref_ical_cav32)
        icav33 = h.Vector().record(dend(loc)._ref_ical_cav33)
        isk = h.Vector().record(dend(loc)._ref_ik_sk)
        ibk = h.Vector().record(dend(loc)._ref_ik_bk)
        # make df
        i_mechs_dend = [t,ipas,ikdr,inaf,ikaf,ikas,ikir,ical12,ical13,ican,icar,icav32,icav33,isk,ibk]
        return i_mechs_dend

# for plotting, returns all branch dendrites
def dend2plot(cell, cell_type='dspn'):
    [branch1_dends, branch2_dends, branch3_dends, branch4_dends, branch5_dends] = branch_selection(cell, cell_type) 
    branch_dends = [branch1_dends] + [branch2_dends] + [branch3_dends] + [branch4_dends] + [branch5_dends]
    branch_dends = [num for sublist in branch_dends for num in sublist]
    return [cell.soma] + branch_dends

def plot1(cell=None, dend=None, t=None, v=None, seclist=None, sparse=False, protocol=''):
    import plotly.graph_objects as go
    v_data = []
    if sparse:
        for group in seclist: # for each nrn dendrite sec, one plot per branch
            if dend in group: # Use if you want sparse plotting
                for sec in group:
                    v_data.append(go.Scatter(x=t, y=v[sec.name()], name='{}:{}'.format(sec.name(), round(h.distance(sec(0.5)), 2))))
                v_data.append(go.Scatter(x=t, y=v['soma[0]'], name='soma'))
    else:
        for sec in seclist: # for each nrn dendrite sec, one plot per branch
            v_data.append(go.Scatter(x=t, y=v[sec.name()], name='{}:{}'.format(sec.name(), round(h.distance(sec(0.5)), 2))))
    
    # Plot vdata
    fig = go.Figure(data=v_data)
    fig.update_layout(
        title="{}".format(protocol),
        title_x=0.5,
        xaxis_title="time (ms)",
        yaxis_title="V (mV)",
        legend_title="section")
    return fig   

def plot1_Ca(cell=None, dend=None, t=None, Ca=None, seclist=None, sparse=False, protocol=''):
    import plotly.graph_objects as go
    Ca_data = []
    if sparse:
        for group in seclist: # for each nrn dendrite sec, one plot per branch
            if dend in group: # Use if you want sparse plotting
                for sec in group:
                    Ca_data.append(go.Scatter(x=t, y=Ca[sec.name()]*1e3, name='{}:{}'.format(sec.name(), round(h.distance(sec(0.5)), 2))))
                Ca_data.append(go.Scatter(x=t, y=Ca['soma[0]']*1e3, name='soma'))
    else:
        for sec in seclist: # for each nrn dendrite sec, one plot per branch
            Ca_data.append(go.Scatter(x=t, y=Ca[sec.name()]*1e3, name='{}:{}'.format(sec.name(), round(h.distance(sec(0.5)), 2))))
    
    # Plot vdata
    fig = go.Figure(data=Ca_data)
    fig.update_layout(
        title="{}".format(protocol),
        title_x=0.5,
        xaxis_title="time (ms)",
        yaxis_title="[Ca] (uM)",
        legend_title="section")
    return fig   


def plot2(soma_v_data, dend_v_data, glut_placement=None, yaxis='V (mV)'):
    import plotly.graph_objects as go
    if yaxis=='V (mV)':
        title1='soma PSP'
        if glut_placement == 'distal':
            title2 = 'distal dendrite PSP'
        elif glut_placement == 'proximal':
            title2 = 'proximal dendrite PSP'
        else:
            title2 = 'dendrite PSP'
    else:
        title1='soma impedance'
        if glut_placement == 'distal':
            title2 = 'distal dendrite impedance'
        elif glut_placement == 'proximal':
            title2 = 'proximal dendrite impedance'
        else:
            title2 = 'dendrite impedance'
            
    fig_soma = go.Figure(data=soma_v_data)
    fig_soma.update_layout(
        title=title1,
        title_x=0.5,
        xaxis_title='time (ms)',
        yaxis_title=yaxis,
        legend_title='sim')
#     fig_soma.show() 

    fig_dend = go.Figure(data=dend_v_data)
    fig_dend.update_layout(
        title=title2,
        title_x=0.5,
        xaxis_title='time (ms)',
        yaxis_title=yaxis,
        legend_title='sim')
#         fig_dend.show() 
    return fig_soma, fig_dend  

 
def plot3(somaV, dendV, glut_placement=None, yaxis='V (mV)', yrange_soma=[-85,-60], yrange_dend=[-85,-30], stim_time = 150, sim_time=400, black_trace=0, gray_trace=None, err_bar=50, baseline=20, dt=0.025, width=500, height=500):

    def hex_palette(n):
        import seaborn as sns
        import matplotlib as mpl
        colors = ['#6A5ACD', '#CD5C5C', '#458B74', '#9932CC', '#FF8247'] # Set your custom color palette
        if n < len(colors):
            colors = colors[0:n]
        else:
            colors = sns.blend_palette(colors,n)
        cols = list(map(mpl.colors.rgb2hex, colors))
        return cols
    
    n = len(somaV)
    cols = hex_palette(n)
    # this routine places black and gray traces if required
    if black_trace == None and gray_trace == None:
        cols = hex_palette(n)
    elif black_trace is not None and gray_trace == None:
        cols = hex_palette(n-1)
        cols.insert(black_trace,'#000000')
    elif black_trace is not None and gray_trace is not None:  
        cols = hex_palette(n-2)
        if (black_trace>gray_trace):
            cols.insert(gray_trace,'#D3D3D3')
            cols.insert(black_trace,'#000000')
        elif (black_trace<gray_trace):
            cols.insert(black_trace,'#000000')
            cols.insert(gray_trace,'#D3D3D3')
            
    def update_layout(fig, main=None, yaxis=None, yrange=[-85,-60], width=500, height=500):
        font = 'Droid Sans'
        font_size = 18
        fig.update_layout(
        autosize=False,
        width=width,
        height=height,
        margin=dict(l=20, r=20, t=30, b=20),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        title=main,
        title_x=0.45,
        title_font_family=font,
        title_font_size=font_size,
        xaxis=dict(showticklabels=False, titlefont=dict(size=font_size, family=font), tickfont=dict(size=font_size, family=font), showgrid=False), # title='time (ms)', 
        yaxis=dict(side='right', tick0 = yrange[0], dtick = (yrange[1]-yrange[0]), tickfont=dict(size=font_size, family=font), showgrid=False),
        legend = dict(title='sim', x=1.1, y=0.95)
    )

    def plot3_(somaV, dendV, glut_placement, yaxis='V (mV)', cols=None, yrange_soma=[-85,-60], yrange_dend=[-85,-30], err_bar=50, bl = 20):
        import plotly.graph_objects as go
        if yaxis=='V (mV)':
            title1='soma PSP'
            if glut_placement == 'distal':
                title2 = 'distal dendritic PSP'
            elif glut_placement == 'proximal':
                title2 = 'proximal dendritic PSP'
            else:
                title2 = 'dendritic PSP'
                
        else:
            title1='soma impedance'
            if glut_placement == 'distal':
                title2 = 'distal dendritic impedance'
            elif glut_placement == 'proximal':
                title2 = 'proximal dendritic impedance'
            else:
                title2 = 'dendritic impedance'

        figSoma = go.Figure()
        
        ind1 = 0
        ind2 = int((sim_time - stim_time + bl)/dt)
        ind3 = int((stim_time - bl)/dt)
        ind4 = int(sim_time/dt)
        
        for ii in range(len(somaV)):
            dat = somaV[ii]
            figSoma.add_trace( go.Scatter(x=dat['x'][ind1:ind2], y=dat['y'][ind3:ind4], mode='lines', line=dict(color=cols[ii])) )
        figSoma.add_hline(y=yrange_soma[0], line_width=2, line_dash="dot", line_color="gray")
        figSoma.add_hline(y=yrange_soma[1], line_width=2, line_dash="dot", line_color="gray")        
        figSoma.add_shape(
                    type='line',
                    x0=ind2*dt-err_bar,
                    y0=yrange_soma[0]+2,
                    x1=ind2*dt,
                    y1=yrange_soma[0]+2,
                    line=dict(color='black'),
                    xref='x',
                    yref='y'
        )
        update_layout(fig=figSoma, main=title1, yaxis=yaxis, yrange=yrange_soma, width=width, height=height)   
        
        figDend = go.Figure()
        for ii in range(len(dendV)):
            dat = dendV[ii]
            figDend.add_trace( go.Scatter(x=dat['x'][ind1:ind2], y=dat['y'][ind3:ind4], mode='lines', line=dict(color=cols[ii])) )
        figDend.add_hline(y=yrange_dend[0], line_width=2, line_dash="dot", line_color="gray")
        figDend.add_hline(y=yrange_dend[1], line_width=2, line_dash="dot", line_color="gray")
        update_layout(fig=figDend, main=title2, yaxis=yaxis, yrange=yrange_dend, width=width, height=height)
        
        return figSoma, figDend


    fig_soma_master, fig_dend_master =  plot3_(somaV=somaV, dendV=dendV, glut_placement=glut_placement, yaxis=yaxis, cols=cols, yrange_soma=yrange_soma, yrange_dend=yrange_dend, err_bar=err_bar, bl=baseline)    

    return fig_soma_master, fig_dend_master 


def save_fig2(soma_fig=None, dend_fig=None, cell_type='dspn', model=None, physiological=True, sim=None, g_name=None):
    import datetime 
    time = datetime.datetime.now()
    path_cell = "{}".format(cell_type)
    if not os.path.exists(path_cell):
        os.mkdir(path_cell)    
    path1 = "{}/model{}".format(path_cell, model)
    if not os.path.exists(path1):
        os.mkdir(path1)
    if physiological: 
        path2 = "{}/physiological".format(path1)
    else:
        path2 = "{}/nonphysiological".format(path1)
    if not os.path.exists(path2):
        os.mkdir(path2)  
    path3 = "{}/images".format(path2)
    if not os.path.exists(path3):
        os.mkdir(path3)

    image_dir = "{}/sim{}".format(path3, sim)
    if not os.path.exists(image_dir):
        os.mkdir(image_dir)    

    if (g_name is None):
        soma_fig.write_image("{}/soma_fig{}.svg".format(image_dir, time))
        dend_fig.write_image("{}/dend_fig{}.svg".format(image_dir, time))
        soma_fig.write_html("{}/soma_fig{}.html".format(image_dir, time))
        dend_fig.write_html("{}/dend_fig{}.html".format(image_dir, time))
    else:
        soma_fig.write_image("{}/{}_soma_fig{}.svg".format(image_dir, g_name, time))
        dend_fig.write_image("{}/{}_dend_fig{}.svg".format(image_dir, g_name, time))  
        soma_fig.write_html("{}/{}_soma_fig{}.html".format(image_dir, g_name, time))
        dend_fig.write_html("{}/{}_dend_fig{}.html".format(image_dir, g_name, time))

    
def convert2df(d, g_name):
    df = pd.DataFrame() 
    lists = sorted(d.items()) # sorted by key, return a list of tuples
    x, y = zip(*lists) # unpack a list of pairs into two tuples
    df['dist'] = x
    df[g_name] = y
    return df

def dist_(cell, g_name):
    if g_name in ['naf', 'kaf', 'kas', 'kdr', 'kir', 'sk', 'bk', 'gaba1', 'gaba2']:
        gbar = 'gbar_{}'.format(g_name)
    else:
        gbar = 'pbar_{}'.format(g_name)        
    d__ = {}
    for sec in cell.dendlist:
         d__[h.distance(sec(0.5))] = eval('sec.{}'.format(gbar))
    return convert2df(d__, g_name)

def plot4(data, g_name): 
    import plotly.graph_objects as go
    y = data[0].y
    num = y.max()
    sig_fig = len(str(num)) - str(num).find('.') - 1
    y1 = 2*round(num, sig_fig)
    if (y1==0):
        y1 = 1e-4
    fig = go.Figure(data=data)
    fig.update_layout(
        title="{}".format(g_name),
        title_x=0.5,
        yaxis=dict(range=[0, y1]),
        xaxis_title="distance (um)",
        yaxis_title="conductance (S/cm2)",
        legend_title="cond")
    return fig

def plot5(X, dt, dists, xaxis_range=[0,150], yaxis_range=[0,8], normalised=True, title='', voltage=True):
    t2 = np.arange(0, len(X[0]), 1) * dt
    import plotly.graph_objects as go
    v_data = []
    if normalised:
        yaxis_title="normalised amplitude"
    else:
        if voltage:
            yaxis_title="V (mV)" 
        else:
            yaxis_title="I (pA)" 
            
    for ii in list(range(len(X))):
        v_data.append(go.Scatter(x=t2, y=X[ii], name='{}'.format(round(dists[ii], 2))))
    # Plot vdata
    fig = go.Figure(data=v_data)
    fig.update_layout(
        title="{}".format(title),
        title_x=0.5,
        xaxis_title="time (ms)",
        yaxis_title=yaxis_title,
        xaxis_range = xaxis_range,
        yaxis_range = yaxis_range,
        legend_title="distance (um)")
    return fig

def plot5a(X, dt, locs, xaxis_range=[0,150], yaxis_range=[0,-30], normalised=False, col=[], title=''):
    t2 = np.arange(0, len(X[0]), 1) * dt
    import plotly.graph_objects as go
    v_data = []
    if normalised:
        yaxis_title="normalised amplitude"
    else:
        yaxis_title="I (pA)"        
    for ii in list(range(len(X))):
        if len(col) == 0:
            v_data.append(go.Scatter(x=t2, y=X[ii], name='{}'.format(locs[ii])))
        else:
            v_data.append(go.Scatter(x=t2, y=X[ii], line=dict(color=col[ii]), name='{}'.format(locs[ii])))
    # Plot vdata
    fig = go.Figure(data=v_data)
    fig.update_layout(
        title="{}".format(title),
        title_x=0.5,
        xaxis_title="time (ms)",
        yaxis_title=yaxis_title,
        xaxis_range = xaxis_range,
        yaxis_range = yaxis_range,
        legend_title="location")
    return fig

def plot5b(X, dt, locs, xaxis_range=[0,150], yaxis_range=[0,-30], normalised=False, dotted=False, col=[], title=''):
    t2 = np.arange(0, len(X[0]), 1) * dt
    import plotly.graph_objects as go
    v_data = []
    if normalised:
        yaxis_title="normalised amplitude"
    else:
        yaxis_title="V (mV)"        
    for ii in list(range(len(X))):
        if dotted:
            v_data.append(go.Scatter(x=t2, y=X[ii], line=dict(dash='dot', color='gray'), showlegend=False))
        else:
            if len(col) == 0: 
                v_data.append(go.Scatter(x=t2, y=X[ii], name='{}'.format(locs[ii])))
            else:
                v_data.append(go.Scatter(x=t2, y=X[ii], line=dict(color=col[ii]), name='{}'.format(locs[ii])))
    # Plot vdata
    fig = go.Figure(data=v_data)
    fig.update_layout(
        title="{}".format(title),
        title_x=0.5,
        xaxis_title="time (ms)",
        yaxis_title=yaxis_title,
        xaxis_range = xaxis_range,
        yaxis_range = yaxis_range,
        legend_title="location")
    return fig

# remove offsets
def normalise(X, stim_time, burn_time, dt):    
    def mean(x):
        n = len(x)
        sum = 0
        for i in x:
            sum = sum + i
        return(sum/n)
    ind1 = int(burn_time/dt)
    ind2 = int(stim_time/dt)
    return(X[ind1:len(X)] - mean(X[ind1:ind2]) )

def plot6(y, x, xaxis_range=[200,0], yaxis_range=[0,1.01], normalised=True):
    import plotly.express as px
    if normalised:
        yaxis_title="normalised amplitude"
    else:
        yaxis_title="V (mV)"        
    fig2 = px.scatter(x=x, y=y)
    fig2.update_layout(
        title="{}".format(''),
        title_x=0.5,
        xaxis_title="distance (um)",
        yaxis_title=yaxis_title,
        xaxis_range = xaxis_range,
        yaxis_range = yaxis_range,
        legend_title="attenuation")
    return fig2

def plot6a(mat, x, xaxis_range=[200,0], yaxis_range=[0,1.01], normalised=True, col=[], current=True):
    import plotly.graph_objects as go
    i_data = []
    if normalised:
        if current:
            yaxis_title="normalised PSC"
        else:
            yaxis_title="normalised PSP"
    else:
        if current:
            yaxis_title="I (pA)" 
        else:
            yaxis_title="V (mV)" 

    rows, columns = mat.shape
    if rows == 3:
        names =['spine', 'dendrite', 'soma']
    else:
        names =['dendrite', 'soma']
    for ii in list(range(rows)):
        if len(col) == 0:
            i_data.append(go.Scatter(x=x, y=mat[ii,:], line=dict(color='gray'), name='{}'.format(names[ii]), showlegend=True))
        else:
            i_data.append(go.Scatter(x=x, y=mat[ii,:], line=dict(color=col[ii]), name='{}'.format(names[ii]), showlegend=True))

    fig2 = go.Figure(data=i_data)
    fig2.update_layout(
        title="{}".format(''),
        title_x=0.5,
        xaxis_title="distance (um)",
        yaxis_title=yaxis_title,
        xaxis_range = xaxis_range,
        yaxis_range = yaxis_range,
        legend_title="attenuation")
    return fig2

def plot6aa(mat, x, xaxis_range=[200,0], yaxis_range=[0,1.01], normalised=True, col=[], current=True):
    import plotly.graph_objects as go
    i_data = []
    if normalised:
        if current:
            yaxis_title="normalised PSC"
        else:
            yaxis_title="normalised PSP"
    else:
        if current:
            yaxis_title="I (pA)" 
        else:
            yaxis_title="V (mV)" 

    rows, columns = mat.shape
    if rows == 3:
        names =['spine', 'dendrite', 'soma']
    else:
        names =['dendrite', 'soma']
    for ii in list(range(rows)):
        if len(col) == 0:
            i_data.append(go.Scatter(x=x, y=mat[ii,:], line=dict(color='gray'), name='{}'.format(names[ii]), showlegend=True))
        else:
            i_data.append(go.Scatter(x=x, y=mat[ii,:], line=dict(color=col[ii]), name='{}'.format(names[ii]), showlegend=True))

    fig2 = go.Figure(data=i_data)
    fig2.update_layout(
        title="{}".format(''),
        title_x=0.5,
        xaxis_type='log',
        xaxis_title="series resistance (MOhm)",
        yaxis_title=yaxis_title,
        xaxis_range = xaxis_range,
        yaxis_range = yaxis_range,
        legend_title="attenuation")
    return fig2

# Find all unique synapse locations in dendritic path 
# if GABA then want all unique locations
# if Glutamate then only where there is a spine
# start from end of particular dendrite and move towards soma

def gaba_idx(dend):
    locs = []
    for seg in dend:
        locs.append(seg.x)
    return locs

def all_synapses_tree(cell, dend_tree, dendrite, glut):
    candidate_list = []
    locs_list = []
    for dend in cell.dendlist:
        if dend.name() == dendrite:
            dendrite = dend
        # get path to soma
    pathlist = path_finder(cell=cell, dend_tree=dend_tree, dend=dendrite)
    if glut:
        # get all unique spine locations to place glut synapses
        for sec in pathlist:
            dend_dist = h.distance(cell.soma(), sec(1)) # distance to end of that dendrite from middle of soma
            locs = []
            if dend_dist < 30:
                candidates = []
            else:
                candidates = spine_idx(cell=cell, spines=spines, dend=sec.name())
                for candidate in candidates:
                    locs.append(candidate.x)

            candidate_list.append(candidates)
            locs_list.append(locs)
    else:
        # get all unique gaba synapse locations
        locs_list = []
        for sec in pathlist:
            locs = gaba_idx(sec)
            locs_list.append(locs)

    return pathlist, locs_list, candidate_list

def space_clamped(cell, spines, Ra = 1.59e-10):
    for sec in cell.allseclist:
        sec.Ra = Ra
    for sec in cell.dendlist:
        sec_spines = list(spines[sec.name()].items())
        for spine_i, spine_obj in sec_spines: 
            spine_obj.head.Ra = Ra
            spine_obj.neck.Ra = Ra

def cap(cell, spines, cm = 1):
    for sec in cell.allseclist:
        sec.cm = cm
    for sec in cell.dendlist:
        sec_spines = list(spines[sec.name()].items())
        for spine_i, spine_obj in sec_spines: 
            spine_obj.head.cm = cm
            spine_obj.neck.cm = cm
            
def find_closest_value(test, target_value):
    import numpy as np
    test = np.array(test)
    distances = np.abs(test - target_value)
    closest_index = np.argmin(distances)
    closest_value = test[closest_index]
    return closest_value, closest_index

def rounded(number, n=10):
    if number >= 0:
        rounded = n*math.ceil(number/n)
    else:
        rounded = n*math.floor(number/n)
    return rounded


def IR(X, step_start, step_end, step):   
    ind1 = int((step_start-5)/dt)
    ind2 = int(step_start/dt)
    ind3 = int((step_end-5)/dt)
    ind4 = int(step_end/dt)

    def mean(x):
        n = len(x)
        sum = 0
        for i in x:
            sum = sum + i
        return(sum/n)
 
    return(1e3 *( mean(X[ind1:ind2]) - mean(X[ind3:ind4]) )  / -step) # MOhm


def whole_cell_capacitance(cell, spines=None, Cm=1):
    # for seg in sec.allseg():
    #     print(seg.area())

    # for sec in cell.dendlist:
    #     print(seg.area())
    areas = []
    for sec in cell.dendlist:
        for i,seg in enumerate(sec):
            areas.append(seg.area())

    if spines is None:
        AREA = sum(areas)
    else:
        spine_areas = []
        for sec in cell.dendlist:
            sec_spines = list(spines[sec.name()].items())
            for spine_i, spine_obj in sec_spines: # area of sphere + cylinder less areas that are connections to head and dendrite
                spine_areas.append(4 * math.pi * ((spine_obj.head.diam/2) ** 2) + 2 * math.pi * spine_obj.neck.L * spine_obj.neck.diam/2  - 2 * math.pi * ((spine_obj.neck.diam/2)**2) ) # diam

        AREA = sum(spine_areas) + sum(areas) +  4 * math.pi * cell.soma.diam/2 ** 2 # in um sq
    # AREA in um sq
    # Cm in uF/cm2
    # convert uF/cm2 to F/um2
    # To convert square centimeters (cm²) to square micrometers (μm²), you can use the following conversion factor:
    # 1 cm² = 1e8 μm²
    # convert uF to pF conversion factor
    # 1uF = 1e6 pF
    
    # AREA * Cm * 1e-8 * 1e6 # (cm2)
    cap = AREA * Cm * 1e-2 # pF
    return cap
 
def sampler(names, n, replacement=True):
    if replacement:
        return [random.choice(names) for _ in range(n)]
    else:
        return random.sample(names, n)
    
def uniform_values(n):
    return [random.uniform(0, 1) for _ in range(n)]

# function to determine what is varying
def variable_detector(xrange):
    vary = False
    differences = [abs(xrange[i] - xrange[i + 1]) for i in range(len(xrange) - 1)]
    sum_of_differences = sum(differences)
    if sum_of_differences > 0:
        vary = True
    return vary

def spine_locator(cell_type, specs, spine_per_length, dend_glut, num_gluts, method=0):
    
    # Create the cell and get the dendrite list
    cell, spines, dend_tree = cell_build(cell_type=cell_type, specs=specs, addSpines=True, spine_per_length=spine_per_length)
    # Get target dendrites for glutamate synapses
    glut_secs = [sec for target_dend in dend_glut for sec in cell.dendlist if sec.name() == target_dend] * num_gluts

    final_spine_locs = []
    glut_id = 0 # Index used for glut_synapse list and printing

    for dend_glut in glut_secs:
        # Get candidate spines from section
        candidate_spines = []
        sec_spines = list(spines[dend_glut.name()].items())

        for spine_i, spine_obj in sec_spines: 
            candidate_spines.append(spine_obj)

        if len(glut_secs) < len(sec_spines):
            if method==1:
                # Reverse order so activate along dendrite towards soma
                spine_idx = 2*len(candidate_spines)//3-1 # Arbitrary start point at 2/3 of spines
                spine = candidate_spines[spine_idx - glut_id] 
            else:
                spine_idx = 2*len(candidate_spines)//3 - num_gluts # Arbitrary start point at 1/3 of spines
                if spine_idx < 0:
                    if len(candidate_spines) >= num_gluts:
                        spine_idx = len(candidate_spines) - num_gluts
                    else:
                        spine_idx = 0        
                spine = candidate_spines[spine_idx + glut_id] 
        else:
            spine = random.choice(candidate_spines)

        spine_loc = spine.x
        final_spine_locs.append(spine_loc) 
        glut_id = glut_id + 1
    return final_spine_locs

# relative version of gaba_onsets
def rel_gaba_onset(n, N):
    if N in [0,1]:
        if (n < 4):
            gaba_onsets = list(range(0, 0 + n)) 
        else:
            if n % 3 == 0:
                gaba_onsets = list(range(0, 0 + int(n/3))) * 3 * N
            else:
                gaba_onsets = list(range(0, 0 + int(n/3)+1)) * 3 * N
        gaba_onsets = gaba_onsets[:n]
    else:
        n1 = round(n / N)
        if n1 % 3 == 0:
            n1 = round(n / N)
            onsets = list(range(0, n1)) 
        else:
            onsets = list(range(0, n1+1)) 
        gaba_onsets = [x for x in onsets for _ in range(N)]
        gaba_onsets = gaba_onsets[:n]
    return gaba_onsets


def plt1(data_dict, sim_time, dt, model, cell_type, stim_time, showplot, save, time, physiological, sim):
    idx1 = 0 
    idx2 = int(sim_time/dt)

    black_trace = 0
    gray_trace = 1

    soma_v_master = []
    dend_v_master = []

    for ii in list(range(len(data_dict['soma_v']))):
        x = data_dict['time'][ii]
        y = data_dict['soma_v'][ii]
        soma_v_master.append(go.Scatter(x=extract2(x)[idx1:idx2], y=extract2(y)[idx1:idx2]))
        y = data_dict['dend_v'][ii]
        dend_v_master.append(go.Scatter(x=extract2(x)[idx1:idx2], y=extract2(y)[idx1:idx2]))
        

    if model == 2 or cell_type == 'ispn':
        yrange_soma1=[-85,-50]
        yrange_dend1=[-85,-20]
    else:
        yrange_soma1=[-85,-60]
        yrange_dend1=[-85,-30]  
        
    fig_soma_master, fig_dend_master =  plot3(somaV=soma_v_master, dendV=dend_v_master, glut_placement='', 
        yaxis='V (mV)', yrange_soma=yrange_soma1, yrange_dend=yrange_dend1, stim_time=stim_time, 
        sim_time=sim_time, black_trace=black_trace, gray_trace=gray_trace, err_bar=50, baseline=20, dt=dt,
        width=800, height=400)    

    if showplot:
        fig_soma_master.show()
        fig_dend_master.show()
        
    if save:
        images_path = '{}/model{}/{}/images/sim{}'.format(cell_type, model, 'physiological' if physiological else 'nonphysiological', sim)
        os.makedirs(images_path, exist_ok=True)
        fig_soma_master.write_image('{}/fig1_soma{}.svg'.format(images_path, time))
        fig_soma_master.write_html('{}/fig1_soma{}.html'.format(images_path, time))
        fig_dend_master.write_image('{}/fig1_dend{}.svg'.format(images_path, time))
        fig_dend_master.write_html('{}/fig1_dend{}.html'.format(images_path, time))

        
        
def plt2(data_dict, sim_time, dt, model, cell_type, stim_time, showplot, save, time, physiological, sim, offset=40
    ):
    idx1, idx2, idx3  = 0, int(sim_time/dt), int(stim_time/dt)
    soma_v_master, dend_v_master, spine_v_master = [], [], []
    peaks_v_soma, peaks_v_dend, peaks_v_spine = [], [], []
    mins_v_soma, mins_v_dend, mins_v_spine = [], [], []
    
    for ii, (soma_v, dend_v, spine_v) in enumerate(zip(data_dict['soma_v'], data_dict['dend_v'], data_dict['spine_v'])):
        x=extract2(data_dict['time'][ii])[idx1:idx2]
        
        y=extract2(soma_v)[idx1:idx2]
        soma_v_master.append(go.Scatter(x=x, y=y))
        peaks_v_soma.append(max(soma_v))
        mins_v_soma.append(soma_v[idx3])
        
        y=extract2(dend_v)[idx1:idx2]
        dend_v_master.append(go.Scatter(x=x, y=y))
        peaks_v_dend.append(max(dend_v))
        mins_v_dend.append(dend_v[idx3])
        
        y=extract2(spine_v)[idx1:idx2]
        spine_v_master.append(go.Scatter(x=x, y=y))
        peaks_v_spine.append(max(spine_v))
        mins_v_spine.append(spine_v[idx3])
    ysoma_range = [0.1*math.floor(min(mins_v_soma)/0.1), 0.1*math.ceil(max(peaks_v_soma)/0.1)]
    ydend_range = [math.floor(min(mins_v_dend)), math.ceil(max(peaks_v_dend))]
    yspine_range = [math.floor(min(mins_v_spine)), math.ceil(max(peaks_v_spine))]

    figs = plot3a(somaV=soma_v_master, dendV=dend_v_master, spineV=spine_v_master, ysoma_range=ysoma_range, 
                  ydend_range=ydend_range, yspine_range=yspine_range, stim_time=stim_time, sim_time=sim_time, 
                  dt=dt, offset=offset)

    if showplot:
        for fig in figs:
            fig.show()
              
    if save:
        path_format = f'{cell_type}/model{model}/{{}}/images/sim{sim}'
        folder = path_format.format('physiological' if physiological else 'nonphysiological')
        os.makedirs(folder, exist_ok=True)

        for idx, name in enumerate(['soma', 'dend']):
            figs[idx].write_image(f'{folder}/fig1_{name}{time}.svg')
            figs[idx].write_html(f'{folder}/fig1_{name}{time}.html')

def hex_palette(n):
    colors = ['#6A5ACD', '#CD5C5C', '#458B74', '#9932CC', '#FF8247'] # Set your custom color palette
    if n < len(colors):
        colors = colors[0:n]
    else:
        colors = sns.blend_palette(colors,n)
    cols = list(map(mpl.colors.rgb2hex, colors))
    return cols

def update_layout(fig, title, yaxis, yrange, width, height):
    font = 'Droid Sans'
    font_size = 18
    fig.update_layout(
        autosize=False,
        width=width,
        height=height,
        margin=dict(l=20, r=20, t=30, b=20),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)',
        title=title,
        title_x=0.45,
        title_font_family=font,
        title_font_size=font_size,
        xaxis=dict(showticklabels=False, titlefont=dict(size=font_size, family=font), tickfont=dict(size=font_size, family=font), showgrid=False),
        yaxis=dict(side='right', tick0=yrange[0], dtick=yrange[1]-yrange[0], tickfont=dict(size=font_size, family=font), showgrid=False),
        legend=dict(title='sim', x=1.1, y=0.95)
    )

# offset corrects baseline
# stim_time - offset gives some baseline dependent on the smallest value of timing_range ie 
# baseline = offset - stim_time + min(timing_range) ie if stim_time = 150, first timing at 120 then offset=40 gives a 10 ms baseline

def plot_trace(fig, data_list, cols, dt, err_bar, yrange, stim_time, sim_time, offset):
    ind1, ind2 = 0, int((sim_time - stim_time + offset)/dt)
    ind3, ind4 = int((stim_time - offset)/dt), int(sim_time/dt)
    dy = (yrange[1] - yrange[0])*9/10
    for dat, color in zip(data_list, cols):
        fig.add_trace(go.Scatter(x=dat['x'][ind1:ind2], y=dat['y'][ind3:ind4], mode='lines', line=dict(color=color)))
    fig.add_hline(y=yrange[0], line_width=2, line_dash="dot", line_color="gray")
    fig.add_hline(y=yrange[1], line_width=2, line_dash="dot", line_color="gray")
    fig.add_shape(type='line', x0=ind2*dt-err_bar, y0=yrange[0]+dy, x1=ind2*dt, y1=yrange[0]+dy, line=dict(color='black'))

def plot3a(somaV, dendV, spineV, ysoma_range, ydend_range, yspine_range, stim_time, sim_time, dt, offset, width=800, height=400):
    n = len(somaV)
    cols = hex_palette(n-2)
    cols.insert(0,'#000000')
    cols.insert(1,'#D3D3D3')
   
    titles = ['spine PSP', 'dendritic PSP', 'soma PSP']
    yranges = [yspine_range, ydend_range, ysoma_range]
    data_list = [spineV, dendV, somaV]
    
    figs = []
    for data, title, yrange in zip(data_list, titles, yranges):
        fig = go.Figure()
        plot_trace(fig, data, cols, dt, 25, yrange, stim_time, sim_time, offset)
        update_layout(fig, title, 'V (mV)', yrange, width, height)
        figs.append(fig)
    return figs


def check_sim(sim, sims):
    """
    Checks if sim starts with or is equivalent to any value in values_to_check.
    
    :param sim: Value to check (can be integer or string)
    :param values_to_check: List of values to check against (mix of integers and strings)
    :return: True if there's a match, False otherwise
    """
    strings_to_check = [str(val) for val in sims]  # Convert all values to strings

    sim_str = str(sim)  # Convert sim to a string regardless of its type

    return any(sim_str.startswith(s) or sim_str == s for s in strings_to_check)


# # Example Usage:
# sims = [10011, 10012, '10011a']
# sim = '101a'
# result = check_sim(sim, sims)
# print(result)  # Outputs True or False based on whether sim matches the pattern

def update_data_dict(data_dict, protocol, v_tree, norm_v_tree, soma_v, dend_v, spine_v, timing, t, dists, i_dend_mechs,
                     ampa_currents, nmda_currents, gaba_currents, gaba_conductances, time, impedance=False, return_currents=False, record_spine=False):
    
    data_dict['v_tree'][protocol] = v_tree
    data_dict['norm_v_tree'][protocol] = norm_v_tree
    data_dict['soma_v'].append(soma_v[0])
    data_dict['dend_v'].append(dend_v[0])
    if record_spine:
        data_dict['spine_v'].append(spine_v[0])
    data_dict['timing'].append(timing)
    data_dict['time'].append(np.asarray(t))
    data_dict['dists'].append(dists)
    data_dict['i_mechs'][protocol] = i_dend_mechs
    
    data_dict['i_ampa'][protocol] = pd.DataFrame(ampa_currents).transpose()
    data_dict['i_nmda'][protocol] = pd.DataFrame(nmda_currents).transpose()
    data_dict['i_gaba'][protocol] = pd.DataFrame(gaba_currents).transpose()
    
    if impedance:
        data_dict['z_input'].append(z_input[0])
        data_dict['z_transfer'].append(z_transfer[0])
        
    if return_currents:
        data_dict['i_ampa'][protocol] = pd.DataFrame(ampa_currents).transpose()
        data_dict['i_nmda'][protocol] = pd.DataFrame(nmda_currents).transpose()
        data_dict['i_gaba'][protocol] = pd.DataFrame(gaba_currents).transpose()
        data_dict['g_gaba'][protocol] = pd.DataFrame(gaba_conductances).transpose()

    data_dict['timestamp'][protocol] = time
    
def plot_cumulative_frequency(cell_type, specs):
    """
    Plots the cumulative frequency of spine distances from the soma.

    Parameters:
    - cell_type: the type of cell
    - specs: cell specifications
    - spines: spine data

    Returns:
    - A plotly figure displaying the cumulative frequency
    """

    import plotly.graph_objects as go

    # Build the cell
    cell, spines, dend_tree = cell_build(cell_type=cell_type, specs=specs, addSpines=True)

    # Compute distances of each spine from the soma
    dist_spine = []
    for dend in cell.dendlist:
        sec_spines = list(spines[dend.name()].items())
        for spine in sec_spines:
            dist_spine.append(h.distance(cell.soma(0.5), dend(spine[1].x)))

    # Sort the data in ascending order
    sorted_data = sorted(dist_spine)

    # Calculate cumulative frequencies
    cumulative_freq = [i / len(sorted_data) for i in range(1, len(sorted_data) + 1)]

    # Create the cumulative frequency plot
    fig = go.Figure(data=go.Scatter(x=sorted_data, y=cumulative_freq, mode='lines'))

    # Set plot labels and title
    fig.update_layout(
        title={
            'text': 'cumulative frequency',
            'y': 0.9,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        xaxis_title='distance (um)',
        yaxis_title='cumulative frequency',
        xaxis_range=[0, 300],  # Set x-axis range
    )

    # Display the plot
    fig.show()
    
    
def plot_spine_distance_histogram(cell_type, specs, bin_size=10):
    """
    Plots the histogram of spine distances from the soma.

    Parameters:
    - cell_type: the type of cell
    - specs: cell specifications
    - bin_size: size of the bins for the histogram

    Returns:
    - A plotly figure displaying the histogram
    """

    import plotly.graph_objects as go

    # Build the cell
    cell, spines, dend_tree = cell_build(cell_type=cell_type, specs=specs, addSpines=True)

    # Compute distances of each spine from the soma
    dist_spine = []
    for dend in cell.dendlist:
        sec_spines = list(spines[dend.name()].items())
        for spine in sec_spines:
            dist_spine.append(h.distance(cell.soma(0.5), dend(spine[1].x)))

    # Create the histogram plot
    fig = go.Figure(data=go.Histogram(x=dist_spine, nbinsx=int(max(dist_spine) / bin_size)))

    # Set plot labels and title
    fig.update_layout(
        title={
            'text': 'Histogram of Spine Distances from Soma',
            'y': 0.9,
            'x': 0.5,
            'xanchor': 'center',
            'yanchor': 'top'
        },
        xaxis_title='distance (um)',
        yaxis_title='frequency',
        xaxis_range=[0, 300],  # Set x-axis range
    )

    # Display the plot
    fig.show()


