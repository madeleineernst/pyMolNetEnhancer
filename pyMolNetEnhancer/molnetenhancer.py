#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 16:19:11 2018

@author: Madeleine Ernst (https://github.com/madeleineernst)
"""
# Standard library imports
import os
import re
import sys
import time

# Third party imports 
import collections
from collections import Counter
from collections import OrderedDict
import csv  
import functools
from functools import reduce
from joblib import Parallel, delayed
import json
import multiprocessing
from networkx import *
import operator
import pandas as pd
import requests
import requests_cache

pd.options.mode.chained_assignment = None 
requests_cache.install_cache('demo_cache')

def Mass2Motif_2_Network(edges,motifs,prob = 0.01,overlap = 0.3, top = 5):
    """Map Mass2Motifs onto a mass spectral molecular network

    :param edges: An edges file downloaded from GNPS 
    :type edges: pandas.core.frame.DataFrame
    :param motifs: A motif summary file downloaded from MS2LDA
    :type motifs: pandas.core.frame.DataFrame
    :param prob: Minimal probability score for a Mass2Motif to be included 
    :type prob: float
    :param overlap: Minimal overlap score for a Mass2Motif to be included
    :type overlap: float
    :param top: Specifies how many most shared motifs per molecular family (network component index) should be shown
    :type top: int
    :return: A dictionary of two dataframes containing network nodes and edges with motifs mapped
    :rtype: dict

    """	
    motifs = motifs[motifs.probability > prob]
    motifs = motifs[motifs.overlap > overlap]
    
    if 'MEH' not in edges.columns and 'OtherScore' not in edges.columns:
        edges['MEH'] = 0.0 
        edges['OtherScore'] = 0.0
    
    motifs_con = motifs.groupby('scans').agg(lambda x: x.tolist())
    df = pd.DataFrame(0.0, index=motifs_con.index, columns=motifs.motif.unique())
    
    for index, row in motifs_con.iterrows():
        for idx, val in enumerate(row["motif"]):
            df.loc[index][val] = row["overlap"][idx]
    
    comb = pd.merge(motifs_con, df, left_index= True, right_index=True)
    
    edges['shared_motifs'] = 'None'

    for index, row in edges.iterrows():
        if row['CLUSTERID1'] in motifs_con.index and row['CLUSTERID2'] in motifs_con.index:
            edges['shared_motifs'].iloc[index] = list(set(motifs_con.loc[row['CLUSTERID1'],'motif']) & set(motifs_con.loc[row['CLUSTERID2'],'motif']))

    # calculate most shared motifs per molecular family
    topmotifs = edges.groupby('ComponentIndex')['shared_motifs'].agg(lambda x: x.tolist()).to_frame(name='topmotifs')

    for index,row in topmotifs.iterrows():
        if row['topmotifs'] != ['None']:
            mcounts = list(filter(lambda a: a != 'None', row['topmotifs']))
            topmotifs.loc[index,'topmotifs'] = [x for sub_list in mcounts for x in sub_list]

    for index,row in topmotifs.iterrows():
        counts = collections.Counter(row['topmotifs'])
        mcs = sorted(row['topmotifs'], key=counts.get, reverse=True)
        row['topmotifs'] = list(OrderedDict.fromkeys(mcs))[0:top]

    topmotifs['ComponentIndex'] = topmotifs.index

    edges['TopSharedMotifs'] = edges['ComponentIndex'].map(topmotifs.set_index('ComponentIndex')['topmotifs'])

    # add separate edge for each shared motif
    edges.insert(loc=1, column='interaction', value= 'cosine')

    cols = edges.columns
    lst = []

    for index, row in edges.iterrows():
        if row['shared_motifs'] != 'None':
            for idx,val in enumerate(row['shared_motifs']):
                if 'EdgeAnnotation' not in edges.columns:
                	lst.append([row['CLUSTERID1'],row['shared_motifs'][idx],row['CLUSTERID2'], 
                        row['DeltaMZ'], row['MEH'], row['Cosine'], row['OtherScore'],row['ComponentIndex'],row['shared_motifs'], row['TopSharedMotifs']])
                else:
                	lst.append([row['CLUSTERID1'],row['shared_motifs'][idx],row['CLUSTERID2'], 
                        row['DeltaMZ'], row['MEH'], row['Cosine'], row['OtherScore'],row['ComponentIndex'],row['EdgeAnnotation'],row['shared_motifs'], row['TopSharedMotifs']])


    motifedges = pd.DataFrame(lst, columns=cols)
    edges = edges.append(motifedges)
    
    return {'nodes':comb,'edges':edges}

def make_inchidic(smilesdic):
    """Convert a dictionary of SMILES to a dictionary of InChIKeys

    :param smilesdic: A dictionary of SMILES
    :type smilesdic: dict
    :return: A dictionary of InChIKeys
    :rtype: InChIKeys

    """
    inchi_dic = smilesdic['dic'].copy()
    d = {k: v for k, v in zip(smilesdic['df'].SMILES, smilesdic['df'].inchikey)}
    for k in inchi_dic:
        inchi_dic[k] = [d.get(k, k) for k in inchi_dic[k] if d.get(k, k) != k]
    
    return inchi_dic

def unique_smiles(matches):
    """Retrieve overall unique SMILES and unique SMILES per molecular feature 

    :param matches: A list of dataframes with feature IDs and corresponding SMILES, where each dataframe can correspond to a different source of chemical structural annotation (e.g. GNPS library matches, in silico strucutral prediction through NAP, Dereplicator or SIRIUS+CSI:FingerID). Feature IDs need to correspond in all dataframes as well as the mass spectral molecular network.
    :type matches: list
    :return: A dictionary containing a dataframe with doverall unique SMILES and a dictionary with unique SMILES per molecular feature
    :rtype: dict

    """
    # combine SMILES for same feature into one string of features
    for index, item in enumerate(matches):
        if 'Scan' in matches[index].columns:
            matches[index] = matches[index].groupby('Scan', as_index=False).agg(lambda x: ','.join(set(x.dropna()))) 
            matches[index] = matches[index].rename(columns = {'Scan':'cluster.index'})
        if '#Scan#' in matches[index].columns:
        	matches[index] = matches[index].groupby('#Scan#', as_index=False).agg(lambda x: ','.join(set(x.dropna())))
        	matches[index] = matches[index].rename(columns = {'#Scan#':'cluster.index'})
            
    comb = reduce(lambda left,right: pd.merge(left,right,on='cluster.index', how = "outer"), matches)
    if 'FusionSMILES' in comb.columns:
        comb = comb.drop(['FusionSMILES', 'ConsensusSMILES'], axis=1)
    
    # concatenate all SMILES
    comb['AllSmiles'] = comb.filter(regex='^.*(Smiles|SMILES).*$').apply(lambda x: ','.join([y for y in x.tolist() if str(y) != 'nan']), axis=1)
    
    # create dictionary of SMILES
    comb_dic = comb.set_index('cluster.index')['AllSmiles'].to_dict()
    
    # split comma separated strings of SMILES 
    for k in comb_dic:
        comb_dic[k] = comb_dic[k].split(',')
    
    # retrieve unique SMILES per feature
    for i in comb_dic:
        comb_dic[i] = list(set(comb_dic[i]))
        
    # remove empty values
    comb_dic = {k: comb_dic[k] for k in comb_dic if not comb_dic[k] == ['']} 
    comb_dic = {k: comb_dic[k] for k in comb_dic if not comb_dic[k] == []}
    comb_dic = {k: comb_dic[k] for k in comb_dic if not comb_dic[k] == [' ']}
    
    # convert dictionary into list of unique SMILES
    l = list(set([item for sublist in list(comb_dic.values()) for item in sublist]))
    l = [x for x in l if x not in ['',' ',None]]
    
    # convert list into dataframe
    df = pd.DataFrame({"SMILES": l})
    
    # remove white space from SMILES
    df.SMILES = df.SMILES.str.replace(' ', '')
    
    return {'df':df, 'dic':comb_dic}

def highestscore(a, chem_dic, score):
    """Retrieve most predominant chemical class per componentindex at a single level of the ClassyFire chemical ontology

    :param a: list of all cluster indexes per componentindex
    :type a: list
    :param chem_dic: A dictionary of chemical classes at a single level of the ClassyFire chemical ontology
    :type chem_dic: dict
    :param score: A list of number of nodes per compontentindex
    :type score: list
    :return: List of componentindexes with name and score of the most predominant chemical class at a signle level of the ClassyFire chemical ontology
    :rtype: list

    """
    # creates a list of lists, where each list item corresponds to one componentindex
    # each sublist item corresponds to the chemical classes of the SMILES matched to one node
    chem_ci = []
    for i in a:
        chem_ci.append([chem_dic[x] for x in i])

    # Creates a list of dictionaries, where each list item corresponds to one componentindex and each
    # dictionary corresponds to the chemical categories found within this node with a corresponding score.
    #
    # This score sums up to 1 per node and corresponds to the number of SMILES associated with this class 
    # divided by the total number of SMILES per node.
    chem_scores = [] 
    for i in chem_ci:
        subl = []
        for j in i:
            dic_a = dict(Counter(j))
            out = {k: v / total for total in (sum(dic_a.values()),) for k, v in dic_a.items()}
            subl.append(out)
        chem_scores.append(subl)

    # List of dictionaries, each list item corresponds to one componentindex and dictionaries contain
    # summed values per chemical category per component index.
    #
    # E.g. Organic acids and derivatives': 2.1246753246753247 means that 2.12 nodes
    # in this molecular family were classified as Organic acids and derivatives
    chem_sums = []
    for i in chem_scores:
        if i:
            sums = functools.reduce(operator.add, map(Counter, i))
        if not i:
            sums = []
        chem_sums.append(sums)

    # List of items, where each list item corresponds to one componentindex and name and score of the most highly abundant
    # chemical class are found in it.
    # 
    # E.g. Lipids and 0.75 means that 0.75 of the molecular family was associated to lipids, 
    # this could be the case if e.g. Lipids score is 2.25 (nodes) and the molecular family consits of a total of 3 nodes.
    chem_finalscore = []
    for index, item in enumerate(chem_sums): 
        if chem_sums[index]:
            char = max(chem_sums[index].items(), key=operator.itemgetter(1))
            num = char[1]/score[index]
        if not chem_sums[index]:
            char = ["no matches","no matches"]
            num = ""
        chem_finalscore.append([char[0],num])
    
    return chem_finalscore

def molfam_classes(net, df, smilesdict):
    """Retrieve most predominant chemical class for each level of the ClassyFire chemical ontology

    :param net: GNPS network data
    :type net: pandas.core.frame.DataFrame
    :param df: A dataframe comprising all unique SMILES, InChIKeys and corresponding chemical classes at each level of the ClassyFire chemical ontology
    :type df: pandas.core.frame.DataFrame
    :param smilesdict: A dictionary of nodes with corresponding unique SMILES
    :type smilesdict: dict
    :return: A dataframe containing most predominant chemical classes per node at each level of the ClassyFire chemical ontology
    :rtype: pandas.core.frame.DataFrame

    """
    # create dictionaries for each hierarchical level of Classyfire
    
    # rename componentindex of selfloops, so they are considered independently of each other
    selfs = list(range(1,len(net.componentindex[net.componentindex == -1])+1))
    net.loc[net.componentindex == -1,['componentindex']] = ["S" + str(s) for s in selfs]


    kingdom_dic = smilesdict.copy()
    d = {k: v for k, v in zip(df.inchikey, df.kingdom)}
    for k in kingdom_dic:
            kingdom_dic[k] = [d.get(k, k) for k in kingdom_dic[k] if d.get(k, k) != k]

    # Superclass
    superclass_dic = smilesdict.copy()
    d = {k: v for k, v in zip(df.inchikey, df.superclass)}
    for k in superclass_dic:
            superclass_dic[k] = [d.get(k, k) for k in superclass_dic[k] if d.get(k, k) != k]
            
    # Class
    class_dic = smilesdict.copy()
    d = {k: v for k, v in zip(df.inchikey, df.CF_class)}
    for k in class_dic:
            class_dic[k] = [d.get(k, k) for k in class_dic[k] if d.get(k, k) != k]
    
    # Sublass
    subclass_dic = smilesdict.copy()
    d = {k: v for k, v in zip(df.inchikey, df.subclass)}
    for k in subclass_dic:
            subclass_dic[k] = [d.get(k, k) for k in subclass_dic[k] if d.get(k, k) != k]
            
    # Direct parent
    Dparent_dic = smilesdict.copy()
    d = {k: v for k, v in zip(df.inchikey, df.direct_parent)}
    for k in Dparent_dic:
            Dparent_dic[k] = [d.get(k, k) for k in Dparent_dic[k] if d.get(k, k) != k]
            
    # Molecular framework
    MFramework_dic = smilesdict.copy()
    d = {k: v for k, v in zip(df.inchikey, df.molecular_framework)}
    for k in MFramework_dic:
            MFramework_dic[k] = [d.get(k, k) for k in MFramework_dic[k] if d.get(k, k) != k]
    
    # iterate over all cluster indexes
    ci = list(set(net['componentindex']))
    
    # get list of all cluster indexes per componentindex
    a = []
    for i in ci:
        a.append(list(net.loc[net['componentindex'] == i, 'cluster index']))
        
    # remove all items that do not have any SMILES (are not contained in keys of smilesdict)
    # score retrieves number of nodes per cluster index
    score = []
    for index, item in enumerate(a): 
        red = list(set(list(smilesdict.keys())).intersection(a[index]))
        score.append(len(a[index]))
        a[index] = red 
    
    kingdom_finalscore = highestscore(a, kingdom_dic, score)
    superclass_finalscore = highestscore(a, superclass_dic, score)
    class_finalscore = highestscore(a, class_dic, score)
    subclass_finalscore = highestscore(a, subclass_dic, score)
    Dparent_finalscore = highestscore(a, Dparent_dic, score)
    MFramework_finalscore = highestscore(a, MFramework_dic, score)
   
    kingdom = [item[0] for item in kingdom_finalscore]
    superclass = [item[0] for item in superclass_finalscore]
    CF_class = [item[0] for item in class_finalscore]
    subclass = [item[0] for item in subclass_finalscore]
    Dparent = [item[0] for item in Dparent_finalscore]
    MFramework = [item[0] for item in MFramework_finalscore]
    kingdom_score = [item[1] for item in kingdom_finalscore]
    superclass_score =[item[1] for item in superclass_finalscore]
    CF_class_score =[item[1] for item in class_finalscore]
    subclass_score =[item[1] for item in subclass_finalscore]
    Dparent_score =[item[1] for item in Dparent_finalscore]
    MFramework_score =[item[1] for item in MFramework_finalscore]
    
    data_tuples = list(zip(ci,score,kingdom,kingdom_score,superclass,superclass_score,CF_class,CF_class_score,subclass,subclass_score,Dparent,Dparent_score,MFramework,MFramework_score))
    sumary = pd.DataFrame(data_tuples, columns=['componentindex','CF_NrNodes','CF_kingdom','CF_kingdom_score','CF_superclass','CF_superclass_score','CF_class','CF_class_score','CF_subclass','CF_subclass_score','CF_Dparent','CF_Dparent_score','CF_MFramework','CF_MFramework_score'])
        
    final = pd.merge(net[['componentindex','cluster index']], sumary, on='componentindex')
    
    # make cluster index first column
    final = final[['cluster index','componentindex','CF_NrNodes','CF_kingdom','CF_kingdom_score','CF_superclass','CF_superclass_score','CF_class','CF_class_score','CF_subclass','CF_subclass_score','CF_Dparent','CF_Dparent_score','CF_MFramework','CF_MFramework_score']]
    final = final.rename(columns = {'componentindex':'CF_componentindex'})
    
    return final
    
def make_classyfire_graphml(graphML,final):
    """Create a network file with chemical classes mapped

    :param graphml: GNPS network file
    :type graphml: networkx.classes.graph.Graph
    :param final: A dataframe containing most predominant chemical classes per node at each level of the ClassyFire chemical ontology
    :type final: pandas.core.frame.DataFrame
    :return: A network file with most predominant chemical classes per node mapped at each level of the ClassyFire chemical ontology
    :rtype: networkx.classes.graph.Graph

    """
    for v in graphML.nodes():  

        graphML.node[v]['CF_componentindex']= str(final[final['cluster index'] == int(v)]['CF_componentindex'].iloc[-1])
        graphML.node[v]['CF_NrNodes']= float(final[final['cluster index'] == int(v)]['CF_NrNodes'].iloc[-1])
        graphML.node[v]['CF_kingdom']= str(final[final['cluster index'] == int(v)]['CF_kingdom'].iloc[-1])

        try:
            graphML.node[v]['CF_kingdom_score']= float(final[final['cluster index'] == int(v)]['CF_kingdom_score'].iloc[-1])
        except:
            graphML.node[v]['CF_kingdom_score']= object()

        graphML.node[v]['CF_superclass']= str(final[final['cluster index'] == int(v)]['CF_superclass'].iloc[-1])

        try:
            graphML.node[v]['CF_superclass_score']= float(final[final['cluster index'] == int(v)]['CF_superclass_score'].iloc[-1])
        except:
            graphML.node[v]['CF_superclass_score']= object()

        graphML.node[v]['CF_class']= str(final[final['cluster index'] == int(v)]['CF_class'].iloc[-1])

        try:
            graphML.node[v]['CF_class_score']= float(final[final['cluster index'] == int(v)]['CF_class_score'].iloc[-1])
        except:
            graphML.node[v]['CF_class_score']= object()

        graphML.node[v]['CF_subclass']= str(final[final['cluster index'] == int(v)]['CF_subclass'].iloc[-1])

        try:
            graphML.node[v]['CF_subclass_score']= float(final[final['cluster index'] == int(v)]['CF_subclass_score'].iloc[-1])
        except:
            graphML.node[v]['CF_subclass_score']= object()

        graphML.node[v]['CF_Dparent']= str(final[final['cluster index'] == int(v)]['CF_Dparent'].iloc[-1])

        try:
            graphML.node[v]['CF_Dparent_score']= float(final[final['cluster index'] == int(v)]['CF_Dparent_score'].iloc[-1])
        except:
            graphML.node[v]['CF_Dparent_score']= object()

        graphML.node[v]['CF_MFramework']= str(final[final['cluster index'] == int(v)]['CF_MFramework'].iloc[-1])

        try:
            graphML.node[v]['CF_MFramework_score']= float(final[final['cluster index'] == int(v)]['CF_MFramework_score'].iloc[-1])
        except ValueError:
            graphML.node[v]['CF_MFramework_score']= object()
        
    return graphML
    
def make_motif_graphml(nodes, edges):
    """Create a network file with Mass2Motifs mapped on nodes and shared Mass2Motifs mapped as multiple edges

    :param nodes: A dataframe showing Mass2Motifs per node
    :type nodes: pandas.core.frame.DataFrame
    :param edges: A dataframe showing shared Mass2Motifs for each network pair
    :type edges: pandas.core.frame.DataFrame
    :return: A network file with Mass2Motifs mapped on nodes and shared Mass2Motifs mapped as multiple edges
    :rtype: networkx.classes.graph.Graph

    """
    # convert lists to strings
    edges['shared_motifs'] = edges['shared_motifs'].replace('None', '')
    edges['TopSharedMotifs'] = edges['TopSharedMotifs'].replace('None', '')
    edges['shared_motifs'] = edges['shared_motifs'].agg(lambda x: ','.join(map(str, x)))
    edges['TopSharedMotifs'] = edges['TopSharedMotifs'].agg(lambda x: ','.join(map(str, x)))
    
    # create motif network with multiple edges
    MG = nx.from_pandas_edgelist(edges, 'CLUSTERID1', 'CLUSTERID2', edge_attr=list(set(list(edges.columns)) - set(['CLUSTERID1','CLUSTERID2'])), 
                             create_using=nx.MultiGraph())
    
    # map node attributes to network
    nodes['precursormass'] = nodes['precursormass'].agg(lambda x: ','.join(map(str, x)))
    nodes['parentrt'] = nodes['parentrt'].agg(lambda x: ','.join(map(str, x)))
    nodes['document'] = nodes['document'].agg(lambda x: ','.join(map(str, x)))
    nodes['motif'] = nodes['motif'].agg(lambda x: ','.join(map(str, x)))
    nodes['probability'] = nodes['probability'].agg(lambda x: ','.join(map(str, x)))
    nodes['overlap'] = nodes['overlap'].agg(lambda x: ','.join(map(str, x)))
    
    for column in nodes:
        nx.set_node_attributes(MG, pd.Series(nodes[column], index=nodes.index).to_dict(), column)
        
    return MG
    
"""
@author: Ricardo Silva (https://github.com/rsilvabioinfo)
"""

def make_classy_table(jsondic):  
    
    dmetadatalist = []
    
    for idx,entry in enumerate(jsondic):
        
        mdict = {}
        
        if (entry != {} and entry != None):
    
            if sum([bool(re.match('smiles', x)) for x in entry.keys()]) > 0 and entry['smiles'] is not None:
                mdict['smiles'] = entry['smiles']
            if sum([bool(re.match('inchikey', x)) for x in entry.keys()]) > 0 and entry['inchikey'] is not None:
                mdict['inchikey'] = entry['inchikey']
            if sum([bool(re.match('kingdom', x)) for x in entry.keys()]) > 0 and entry['kingdom'] is not None:
                mdict['kingdom'] = entry['kingdom']['name']
            if sum([bool(re.match('superclass', x)) for x in entry.keys()]) > 0 and entry['superclass'] is not None :
                mdict['superclass'] = entry['superclass']['name']
            if sum([bool(re.match('class', x)) for x in entry.keys()]) > 0 and entry['class'] is not None :
                mdict['class'] = entry['class']['name']
            if sum([bool(re.match('subclass', x)) for x in entry.keys()]) > 0 and entry['subclass'] is not None:
                mdict['subclass'] = entry['subclass']['name']
            if sum([bool(re.match('direct_parent', x)) for x in entry.keys()]) > 0 and entry['direct_parent'] is not None:
                mdict['direct_parent'] = entry['direct_parent']['name']
            if sum([bool(re.match('molecular_framework', x)) for x in entry.keys()]) > 0 and entry['molecular_framework'] is not None:
                mdict['molecular_framework'] = entry['molecular_framework']
        else:
            
            mdict['smiles'] = "None"
            mdict['inchikey'] = "None"
            mdict['kingdom'] = "None"
            mdict['superclass'] = "None"
            mdict['class'] = "None"
            mdict['subclass'] = "None"
            mdict['direct_parent'] = "None"
            mdict['molecular_framework'] = "None"
            
        dmetadatalist.append(mdict)
    
    df_metares = pd.DataFrame.from_dict(dmetadatalist)  
    return(df_metares)
    
"""
@author: Ming Wang (https://github.com/mwang87)
"""

def get_structure_class_entity(inchikey):
    print(inchikey)
    entity = None
    try:
        entity = json.loads(get_entity(inchikey))
    except KeyboardInterrupt:
        raise
    except:
        new_inchi_no_stereo = inchikey.split("-")[0] + "-UHFFFAOYSA-N"
        try:
            entity = json.loads(get_entity(new_inchi_no_stereo))
        except KeyboardInterrupt:
            raise
        except:
            return entity

    return entity


def get_structure_class(inchikey):
    return_dict = {}
    print(inchikey)
    try:
        entity = json.loads(get_entity(inchikey))
    except KeyboardInterrupt:
        raise
    except:
        new_inchi_no_stereo = inchikey.split("-")[0] + "-UHFFFAOYSA-N"
        try:
            entity = json.loads(get_entity(new_inchi_no_stereo))
        except:
            print(inchikey, "Not Cached")
            return_dict["inchikey"] = inchikey
            return_dict["superclass"] = "None"
            return_dict["class"] = "None"
            return_dict["subclass"] = "None"
            return return_dict

    try:
        print(json.dumps(entity))
        return_dict["inchikey"] = inchikey
        if "superclass" in entity:
            return_dict["superclass"] = entity["superclass"]["name"]
        else:
            return_dict["superclass"] = "None"

        if "class" in entity:
            return_dict["class"] = entity["class"]["name"]
        else:
            return_dict["class"] = "None"

        if "subclass" in entity and entity["subclass"] != None:
            return_dict["subclass"] = entity["subclass"]["name"]
        else:
            return_dict["subclass"] = "None"

        return return_dict
    except:
        return_dict["inchikey"] = inchikey
        return_dict["superclass"] = "None"
        return_dict["class"] = "None"
        return_dict["subclass"] = "None"
        return return_dict

"""A client for the ClassyFire API which enables efficient querying with
 chemical database files"""

url = "http://classyfire.wishartlab.com"
chunk_size = 1000
sleep_interval = 60


def structure_query(compound, label='pyclassyfire'):
    """Submit a compound information to the ClassyFire service for evaluation
    and receive a id which can be used to used to collect results

    :param compound: The compound structures as line delimited inchikey or
         smiles. Optionally a tab-separated id may be prepended for each
         structure.
    :type compound: str
    :param label: A label for the query
    :type label:
    :return: A query ID number
    :rtype: int

    >>> structure_query('CCC', 'smiles_test')
    >>> structure_query('InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)')

    """
    r = requests.post(url + '/queries.json', data='{"label": "%s", '
                      '"query_input": "%s", "query_type": "STRUCTURE"}'
                                                  % (label, compound),
                      headers={"Content-Type": "application/json"})
    r.raise_for_status()
    return r.json()['id']


def iupac_query(compound, label='pyclassyfire'):
    """Submit a IUPAC compound name to the ClassyFire service for evaluation
     and receive a id which can be used to used to collect results.

    :param compound: The line delimited compound names. Optionally a
         tab-separated id may be prepended for each compound.
    :type compound: str
    :param label: A label for the query
    :type label:
    :return: A query ID number
    :rtype: int

    >>> iupac_query('ethane', 'iupac_test')
    >>> iupac_query('C001\\tethane\\nC002\\tethanol', 'iupac_test')

    """
    r = requests.post(url + '/queries.json', data='{"label": "%s", '
                      '"query_input": "%s", "query_type": "IUPAC_NAME"}'
                                                  % (label, compound),
                      headers={"Content-Type": "application/json"})
    r.raise_for_status()
    return r.json()['id']


def get_results(query_id, return_format="json", blocking=False):
    """Given a query_id, fetch the classification results.

    :param query_id: A numeric query id returned at time of query submission
    :type query_id: str
    :param return_format: desired return format. valid types are json, csv or sdf
    :type return_format: str
    :return: query information
    :rtype: str

    >>> get_results('595535', 'csv')
    >>> get_results('595535', 'json')
    >>> get_results('595535', 'sdf')

    """
    if blocking == False:
        r = requests.get('%s/queries/%s.%s' % (url, query_id, return_format),
                         headers={"Content-Type": "application/%s" % return_format})
        r.raise_for_status()
        return r.text
    else:
        while True:
            r = requests.get('%s/queries/%s.%s' % (url, query_id, return_format),
                             headers={"Content-Type": "application/%s" % return_format})
            r.raise_for_status()
            result_json = r.json()
            if result_json["classification_status"] != "In Queue":
                return r.text
            else:
                print("WAITING")
                time.sleep(10)


def get_entity(inchikey, return_format="json"):
    """Given a InChIKey for a previously queried structure, fetch the
     classification results.

    :param inchikey: An InChIKey for a previously calculated chemical structure
    :type inchikey: str
    :param return_format: desired return format. valid types are json, csv or sdf
    :type return_format: str
    :return: query information
    :rtype: str

    >>> get_entity("ATUOYWHBWRKTHZ-UHFFFAOYSA-N", 'csv')
    >>> get_entity("ATUOYWHBWRKTHZ-UHFFFAOYSA-N", 'json')
    >>> get_entity("ATUOYWHBWRKTHZ-UHFFFAOYSA-N", 'sdf')

    """
    inchikey = inchikey.replace('InChIKey=', '')
    r = requests.get('%s/entities/%s.%s' % (url, inchikey, return_format),
                     headers={
                         "Content-Type": "application/%s" % return_format})
    r.raise_for_status()
    return r.text


def get_chemont_node(chemontid):
    """Return data for the TaxNode with ID chemontid.

    :param chemontid: the ChemOnt ID of the entity.
    :type chemontid: str
    :return: The classification results for the entity as json.
    :rtype: str

    >>> get_chemont_node('CHEMONTID:0004253')

    """
    chemontid = chemontid.replace("CHEMONTID:", "C")
    r = requests.get('%s/tax_nodes/%s.json' % (url, chemontid),
                     headers={"Content-Type": "application/json" })
    r.raise_for_status()
    return r.text


def tabular_query(inpath, structure_key, dialect='excel', outpath=None,
                  outfields=('taxonomy', 'description', 'substituents')):
    """Given a path to a compound set in tabular form (comma or tab delimited)
     annotate all compounds and write results to an expanded table.

    :param inpath: path to compound file to be annotated
    :type inpath: str
    :param structure_key: column heading which contains the compounds InChIKey
         or SMILES
    :type structure_key: str
    :param dialect: dialect for parsing table (generally 'excel' for csv,
         'excel-tab' for tsv)
    :type dialect: str
    :param outpath: Path to desired output location
    :type outpath: str
    :param outfields: Fields to append to table from ClassyFire output
    :type outfields: tuple(string)

    >>> tabular_query('/tabulated_data.tsv', 'structure', 'excel-tab')

    """
    tax_fields = ('kingdom', 'superclass', 'class', 'subclass')
    query_ids = []
    infile = open(inpath, 'rU')
    if not outpath:
        outpath = _prevent_overwrite(inpath)
    comps = []
    for line in csv.DictReader(infile, dialect=dialect):
        comps.append(line[structure_key])
        if not len(comps) % chunk_size:
            query_ids.append(structure_query('\\n'.join(comps)))
            comps = []
    if comps:
        query_ids.append(structure_query('\\n'.join(comps)))
    print('%s queries submitted to ClassyFire API' % len(query_ids))
    i = 0
    infile.seek(0)
    with open(outpath, 'w') as outfile:
        reader = csv.DictReader(infile, dialect=dialect)
        writer = csv.DictWriter(outfile, reader.fieldnames+list(outfields),
                                dialect=dialect)
        writer.writeheader()
        while i < len(query_ids):
            result = json.loads(get_results(query_ids[i]))
            if result["classification_status"] == "Done":
                for j, line in enumerate(reader):
                    if result['entities'] and str(j+1) == result['entities'][0]['identifier'].split('-')[1]:
                        hit = result['entities'].pop(0)
                        if 'taxonomy' in outfields:
                            hit['taxonomy'] = ";".join(
                                ['%s:%s' % (hit[x]['name'], hit[x]['chemont_id'])
                                 for x in tax_fields if hit[x]])
                        for field in outfields:
                            if isinstance(hit[field], list):
                                line[field] = ';'.join(hit[field])
                            else:
                                line[field] = hit[field]
                    writer.writerow(line)
                i += 1
            else:
                print("%s percent complete" % round(i/len(query_ids)*100))
                time.sleep(sleep_interval)
    infile.close()


def sdf_query(inpath, outpath=None):
    """Given a path to a compound set in a sdf file, annotate all compounds
     and write results as attributes in a sdf file.

    :param inpath: path to compound file to be annotated
    :type inpath: str
    :param outpath: Path to desired output location
    :type outpath: str

    >>> sdf_query('/sdf_data.sdf')

    """
    from rdkit.Chem import AllChem
    query_ids = []
    if not outpath:
        outpath = _prevent_overwrite(inpath)
    comps = []
    for mol in AllChem.SDMolSupplier(inpath):
        if mol:
            comps.append(AllChem.MolToSmiles(mol))
        if not len(comps) % chunk_size:
            query_ids.append(structure_query('/n'.join(comps)))
            comps = []
    if comps:
        query_ids.append(structure_query('\\n'.join(comps)))
    print('%s queries submitted to ClassyFire API' % len(query_ids))
    i = 0
    with open(outpath, 'w') as outfile:
        while i < len(query_ids):
            result = json.loads(get_results(query_ids[i]))
            if result["classification_status"] == "Done":
                outfile.write(get_results(query_ids[i], return_format='sdf'))
                i += 1
            else:
                print("%s percent complete" % round(i / len(query_ids) * 100))
                time.sleep(sleep_interval)


def _prevent_overwrite(write_path, suffix='_annotated'):
    """Prevents overwrite of existing output files by appending a suffix when
    needed

    :param write_path: potential write path
    :type write_path: string
    :return:
    :rtype:
    """
    while os.path.exists(write_path):
        sp = write_path.split('.')
        if len(sp) > 1:
            sp[-2] += suffix
            write_path = '.'.join(sp)
        else:
            write_path += suffix
    return write_path

def run_shell_command(script_to_run):
	os.system(script_to_run)
	return "DONE"

# Wraps running in parallel a set of shell scripts
def run_parallel_shellcommands(input_shell_commands, parallelism_level):
	return run_parallel_job(run_shell_command, input_shell_commands, parallelism_level)

# Wraps the parallel job running, simplifying code
def run_parallel_job(input_function, input_parameters_list, parallelism_level):
	if parallelism_level == 1:
		output_results_list = []
		for input_param in input_parameters_list:
			result_object = input_function(input_param)
			output_results_list.append(result_object)
		return output_results_list
	else:
		results = Parallel(n_jobs = parallelism_level)(delayed(input_function)(input_object) for input_object in input_parameters_list)
		return results
		
def get_classifications(inchifile):
    with open(inchifile) as csvfile:
        all_inchi_keys = []
    
        reader = csv.DictReader(csvfile)
        row_count = 0
        for row in reader:
            row_count += 1
    
            if row_count % 1000 == 0:
                print(row_count)
    
            all_inchi_keys.append(row["InChIKey"].split("=")[1])
    
            continue
    
        #all_inchi_keys = all_inchi_keys[-1000:]
        all_json = run_parallel_job(get_structure_class_entity, all_inchi_keys, 50)
    
        open("all_json.json", "w").write(json.dumps(all_json))
