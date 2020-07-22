#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 16:01:33 2018

@author: mgerard
"""

import re
import time
import codecs
import argparse
import os


class BioCyc(object):
    '''
    '''
    
    #=================================================
    def __init__(self, DB_cpd, DB_rxn):
        '''
        '''
        
        path,filename = os.path.split(DB_rxn)
        
        self._path = [path, filename]
        
        #------------------------------------
        self._id = ''
        self._version = ''
        self._cpd = dict()
        self._rxn = dict()
        
        #------------------------------------
        
        self._cpd = self.load_compounds_from_file(DB_cpd)
        
        species,version,REACTIONS = self.load_reactions_from_file(DB_rxn)
        
        self._id = species
        self._version = version
        self._rxn = REACTIONS
        
    
    #=================================================
    def load_compounds_from_file(self, DB_cpd):
        '''
        '''
        with codecs.open(DB_cpd,'r', "ISO-8859-1") as fp:
            DATA = fp.read()
        
        #===============================
        # EXTRACT METADATA AND DATA
        #===============================
        metadata = DATA[:DATA.find('UNIQUE-ID -')-1]
        data = DATA[DATA.find('UNIQUE-ID -'):]
        
        #==============================================
        # EXTRACT GENERAL INFORMATION AND ATTRIBUTES
        #==============================================
        Info,Attributes = metadata.split('Attributes:\n')
                
        attributes = Attributes.replace(' ','').replace('#','').split('\n')[:-1] # UNIQUE-ID, EC-NUMBER, SPONTANEOUS?, etc
        
        #==============================================
        # EXCLUDING SOME ATTRIBUTES
        #==============================================
        attributes.remove('COMMON-NAME')
        attributes.remove('DBLINKS')
        attributes.remove('SYNONYMS')
        
        
        compounds = data.strip().split('//\n')
        
        COMPOUNDS = dict()
        for compound in compounds[:-1]:
            
            unique_id = re.findall('UNIQUE-ID - (.{1,100}?)\n',compound)[0]
            
            COMPOUNDS[unique_id] = self.initialize_compound()
            
            #-------------------------
            # EXTRACTING COMMON-NAME
            #-------------------------
            COMMON_NAME = re.findall('COMMON-NAME - (.{1,100}?)\n',compound)
            if COMMON_NAME:
                COMPOUNDS[unique_id]['name'] = COMMON_NAME
            
            
            #-------------------------
            # EXTRACTING SYNONIMS
            #-------------------------
            SYNONYMS = re.findall('SYNONYMS - (.{1,100}?)\n',compound)
            if SYNONYMS:
                COMPOUNDS[unique_id]['name'].extend(SYNONYMS)
            
            
            #-------------------------------------------
            # FIXING NAMES
            #-------------------------------------------
            T = ['<i>', '</i>', '<sub>', '</sub>', '<SUB>', '</SUB>', '<sup>', '</sup>', '<SUP>', '</SUP>']
            
            for idx in range(len(COMPOUNDS[unique_id]['name'])):
                for t in T:
                    COMPOUNDS[unique_id]['name'][idx] = COMPOUNDS[unique_id]['name'][idx].replace(t,'')
                    
        
        return COMPOUNDS
    
    
    #=================================================
    def load_reactions_from_file(self, DB_rxn):
        '''
        '''
        
        with codecs.open(DB_rxn,'r', "ISO-8859-1") as fp:
            DATA = fp.read()
        
        
        #-----------------------------
        # EXTRACT METADATA AND DATA
        #-----------------------------
        metadata = DATA[:DATA.find('UNIQUE-ID -')-1]
        data = DATA[DATA.find('UNIQUE-ID -'):]
        
        
        #------------------------------------------------
        # EXTRACT GENERAL INFORMATION AND ATTRIBUTES
        #------------------------------------------------
        Info,Attributes = metadata.split('Attributes:\n')
                
        species = re.findall('# Species: (.{1,50})\n',metadata)[0]
        version = re.findall('# Version: (.{1,20})\n',metadata)[0]
        
        attributes = Attributes.replace(' ','').replace('#','').split('\n')[:-1] # UNIQUE-ID, EC-NUMBER, SPONTANEOUS?, etc
        
        
        #--------------------------------
        # EXCLUDING SOME ATTRIBUTES
        #----------------------------
        attributes.remove('LEFT')
        attributes.remove('RIGHT')
        attributes.remove('EC-NUMBER')
        
        
        #------------------------------------------
        # SPLIT data into indivudual REACTIONS
        #---------------------------------------
        reactions = data.strip().split('//\n')
        
        REACTIONS = dict()
        
#        REACTIONS['id'] = species
#        REACTIONS['version'] = version
        #-------------------------------------------------------
        for reaction in reactions:
            
            unique_id = re.findall('UNIQUE-ID - (.{1,100}?)\n',reaction)[0]
            
            REACTIONS[unique_id] = self.initialize_reaction()
            
            REACTIONS[unique_id]['enzymes'] = re.findall('EC-NUMBER - EC-(\d.\d{1,3}.\d{1,3}.{0,1}\d{0,3})\n',reaction)
            
            
            #-----------------------------------------
            # EXTRACTING SUBSTRATES AND PRODUCTS
            #--------------------------------------
            LEFT = re.findall('LEFT - (.{1,100}?)\n\^COEFFICIENT - (\d{0,3})\n|LEFT - (.{1,100}?)\n',reaction)
            
            RIGHT = re.findall('RIGHT - (.{1,100}?)\n\^COEFFICIENT - (\d{0,3})\n|RIGHT - (.{1,100}?)\n',reaction)
            
            
            #------------------------------------
            # FILETRING EMPTY REACTIONS
            #----------------------------
            if LEFT:
                
                for left in LEFT:
                    
                    if left[2] != '':
                        REACTIONS[unique_id]['Scoef'].append(1)
                        REACTIONS[unique_id]['substrates'].append(left[2])
                    else:
                        REACTIONS[unique_id]['substrates'].append(left[0])
                        REACTIONS[unique_id]['Scoef'].append(int(left[1]))
                
                # SUBSTRATES
                for coef,S in zip(REACTIONS[unique_id]['Scoef'],REACTIONS[unique_id]['substrates']):
                    
                    REACTIONS[unique_id]['equation'] += str(coef) + ' ' + S + ' + '
    
                
            
            if RIGHT:
                
                for right in RIGHT:
                    if right[2] != '':
                        REACTIONS[unique_id]['Pcoef'].append(1)
                        REACTIONS[unique_id]['products'].append(right[2])
                    else:
                        REACTIONS[unique_id]['Pcoef'].append(int(right[1]))
                        REACTIONS[unique_id]['products'].append(right[0])
                
                
                if LEFT:
                    REACTIONS[unique_id]['equation'] = REACTIONS[unique_id]['equation'][:-3] + ' --> '
                else:
                    REACTIONS[unique_id]['equation'] = ' --> '
                
                # PRODUCTS
                for coef,P in zip(REACTIONS[unique_id]['Pcoef'],REACTIONS[unique_id]['products']):
                    
                    REACTIONS[unique_id]['equation'] += str(coef) + ' ' + P + ' + '
                
            
            
            REACTIONS[unique_id]['equation'] = REACTIONS[unique_id]['equation'][:-3]
            #======================================
            
            
            #---------------------------
            # IS A TRANSPORT REACTION?
            #---------------------------
            if 'TYPES - Transport-Reactions' in reaction:
                REACTIONS[unique_id]['transport'] = True
                
            if not REACTIONS[unique_id]['substrates'] or not REACTIONS[unique_id]['products']:
                REACTIONS[unique_id]['transport'] = True
            
            if set(REACTIONS[unique_id]['substrates']) == set(REACTIONS[unique_id]['products']):
                REACTIONS[unique_id]['transport'] = True
            
                
            
            #========================
            # REVERSIBILITY
            #==================
            
            REVERSIBILITY = re.findall('REACTION-DIRECTION - (PHYSIOL-LEFT-TO-RIGHT|PHYSIOL-RIGHT-TO-LEFT|LEFT-TO-RIGHT|RIGHT-TO-LEFT|REVERSIBLE)\n',reaction)
            
            if REVERSIBILITY:
                
                REVERSIBILITY = REVERSIBILITY[0]
                
                if (REVERSIBILITY == 'PHYSIOL-LEFT-TO-RIGHT') or (REVERSIBILITY == 'LEFT-TO-RIGHT'):
                    REACTIONS[unique_id]['reversibility'] = 1
                    
                elif (REVERSIBILITY == 'PHYSIOL-RIGHT-TO-LEFT') or (REVERSIBILITY == 'RIGHT-TO-LEFT'):
                    REACTIONS[unique_id]['reversibility'] = -1
                    
                else:
                    REACTIONS[unique_id]['reversibility'] = 0
            
        
        
        return (species,version,REACTIONS)
    
    
    #=================================================
    def initialize_reaction(self):
        
        #if (rxn_code not in self._rxn.keys() ) or (overwrite):
         reaction = {'name': [],
                     'identifiers': {'KEGG': []},
                     'enzymes': [],
                     'equation': '',
                     'substrates': [],
                     'products': [],
                     'Scoef': [],
                     'Pcoef': [],
                     'transport': False,
                     'reversibility': 0
                     }
         
         return reaction
    
    
    #=================================================
    def initialize_compound(self):
        
        compound = {'name': [],
                    'identifiers': {'KEGG': []},
                   }
        
        return compound
    
    #==========================================
    def save_reactions_to_file(self, reversibility=False):
        '''
        '''
        filename_reactions = 'REACTIONS_' + self._id + '-' + self._version + '_' + time.strftime("%Y%m%d") + '.txt'
        
        REACTIONS = ''
        
        #====================
        # SAVE REACTIONS
        #=================
        KEYS = list(self._rxn.keys())
        KEYS.sort()
        
        for R in KEYS:
            
            if not self._rxn[R]['transport']:
            
                
                REACTIONS += R + ': '
                
                
                SUBSTRATES = ''
                for coef, S in zip(self._rxn[R]['Scoef'],self._rxn[R]['substrates']):
                    SUBSTRATES += str(int(coef)) + ' ' + S + ' + '
                
                
                PRODUCTS = ''
                for coef, P in zip(self._rxn[R]['Pcoef'],self._rxn[R]['products']):
                    PRODUCTS += str(int(coef)) + ' ' + P + ' + '
                
                
                
                if not reversibility:
                    
                    if self._rxn[R]['reversibility'] == 1:
                        
                        REACTIONS += SUBSTRATES[:-3] + ' --> ' + PRODUCTS[:-3] + '\n'
                    
                    elif self._rxn[R]['reversibility'] == -1:
                        
                        REACTIONS += PRODUCTS[:-3] + ' --> ' + SUBSTRATES[:-3] + '\n'
                        
                    else:
                        
                        REACTIONS += SUBSTRATES[:-3] + ' <=> ' + PRODUCTS[:-3] + '\n'
                
                else:
                    
                    REACTIONS += SUBSTRATES[:-3] + ' <=> ' + PRODUCTS[:-3] + '\n'
                    
        
        
        REACTIONS = REACTIONS.strip().split('\n')
        
        
        FIXED_REACTIONS = ''
        
        for R in REACTIONS:
            
            if 'NAD-P-OR-NOP' in R:
                
                R_NAD_NADH = R + ' '
                R_NAD_NADH = R_NAD_NADH.replace('NAD-P-OR-NOP','NAD').replace('NADH-P-OR-NOP','NADH')
                FIXED_REACTIONS += R_NAD_NADH.strip() + '\n'
                
                R_NADP_NADPH = R + '  '
                R_NADP_NADPH = R_NADP_NADPH.replace('NAD-P-OR-NOP','NADP').replace('NADH-P-OR-NOP','NADPH')
                FIXED_REACTIONS += R_NADP_NADPH.strip() + '\n'
                
            else:
                
                FIXED_REACTIONS += R + '\n'
        
        
        with open(os.path.join(self._path[0],filename_reactions), 'w') as fp:
            fp.writelines(FIXED_REACTIONS.strip())
    
    #==========================================
    def save_enzymes_to_file(self):
        '''
        '''
        
        filename_enzymes = 'ENZYMES_' + self._id + '-' + self._version + '_' + time.strftime("%Y%m%d") + '.txt'
        
        REACTIONS = ''
        
        #======================
        # SAVE ENZYMES
        #===================
        
        KEYS = list(self._rxn.keys())
        KEYS.sort()
        
        for R in KEYS:
            
            if not self._rxn[R]['transport']:
                
                if self._rxn[R]['enzymes']:
                    
                    enzymes = ', '.join(self._rxn[R]['enzymes'])
                    
                    REACTIONS += R + ': ' + enzymes + '\n'
    
        
        with open(os.path.join(self._path[0],filename_enzymes), 'w') as fp:
            fp.writelines(REACTIONS)
    
    
    #==========================================
    def save_cpd_names_to_file(self):
        '''
        '''
        
        filename_cpd_names = 'NAMES_' + self._id + '-' + self._version + '_' + time.strftime("%Y%m%d") + '.txt'
        
        #========================
        # SAVE COMPOUND NAMES
        #=====================
        KEYS = list(self._cpd.keys())
        KEYS.sort()
        
        COMPOUNDS = ''
        for C in KEYS:
            
            if self._cpd[C]['name']:
                
#                names = ', '.join(self._cpd[C]['name'])
                names = self._cpd[C]['name'][0]
                
                COMPOUNDS += C + ': ' + names + '\n'
        
        with open(os.path.join(self._path[0],filename_cpd_names), 'w') as fp:
            fp.writelines(COMPOUNDS)
    
    
    #==========================================
    def save_DB_to_file(self, reversibility=False):
        '''
        '''
        
        self.save_reactions_to_file(reversibility=reversibility)
        
        self.save_cpd_names_to_file()
        
        self.save_enzymes_to_file()


#======================================================================================================================






parser = argparse.ArgumentParser(description='Build PhDSeeker files from BioCyc database files.')
parser.add_argument('--cpdDB', help='"compounds.dat" file for a given BioCyc database. i.e EcoCyc.')
parser.add_argument('--rxnDB', help='"reactions.dat" file for a given BioCyc database. i.e EcoCyc.')
parser.add_argument('--allreversible', default='False', type=(lambda x: x.lower() in ("yes", "true", "t", "1")), help='Boolean value indicating if reactions should be taken always as reversible.')


args = vars(parser.parse_args())


if __name__ == '__main__':
    
    #==========================================
    #path,filename = os.path.split(args['file'])

    DATA = BioCyc(args['cpdDB'], args['rxnDB'])

    DATA.save_DB_to_file(reversibility=args['allreversible'])

    #==========================================
