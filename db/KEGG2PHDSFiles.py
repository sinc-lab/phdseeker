#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 09:30:44 2018

@author: mgerard
"""
import os
import json
import argparse

import urllib3, re, time


class KEGG(object):
    '''
    Clase para manejar los datos de KEGG.
    '''
    
    #===================================================
    def __init__(self, usedict=None):
        
        self._query = urllib3.PoolManager()
        
        self._rxn = dict()
        self._cpd = dict()
        self._gly = dict()
        
        if not os.path.exists('KEGG'):
            os.makedirs('KEGG')
            
        self.destination_folder = 'KEGG'
        
        self._organism = dict()
        
        self._translator = dict()
        
        print('Downloading information of compounds...\n')
        self._download_relation_cpd_name()
        print('Done!!\n')
        
        print('Downloading information of glycans...\n')
        self._download_relation_glycan_cpd()
        print('Done!!\n')
        
        print('Downloading information of reactions...\n')
        self._download_reactions()
        print('Done!!\n')
        
        print('Downloading information of enzymes...\n')
        self._download_enzymes()
        print('Done!!\n')
        
        
        if usedict is None:
            
            # BUILD TRANSLATION DICTIONARY FOR COMPOUNDS
            print('Building dictionary to translate compounds...\n')
            self.build_translator()
            
            translator_filename = 'CompDict_rn_' + time.strftime("%Y%m%d") + '.json'
            
            with open(self.destination_folder + '/' + translator_filename, 'w') as fp:
                json.dump(self._translator, fp)
            
        else:
            # LOAD DICTIONARY FOR COMPOUNDS
            print('Loading dictionary to translate compounds...\n')
            with open(usedict, 'r') as fp:
                self._translator = json.load(fp)
            
        print('Done!!\n')
    
    
    #===================================================
    def _download_relation_cpd_name(self):
        
        compounds = self._query.request('GET','http://rest.kegg.jp/list/compound')
        compounds = compounds.data.decode('UTF-8')
        compounds = compounds.strip().split('\n')
        
        for compound in compounds:
            
            code,names = compound.replace('cpd:','').split('\t')
            
            names = names.split(';')
            
            self._cpd[code] = self.initialize_compound()
            self._cpd[code]['cpd_code'] = code
            self._cpd[code]['identifiers']['KEGG'].append(code)
            
            for name in names:
                
                self._cpd[code]['name'].append(name.strip())
    
    
    
    #===================================================
    def _download_relation_glycan_cpd(self):
        
        #----------------------------------------------------------------------
        GL2CPD = self._query.request('GET','http://rest.kegg.jp/link/compound/glycan')
        GL2CPD = GL2CPD.data.decode('UTF-8')
        GL2CPD = GL2CPD.strip().split('\n')
        
        for gl2cpd in GL2CPD:
            
            gl,cpd = re.findall('gl:(G\d{5})\tcpd:(C\d{5})', gl2cpd)[0]
            
            
            if cpd not in self._cpd.keys():
                
                self._cpd[cpd] = self.initialize_compound()
                self._cpd[cpd]['identifiers']['KEGG'].append(cpd)
            
            if gl not in self._cpd.keys():
                
                self._cpd[gl] = self.initialize_compound()
                self._cpd[gl]['identifiers']['KEGG'].append(gl)
            
            
            self._cpd[cpd]['cpd_code'] = cpd
            self._cpd[cpd]['glycan_code'] = gl
            
            self._cpd[gl]['cpd_code'] = cpd
            self._cpd[gl]['glycan_code'] = gl
            
            
        #----------------------------------------------------------------------
        
        #----------------------------------------------------------------------
        GL2NAME = self._query.request('GET','http://rest.kegg.jp/list/glycan')
        GL2NAME = GL2NAME.data.decode('UTF-8')
        GL2NAME = GL2NAME.strip().split('\n')
        
        for gl2name in GL2NAME:
            
            gl,names = gl2name.replace('gl:','').split('\t')
            
            if gl not in self._cpd.keys():
                
                self._cpd[gl] = self.initialize_compound()
                self._cpd[gl]['glycan_code'] = gl
                
                self._cpd[gl]['identifiers']['KEGG'].append(gl)
                
            
            names = names.split(';')
            
            for name in names:
                
                if name not in self._cpd[gl]['name']:
                    self._cpd[gl]['name'].append(name.strip())
                
                cpd_code = self._cpd[gl]['cpd_code']
                
                if (cpd_code != '') and self._cpd[cpd_code]['name']:
                    
                    if name not in self._cpd[cpd_code]['name']:
                        
                        self._cpd[cpd_code]['name'].append(name.strip())
    
    
    
    #==========================================
    def _download_reactions(self):
        '''
        '''
        
        COMPOUNDS = set()
        
        reactions = self._query.request('GET','http://rest.kegg.jp/list/reaction')
        reactions = reactions.data.decode('UTF-8')
        
        #=======================================================
        # CORRECCION DE ERRORES EN LA BASE DE DATOS
        #===============================================
        # Ferrocytochrome b51 --> Ferrocytochrome b5
        #--------------------------------------------
        reactions = reactions.replace('Ferrocytochrome b51','Ferrocytochrome b5')
        #=======================================================
        
        reactions = reactions.replace('rn:','').strip()
        reactions = reactions.split('\n')
        
        
        for reaction in reactions:
            
            R,name_equation = reaction.split('\t')
            rxn_code = R.strip()
            
            data = name_equation.split('; ')
            equation = data[-1]
            if len(data) == 2:
                name = data[0]
            elif len(data) == 1:
                name = ['']
            else:
                name = []
                for d in data[:-1]:
                    name.append(d.strip())
            
            
            self._rxn[rxn_code] = self.initialize_reaction()
            self._rxn[rxn_code]['identifiers']['KEGG'].append(rxn_code)
            
            self._rxn[rxn_code]['name'] = name[:]
            self._rxn[rxn_code]['equation'] = equation.strip()
            
            S,P = equation.rstrip().split(' <=> ')
            
            
            # PROCESAR Y PARSEAR SUBSTRATES
            S = S.split(' + ')
            for s in S[:]:
                
                coef = re.findall('^\d+\ |^[1-9nm]+\ |^\([1-9nm]+\)\ |\([1-9nmx]+\)$|^\([nm][+-][nm\d]\)\ |^[nm][+-][nm1-3]\ |\([nm][+-][nmx\d]\)$',s)
                
                # REACOMODO COEFICIENTES
                if len(coef) > 0:
                    s = s.strip().replace(coef[0],'',1)
                    
                else:
                    coef = ['1']
                
                
                #----------------------------------------------------
                # CORRECCIONES
                #----------------------
                coef = coef[0].strip()
                coef = coef.replace('(n-2)' ,'98')
                coef = coef.replace('(n-1)' ,'99')
                coef = coef.replace('(n)'   ,'100')
                coef = coef.replace('(n+1)' ,'101')
                coef = coef.replace('(n+2)' ,'102')
                coef = coef.replace('(m-2)' ,'98')
                coef = coef.replace('(m-1)' ,'99')
                coef = coef.replace('(m)'   ,'100')
                coef = coef.replace('(m+1)'   ,'101')
                coef = coef.replace('(m+2)'   ,'102')
                coef = coef.replace('(n+m)' ,'200')
                coef = coef.replace('(m+n)' ,'200')
                coef = coef.replace('(n-x)' ,'90')
                coef = coef.replace('n-1'   ,'99')
                coef = coef.replace('n+1'   ,'101')
                coef = coef.replace('n'     ,'100')
                coef = coef.replace('(x)'     ,'10')
                #----------------------------------------------------
                s = s.replace('(in)','')
                s = s.replace('(out)','')
                #----------------------------------------------------
                
                COMPOUNDS.add(s.strip())
                
                self._rxn[rxn_code]['substrates'].append(s.strip())
                self._rxn[rxn_code]['Scoef'].append(float(coef[0]))
                
            
            # PROCESAR Y PARSEAR PRODUCTS
            P = P.split(' + ')
            for p in P[:]:
                
                coef = re.findall('^\d+\ |^[1-9nm]+\ |^\([1-9nm]+\)\ |\([1-9nmx]+\)$|^\([nm][+-][nm\d]\)\ |^[nm][+-][nm1-3]\ |\([nm][+-][nmx\d]\)$',p)
                
                # REACOMODO COEFICIENTES
                if len(coef) > 0:
                    p = p.replace(coef[0],'',1)
                else:
                    coef = ['1']
                
                #----------------------------------------------------
                # CORRECCIONES
                #----------------------
                coef = coef[0].strip()
                coef = coef.replace('(n-2)' ,'98')
                coef = coef.replace('(n-1)' ,'99')
                coef = coef.replace('(n)'   ,'100')
                coef = coef.replace('(n+1)' ,'101')
                coef = coef.replace('(n+2)' ,'102')
                coef = coef.replace('(m-2)' ,'98')
                coef = coef.replace('(m-1)' ,'99')
                coef = coef.replace('(m)'   ,'100')
                coef = coef.replace('(m+1)'   ,'101')
                coef = coef.replace('(m+2)'   ,'102')
                coef = coef.replace('(n+m)' ,'200')
                coef = coef.replace('(m+n)' ,'200')
                coef = coef.replace('(n-x)' ,'90')
                coef = coef.replace('n-1'   ,'99')
                coef = coef.replace('n+1'   ,'101')
                coef = coef.replace('n'     ,'100')
                coef = coef.replace('(x)'     ,'10')
                #----------------------------------------------------
                p = p.replace('(in)','')
                p = p.replace('(out)','')
                #----------------------------------------------------
                
                COMPOUNDS.add(p.strip())
                
                self._rxn[rxn_code]['products'].append(p.strip())
                self._rxn[rxn_code]['Pcoef'].append(float(coef[0]))
        
        
    #==========================================
    def _download_enzymes(self):
        '''
        '''
        
        #------------------
        # REACTIONS
        #------------
        enzymes = self._query.request('GET','http://rest.kegg.jp/link/reaction/enzyme')
        enzymes = enzymes.data.decode('UTF-8')
        enzymes = enzymes.strip().split('\n')
        
        for enzyme in enzymes:
            
            ec,rxn_code = re.findall('ec:(.{4,14})\trn:(R\d{5})', enzyme)[0]
            
            if rxn_code not in self._rxn.keys():
                
                self._rxn[rxn_code] = self.initialize_reaction()
                
                self._rxn[rxn_code]['identifiers']['KEGG'].append(rxn_code)
            
            
            self._rxn[rxn_code]['enzymes'].append(ec)
    
    
    #==========================================
    def download_organism_data(self, organism='rn', pathways='all', reversibility=False):
        '''
        reversibility=False --> Download information of reversibility
        '''
        
        self.reversibility = reversibility
        
        print('Initializing organism...\n')
        self._initialize_organism()
        self._organism['id'] = organism
        print('Done!!\n')
        
        #------------------------------------------------------
        print('Downloading information of enzymes for {}...\n'.format(organism))
        enzymes = self._query.request('GET', 'http://rest.kegg.jp/link/enzyme/%s' % organism)
        enzymes = enzymes.data.decode('UTF-8')
        enzymes = enzymes.strip().split('\n')
        
        ENZYMES = []
        for enzyme in enzymes:
            ENZYMES.append(enzyme.split(':')[-1])

        
        for R in self._rxn.keys():
            
            if ( R not in self._organism['reactions'].keys() ):
                
                insert = False
                
                for enzyme in self._rxn[R]['enzymes']:
                    
                    if enzyme in ENZYMES:
                        insert = True
                
                if insert:
                    self._organism['reactions'][R] = self._rxn[R].copy()
                    self._organism['reactions'][R]['use'] = False
        
        print('Done!!\n')
        #------------------------------------------------------
        
        #------------------------------------------------------
        # REVERSIBILITY
        #----------------
        REVERSIBILITY = dict()
        
        if pathways == 'all':
            self._Npathways = 0
            
        else:
            if ',' in pathways:
                pathways = pathways.replace(' ','').split(',')
            else:
                pathways = [pathways]
            
            self._Npathways = len(pathways)
            
        
        if (reversibility):
            
            print('Downloading information of the reversibility of reactions in the organism...\n')
            _pathways = self._query.request('GET', 'http://rest.kegg.jp/list/pathway/%s' % organism)
            _pathways = _pathways.data.decode('UTF-8').strip().split('\n')
            
                
            _pathways = [[_pathway for _pathway in _pathways if pathway_code in _pathway][0] for pathway_code in pathways]
            
            
            
            REACTIONS_org = ''
            for pathway in _pathways:
                
                pw = pathway.split('\t')[0][5:]
                
                kgml = self._query.request('GET', 'http://rest.kegg.jp/get/' + pw + '/kgml')
                
                REACTIONS_org += kgml.data.decode('UTF-8').strip()
                
            REACTIONS_org = re.findall('<reaction id=.\d*. name=.rn:R\d{5}. type="[ir]{1,3}eversible">[\s\W\w]{1,1000}?</reaction>\n|<reaction id=.\d*. name=.rn:R\d{5} rn:R\d{5}. type="[ir]{1,3}eversible">[\s\W\w]{1,1000}?</reaction>\n|<reaction id=.\d*. name=.rn:R\d{5} rn:R\d{5} rn:R\d{5}. type="[ir]{1,3}eversible">[\s\W\w]{1,1000}?</reaction>\n|<reaction id=.\d*. name=.rn:R\d{5} rn:R\d{5} rn:R\d{5} rn:R\d{5}. type="[ir]{1,3}eversible">[\s\W\w]{1,1000}?</reaction>\n|<reaction id=.\d*. name=.rn:R\d{5} rn:R\d{5} rn:R\d{5} rn:R\d{5} rn:R\d{5}. type="[ir]{1,3}eversible">[\s\W\w]{1,1000}?</reaction>\n|<reaction id=.\d*. name=.rn:R\d{5} rn:R\d{5} rn:R\d{5} rn:R\d{5} rn:R\d{5} rn:R\d{5}. type="[ir]{1,3}eversible">[\s\W\w]{1,1000}?</reaction>\n',REACTIONS_org)
            print('Done!!\n')
            #-----------------------------------------------------
            
            #-----------------------------------------------------
            print('Extracting reversibility information...\n')
            for REACTION in REACTIONS_org:
                
                reactions,rev = re.findall('<reaction id=.\d*. name=.(rn:R\d{5}.*). type="(reversible|irreversible)">',REACTION)[0]
                reactions = re.findall('(R\d{5})', reactions)
                
                
                S = set(re.findall('<substrate id=.\d*. name=.cpd:([CG]\d{5})./>\n',REACTION))
                P = set(re.findall('<product id=.\d*. name=.cpd:([CG]\d{5})./>\n',REACTION))
                
                for reaction in reactions:
                    
                    if reaction in self._organism['reactions'].keys():
                        self._organism['reactions'][reaction]['use'] = True
                    
                    
                    if REVERSIBILITY.get(reaction,0) != 0:
                        
                        if rev == 'reversible':
                            REVERSIBILITY[reaction]['reversibility'] = 'reversible'
                            REVERSIBILITY[reaction]['S'] = S
                            REVERSIBILITY[reaction]['P'] = P
                        
                        
                        else:
                            
                            if REVERSIBILITY[reaction]['reversibility'] != 'reversible':
                                
                                if (len(REVERSIBILITY[reaction]['S'].intersection(S)) == 0) and len(REVERSIBILITY[reaction]['P'].intersection(P)) == 0:
                                    
                                    REVERSIBILITY[reaction]['reversibility'] = 'reversible'
                                
                                else:
                                    
                                    REVERSIBILITY[reaction]['S'].update(S)
                                    REVERSIBILITY[reaction]['P'].update(P)
                    else:
                        
                        REVERSIBILITY[reaction] = {'reversibility':rev, 'S': S, 'P': P}
            
            print('Done!!\n')
            #-------------------------------------------------------
            
            
        #-------------------------------------------------------
        print('Building organism...\n')
        
        for R in self._organism['reactions'].keys():
            
            if R not in REVERSIBILITY.keys(): # SUPONGO QUE ES REVERSIBLE
                
                self._organism['reactions'][R]['reversibility'] = 0
                
            else:
                
                if REVERSIBILITY[R]['reversibility'] == 'reversible':
                    
                    self._organism['reactions'][R]['reversibility'] = 0
                    
                else:
                    
                    substrates = self._organism['reactions'][R]['substrates']
                    
                    if REVERSIBILITY[R]['S'].issubset(set(substrates)):
                        
                        self._organism['reactions'][R]['reversibility'] = -1
                        
                    else:
                        
                        self._organism['reactions'][R]['reversibility'] = 1
    
        print('Done!!\n')
        #-------------------------------------------------------

    
    #==========================================
    def save_reactions_to_file(self, translate=True, organism=False):
        '''
        '''
        
        if (self.destination_folder != ''):
            filename_reactions = self.destination_folder + '/' + 'REACTIONS_rn_' + time.strftime("%Y%m%d") + '.txt'
        else:
            filename_reactions = 'REACTIONS_rn_' + time.strftime("%Y%m%d") + '.txt'
        
        #====================
        # SAVE REACTIONS
        #=================
        REACTIONS = ''
        
        if (organism != False):
            
            if self._Npathways == 0:
                
                filename_reactions = filename_reactions.replace('_rn_', '_' + self._organism['id'] + '_')
                
            else:
                
                Npathways = 'pathways_' if (self._Npathways > 1) else 'pathway_'
                
                filename_reactions = filename_reactions.replace('_rn_', '_' + self._organism['id'] + '-' + str(self._Npathways) + Npathways)
            
            
            KEYS = list(self._organism['reactions'].keys())
            KEYS.sort()
            
            for R in KEYS:
                
                if (self._organism['reactions'][R]['use'] == True) or (not self.reversibility == True):
                
                    REACTIONS += R + ': '
                    
                    reversibility = self._organism['reactions'][R]['reversibility']
                    
                    #········································································
                    if reversibility == 0:
                        
                        for coef, S in zip(self._organism['reactions'][R]['Scoef'],self._organism['reactions'][R]['substrates']):
                        
                            if translate:
                                S = self._translator[S]
                            
                            REACTIONS += str(int(coef)) + ' ' + S + ' + '
                        
                        
                        REACTIONS = REACTIONS[:-3] + ' <=> '
                        
                        for coef, P in zip(self._organism['reactions'][R]['Pcoef'],self._organism['reactions'][R]['products']):
                        
                            if translate:
                                P = self._translator[P]
                                
                            REACTIONS += str(int(coef)) + ' ' + P + ' + '
                    #········································································
                    
                    
                    #········································································
                    if reversibility == 1:
                        
                        for coef, S in zip(self._organism['reactions'][R]['Scoef'],self._organism['reactions'][R]['substrates']):
                        
                            if translate:
                                S = self._translator[S]
                            
                            REACTIONS += str(int(coef)) + ' ' + S + ' + '
                        
                        
                        REACTIONS = REACTIONS[:-3] + ' --> '
                        
                        
                        for coef, P in zip(self._organism['reactions'][R]['Pcoef'],self._organism['reactions'][R]['products']):
                        
                            if translate:
                                P = self._translator[P]
                                
                            REACTIONS += str(int(coef)) + ' ' + P + ' + '
                    #········································································    
                    
                    
                    #········································································
                    if reversibility == -1:
                        
                        for coef, P in zip(self._organism['reactions'][R]['Pcoef'],self._organism['reactions'][R]['products']):
                        
                            if translate:
                                P = self._translator[P]
                                
                            REACTIONS += str(int(coef)) + ' ' + P + ' + '
                        
                        
                        REACTIONS = REACTIONS[:-3] + ' --> '
                        
                        
                        for coef, S in zip(self._organism['reactions'][R]['Scoef'],self._organism['reactions'][R]['substrates']):
                        
                            if translate:
                                S = self._translator[S]
                            
                            REACTIONS += str(int(coef)) + ' ' + S + ' + '
                    
                    
                    REACTIONS = REACTIONS[:-3] + '\n'
                    #········································································
            
            
        else:

            KEYS = list(self._rxn.keys())
            KEYS.sort()

            for R in KEYS:
                REACTIONS += R + ': '
            
            for coef, S in zip(self._rxn['reactions'][R]['Scoef'],self._rxn['reactions'][R]['substrates']):
                
                if translate:
                    S = self._translator[S]
                    
                REACTIONS += str(int(coef)) + ' ' + S + ' + '
                
            REACTIONS = REACTIONS[:-3] + ' <=> '
            
            for coef, P in zip(self._rxn['reactions'][R]['Pcoef'],self._rxn['reactions'][R]['products']):
                
                if translate:
                    P = self._translator[P]
                    
                REACTIONS += str(int(coef)) + ' ' + P + ' + '
            
            
            REACTIONS = REACTIONS[:-3] + '\n'
            
        
        with open(filename_reactions, 'w') as fp:
            fp.writelines(REACTIONS)
    
    
    #==========================================
    def save_enzymes_to_file(self, organism=False):
        '''
        '''
        
        if (self.destination_folder != ''):
            filename_enzymes = self.destination_folder + '/' + 'ENZYMES_rn_' + time.strftime("%Y%m%d") + '.txt'
        else:
            filename_enzymes = 'ENZYMES_rn_' + time.strftime("%Y%m%d") + '.txt'
        
        
        
        REACTIONS = ''
        #======================
        # SAVE ENZYMES
        #===================
        if organism:
            
            filename_enzymes = filename_enzymes.replace('_rn_', '_' + self._organism['id'] + '_')
            
            KEYS = list(self._organism['reactions'].keys())
            KEYS.sort()
            
            for R in KEYS:
                
                if self._organism['reactions'][R]['enzymes']:
                    enzymes = ', '.join(self._organism['reactions'][R]['enzymes'])
                    
                    REACTIONS += R + ': ' + enzymes + '\n'
            
        else:
            
            KEYS = list(self._rxn.keys())
            KEYS.sort()
            
            for R in KEYS:
                
                if self._rxn[R]['enzymes']:
                    enzymes = ', '.join(self._rxn[R]['enzymes'])
                    
                    REACTIONS += R + ': ' + enzymes + '\n'
        
        with open(filename_enzymes, 'w') as fp:
            fp.writelines(REACTIONS)
    
    
    #==========================================
    def save_cpd_names_to_file(self, organism=True):
        '''
        '''
        
        if (self.destination_folder != ''):
            filename_cpd_names = self.destination_folder + '/' + 'NAMES_rn_' + time.strftime("%Y%m%d") + '.txt'
        else:
            filename_cpd_names = 'NAMES_rn_' + time.strftime("%Y%m%d") + '.txt'
        
        
        #========================
        # SAVE COMPOUND NAMES
        #=====================
        if organism:
            
            filename_cpd_names = filename_cpd_names.replace('_rn_', '_' + self._organism['id'] + '_')
            
        
        KEYS = list(self._cpd.keys())
        KEYS.sort()
        
        COMPOUNDS = ''
        for C in KEYS:
            
            if self._cpd[C]['name']:
                
                #names = ', '.join(self._cpd[C]['name'])  # <--- ESCRIBE TODOS LOS NOMBRES
                names = self._cpd[C]['name'][0]           # <--- ESCRIBE EL PRIMER NOMBRE
                
                COMPOUNDS += C + ': ' + names + '\n'
        
        with open(filename_cpd_names, 'w') as fp:
            fp.writelines(COMPOUNDS)
    
    
    
    #==========================================
    def save_organism_to_file(self, translate=True):
        '''
        '''
        
        self.save_reactions_to_file(translate=translate, organism=True)
        
        self.save_enzymes_to_file(organism=True)
        
        self.save_cpd_names_to_file(organism=True)
    
    #=================================================
    def _translate(self, compound):
        '''
        '''
        
        translated_cpd = compound.replace(' ', '_')
        
        
        #----------------------------------------
        # VERIFICO SI EL COMPUESTO ES UN CODIGO
        #----------------------------------------
        if compound in self._cpd.keys():
            
            cpd = self._cpd[compound]
            
            if cpd['glycan_code'] != '':
                
                translated_cpd = cpd['glycan_code']
                
            if cpd['cpd_code'] != '':
                
                translated_cpd = cpd['cpd_code']
        #-------------------------------------
        
        
        #-------------------------------------
        else:
            
            for key in self._cpd.keys():
                
                cpd = self._cpd[key]
                
                if compound in cpd['name']:
                    
                    if cpd['cpd_code'] != '':
                        translated_cpd = cpd['cpd_code']
                    else:
                        translated_cpd = cpd['glycan_code']
        #-------------------------------------
        
        
        return translated_cpd
    
    
    #====================================================
    def build_translator(self):
        '''
        '''
        
        self._translator = dict()
        
        COMPOUNDS = set()
        for R in self._rxn.keys():
            
            COMPOUNDS.update(self._rxn[R]['substrates'])
            COMPOUNDS.update(self._rxn[R]['products'])
            
        COMPOUNDS = list(COMPOUNDS)
        for compound in COMPOUNDS:
            
            self._translator[compound] = self._translate(compound)
    
    
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
                     'reversibility': 0,
                     'use': True
                     }
         
         return reaction
    
    
    #=================================================
    def initialize_compound(self):
        
        compound = {'name': [],
                    'identifiers': {'KEGG': []},
                    'cpd_code': '',
                    'glycan_code': ''
                   }
        
        return compound
    
    
    
    #=================================================
    def _initialize_organism(self):
        
        self._organism = {'id': '',
                          'compounds': dict(),
                          'reactions': dict()
                         }
        
        
###############################################################################








###############################################################################
# MAIN
########

parser = argparse.ArgumentParser(description='Build PhDSeeker files from KEGG rest service (requires internet connection).')
parser.add_argument('-o', '--organism', default='rn', help='Organism for which it is wanted to build the files.')
parser.add_argument('-p', '--pathways', default='all', help='List of pathways to be downloaded. By default, all available pathways are used to build the dataset.')
parser.add_argument('-r', '--allreversible', default='False', type=(lambda x: x.lower() in ("yes", "true", "t", "1")), help='Boolean value indicating if reactions should be taken always as reversible.')
parser.add_argument('-d', '--usedict', default=None, help='Specify a dictionary to translate KEGG codes into compound names.')

args = vars(parser.parse_args())


if __name__ == '__main__':

    #==========================================
    DATA = KEGG(usedict=args['usedict'])

    DATA.download_organism_data(organism=args['organism'], pathways=args['pathways'], reversibility=args['allreversible'])

    DATA.save_organism_to_file()

    #==========================================
