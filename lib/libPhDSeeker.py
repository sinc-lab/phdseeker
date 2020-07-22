#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 20:31:15 2014

@author: matias
"""

import yaml

import copy

import numpy as np
from scipy.sparse import lil_matrix #, csc_matrix
import sys
import json


#==============================================================
# FUNCTIONS FOR PARALELLIZATION
#==============================================================


#----------------
# INITIALIZATION
#----------------
def initializing(ARGUMENTS):
    '''
    Executes the "Initialize" method implemented in ants.
    '''
    
    ant, pool, idx_initial, idx_produce, compounds, Rinitial = ARGUMENTS
    
    ant.Initialize(pool, idx_initial, idx_produce, compounds, Rinitial)
    
    return ant

#---------------------
# SEARCH OF SOLUTIONS
#---------------------
def searching(ARGUMENTS):
    
    '''
    Executes the "Search" method implemented in ants.
    '''
    
    ant, S, P, pheromones, Cnames, Rnames = ARGUMENTS
    
    ant.Search(S, P, pheromones, Cnames, Rnames)
    return ant
#==============================================================


#==============================================================================
class Ant(object):
    '''
    Esta clase implementa el agente "hormiga", definiendo los atributos y
    métodos que tendrá disponible.
    '''
    
    def __init__(self, Rmax=100, idx=0, strict_initialization=False):
        
        self.idx = idx
        self.cost = np.inf
        self.found = False
        
        # REACTIONS
        self.semireactions = list() # i.e. 1, 2, etc
        self.filtered_semireactions = list() # REACTIONS EFFECTIVELY USED 
                                             # TO RELATE THE SPECIFIED
                                             # COMPOUNDS
        
        self.strict_initialization = strict_initialization
        

    #===========================================
    def Cost(self, S, P, Rnames):
        """EXPLICACION"""
#        import timeit
#        S = np.array(S.todense()).astype(np.bool)
#        P = np.array(P.todense()).astype(np.bool)
        
        
        #------------------------------------------------
        # EVALUACION DE LA CONECTIVIDAD
        #------------------------------------------------
        Pcum = np.zeros(S.shape[0]).astype(np.bool)
        
        Z = np.zeros(S.shape[0]).astype(np.bool)
        Z[self.idx_source] = True
        
        
        NCG = 0
        for r in self.filtered_semireactions:
            
            s = S[:,r].flatten()
            p = P[:,r].flatten()
            
            if Z[s].sum() > 0:
                Z[self.nonCa] += p[self.nonCa]
            
            
            if np.not_equal(Pcum[self.nonCa], Pcum[self.nonCa]+p[self.nonCa]).sum() > 0:
                Pcum[self.nonCa] += p[self.nonCa]
                NCG += 1
        
        
        if ( Z[self.idx_target].sum() == len(self.idx_target)):
            Conectividad = np.float64(0.1)
        else:
            Conectividad = np.float64(2.0)
        
        
        #------------------------------------------------
        # NUMERO TOTAL DE REACCIONES EN LA VIA METABOLICA
        #------------------------------------------------
        NRT = np.float64(len(self.filtered_semireactions))
        
        #------------------------------------------------
        # NUMERO DE REACCIONES UNICAS EN LA VIA METABOLICA
        #------------------------------------------------
        NRU = np.float64(len(set([Rnames[r] for r in self.filtered_semireactions])))
        
        
        #===========================
        # COST FUNCTION
        #===========================
        Costo = (NRT/NCG) * NRU * Conectividad
        
        
        SALIDA = dict()
        SALIDA['COSTO'] = Costo
        SALIDA['NRT'] = NRT
        SALIDA['NRU'] = NRU
        SALIDA['CONECTIVIDAD'] = Conectividad
        
        
        return SALIDA
    
    
    
    #===========================================
    
    
    def ExportPathway(self, pool):
        """DEVUELVE EL PATHWAY ENCONTRADO (CON TODOS LOS DATOS NECESARIOS)"""
        
        pathway = Pathway(pool, self.filtered_semireactions)
        
        pathway.compounds['initials'] = [pool.Cnames[self.idx_source]]
        pathway.compounds['abundant'] = [pool.Cnames[idx] for idx in np.argwhere(self.compounds['abundant'] > 0).flatten()]
        pathway.compounds['external'] = [pool.Cnames[idx] for idx in np.argwhere(self.compounds['external'] > 0).flatten()]
        pathway.compounds['produce'] = [pool.Cnames[idx] for idx in self.idx_target]
        
        pathway.cost = self.Cost(pool.S, pool.P, pool.Rnames)
        
        return pathway
    
    
    #===========================================
    def Initialize(self, pool, idx_initial, idx_produce, compounds, Rinitial):
        '''
         def Initialize(pool, relate_compounds, initial_compounds, abundant_compounds, external_compounds, pheromones, memory)
        
        Este método ejecuta todas las tareas necesarias para inicializar la hormiga.
        - Limpia los conjuntos de compuestos de la hormiga.
        - Establece los diferentes conjuntos de compuestos.
        - Selecciona una reacción, de forma aleatoria, para cada uno de los posibles compuestos iniciales y la inserta en el camino.
        - Actualiza los distintos conjuntos de compuestos luego de insertar cada reacción.
        
        '''
        
        self.cost = np.inf
        self.found = False
        
        self.Rinitials = Rinitial
        self.semireactions.append(self.Rinitials)
        
        self.filtered_semireactions = list()
        
        
        self.idx_source = idx_initial
        
        self.idx_target = idx_produce
        
#        self.nonCa = compounds['non abundant'].copy()
        self.nonCa = np.array(compounds['non abundant']).astype(np.bool)
        
        self.Co = (compounds['abundant'] + compounds['external']).astype(np.bool)
        
        self.compounds = dict()
        self.compounds['initials'] = [pool.Cnames[idx_initial]]
        self.compounds['produce'] = [pool.Cnames[idx] for idx in idx_produce]
        self.compounds['abundant'] = compounds['abundant'].copy()
        self.compounds['external'] = compounds['external'].copy()
        
        
        
        if self.strict_initialization:
            
            self.Co[idx_initial] = True
        
        else:
            
            self.Co[pool.S[:,Rinitial]] = True
            
#            self.Co = ((self.Co + pool.S[:,Rinitial].flatten()) > 0).astype(np.float32)
        
#        self.Co = (self.Co).astype(np.bool)
        
        
        return self
    
    
    #===========================================
    def Search(self, S, P, FEROMONAS, Cnames, Rnames):
        """EXPLICACION"""
#        import timeit
        
        
        np.random.seed()
        
        #-------------------------
        # ENTRADAS:
        #-------------------------
        # Ri
        # Co
        # S, P
        # FEROMONAS
        # idx_source --> GENERO "Z", idx_target
        self.found = False
        FOUND = False
        
        Nr = S.shape[1]
        
        # INICIALIZO PATHWAY
        pathway = [self.Rinitials]
        binary_pathway = np.ones(Nr).astype(np.float32)
        binary_pathway[self.semireactions[0]] = 0
        
        reactions_binary = np.ones(Nr).astype(np.float32)
        reactions_binary[[idx for idx in range(len(Rnames)) if Rnames[idx] == Rnames[self.Rinitials]]] = 0
        
        
        s = S[:,self.semireactions[0]].flatten()
        p = P[:,self.semireactions[0]].flatten()
        
        # ACTUALIZO COMPUESTOS DISPONIBLES
        self.Co[s] = True
        
        Z = np.zeros(S.shape[0]).astype(np.bool)
        Z[self.idx_source] = True
        Z[self.nonCa] += p[self.nonCa]
        
        STOP = False
        #-------------------------------------------
        
        # INDICE DE LAS REACCIONES QUE AUN NO SON FACTIBLES
        idx_infact = np.ones(Nr).astype(np.bool) # <--- Puedo determinar luego las factibles
        
        # NUMERO DE SUSTRATOS QUE REQUIERE CADA REACCION
        NS = S.sum(axis=0)
        
        
        while not STOP:
            
            feromonas = FEROMONAS[pathway[-1],:]
            
            #------------------------------------
            # AGREGO PRODUCTOS
            #-------------------------
            old_Co = self.Co.copy()
            
            self.Co[p] = True
            
            if np.all(np.equal(old_Co, self.Co)):
                Co_updated = False
            else:
                Co_updated = True
            #------------------------------------
            
            if (idx_infact.sum() > 0) and Co_updated:
                idx_infact_upd = np.not_equal(S[self.Co,:][:,idx_infact].sum(axis=0),NS[idx_infact])
                
                if (idx_infact_upd.sum() > 0):
                    idx_infact[np.argwhere(idx_infact == True).flatten()] = idx_infact_upd
                
            
            # ELIJO NUEVA REACCION
#            search = (binary_pathway[idx_infact == False] * feromonas[idx_infact == False]).cumsum()
            search = ((binary_pathway * reactions_binary)[idx_infact == False] * feromonas[idx_infact == False]).cumsum()
            
            if (search.sum() == 0):
                STOP = True
                
            else:
                
                r_selected = np.argmin((search/search.max()) < np.random.rand())
                
                idx = np.argwhere(idx_infact == False).flatten()[r_selected]
                
                s = S[:,idx].flatten()
                p = P[:,idx].flatten()
                
                
                if Z[s].sum() > 0:
                    Z[self.nonCa] += p[self.nonCa]
                
                
                if Z[self.idx_target].sum() == len(self.idx_target):
                    STOP = True
                    FOUND = True
                
                
                
                pathway.append(idx)
                binary_pathway[idx] = 0
                reactions_binary[[i for i in range(len(Rnames)) if Rnames[i] == Rnames[idx]]] = 0
        
        
        #----------------------
        # PATHWAY FILTERING
        #----------------------
        if FOUND:
            
            self.filtered_semireactions, C = self._filtering(pathway, S, P, self.nonCa, self.idx_target)
            
            self.found = True
            
            self.Co = [Cnames[idx] for idx in range(self.Co.shape[0]) if self.Co[idx]]
        
        else:
            
            self.semireactions = pathway[:]
            self.filtered_semireactions = []
        
        return self
        
    #===========================================================================================
    def _filtering(self, pathway, S, P, nonCa, idx_target):
        '''
        '''
        
        C = []
        filtered_pathway = []
        
        # FILTRADO DEL PATHWAY
        Xs = S[:,pathway]
        Xp = P[:,pathway]
        
        
        for idx in range(len(pathway)-1,-1,-1):
            
            used_products = 0
            if idx < len(pathway)-1:
                
                used_products = Xs[nonCa,idx+1:]
                used_products = used_products[Xp[nonCa,idx],:].sum()
                
                
            if (Xp[idx_target,idx].sum() > 0) or (used_products > 0):
                filtered_pathway.insert(0, pathway[idx])
                
            else:
                # ELIMINO COLUMNA DE Xs y Xp
                Xs = np.delete(Xs,idx,1)
                Xp = np.delete(Xp,idx,1)
        
        C = np.any(Xs + Xp, axis=1)
            
        return (filtered_pathway, C)
    #===========================================================================================
    
        

    #===========================================
    def Size(self):
        """EXPLICACION"""
        return len(self.semireactions)
    
    
    #===========================================
    def ExportData(self, pool):
        '''
        Exporta los datos completos de la hormiga a un diccionario.
        '''
        data = dict()
        
        data['semireactions'] = self.semireactions
        data['filtered_semireactions'] = self.filtered_semireactions
        
        data['reactions'] = [pool.Rnames[r] for r in self.semireactions]
        data['filtered_reactions'] = [pool.Rnames[r] for r in self.filtered_semireactions]
        
        data['Rinitials'] = pool.Rnames[self.Rinitials]
        
        
        # ACA!!!
        data['enzymes'] = [pool.Enzyme(r, onlyfirst=False) for r in data['filtered_semireactions']]
        
        Compounds = dict()
        for key,value in self.compounds.items():
            Compounds[key] = list(value)
        
        
        data['compounds'] = Compounds
        
        data['cost'] = self.Cost(pool.S, pool.P, pool.Rnames)
        
        return data
    

#==============================================================================
class Anthill(object):
    """EXPLICACION"""
    
    
    #===========================================
    def __init__(self, pool, Nants=100, Rmax=100, compounds=dict(), strict_initialization=False):
        
        self.size = Nants
        
        self.ants = [Ant(Rmax=Rmax, idx=idx) for idx in range(0,self.size)]
        
        
        
        #======================================================================
        # COMPOUNDS - v2
        #======================================================================
        self.compounds = dict()
        
        binary_compounds = np.zeros(pool.S.shape[0], dtype=np.float32)
        
#        print([pool.Cnames.index(c) for c in compounds['relate']])
#        binary_compounds[[pool.Cnames.index(c) for c in compounds['relate']]] = 1
        
        self.compounds['relate'] = binary_compounds
        self.compounds['relate'][[pool.Cnames.index(c) for c in compounds['relate']]] = 1
        
        self.compounds['abundant'] = np.zeros(pool.S.shape[0], dtype=np.float32)
        self.compounds['abundant'][[pool.Cnames.index(c) for c in compounds['abundant'] if (c not in compounds['relate']) and (c in pool.Cnames)]] = 1
        
        self.compounds['external'] = np.zeros(pool.S.shape[0], dtype=np.float32)
        self.compounds['external'][[pool.Cnames.index(c) for c in compounds['external'] if (c not in compounds['relate']) and (c in pool.Cnames)]] = 1
        
        self.compounds['non abundant'] = 1.0 - (self.compounds['abundant'] + self.compounds['external']) > 0
        
        self.Cinitials = []
        self.Rinitials = []
        self.pRinitials = []
        
        
        # NON FORCED INITIAL COMPOUNDS
        for cinitial in compounds['initials']:
            
            idx_cinitial = pool.Cnames.index(cinitial)
            
            R = [int(idx) for idx in np.argwhere(pool.S[idx_cinitial,:] == True)]
            
            #------------------------------------------------
            if strict_initialization:
                
                Co = np.zeros(pool.S.shape[0], dtype=np.bool)
                Co[idx_cinitial] = True
                Co[self.compounds['abundant'].astype(np.bool)] = True
                Co[self.compounds['external'].astype(np.bool)] = True
                
                Rf = pool.FeasibleReactions(Co)
                
                R = [int(r) for r in R if r in Rf]
            #------------------------------------------------
            if len(R) > 0:
                
                self.Cinitials.extend([cinitial] * len(R))
                self.Rinitials.extend(R)        
                
        self.pRinitials = ((1.0/len(self.Rinitials)) * np.ones(len(self.Rinitials))).astype(np.float32)
        #======================================================================
        
        
        # PHEROMONES INITIALIZATION
        self.pheromones = np.ones((pool.size,pool.size), dtype=np.float32)
        
        
        
        # PATHWAYS
        self.bestpathway = Pathway()
        self.bestcost = np.inf
        
        self.convergence = 0.0
    
    
    #===========================================
    def BestPathway(self, pathway=None):
        
        if pathway == None:
            return self.bestpathway
        else:
            self.bestpathway = copy.deepcopy(pathway)
    
    
    #===========================================
    def Initialize(self, pool):
        '''
        Inicializa el hormiguero
        '''
        
        pRinitials = self.pRinitials.cumsum()
        
        idxs = [np.argwhere(pRinitials >= np.random.rand()).min() for n in range(self.size)]
        
        Rinitials = [self.Rinitials[idx] for idx in idxs]
        
        idx_relate = np.argwhere(self.compounds['relate'] > 0)
        
        idxs_initial = [pool.Cnames.index(self.Cinitials[idx]) for idx in idxs]
        
        idxs_produce = [ [int(idx) for idx in idx_relate if idx != idx_initial] for idx_initial in idxs_initial]
        
        for idx in range(self.size):
            
            self.ants[idx] = self.ants[idx].Initialize(pool, idxs_initial[idx], idxs_produce[idx], self.compounds, Rinitials[idx])

    
    
    #===========================================
    def Search(self, pool, jobs_server=''):
        '''
        Realiza una iteración de búsqueda (cada hormiga realiza la búsqueda una
        vez.
        
        --- MULTIPROCESSING ---
        jobs_server = multiprocessing.Pool(ncpus)
        '''
        
        if jobs_server == '':
            
            self.ants = [ant.Search(pool.S, pool.P, self.pheromones, pool.Cnames, pool.Rnames) for ant in self.ants]
            
            
        
        else:
            
            # MULTIPROCESSING
            self.ants = jobs_server.map(searching,[(ant, pool.S, pool.P, self.pheromones, pool.Cnames, pool.Rnames) for ant in self.ants])
            
    
    
    #===========================================
    def Size(self):
        '''
        Returns the anthill size
        '''
        return len(self.ants)
    
    
    
    #===========================================
    def UpdatePheromones(self):
        """EXPLICACION"""
        
        Values = [[idx,self.delta_pheromones.Get([idx[0]],[idx[1]])] for idx in self.delta_pheromones.StoredIdxs()]

        self.pheromones.Update(Values)
        
        #self.pheromones.Normalize()
    
    
    #===========================================
    def Update_pRinitials(self, pR_unforced=[], umbral=1E-6):
        """EXPLICACION"""
        
        # ACTUALIZO LOS COMPUESTOS FORZADOS
        #-----
        
        # ACTUALIZO LOS COMPUESTOS NO FORZADOS
        if self.compounds['initials_unforced']['p'].shape[0] > 1:
            
            self.compounds['initials_unforced']['p'] += pR_unforced
    
            idxs = np.argwhere(self.compounds['initials_unforced']['p'] < umbral).flatten().tolist()
            
            if idxs:
                
                Rinitials = [self.compounds['initials_unforced']['R'][idx] for idx in idxs]
                
                
                for Rinitial in Rinitials:
                    
                    idx = self.compounds['initials_unforced']['R'].index(Rinitial)
                    
                    # ACTUALIZO Rinitials
                    self.compounds['initials_unforced']['R'].remove(Rinitial)
                    
                    # ACTUALIZO Cinitials
                    self.compounds['initials_unforced']['C'].remove(self.compounds['initials_unforced']['C'][idx])
                    
                    
                    # ACTUALIZO pRinitials
                    self.compounds['initials_unforced']['p'] = np.delete(self.compounds['initials_unforced']['p'],idx)
                
            
            self.compounds['initials_unforced']['p'] /= np.sum(self.compounds['initials_unforced']['p'])
        
    
    
    #===========================================
    def Update_pCinitials(self, pCinitials=[], umbral=1E-6):
        """EXPLICACION"""
        
        if len(pCinitials) > 1:
            
            self.pCinitials *= ( 1.0 - pCinitials )
            
            idx = np.argwhere(self.pCinitials < umbral).flatten().tolist()
            
            if idx:
                self.pCinitials = np.delete(self.pCinitials, idx)
                self.Cinitials = [self.Cinitials[ii] for ii in range(len(self.Cinitials)) if ii not in idx]
                
                self.compounds['initials'] = set(self.Cinitials)
                self.compounds['produce'] = self.compounds['relate'] - self.compounds['initials']
            
            self.pCinitials /= np.sum(self.pCinitials)
    
    
    
    #===========================================
    def ExportData(self, pool, idx):
        '''
        '''
        
        data = self.ants[idx].ExportData(pool)
        
        data['pathway'] = pool.ShowReactions(self.ants[idx].filtered_semireactions, export=True)
        
        
        return data
    


#==============================================================================
class Settings(dict):
    """EXPLICACION"""
    
    def __init__(self, filename):
        """EXPLICACION"""
        
        #-----------------------------
        # OPEN SETTINGS FILE
        #-----------------------------
        with open(filename, 'r') as f:
            if (sys.version_info.major < 3) or (sys.version_info.minor < 6):
                settings = yaml.load(f)  # 2.7 and < 3.6
            else:
                settings = yaml.load(f, Loader=yaml.FullLoader)
        
        
        for key,value in settings.items():
            self[key] = value



#==============================================================================





#==============================================================================
class Pool(object):
    """
    EXPLICACION
    """
    
    
    def __init__(self, RXNfile, ENZfile=None, COMPNAMEfile=None, DTYPE=np.float32):
        
        
        self._code2compname = dict()
        self._enzymes = dict()
        
        self.Cnames = set()
        self.Rnames = list()
        
        #------------------------------------------------------------------
        # ENZYMES
        #----------------------
        ENZYMES = dict()
        
        if ENZfile is not None:
            
            with open(ENZfile, 'r') as fp:
                enzymes = fp.read()
                enzymes = enzymes.strip().split('\n')
                
            for enzyme in enzymes:
                R, ec_codes = enzyme.split(': ')
                self._enzymes[R] = ec_codes.split(', ')
        #------------------------------------------------------------------        
        
        
        #------------------------------------------------------------------
        # COMPOUND NAMES
        #----------------------
        if COMPNAMEfile is not None:
            
            with open(COMPNAMEfile, 'r') as fp:
                names = fp.read()
                names = names.strip().split('\n')
                
            for name in names:
                code, translations = name.split(': ')
                self._code2compname[code] = translations.split(', ')
        #------------------------------------------------------------------        
        
        
        
        #------------------------------------------------------------------
        # REACTIONS
        #----------------------
        SEMIREACTIONS = []
        
        with open(RXNfile,'r') as fp:
            reactions = fp.readlines()
        
        
        for reaction in reactions:
            
            R,equation = reaction.split(': ')
            
            if ' <=> ' in equation:
                
                left,right = equation.split(' <=> ')
                
                reverse_reaction = R + ': ' + right.strip() + ' --> ' + left.strip() + '\n'
                reactions.append(reverse_reaction)
                
                
                left = left.replace(' + ', ' ').replace('\n', '')
                right = right.replace(' + ', ' ').replace('\n', '')
            
            else:
                equation = equation.replace(' + ', ' ')
                equation = equation.replace('\n', '')
                left,right = equation.split(' --> ')
            
            
            
            Substrates = left.split(' ') # [Coeficiente - Compuesto - Coeficiente - Compuesto]
            Products = right.split(' ')  # [Coeficiente - Compuesto - Coeficiente - Compuesto]
            
            
            
            Coefficients = list()
            Compounds = list()
            S = list()
            P = list()
            
            for idx in range(0,len(Substrates),2):
                Coefficients.append(-int(Substrates[idx]))
                Compounds.append(Substrates[idx+1])
                S.append(Substrates[idx+1])
                
                self.Cnames.add(Substrates[idx+1])
                
                
            for idx in range(0,len(Products),2):
                Coefficients.append(int(Products[idx]))
                Compounds.append(Products[idx+1])
                P.append(Products[idx+1])
                
                self.Cnames.add(Products[idx+1])
            
            self.Rnames.append(R)
            
            SEMIREACTIONS.append({"reaction": R, "substrates": S, "products": P, "compounds": Compounds, "coefficients": Coefficients, "enzyme": ENZYMES[R] if R in ENZYMES.keys() else ['-.-.-.-']})
            
        
        #==================================
        # BUILD STOICHIOMETRIC MATRIX
        #==================================
        # COMPOUNDS
        #=============
        self.Cnames = list(self.Cnames)
        self.Cnames.sort()
        Nc = len(self.Cnames)
        #==================================
        
        
        #==================================
        # REACTIONS
        Nr = len(SEMIREACTIONS)
        
        self.size = Nr
        #============
        
        self.S = np.zeros((Nc,Nr), dtype=np.bool)
        self.P = np.zeros((Nc,Nr), dtype=np.bool)
        self.M = lil_matrix((Nc,Nr), dtype=np.float32)
        
        #
        
        for idx_col in range(len(SEMIREACTIONS)):
            
            if idx_col not in self._enzymes.keys():
                self._enzymes[R] = ['-.-.-.-']
            
            r = SEMIREACTIONS[idx_col]
            
            for idx_s,s in enumerate(r['substrates']):
                
                idx_row = self.Cnames.index(s)
                self.S[idx_row, idx_col] = True
                self.M[idx_row, idx_col] = r['coefficients'][idx_s]
                
            
            for idx_p,p in enumerate(r['products']):
                
                idx_row = self.Cnames.index(p)
                self.P[idx_row, idx_col] = True
                self.M[idx_row, idx_col] = r['coefficients'][idx_s+idx_p+1]
                
        
        del SEMIREACTIONS
        
    
    
    #===========================================
    def Enzyme(self, semireaction, onlyfirst=True):
        '''
        '''
        reaction = self.Rnames[semireaction]
        
        if (reaction not in self._enzymes):
            
            return ['-.-.-.-']
        
        else:
            
            if onlyfirst:
                return self._enzymes[reaction][0]
            else:
                return self._enzymes[reaction]
    
    #===========================================
    def ExternalCompounds(self):
        """EXPLICACION"""
        
        #idx_E = np.argwhere((self.P.sum(axis=1) > 0) != (self.S.sum(axis=1) > 0) )[:,0]
        
        idx_E = [idx for idx in range(self.P.shape[0]) if (self.P[idx,:].sum() == 0) and (self.S[idx,:].sum() > 0)]
        
        return [self.Cnames[c] for c in idx_E]
    
    
    
    #====================================================
    def FeasibleReactions(self, Co):
        '''
        '''
        
        R = np.argwhere(np.equal(self.S[Co,:].sum(axis=0), self.S.sum(axis=0)) == True)
        
#        R = np.argwhere(np.array((Co @ self.S) == self.S.sum(axis=0), dtype=self.S.dtype) == True)
        
        return R
#        idx = np.argwhere(((1.0/len(R)) * np.ones(len(R))).cumsum() >= np.random.rand()).min()
    
    
    
    
    #====================================================
    def ShowReactions(self, semireactions = list(), export=False, translate=True):
        """
        EXPLICACION
        """
        
        if not semireactions:
            semireactions = range(self.S.shape[0])
        
        linea = str()
        for semireaction in semireactions:
            
            reaction = self.Rnames[semireaction]
            
            enzymes = self.Enzyme(semireaction, onlyfirst=False)
            if len(enzymes) > 1:
                enzyme = enzymes[0] + ' +' + str(len(enzymes)-1)
            else:
                enzyme = enzymes[0]
            
            
            
            linea += reaction + ' (' + enzyme + "): "
            
            
            for idx_cpd in np.argwhere(self.S[:,semireaction] > 0)[:,0]:
                
                cname = self.Cnames[idx_cpd]
                
                if translate and (len(self._code2compname) > 0) and (cname in self._code2compname.keys()):
                    cname = self._code2compname[cname][0]
                
                coef = int(np.abs(self.M[idx_cpd,semireaction]))
                
                linea += str(coef) + ' ' + cname + " + "
        
        
            linea = linea[:-3] + " --> "
            
            
            for idx_cpd in np.argwhere(self.P[:,semireaction] > 0)[:,0]:
                
                cname = self.Cnames[idx_cpd]
                
                if translate and (len(self._code2compname) > 0) and (cname in self._code2compname.keys()):
                    cname = self._code2compname[cname][0]
                
                coef = int(np.abs(self.M[idx_cpd,semireaction]))
                
                linea += str(coef) + ' ' + cname + " + "
                    
                    
            linea = linea[:-3] + "\n"
            
            
            
        if not export:
            print(linea)
        else:
            return linea


#==============================================================================





#==============================================================================
class Pathway(object):
    """EXPLICACION"""
    
    def __init__(self, pool=[], semireactions = [], compounds=None):
        
        self.compounds = {'initials': set(),
                          'produce': set(),
                          'abundant': set(),
                          'external': set(),
                          'availables': set()
                         }
        
        
        self.cost = np.inf
        
        self.semireactions = list()
        self.reactions = list()
        
        self.enzymes = list()
        
        self.substrates = list()
        self.products = list()
        
        self.compounds['all'] = list()
        self.coefficients = list()
        
        
        if semireactions:
            
            self.cost = np.inf
            
            for r in semireactions:
                
                self.semireactions.append(r)
                self.reactions.append(pool.Rnames[r])
                self.enzymes.append(pool.Enzyme(r, onlyfirst=False))
                
                coefficients = []
                comp_all = []
                substrates = []
                products = []
                #--------------------
                # SUBSTRATES
                #--------------------
                idxS = np.argwhere(pool.S[:,r] > 0)
                
                for idxs in idxS:
                    
                    idxs = idxs[0]
                    
                    s = pool.Cnames[idxs]
                    
                    substrates.append(s)
                    comp_all.append(s)
                    coefficients.append(int(pool.M[idxs,r]))
                
                
                
                #--------------------
                # PRODUCTS
                #--------------------
                idxP = np.argwhere(pool.P[:,r] > 0)
                
                for idxp in idxP:
                    
                    idxp = idxp[0]
                    
                    p = pool.Cnames[idxp]
                    
                    products.append(p)
                    comp_all.append(p)
                    coefficients.append(int(pool.M[idxp,r]))
                    
                    
                self.substrates.append(substrates)
                self.products.append(products)
                self.coefficients.append(coefficients[:])
                self.compounds['all'].append(comp_all[:])
                
            
            
            if compounds is not None:
                
                self.compounds['initials'] = set()
                self.compounds['produce'] = set()
                self.compounds['abundant'] = compounds['abundant']
                self.compounds['external'] = compounds['external']
                self.compounds['availables'] = set()
        
   
    
    #===========================================
    def Size(self):
        """EXPLICACION"""
        
        return len(self.semireactions)


#==============================================================================



#==============================================================================
class Measures(dict):
    '''
    '''
    
    #===========================================
    def __init__(self, maxIter=1000, *keys):
        
        self.iteration = 0
        
        for key in keys:
            
            self[key] = list()
    
    
    #===========================================
    def summary(self, iteration=None, *keys):
        '''
        '''
        if iteration== None:
            Iteration = self.iteration
            
        else:
            Iteration = iteration
        
        
        if len(keys) == 0:
            keys = self.keys()
        
            
        M = np.zeros((len(keys),Iteration),dtype=np.float64)
        
        S = np.zeros((len(keys),4),dtype=np.float64)
        
        k = 0
        for key in keys:
            
            M[k,:] = self[key][:Iteration]
            k += 1
        
        
        # MEAN
        S[:,0] = np.mean(M,axis=1)
        
        # STD
        S[:,1] = np.std(M,axis=1)
        
        # MEDIAN
        S[:,2] = np.median(M,axis=1)
        
        # MAD
        S[:,3] = np.mean(np.absolute(M - np.reshape(np.repeat(np.mean(M,axis=1),Iteration), (len(keys),Iteration))), axis=1)
        
        print('%s \t %s \t %s \t %s \t %s') % ('  ','MEAN',' STD','MEDIAN','MAD')
        k = 0
        for key in keys:
            print('%s \t %0.4f \t %0.4f \t %0.4f \t %0.4f') % (key,S[k,0],S[k,1],S[k,2],S[k,3])
            k += 1

#==============================================================================




#==============================================================================
class Compounds(dict):
    """EXPLICACION"""
    
    def __init__(self, filename):
        """EXPLICACION"""
        
        #-----------------------------
        # OPEN COMPOUNDS FILE
        #-----------------------------
        with open(filename, 'r') as f:
            compounds = yaml.load(f, Loader=yaml.FullLoader)
        f.close()
        
        # CREATE COMPOUNDS SETS
        self['abundant'] = set(compounds['abundant'])
        
        self['relate'] = set()
        self['initials'] = set()
        
        
        for compound in compounds['relate']:
            self['relate'].add( compound['compound'] )
            
            if compound['initial']:
                
                self['initials'].add( compound['compound'] )




#==============================================================================





#==============================================================================
class SOLUTIONS(object):
    '''
    Permite guardar el historial de las soluciones que fueron obtenidas
    en cada iteración.
    '''
    
    
    #===========================================
    def __init__(self):
        '''
        '''
        
        self.N = 0.0                        # Numero de soluciones
        self.pathways = dict()              # [VIA METABOLICA (lista reacciones), size]
        self.history = dict()
        
    
    #===========================================
    def Update(self,key,N,iteration):
        '''
        key: pathway
        N: numero de reacciones
        iteration: iteracion evaluada
        '''
        
        if key in self.pathways.keys():
            
            if iteration in self.history[key]:
                
                self.history[key][iteration] += 1.0
            
            else:
                self.history[key] = {iteration: 1.0}
        
        else:
            self.N += 1.0
            
            self.pathways[key] = N
            
            self.history[key] = {iteration: 1.0}
        
    
    #===========================================
    def get_values(self, iteration):
        '''
        '''
        
        L = []   # L: Pathway size
        V = []   # V: Número de veces encontrada
        
        for key in self.pathways.keys():
            
            if iteration in self.history[key]:
                
                L.append(self.pathways[key])
                V.append(self.history[key][iteration])
        
        L = np.array(L)
        V = np.array(V)
        
        idxsort = L.argsort()
        
        L = L[idxsort].tolist()   # L: Pathway size
        V = V[idxsort].tolist()   # V: Número de veces encontrada        
        
        return (L,V)
    
    
    #===========================================
    def export(self, iteration=None):
        '''
        
        SOLUTIONS['pathways']   ---> KEY: SIZE
        SOLUTIONS['history']    ---> KEY: ITERATION
        
        '''
        
        SOLUTIONS = dict()
        
        if iteration != None:
            
            K = []
            L = []
            V = []
            
            for key in self.pathways.keys():
                
                if iteration in self.history[key]:
                    
                    K.append(key)                          # PATHWAY KEY
                    L.append(self.pathways[key])           # PATHWAY SIZE
                    V.append(self.history[key][iteration]) # ITERATION
            
            
            for idx in range(len(K)):
                SOLUTIONS['pathways'].update({K[idx]: L[idx]})
                SOLUTIONS['history'].update({K[idx]: V[idx]})
        
        else:
            SOLUTIONS['pathways'] = self.pathways   # --> 
            SOLUTIONS['history'] = self.history
        
        
        
        #SOLUTIONS['pathways'] = sorted(self.pathways.keys())
        #SOLUTIONS['history'] = sorted(self.history.keys())
        
        
        
        return SOLUTIONS
        
        
        
        
#==============================================================================






#==============================================================================
class Pheromones(object):
    '''
    '''
    
    #====================================================
    def __init__(self, size, default=1.0):
        '''
        "memory" indicates if must be constructed a pheromone matrix or 
        pheromone vector.
        
        "default" indicates the default value stored into the matrix. Since only
        a few values are different to the default one, just those values are
        stored in the matrix for an efficient usage of memory.
        '''
        
        self.default = default
        
        self.M = lil_matrix(size,dtype=np.float64)
    
    
    #====================================================
    def Clean(self, default=1.0):
        '''
        '''
        self.default = default
        idxs = self.M.nonzero()
        self.M[idxs] = 0.0
    
    
    #====================================================
    def Get(self, rows, cols):
        '''
        '''
        Values = []
        for row in rows:
            
            for col in cols:
                
                if self.M[row,col] == 0.0:
                    
                    Values.append(self.default)
                
                else:
                    Values.append(self.M[row,col])
        
        return np.array(Values)
    
    #====================================================
    def Set(self, row, col, value):
        '''
        Replace pheromone value with specified one.
        '''
        if self.M[row,col] == 0.0:
            self.M[row,col] = self.default + value
            
        else:
            self.M[row,col] = value
    
    
    #====================================================
    def Update(self, values):
        '''
        Add to the specified pheromone value the new one.
        values: Structure --> [(row,col),value]
        '''
        for idx,value in values:
            
            if self.M[idx] == 0.0:
                self.M[idx] = self.default + value
            
            else:
                self.M[idx] += value
    
    
    
    #===========================================
    def Evaporate(self, rho):
        '''
        EXPLICACION
        '''
        value = 1.0 - np.float64(rho)
        
        self.default *= value
        
        idxs = self.M.nonzero()
        self.M[idxs] *= value
    
    
    #===========================================
    def Normalize(self):
        '''
        Normalize pheromone matrix dividing by the maximum value.
        '''
        
        # MAXIMUM VALUE IN THE PHEROMONE MATRIX
        MAX = np.max(self.M.data.max())
        
        # DEFAULT VALUE NORMALIZATION
        self.default /= MAX
        
        # MATRIX NORMALIZATION
        idxs = self.M.nonzero()
        self.M[idxs] /= MAX
    
    
    #===========================================
    def StoredIdxs(self):
        '''
        Devuelve los índices de las celdas NO vacías.
        '''
        return zip(self.M.nonzero()[0], self.M.nonzero()[1])
    
    
    #===========================================
    def ExportData(self):
        '''
        Devuelve los índices de las celdas NO vacías.
        '''
        return self.M.todense().tolist()
 #==============================================================================   




#==============================================================================
class DynamicPathwayPlotter(object):
    '''
    Genera un gráfico interactivo de la vía metabólica.
    '''
    
    def __init__(self, pool, pathway):
        
        
        #---------------------------
        # NODES
        #---------------------------
        #{name = '', label = '', typenode = ''}

        #---------------------------
        # LINKS
        #---------------------------
        #{source = '', target = '', typelink = ''}
        
        
        self.Nnodes = 0
        self.Nlinks = 0
        
        self.nodes = []
        self.links = []
        
        
        pathway.compounds['source'] = pathway.compounds.pop('initials')
        
        pathway.compounds['target'] = pathway.compounds.pop('produce')
        
        pathway.compounds['abundant'] = list(set(pathway.compounds['abundant']) - set(pathway.compounds['source']))
        
        reactions = pathway.reactions[:]                                           # [pool.Reactions(r) for r in pathway.semireactions]
        enzymes = [pool.Enzyme(r, onlyfirst=False) for r in pathway.semireactions] # pathway.enzymes[:] #
        substrates = pathway.substrates[:]                                         # [pool.Substrates(r) for r in pathway.semireactions]
        products = pathway.products[:]                                             # [pool.Products(r) for r in pathway.semireactions]
        
        coefficients = pathway.coefficients[:]
        
        
        #--------------
        self.Nreactions = len(reactions)
        #--------------
        
        
        #-----------------------------
        S = set()
        for s in substrates:
            S.update(s)
        
        P = set()
        for p in products:
            P.update(p)
        
        # COMPOUNDS FILTERING
        Sf = S - set(pathway.compounds['abundant']) - set(pathway.compounds['external'])
        Pf = P - set(pathway.compounds['abundant']) - set(pathway.compounds['external'])
        
        # IDENTIFYING BACKBONE
        backbone = Pf.intersection(Sf)
        backbone.update(set(pathway.compounds['source']))
        backbone.update(set(pathway.compounds['target']))
        
        
        
#        C = P - S - set(pathway.compounds['target'])
        
        ##################################################################
        KEYS = list(pathway.compounds.keys())
        if 'Co' in KEYS:
            KEYS.remove('Co')
        KEYS.remove('availables')
        NODES = []
        #-----------------------------
        
        for idx in range(len(reactions)):
            
            # NODE "REACTION"
            NODES.append({'id': reactions[idx],
                          'code': reactions[idx],
                          'name': "",
                          "enzyme": enzymes[idx],
                          'type': 'reaction',
                          'function': '',
                          'backbone': 'yes',
                          'location': '',
                          'isbase': 'yes',
                          'image': '',
                          'mol_struct': '',
                          'pathways': '',
                          'url': ''#'http://www.brenda-enzymes.org/enzyme.php?ecno=' + enzymes[idx] if enzymes[idx] != '-.-.-.-' else ''
                          })
            
            
            
            #//////////////////////////////////////////////////////////////////
            for ii,substrate in enumerate(substrates[idx]):
                
                NAME = pool._code2compname[substrate] if substrate in pool._code2compname else substrate
                if isinstance(NAME,list) or isinstance(NAME,tuple):
                    NAME = NAME[0]
                
                node = {'id': substrate if substrate not in pathway.compounds['abundant'] else substrate + '_' + reactions[idx],
                        'code': substrate,
                        'name': NAME,
                        'type': 'compound',
                        'function': 'intermediate',
                        'backbone': 'yes' if substrate in backbone else 'no',
                        'substrate': 'yes',
                        'location': '',
                        'reaction': reactions[idx],
                        'image': '',
                        'mol_struct': '',
                        'pathways': '',
                        'stoichiometry': int(np.abs(coefficients[idx][ii])),  #int(np.abs(pool.M[pool.Cnames.index(substrate),pool.Rnames.index(reactions[idx])])),
                        'url': ''}
                
                #==============================================================
                for key in KEYS:
                    
                    if substrate in pathway.compounds[key]:
                        
                        node['function'] = key
                
                
                NODES.append(node)
            #//////////////////////////////////////////////////////////////////
            
            
            #//////////////////////////////////////////////////////////////////
            for ii,product in enumerate(products[idx]):
                
                NAME = NAME = pool._code2compname[product] if product in pool._code2compname else product
                if isinstance(NAME,list) or isinstance(NAME,tuple):
                    NAME = NAME[0]
                
                node = {'id': product if product not in pathway.compounds['abundant'] else reactions[idx] + '_' + product,
                        'code': product,
                        'name': NAME,
                        'type': 'compound',
                        'function': 'intermediate',
                        'backbone': 'yes' if product in backbone else 'no',
                        'substrate': 'no',
                        'location': '',
                        'reaction': reactions[idx],
                        'image': '',
                        'mol_struct': '',
                        'pathways': '',
                        'stoichiometry': int(np.abs(coefficients[idx][len(substrates[idx])+ii])),  #int(np.abs(pool.M[pool.Cnames.index(product),pool.Rnames.index(reactions[idx])])),
                        'url': ''}
                
                #==============================================================
                for key in KEYS:
                    
                    if product in pathway.compounds[key]:
                        
                        node['function'] = key
                
                
                NODES.append(node)
        #//////////////////////////////////////////////////////////////////
            
            
        #//////////////////////////////////////////////////////////////////
        IDs = []
        for idx in range(len(NODES)):
            
            node = NODES[idx].copy()
            
            #---------------------------------
            if node['id'] not in IDs:
                self.nodes.append(node)
                IDs.append(node['id'])
            #---------------------------------
            
            if (node['type'] == 'compound'):
                
                if (node['substrate'] == 'yes'):
                    
                    link = {
                            'backbone': node['backbone'],
                            'function': node['function'],
                            'source': '',
                            'source_name': node['id'],
                            'stoichiometry': node['stoichiometry'],
                            'target': '',
                            'target_name': node['reaction'],
                            'type': 'substrate'
                           }
                    
                
                else:# node['substrate'] == 'no'
                    
                    link = {
                            'backbone': node['backbone'],
                            'function': node['function'],
                            'source': '',
                            'source_name': node['reaction'],
                            'stoichiometry': node['stoichiometry'],
                            'target': '',
                            'target_name': node['id'],
                            'type': 'product'
                           }
                    
                
                self.links.append(link)
        #//////////////////////////////////////////////////////////////////
        self.Nnodes = len(self.nodes)
        self.Nlinks = len(self.links)
        
        
        for idx in range(self.Nlinks):
            
            self.links[idx]['source'] = [idx_node for idx_node in range(self.Nnodes) if self.nodes[idx_node]['id'] == self.links[idx]['source_name']][0]
            self.links[idx]['target'] = [idx_node for idx_node in range(self.Nnodes) if self.nodes[idx_node]['id'] == self.links[idx]['target_name']][0]
        
        
        
#        if added_structure:
#            with open('db/molecular_structures.json','w') as fp:
#                json.dump(molecular_structures,fp)
    
    
    #======================================
    def export_draw(self, FILENAME):
        '''
        '''
        pathway_filename = FILENAME + '.html'
        
        # OPEN THE HTML FILE
        with open(pathway_filename,'w') as fp:
            
            #==========================================================================
            fp.write('<!DOCTYPE html>\n')
            fp.write('<html lang="en">\n')
            fp.write('<meta charset="utf-8">\n')
            fp.write('    <head>\n')
            
            fp.write('        <title>PATHWAY</title>\n')
                    
            fp.write('        <meta charset="utf-8">\n')
#            fp.write('        <script src="https://d3js.org/d3.v3.min.js"></script>\n')
            
            
            #.........................
            # WRITE CSS
            #.........................
            fp.write('<style>\n')
            
            with open('lib/pathway.css','r') as css:
                Style = css.read()
            
            fp.write('%s' % Style)
            
            fp.write('</style>\n')
            
            #·····················································
            
            fp.write('<script type="text/javascript">\n')
            
            #-------------------------------
            with open('lib/d3.v3.min.js','r') as js:
                D3functions = js.read()
            
            fp.write('%s' % D3functions)
            #-------------------------------
            
            fp.write('</script>\n\n')
            
            #·····················································
            
            fp.write('    </head>\n\n')
            
            fp.write('    <body>\n')
                
            fp.write('    <noscript>\n')
            fp.write('    Para utilizar las funcionalidades completas de este sitio es necesario tener\n')
            fp.write('    JavaScript habilitado. Aquí están las <a href="http://www.enable-javascript.com/es/"\n')
            fp.write('    target="_blank"> instrucciones para habilitar JavaScript en tu navegador web</a>.\n')
            fp.write('    </noscript>\n\n')
            
            
            fp.write('    <div id="contenedor">\n')
            
            
            fp.write('        <div id="help_text" style="opacity: 1.0">\n')
            fp.write('            <p><b>Press "H" for help...</b></p>')
            fp.write('        </div>\n')
            
            fp.write('        <div id="help" style="visibility: hidden">\n')
            fp.write('        <p><b>HOTKEYS MENU</b></p>\n')
            fp.write('            <ul>\n')
            fp.write('                <li type="disc"> <i>h</i>: Displays or hides this menu.</li>\n')
            fp.write('                <li type="disc"> <i>a</i>: Displays or hides abundant compounds.</li>\n')
            fp.write('                <li type="disc"> <i>e</i>: Displays or hides external compounds.</li>\n')
            fp.write('                <li type="disc"> <i>u</i>: Displays or hides all compounds not in the main chain of reactions.</li>\n')
            fp.write('                <li type="disc"> <i>Shift + Click</i>: Hides links / compounds when mouse is over the element.</li>\n')
            fp.write('            </ul>\n')
            
            fp.write('            <hr width="90%">\n')
            fp.write('                <p><b>COLOR KEYS</b></p>\n')
            fp.write('                <ul>\n')
            fp.write('                    <li type="disc" style="color: blue;"> Enzymes or Reactions </li>\n')
            fp.write('                    <li type="disc" style="color: red;"> <i>Initial substrate</i></li>\n')
            fp.write('                    <li type="disc" style="color: #FFED1A;"> <i>Final products</i></li>\n')
            fp.write('                    <li type="disc" style="color: #35FFFF;"> <i>Intermediate products</i></li>\n')

            fp.write('                </ul>\n')
            fp.write('        </div>\n')
            
            fp.write('    </div>\n')
            
              
            fp.write('<script type="text/javascript">\n')
            
            #//////////////////////////////////////////////////////////////////
            # WRITE DATA IN JSON FORMAT
            #/////////////////////////////
            fp.write('    json = {\n')
            
            #----------------------------
            # NODES
            #----------------------------
            fp.write('        "nodes": [\n')
            
            NODES = self.nodes[:]
            
            while NODES:
                
                node = NODES.pop(0)
                
                fp.write('        {\n')
                #------------------------------------------------------------
                KEYS = sorted(node.keys())
                for key in KEYS:
                    
                    if (len(NODES) == 0) and (key == KEYS[-1]): # LAST VALUE
                        if isinstance(node[key], dict):
                            fp.write('            "%s": %s\n' % (key,json.dumps(node[key])) )
                        elif (key == 'enzymes'):
                            fp.write('            "%s": %s\n' % (key,node[key]) )
                        else:
                            fp.write('            "%s": "%s"\n' % (key,node[key]) )
                    else:
                        if isinstance(node[key], dict):
                            fp.write('            "%s": %s,\n' % (key,json.dumps(node[key])) )
                        elif (key == 'enzyme'):
                            fp.write('            "%s": %s,\n' % (key,node[key]) )
                        else:
                            fp.write('            "%s": "%s",\n' % (key,node[key]) )
                #------------------------------------------------------------
                fp.write('        },\n')
            
            fp.write('	],\n')
            
            
            
            #----------------------------
            # LINKS
            #----------------------------
            
            fp.write('    "links": [\n')
            
            LINKS = self.links[:]
            
            while LINKS:
                
                link = LINKS.pop(0)
                
                fp.write('        {\n')
                #------------------------------------------------------------
                KEYS = sorted(link.keys())
                for key in KEYS:
                    
                    if (len(LINKS) == 0) and (key == KEYS[-1]): # LAST VALUE
                        if isinstance(link[key], str):
                            fp.write('            "%s": "%s"\n' % (key,link[key]) )
                        else:
                            fp.write('            "%s": %s\n' % (key,link[key]) )
                    else:
                        if isinstance(link[key], str):
                            fp.write('            "%s": "%s",\n' % (key,link[key]) )
                        else:
                            fp.write('            "%s": %s,\n' % (key,link[key]) )
                #------------------------------------------------------------
                fp.write('        },\n')
            
            fp.write('	],\n')
            
            
            fp.write('    "Nnodes": %d,\n' % self.Nnodes)
            fp.write('    "Nlinks": %d\n' % self.Nlinks)
            
            fp.write('}\n')
            #//////////////////////////////////////////////////////////////////
            
            
            fp.write('    </script>\n')
            
            #..........................................................................
            
            fp.write('    <script type="text/javascript">\n')
            
            #-------------------------------
            with open('lib/pathway.js','r') as js:
                JSfunctions = js.read()
            
            fp.write('%s' % JSfunctions)
            #-------------------------------
            
            fp.write('    </script>\n\n')
            
            #..........................................................................
            
            fp.write('</body>\n')
            fp.write('</html>')
            #==========================================================================
        
        return pathway_filename


    #======================================
    def export_data(self):
        '''
        '''
        data = {
                "nodes": [],
                "links": [],
                "Nnodes": self.Nnodes,
                "Nlinks":self.Nlinks,
                "Size": self.Nreactions
               }
        
        #----------------------------
        # NODES
        #----------------------------
        NODES = self.nodes[:]
        
        while NODES:
            
            node = NODES.pop(0)
            
            #------------------------------------------------------------
            KEYS = sorted(node.keys())
            
            node2append = dict()
            
            for key in KEYS:
                
                if isinstance(node[key], dict):
                    node2append[key] = json.dumps(node[key])
                
                else:
                    node2append[key] = node[key]
            
            data["nodes"].append(node2append)
            #------------------------------------------------------------
        
        
        
        #----------------------------
        # LINKS
        #----------------------------
        LINKS = self.links[:]
        
        while LINKS:
            
            link = LINKS.pop(0)
            
            #------------------------------------------------------------
            KEYS = sorted(link.keys())
            
            link2append = dict()
            for key in KEYS:
                link2append[key] = link[key]
                
            data["links"].append(link2append)
        #==========================================================================
        
        return data
        
        
        

#==============================================================================

###############################################################################
