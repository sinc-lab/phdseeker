# -*- coding: utf-8 -*-
"""
Created on Thu Feb  5 15:26:02 2015

@author: mgerard
"""

#-------------------
# IMPORTING MODULES
#-------------------
import numpy

try:
    import multiprocessing
except:
    pass

import timeit
import datetime

import sys
sys.path.append('lib/')
import libPhDSeeker as phds

import json


##===========================================

#def PhDSeeker(reactions_file,compounds_file,settings_file):
def PhDSeeker(settings_file):    
    
    #--------------
    # INITIALIZING
    #--------------
    tic = timeit.default_timer
    start_initialization = tic()
    
    
    #----------------
    # LOAD SETTINGS
    #----------------
    print("\nLoading SETTINGS...")
    settings = phds.Settings(settings_file)
    
    
    
    #----------------
    # LOAD SETTINGS
    #----------------
    #print("\nLoading SETTINGS...")
#    settings = phds.Settings(settings_file)
    
#    settings['NCORES'] = 1 # NUMBER OF CORES
#    settings['AntsInTheSamePath'] = 1.10
#    
#    if settings['Nants'] > 10:
#        settings['Nants'] = 10
#    
#    if settings['maxIterations'] > 500:
#        settings['maxIterations'] = 500
#    
#    if settings['IterationsWithoutChanges'] > settings['maxIterations']:
#        settings['IterationsWithoutChanges'] = settings['maxIterations']
#    
#    if settings['IterationsWithAlignedAnts'] > settings['maxIterations']:
#        settings['IterationsWithAlignedAnts'] = settings['maxIterations']
#    
#    settings['REACTIONS'] = reactions_file
#    settings['COMPOUNDS'] = compounds_file
    
    
    #-----------------------------
    # INSTANTANEOUS MEASURES
    #-----------------------------
    MEASURES = dict()
    MEASURES['Cost'] = list()
    MEASURES['NRT'] = list()
    MEASURES['NRU'] = list()
    MEASURES['Conectivity'] = list()
    MEASURES['History'] = list()
    
    
    
    #-------------------------
    #STARTING PARALLEL SERVER
    #-------------------------
    if settings['NCORES'] == 1:
        jobs_server = ''
    else:
        jobs_server = multiprocessing.Pool(settings['NCORES'])

    time_settings = tic()
    loadig_parameters = time_settings - start_initialization
    
    print("Elapsed time: {:.4} sec\n".format(str(loadig_parameters) if loadig_parameters > 0.01 else '0.01'))
    
    
    #----------------
    # LOAD COMPOUNDS
    #----------------
    print("Loading COMPOUNDS...")
    
    COMPOUNDS = settings['COMPOUNDS']
    
    compounds = {'abundant': set(COMPOUNDS['abundant']),
                 'external': set(),
                 'relate': set(),
                 'initials': set()
                }
    
    
    for c in COMPOUNDS['relate']:
        compounds['relate'].add( c['compound'] )
        
        if c['initial']:
            compounds['initials'].add( c['compound'] )
    
    
    time_compounds = tic()
    loadig_compounds = time_compounds - time_settings
    print("Elapsed time: {:.4} sec\n".format(str(loadig_compounds) if loadig_compounds > 0.01 else '0.01'))
    

    
    
    #-----------------
    # INITIALIZE POOL
    #-----------------
    print("Initializing POOL...")
    pool = phds.Pool(settings['REACTIONS'], settings['ENZYMES'], settings['COMPNAMES'])
    time_pool = tic()
    print(str(pool.size) + " reactions loaded")
    print("Elapsed time: {:.4} sec\n".format(str(time_pool - time_compounds) if (time_pool - time_compounds) > 0.01 else '0.01'))
    
    
    
    #-----------------------------
    # IDENTIFY EXTERNAL COMPOUNDS
    #-----------------------------
    if settings['AllowExternalCompounds']:
        compounds['external'] = set(pool.ExternalCompounds())
    else:
        compounds['external'] = set([])
    
    
    #------------------
    # BUILDING ANTHILL
    #------------------
    anthill = phds.Anthill(pool=pool,
                            Nants=settings['Nants'],
                            Rmax=pool.size,
                            compounds=compounds,
                            strict_initialization=settings['StrictInitialization']
                          )
    
    
    
    #--------------------------
    start_search = tic()
    
    delta_pheromones = numpy.zeros_like(anthill.pheromones)
    
    SOLUTIONS = phds.SOLUTIONS()
    
    
    iteration = 0
    STOP = False
    MAX_SOLUTION_FOUND = 9
    MIN_SOLUTION_FOUND = 10
    MaxL = 0
    ITER_MAX = 0
    
    results = list()
    
    PATHWAYS_HISTORY = []
    
    ###########################################################################
    while (iteration <= settings['maxIterations']) and (not STOP):
        
        if settings['Verbose']:
            print("Iteración: " + str(iteration) + "\n")
        
        #---------------------------------------------------
        # INITIALIZING ANTHILL - N CPUs ---> MULTIPROCESSING
        anthill.Initialize(pool) #,jobs_server)
        #---------------------------------------------------
        
        delta_pheromones.fill(0)
        
        #-----------------------------------------
        # SEARCHING - N CPUs ---> MULTIPROCESSING
        #-----------------------------------------
        anthill.Search(pool,jobs_server)

        
        N = 0
        UPDATED = False
        IDX = None
        
        for ant in anthill.ants:
            
            if (ant.found == True):
                
                # --------------------------
                semireactions_filtered = ant.filtered_semireactions[:]
                semireactions = list(set(semireactions_filtered))
                semireactions.sort()
                solution = str(semireactions)
                
                SOLUTIONS.Update(solution,len(semireactions_filtered),iteration)
                
                # --------------------------
                
                N += 1
                
                measures = ant.Cost(pool.S, pool.P, pool.Rnames)
                
                
                # --------------------------
                # UPDATE DELTA PHEROMONES
                # --------------------------
                
                value = numpy.float64(1.0/measures['COSTO'])
                
                delta_pheromones[semireactions_filtered[:-1],semireactions_filtered[1:]] += value
                
                
                # --------------------------
                # UPDATE BEST SOLUTION
                # --------------------------
                if  measures['COSTO'] < anthill.bestcost:
                    
                    anthill.bestcost = numpy.float64(measures['COSTO'])
                    anthill.bestpathway = ant.ExportPathway(pool)
                    
                    IDX = ant.idx
                    
                    UPDATED = True
                    used_time = tic() - start_search
                
                
                # -------------------------------------------------------------
                # UPDATING INITIAL REACTIONS PROBABILITIES
                # -------------------------------------------------------------
                anthill.pRinitials[anthill.Rinitials.index(ant.Rinitials)] += ( 1.0 / ( measures['COSTO'] ) )
                
        
        
        if UPDATED:
            MEASURES['History'].append((iteration,used_time,anthill.bestcost,anthill.bestpathway.Size()))
            
            PATHWAYS_HISTORY.append((iteration, used_time, anthill.bestcost, anthill.bestpathway.semireactions, anthill.bestpathway))
        
        
        
        
        
        #--------------------------
        # ANALISYS OF SOLUTIONS
        #--------------------------
        L, heights = SOLUTIONS.get_values(iteration)
        Ncol = MAX_SOLUTION_FOUND-MIN_SOLUTION_FOUND + 1
        solutions = list()
        if L:
            if numpy.max(L) > MAX_SOLUTION_FOUND:
                MAX_SOLUTION_FOUND = int(numpy.max(L))
                
            if numpy.min(L) < MIN_SOLUTION_FOUND:
                MIN_SOLUTION_FOUND = int(numpy.min(L))
            
            Ncol = MAX_SOLUTION_FOUND - MIN_SOLUTION_FOUND + 1
            
            
            # GENERO UNA MATRIZ CON TANTAS COLUMNAS COMO SOLUCIONES TENGA,
            # Y TANTAS FILAS COMO VARIANTES DE UNA MISMA SOLUCION TENGA
            solutions = numpy.zeros((Ncol,2),dtype=numpy.float64)
            
            if L:
                
                l = numpy.array(L)
                h = numpy.array(heights)
                
                for n in numpy.arange(Ncol):
                    
                    solutions[n,0] = n+MIN_SOLUTION_FOUND
                    solutions[n,1] = numpy.sum(h[l == n+MIN_SOLUTION_FOUND])
            
            
            idx = numpy.argmax(solutions[:,1])
            maxL = solutions[idx,0]
            maxHeight = solutions[idx,1]
        
        
        # --------------------------
        # STOP CRITERIA
        # --------------------------
        if ((SOLUTIONS.N > 0) and (len(SOLUTIONS.get_values(iteration)[1]) > 0) and (settings['Nants'] == max(SOLUTIONS.get_values(iteration)[1]))):
            
            if (maxL,maxHeight) == (MaxL,settings['Nants']):
                ITER_MAX += 1
            else:
                MaxL = maxL
                ITER_MAX = 0
            
            if ITER_MAX == settings['IterationsWithAlignedAnts']:
                STOP = True
                print('It was reached the maximum number of iterations with aligned ants.')
            
        elif (MEASURES['History'] and iteration-MEASURES['History'][-1][0] >= float(settings['IterationsWithoutChanges'])):
            STOP = True
            print('It was reached the allowed maximum number of iterations without changes on the best solution.')
        
        #======================================================================
        
        
        
        
        
        #----------------------------------------------
        # UPDATING PHEROMONES (EVAPORATION + ADDITION)
        #----------------------------------------------
        anthill.pheromones *= (1.0 - settings['rho'])
        
        if int(N) > 0:
            
            anthill.pheromones += delta_pheromones
        
        
        
        #--------------------------
        # EXPORT NEW BEST SOLUTION
        #--------------------------
        if IDX != None:
            data = anthill.ExportData(pool,IDX)
            data['iteration'] = iteration
            results.append(data)
        
        #==========================================================
        
        
        
        
        if settings['Verbose']:
            
            if anthill.bestpathway.semireactions:
                pool.ShowReactions(anthill.bestpathway.semireactions)
            else:
                print('Pathway not found.\n')
            
            print("----------------------------")
            print("Iteration: {}".format(str(iteration)))
            print("Best Cost: {:.4}".format(str(anthill.bestcost)))
            print("Size: {} reactions".format(str(anthill.bestpathway.Size())))
            print("============================\n\n")
            
            timeit.time.sleep(0.01)
        
        
        iteration += 1
    ###########################################################################
    
    
    #----------------------
    # SAVE SEARCHING TIME
    #----------------------
    end_search = tic()
    execution_time = end_search - start_search
    
    if execution_time < 60:
        print('Execution time: {:.5} sec\n'.format(str(execution_time) if (execution_time) > 0.01 else '0.01'))
        
    elif execution_time < 3600:
        print('Execution time: {} [mm:ss]\n'.format(str(datetime.timedelta(seconds=numpy.round(execution_time)))))
                                                        
    else:
        print('Execution time: {} [hh:mm:ss]\n'.format(str(datetime.timedelta(seconds=numpy.round(execution_time)))))
    
    #print("Execution time: {} [hh:mm:ss]".format(datetime.timedelta(seconds=numpy.round(execution_time))))
    #print("Execution time: {:.4} sec\n".format(str(execution_time) if execution_time > 0.01 else '0.01'))
    
    
    #----------------------
    # SAVE TOTAL TIME
    #----------------------
    total_time = end_search-start_initialization
    
    
    
        
    print("----------------------------")
    print("\nBEST SOLUTION FOUND\n")
    print("----------------------------\n")
    
    if anthill.bestpathway.semireactions:
        pool.ShowReactions(anthill.bestpathway.semireactions)
    else:
        print('Pathway not found.\n')
    
    print("----------------------------")
    print("Iteration: {}".format(str(iteration)))
    print("Best Cost: {:.4}".format(str(anthill.bestcost)))
    print("Best Size: {} reactions".format(str(anthill.bestpathway.Size())))
    end_search = tic()
    
    _time = end_search - start_search
    if _time < 60:
        print('Elapsed time: {:.5} sec'.format(str(_time) if (_time) > 0.01 else '0.01'))
        
    elif _time < 3600:
        print('Elapsed time: {} [mm:ss]'.format(str(datetime.timedelta(seconds=numpy.round(_time)))))
                                                        
    else:
        print('Elapsed time: {} [hh:mm:ss]'.format(str(datetime.timedelta(seconds=numpy.round(_time)))))
    
    #print("Elapsed time: {} [hh:mm:ss]".format(datetime.timedelta(seconds=numpy.round(end_search - start_search))))
    #print('Elapsed time: {:.4} sec\n'.format(str(end_search - start_search) if (end_search - start_search) > 0.01 else '0.01'))
    
    print("============================\n\n\n\n")
    
    
    
    print("========================================")
    print("SUMMARY")
    print("=========\n")
    
    
    print("ITERATIONS: {}".format(str(iteration)))
    
    if total_time < 60:
        print('TOTAL TIME: {:.5} sec'.format(str(total_time) if (total_time) > 0.01 else '0.01'))
              
    elif total_time < 3600:
        print('TOTAL TIME: {} [mm:ss]'.format(str(datetime.timedelta(seconds=numpy.round(total_time)))))
                                                  
    else:
        print('TOTAL TIME: {} [hh:mm:ss]'.format(str(datetime.timedelta(seconds=numpy.round(total_time)))))
    
    #print("TOTAL TIME: {} [hh:mm:ss]".format(numpy.round(total_time,2) if total_time < 3600 else datetime.timedelta(seconds=numpy.round(total_time))))
    #print("TOTAL TIME: {:.5} sec".format(str(total_time) if total_time > 0.01 else '0.010'))
    
    print("BEST COST: {:.4}".format(str(MEASURES['History'][-1][2])))
    print("BEST PATHWAY SIZE: {}".format(str(MEASURES['History'][-1][3])))
    print("\n========================================")
    print("BEST PATHWAY FOUND")
    print("===================\n")
    
    if anthill.bestpathway.semireactions:
        pool.ShowReactions(anthill.bestpathway.semireactions)
    else:
        print('Pathway not found.\n')
    
    print("========================================\n\n")
    
    
    
    #-----------------
    # SAVE RESULTS
    #-----------------
    DBname = settings['REACTIONS'].split('/')[-1].split('.txt')[0].split('_')[1]
    
    FILENAME = timeit.time.strftime("%Y%m%d-%H%M%S__") + DBname + '__'
    
    
    #------------------------
    # BUILDING REPORT
    #------------------------
    output = ''
    
    # 
    #------------------------
    output += '##################################\n'
    output += 'EXPERIMENTAL SETTINGS:\n'
    output += '######################\n\n'
    
    output += '\tNumber of cores: ' + str(settings['NCORES']) + '\n'
    output += '\tNumber of ants: ' + str(settings['Nants']) + '\n'
    output += '\tPheromone evaporation rate: ' + str(settings['rho']) + '\n'
    output += '\tMaximum number of iterations: ' + str(settings['maxIterations']) + '\n'
    output += '\tIterations without changes: ' + str(settings['IterationsWithoutChanges']) + '\n'
    output += '\tIterations with aligned ants: ' + str(settings['IterationsWithAlignedAnts']) + '\n'
    output += '\tStrict initialization: ' + str(settings['StrictInitialization']) + '\n'
    output += '\tAllow external compounds: ' + str(settings['AllowExternalCompounds']) + '\n\n\n'
    
        
    #------------------------
    output += '#################\n'
    output += 'SUMMARY:\n'
    output += '###########\n\n'
    
    output += '\tIterations performed: ' + str(iteration - 1) + '\n'
    
    
    if total_time < 60:
        output += '\tSearching time: {:.4} sec\n'.format(str(total_time) if (total_time) > 0.01 else '0.01')
        
    elif total_time < 3600:
        
        output += '\tSearching time: {} [mm:ss]\n'.format(str(datetime.timedelta(seconds=numpy.round(total_time))))
    else:
        output += '\tSearching time: {} [hh:mm:ss]\n'.format(str(datetime.timedelta(seconds=numpy.round(total_time))))
    
    
    output += '\tPathway size: ' + str(len(anthill.bestpathway.semireactions)) + ' reactions\n\n\n'
    
    
    #------------------------
    output += '#################\n'
    output += 'COMPOUNDS:\n'
    output += '###########\n\n'
    
    
    
    #-----------------------------------------------------------------
    Cinitial = anthill.bestpathway.compounds['initials'][0]
    
    comp = ''
    if Cinitial in pool._code2compname.keys():
        comp += Cinitial + '('+ pool._code2compname[Cinitial][0] + ') '
    else:
        comp += Cinitial + ' '
    
    output += '\tInitial: ' + comp.strip() + '\n\n'
    #-----------------------------------------------------------------
    
    
    #-----------------------------------------------------------------
    Cproduce = anthill.bestpathway.compounds['produce']
    
    comp = ''
    for c in Cproduce:
        if c in pool._code2compname.keys():
            comp += c + '('+ pool._code2compname[c][0] + ') '
        else:
            comp += c + ' '
            
    output += '\tProduce: ' + comp.strip() + '\n\n'
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    compounds = sorted(anthill.bestpathway.compounds['abundant'])
    
    if Cinitial in compounds:
        compounds.remove(Cinitial)
    
    comp = ''
    for c in compounds:
        
        if c in pool._code2compname.keys():
            comp += c + '('+ pool._code2compname[c][0] + ') '
        else:
            comp += c + ' '
    
    output += '\tAbundant: ' + comp.strip() + '\n\n'
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    comp = ''
    compounds = sorted(anthill.bestpathway.compounds['external'])
    
    for c in compounds:
        if c in pool._code2compname.keys():
            comp += c + '('+ pool._code2compname[c][0] + ') '
        else:
            comp += c + ' '
    
    output += '\tExternal: ' + comp.strip() + '\n\n\n'
    #-----------------------------------------------------------------
    
    output += '########################################################\n'
    output += 'PATHWAY:\n'
    output += '###########\n\n'
    
    if list(pool._code2compname.keys()):
        output += '============================================\n'
        output += 'COMMON NAME NOTATION\n'
        output += '--------------------------------------------\n'
        output += 'Reaction (Enzyme): Substrates --> Products\n'
        output += '============================================\n\n'
        pathway = pool.ShowReactions(semireactions=anthill.bestpathway.semireactions, export=True)
        output += pathway  + '\n'
        
    else:
        
        output += '============================================\n'
        output += 'CODE NOTATION\n'
        output += '--------------------------------------------\n'
        output += 'Reaction (Enzyme): Substrates --> Products\n'
        output += '============================================\n\n'
        pathway = pool.ShowReactions(semireactions=anthill.bestpathway.semireactions, export=True)
        output += pathway  + '\n'
    
    
    
    filename_summary = 'out/' + FILENAME + 'SUMMARY.txt'
    filename_history = 'out/' + FILENAME + 'HISTORY.txt'
    filename_html = 'out/' + FILENAME + 'PATHWAY.html'
    
    
    #------------------------
    # SAVING REPORT
    #------------------------
    with open( filename_summary, 'w' ) as fp:
        fp.write(output)
        
    
    #------------------------
    # SAVING HISTORY
    #------------------------
    with open( filename_history, 'w' ) as fp:
        
        output = ''
        
        if PATHWAYS_HISTORY:
            K = 1
            for data in PATHWAYS_HISTORY[::-1]:
                
                iteration = data[0]
                elapsed_time = data[1]
                cost = data[2]
                semireactions = data[3]
                
                output += '########################################################\n'
                output += 'PATHWAY ' + str(K) + ':\n'
                output += '###########\n\n'
                
                output += '============================================\n'
            
                if list(pool._code2compname.keys()):
                    output += 'COMMON NAME NOTATION\n'
                else:
                    output += 'CODE NOTATION\n'
                
                
                output += '--------------------------------------------\n'
                output += 'Reaction (Enzyme): Substrates --> Products\n'
                output += '============================================\n\n'
                pathway = pool.ShowReactions(semireactions=semireactions, export=True)
                output += pathway  + '\n\n'
                
                
                
                output += "----------------------------\n"
                output += "Iteration: {}\n".format(str(iteration))
                output += "Cost: {:.4}\n".format(str(cost))
                output += "Size: {} reactions\n".format(str(len(semireactions)))
                
                if elapsed_time < 60:
                    output += 'Elapsed time: {:.4} sec\n'.format(str(elapsed_time) if (elapsed_time) > 0.01 else '0.01')
                elif elapsed_time < 3600:
                    output += 'Elapsed time: {} [mm:ss]\n'.format(str(datetime.timedelta(seconds=numpy.round(elapsed_time))))
                else:
                    output += 'Elapsed time: {} [hh:mm:ss]\n'.format(str(datetime.timedelta(seconds=numpy.round(elapsed_time))))
                
                output += "============================\n\n\n\n"
                
                K += 1
            
        fp.write(output)
        
    
    
    #------------------------
    # SAVE PATHWAY PLOT
    #------------------------
    #DPP = phds.DynamicPathwayPlotter(pool,anthill.bestpathway)
    #DPP.export_draw(filename_html[:-5])
    
    data = []
    options = '<option value="0" selected>Select pathway...</option>\n'
    K = 0
    for _solution in PATHWAYS_HISTORY[::-1]:
        
        K += 1
        
        _pathway = _solution[-1]
        semireactions = _solution[3]
        
        DPP = phds.DynamicPathwayPlotter(pool, _pathway)
        
        p = DPP.export_data()
        
        data.append(p)
        
        options += '<option value="{}">Pathway {} ({} reactions)</option>\n'.format(K, K, len(semireactions))
    
    
    with open('lib/template.html', 'r') as fp:
        HTML = fp.read()
    
    
    #---------------------
    # ADD OPTIONS IN MENU
    HTML = HTML.replace('<!--INSERT_PATHWAY_MENU_HERE-->', options)
    
    # ADD JSON DATA
    HTML = HTML.replace('/*INSERT_DATA_HERE*/', 'data={}'.format(json.dumps(data)))
    
    
    with open(filename_html,'w') as fp:
        fp.write(HTML)
    
    #--------------------------------------------------------------------------
    
    
    return(filename_summary,filename_history,filename_html)




#==============================================================================
if __name__ == "__main__":
#==============================================================================

    import sys
    
    settings = ''
    
    for ii in range(1,len(sys.argv),2):
        
        if sys.argv[ii] == '--settings':
            settings = sys.argv[ii+1]
        
        else:
            print('Parámetro desconocido.')
    
    PhDSeeker(settings)
