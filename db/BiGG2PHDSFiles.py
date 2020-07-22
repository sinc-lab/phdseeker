#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  5 09:30:44 2018

@author: mgerard
"""
import os
import json
import argparse


parser = argparse.ArgumentParser(description='Build PhDSeeker files from BiGG JSON model.')
parser.add_argument('-f', '--file', help='BiGG JSON model file.')
parser.add_argument('-r', '--allreversible', default='False', type=(lambda x: x.lower() in ("yes", "true", "t", "1")), help='Boolean value indicating if reactions should be taken always as reversible.')


args = vars(parser.parse_args())

#====================================================
def BiGG2PHDSFiles(bigg_filename, allreversible):
    '''
    '''
    
    #-----------------
    # READING FILE
    #-----------------
    with open(bigg_filename, 'r') as fp:
        data = json.load(fp)


    #------------------------
    # COMPOUNDS DICTIONARY
    #------------------------
    NAMES = ''
    for cpd in data['metabolites']:
        NAMES += cpd['id'] + ': ' + cpd['name'] + '\n'



    #--------------------------------------
    # REACTIONS, ENZYMES and REVERSIBILITY
    #--------------------------------------
    REACTIONS = ''
    ENZYMES = ''

    for R in data['reactions']:
        
        if ('BIOMASS' not in R['id'].upper()) and ('biomass' not in ' '.join(R['metabolites'].keys())):
            
            reaction = None
            reversibility = None
            ec_codes = None
            
            #--------------
            # ENZYMES
            #--------------
            if ('annotation' in R.keys()) and R['annotation']:
                
                for annotation in R['annotation']:
                    
                    key = annotation[0]
                    link = annotation[1]
                    
                    if 'EC Number' in key:
                        
                        if ec_codes is None:
                            ec_codes = R['id'] + ': ' + link.split('/')[-1]
                        else:
                            ec_codes = ec_codes + ', ' + link.split('/')[-1]
                    
            
            
            #--------------
            # REACTION
            #--------------
            S = ''
            P = ''
            
            
            for comp, coef in R['metabolites'].items():
                
                if coef < 0:
                    
                    S += str(int(abs(coef)) ) + ' ' + comp + ' + '
                
                if coef > 0:
                    
                    P += str(int(abs(coef)) ) + ' ' + comp + ' + '
            
            
            #----------------
            # REVERSIBILITY
            #----------------
            forward = False
            reverse = False
            
            if R['upper_bound'] > 0:
                forward = True
            
            if R['lower_bound'] < 0:
                reverse = True
                

            if len(S) > 0 and (len(P) > 0):
                
                if allreversible:
                    REACTIONS += R['id'] + ': ' + S[:-3] + ' <=> ' + P[:-3] + '\n'
                    
                else:
                    
                    if forward and reverse:
                        REACTIONS += R['id'] + ': ' + S[:-3] + ' <=> ' + P[:-3] + '\n'
                    
                    elif forward and not reverse:
                        REACTIONS += R['id'] + ': ' + S[:-3] + ' --> ' + P[:-3] + '\n'
                        
                    else:
                        REACTIONS += R['id'] + ': ' + P[:-3] + ' --> ' + S[:-3] + '\n'
                
                
                if ec_codes is not None:
                    ENZYMES += ec_codes + '\n'



    #---------------------

    path,filename = os.path.split(bigg_filename)

    filename_Names = 'NAMES_' + filename
    filename_Names = filename_Names.replace('.json','.txt')

    filename_Reacions = 'REACTIONS_' + filename
    filename_Reacions = filename_Reacions.replace('.json','.txt')

    filename_Enzymes = 'ENZYMES_' + filename
    filename_Enzymes = filename_Enzymes.replace('.json','.txt')


    with open(os.path.join(path, filename_Names), 'w') as fp:
        fp.writelines(NAMES)


    with open(os.path.join(path, filename_Reacions), 'w') as fp:
        fp.writelines(REACTIONS)


    with open(os.path.join(path, filename_Enzymes), 'w') as fp:
        fp.writelines(ENZYMES)



#==============================================================================
if __name__ == "__main__":
#==============================================================================

    BiGG2PHDSFiles(args['file'], args['allreversible'])
