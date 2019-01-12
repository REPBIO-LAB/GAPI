'''
Module 'unix' - Contains wrappers to unix commands
'''

## DEPENDENCIES ##
# External
import os
import sys

# Internal
import log

## FUNCTIONS ##

def mkdir(path):
    '''
    Create directory

    Input:
        1. path: directory to be created
    '''

    exist = os.path.isdir(path)

    # Only attempt to create directory if it does not exists
    if not exist: 
        try:  
            os.mkdir(path)

        except OSError:  
            step = 'ERROR'
            msg = "Creation of the directory %s failed" % path
            log.step(step, msg)
