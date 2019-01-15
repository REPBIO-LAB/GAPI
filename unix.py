'''
Module 'unix' - Contains wrappers to unix commands
'''

## DEPENDENCIES ##
# External
import os
import sys

# Internal
import log
import unix

## FUNCTIONS ##

def mkdir(path):
    '''
    Create directory

    Note: improve function to be able to create lists of directories

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

def rm(path):
    '''
    Delete file

    Note: improve function to be able to delete lists of files and directories

    Input:
        1. path: path to file to be deleted
    '''

    exist = os.path.isfile(path)

    # Only attempt to delete file if it does exists
    if exist: 
        try:  
            os.remove(path)

        except OSError:  
            step = 'ERROR'
            msg = "Deletion of the file %s failed" % path
            log.step(step, msg)