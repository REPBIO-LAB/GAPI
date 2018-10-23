'''
Module 'log' - Contains functions to report log information
'''

## DEPENDENCIES ##
import time

## FUNCTIONS ##
def log(label, string):
    '''
        Display labelled information
    '''
    print('[' + label + ']', string)


def subHeader(string):
    '''
        Display  subheader
    '''
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print(timeInfo, '**', string, '**')


def info(string):
    '''
        Display basic information
    '''
    timeInfo = time.strftime("%Y-%m-%d %H:%M")
    print(timeInfo, string)

