'''
Module 'virus' - for dealing with virus specific needs
'''

## Internal
import bamtools

def is_virusSR(eventsDict, tumourBam, normalBam):

    ## 1. Collect mate sequence of discordant events ##
    bamtools.collectMatesSeq(eventsDict, tumourBam, normalBam)