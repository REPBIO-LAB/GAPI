'''
Module 'virus' - for dealing with virus specific needs
'''
## External
import subprocess

## Internal
import bamtools
import formats

def is_virusSR(events, tumourBam, normalBam, outDir):

    ## 1. Collect mate sequence of discordant events ##
    bamtools.collectMatesSeq(events, tumourBam, normalBam, True, 20)