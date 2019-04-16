'''
Module 'databases' - Contains functions and classes to deal with different types of databases
'''

## DEPENDENCIES ##
# External
import subprocess

# Internal
import unix 
# [SR CHANGE]
import log

## FUNCTIONS ##
def buildRetrotransposonDb(fastaDir, includeTransduced, outDir):
    '''
    Build database containing retrotransposon related sequences (consensus sequences, transduced regions, ...)

    Input:
        1. fastaDir: Directory containing reference fasta files (retrotransposon consensus sequences, transduced regions...)
        2. includeTransduced: Boolean to specify if include (True) or not (False) transduced regions from known source elements in the database
        3. outDir: Output directory

    Output:
        1. retrotransposonDb: Fasta file containing retrotransposon related sequences 
    '''

    ## 0. Create logs directory
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Create database fasta file ##
    retrotransposonDb = outDir + '/retrotransposonDb.fa'

    # A) Include retrotransposon consensus sequences + transduced regions
    if includeTransduced:

        consensusDb = fastaDir + 'consensusDb.fa'
        transducedDb = fastaDir + 'transducedDb.masked.fa'
        files = [consensusDb, transducedDb]

        with open(retrotransposonDb, 'w') as outFile:
            # Iterate over each fasta and write fasta into output database
            for f in files:
                with open(f) as inFile:
                    outFile.write(inFile.read())
        
    # B) Only include retrotransposon consensus sequences
    else:
        consensusDb = fastaDir + 'consensusDb.fa'

        with open(retrotransposonDb, 'w') as outFile:
            with open(consensusDb) as inFile:
                outFile.write(inFile.read())

    ## 2. Index retrotransposon database fasta file ##
    index = outDir + '/retrotransposonDb.mmi'
    err = open(logDir + '/index.err', 'w') 
    command = 'minimap2 -k 10 -w 1 -d ' + index + ' ' + retrotransposonDb 
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'BUILD-DATABASE'
        msg = 'Database indexing failed' 
        log.step(step, msg)

    return retrotransposonDb, index

def buildVirusDb(fastaDir, outDir):
    '''
    Build database containing retrotransposon related sequences (consensus sequences, transduced regions, ...)

    Input:
        1. fastaDir: Directory containing reference fasta files (retrotransposon consensus sequences, transduced regions...)
        2. includeTransduced: Boolean to specify if include (True) or not (False) transduced regions from known source elements in the database
        3. outDir: Output directory

    Output:
        1. retrotransposonDb: Fasta file containing retrotransposon related sequences 
    '''

    ## 0. Create logs directory
    logDir = outDir + '/Logs'
    unix.mkdir(logDir)

    ## 1. Create database fasta file ##
    virusDb = outDir + '/virusDb.fa'

    ## DESILENCIAAAAR
    
    consensusDb = fastaDir + 'consensusViralDb.fa'

    with open(virusDb, 'w') as outFile:
        with open(consensusDb) as inFile:
            outFile.write(inFile.read())

    ## 2. Index retrotransposon database fasta file ##
    

    index = outDir + '/virusDb.mmi'
    ## DESILENCIAAAAR
    
    err = open(logDir + '/index.err', 'w') 
    command = 'minimap2 -k 10 -w 1 -d ' + index + ' ' + virusDb 
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'BUILD-VIRUS-DATABASE'
        msg = 'Database indexing failed' 
        log.step(step, msg)
        

    return virusDb, index


def buildRefIdentityDb(ref, beg, end, identity, specificIdentity, identityDbFasta, reference, side, outDir):

    region = ref +':'+ str(beg) +'-'+ str(end)

    # TODO: poner bien el outfile:
    outFileRefRegion = outDir + '/' + str(beg) + '_' + side + '_reference_region.fa' 
    # 1. Cojo el trozo de la ref que me interesa:
    # TODO poner bien el status y todo eso
    #err = open(logDir + '/index.err', 'w') 
    command = 'samtools faidx  ' + reference + ' ' + region + ' -o ' + outFileRefRegion
    status = subprocess.call(command, shell=True)
    outFilespecificIdentity = outDir + '/' + str(beg) + '_' + side + '_specificIdentity.fa' 
    specificHeader = '"consensus|' + specificIdentity + '|' + identity + '"'
    # 2. Cojo de la db de identities la secuencia que fue asignada como identity:
    # TODO poner bien el status y todo eso
    #err = open(logDir + '/index.err', 'w') 
    command = 'samtools faidx  ' + identityDbFasta + ' ' + specificHeader + ' -o ' + outFilespecificIdentity
    status = subprocess.call(command, shell=True)
    refIdentityDb = outDir + '/' + str(beg) + '_' + side + '_refIdentityDb_region.fa' 
    command = 'cat  ' + outFileRefRegion + ' ' + outFilespecificIdentity + ' > ' + refIdentityDb
    status = subprocess.call(command, shell=True)
    #if status != 0:
        #step = 'BUILD-VIRUS-DATABASE'
        #msg = 'Database indexing failed' 
        #log.step(step, msg)
    
    # indexo
    refIdentityDbIndex = outDir + '/' + str(beg) + '_' + side + '_refIdentityDb_region.mmi' 
    ## DESILENCIAAAAR
    # TODO
    # ponerlo bien
    #err = open(logDir + '/index.err', 'w')
    command = 'minimap2 -k 10 -w 1 -d ' + refIdentityDbIndex + ' ' + refIdentityDb 
    status = subprocess.call(command, shell=True)

    
    if status != 0:
        step = 'BUILD-VIRUS-DATABASE'
        msg = 'Database indexing failed' 
        log.step(step, msg)
        

    return refIdentityDbIndex