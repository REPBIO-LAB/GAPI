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
    '''
    consensusDb = fastaDir + 'consensusViralDb.fa'

    with open(virusDb, 'w') as outFile:
        with open(consensusDb) as inFile:
            outFile.write(inFile.read())
            '''
            

    ## 2. Index retrotransposon database fasta file ##
    

    index = outDir + '/virusDb.mmi'
    ## DESILENCIAAAAR
    '''
    err = open(logDir + '/index.err', 'w') 
    command = 'minimap2 -k 10 -w 1 -d ' + index + ' ' + virusDb 
    status = subprocess.call(command, stderr=err, shell=True)

    if status != 0:
        step = 'BUILD-VIRUS-DATABASE'
        msg = 'Database indexing failed' 
        log.step(step, msg)
        '''
        

    return virusDb, index


def buildIdentityDb(metacluster, db, outDir):

    # Coger la identity del primer evento DISCORDANT que aparezca (pq los clipping no tienen identity)
    identity = next(event.identity for event in metacluster.events if event.type == "DISCORDANT")
    
    specificIdentity = max(set([event.specificIdentity for event in metacluster.events if event.type == "DISCORDANT"]), key=[event.specificIdentity for event in metacluster.events if event.type == "DISCORDANT"].count)

    specificHeader = '"consensus|' + specificIdentity + '|' + identity + '"'

    #dbSpecificIdentity = outDir + '/' + str(CLIPPING_clusterID) + '_specificIdentity.fa'
    dbSpecificIdentity = outDir + '/' + str(metacluster.id) + '_specificIdentity.fa' 

    # Build database con identity + ref
    # TODO coger la identity que ya hemos reconocido, en lugar de toda la db!
    # 2. Cojo de la db de identities la secuencia que fue asignada como identity:
    # TODO poner bien el status y todo eso
    #err = open(logDir + '/index.err', 'w') 
    command = 'samtools faidx  ' + db + ' ' + specificHeader + ' -o ' + dbSpecificIdentity
    status = subprocess.call(command, shell=True)
    # indexo
    indexDbSpecificIdentity = dbSpecificIdentity.replace("_specificIdentity.fa", "_refIdentityDb.mmi")

    ## DESILENCIAAAAR
    # TODO
    # ponerlo bien
    #err = open(logDir + '/index.err', 'w')
    command = 'minimap2 -k 10 -w 1 -d ' + indexDbSpecificIdentity + ' ' + dbSpecificIdentity 
    status = subprocess.call(command, shell=True)

    if status != 0:
        step = 'BUILD-VIRUS-DATABASE'
        msg = 'Database indexing failed' 
        log.step(step, msg)

    return indexDbSpecificIdentity