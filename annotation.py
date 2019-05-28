'''
Module 'annotation' - Contains functions for the annotation of genomic intervals according to different annotation resources
'''

def annotate_repeats(events, repeatsDb, buffer):
    '''
    For each input event assess if overlaps with an annotated repeat in the reference genome

    Input: 
        1. events: list containing input events to be annotated. Events should be objects containing ref, beg and end attributes.
        2. repeatsDb: dictionary containing annotated repeats organized per chromosome (keys) into genomic bins (values)
        3. buffer: number of base pairs to extend begin and end coordinates for each event prior assessing overlap

    Output:
        1) events: list of input events with a new 'repeats' attribute set for each event. Attribute contains a not redundant list of overlapping repeats
    '''
    ## Assess for each input event if it overlaps with an annotated repeat
    for event in events:

        # A) Annotated repeat in the same ref where the event is located
        if event.ref in repeatsDb:
            
            ## Select repeats bin database for the corresponding reference 
            repeatsBinDb = repeatsDb[event.ref]        

            ## Retrieve all the annotated repeats overlapping with the event interval
            overlaps = repeatsBinDb.collect_interval(event.beg - buffer, event.end + buffer, 'ALL')    

            ## Make list of overlapping repeats
            repeats = [overlap[0].name for overlap in overlaps]

        # B) No repeat in the same ref as the event
        else:
            repeats = []
        
        ## Generate repeats string and add to the input event
        # A) Event overlapping repeat
        if repeats:
            event.repeats = ','.join(set(repeats))

        # B) Event NOT overlapping repeat
        else:    
            event.repeats = None
