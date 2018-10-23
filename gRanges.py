'''
Module 'gRanges' - Contains functions to report log information
'''

## DEPENDENCIES ##
import pysam

## FUNCTIONS ##
def overlap(begA, endA, begB, endB):
    '''
    Check if two ranges overlap. 2 criteria for defining overlap: 
    A) Begin of the range A within the range B         
           *beg* <---------range_A---------->                         
       <---------range_B----------> 
                
          *beg* <-------range_A----->
      <-------------range_B------------------>
    
    B) Begin of the range B within the range A     
          <---------range_A----------> 
                  *beg* <---------range_B---------->
            
          <-------------range_A----------------->
            *beg* <-------range_B------>
    '''    
       
    # a) Begin of the range A within the range B   
    if ((begA >= begB) and (begA <= endB)):
        overlap = True
        
    # b) Begin of the range B within the range A            
    elif ((begB >= begA) and (begB <= endA)):
        overlap = True

    # c) Ranges do not overlapping
    else:
        overlap = False

    return overlap
        

