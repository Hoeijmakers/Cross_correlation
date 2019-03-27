def construct_outlier_mask(order,w,c_thresh):
    print('THIS DOESNT EXIST YET')
    return(order*0.0+1.0)



def clean_order(order,w,thresh,mask=None):
    import numpy as np
    import lib.utils as ut
    ut.typetest('order',order,np.ndarray)
    ut.typetest('w',w,float)
    ut.typetest('thresh',thresh,float)
    shape=np.shape(order)

    #THIS IS A WORK IN PROGRESS. MAY NOT NEED TO CLEAN THE SPECTRUM BEFORE DOING THE CCF...
    #ONLY TO CONSTRUCT AN OUTLIER MASK PERHAPS?
