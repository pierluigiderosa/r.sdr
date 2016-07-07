#!/usr/bin/env python
#
############################################################################
#
# MODULE:      r.sdr
# AUTHOR(S):   Pierluigi De Rosa
#
# PURPOSE:     compure the SDR from DEM using several equations
# COPYRIGHT:   (C) 2016, pierluigi de rosa
#
#              This program is free software under the GNU General Public
#              License (>=v2). Read the file COPYING that comes with GRASS
#              for details.
#
#############################################################################

#%Module
#% description: Calculate the SDR raster map from a DEM
#% keyword: display
#% keyword: raster
#%End
#%option
#% key: demraster
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of DEM map
#% required : yes
#%end
#%option
#% key: weightmap
#% type: string
#% gisprompt: old,cell,raster
#% description: Name of weighted raster map
#% required : yes
#%end
# %option G_OPT_R_OUTPUT
# % description: Name for output raster map
# % guisection: Output
# %end

import sys
import grass.script as grass
import grass.script.array as garray
import numpy
import math


import time

def main():
    dem = options['demraster']
    dem = 'dem_tinitaly_rocchetta'
    #SDR = options['output']
    weightmap = options['weightmap']


    #region setting
    gregion=grass.region()
    cell_s=gregion['nsres']
    cell_s=float(cell_s)



    # r.slope.aspect
    # elevation = dem_tinitaly_rocchetta @ SDR
    # slope = slope
    rasterTemp = []
    vectTemp = []
    grass.run_command('r.slope.aspect', elevation=dem, slope="slope1", overwrite=True)
    rasterTemp.append('slope1')

    grass.run_command('r.watershed', flags='s', elevation=dem, accumulation='accD8', drainage='drainD8', overwrite=True)
    rasterTemp.append('accD8')
    rasterTemp.append('drainD8')

    grass.run_command('r.slope.aspect',elevation = dem,slope = 'slope',format ='percent',overwrite=True)
    rasterTemp.append('slope')

    # tif_fdir8 coincide con drainD8

    # read drainage direction map
    tif_fdir8_ar = garray.array()
    tif_fdir8_ar.read('drainD8')
    # converto il float
    tif_fdir8_ar = tif_fdir8_ar.astype(numpy.float)  # otherwise overflow in future operations
    # r.watershead: Negative numbers indicate that those cells possibly have surface runoff from outside of the current geographic region.
    tif_fdir8_ar[(tif_fdir8_ar <= 0)] = 0
    ndv = numpy.min(tif_fdir8_ar)
    tif_fdir8_ar[tif_fdir8_ar == ndv] = numpy.NaN

    # create constant array to trasform into raster
    const_ar = tif_fdir8_ar * 0 + cell_s

    ### zero matrix bigger than F_dir8, to avoid border indexing problems
    # sorrounding tif_fdir8_ar with one width zeros cells
    Fd8 = numpy.zeros(shape=((tif_fdir8_ar.shape[0]) + 1, (tif_fdir8_ar.shape[1]) + 1), dtype=numpy.float32)
    # popolo la matrice
    Fd8[1:Fd8.shape[0], 1:Fd8.shape[1]] = Fd8[1:Fd8.shape[0], 1:Fd8.shape[1]] + tif_fdir8_ar
    # adding bottom row and right y axis with zeros
    Fdir8 = numpy.zeros(shape=((Fd8.shape[0]) + 1, (Fd8.shape[1]) + 1), dtype=numpy.float32)
    Fdir8[:Fdir8.shape[0] - 1, :Fdir8.shape[1] - 1] = Fd8
    ##------------
    # read weight map an slope
    tif_wgt_ar = garray.array()
    #TODO controllare la mappa weight che va presa da input
    tif_wgt_ar.read('weight')
    tif_slope = garray.array()
    tif_slope.read('slope')

    tif_slope=tif_slope/100. #converting percentage from r.slope.aspect to value in range 0 - 1

    # imposing upper and lower limits to slope, no data here are -1
    tif_slope[(tif_slope >= 0) & (tif_slope < 0.005)] = 0.005
    tif_slope[(tif_slope > 1)] = 1
    tif_slope[(tif_slope < 0)] = -1

    #imposing a value bigger than zero in weight map
    tif_wgt_ar[tif_wgt_ar==0]=1e-10

    Ws_1 = 1 / (tif_wgt_ar * tif_slope)
    # converto il float
    Ws_1 = Ws_1.astype(numpy.float)  # otherwise overflow in future operations
    # r.watershead: Negative numbers indicate that those cells possibly have surface runoff from outside of the current geographic region.
    # tif_fdir8_ar[(tif_fdir8_ar <= 0)] = 0
    ndv = numpy.min(Ws_1)
    Ws_1[Ws_1 == ndv] = numpy.NaN
    #
    # zero matrix bigger than weight, to avoid border indexing problems, and have same indexing as Fdir8
    Wg = numpy.zeros(shape=((tif_wgt_ar.shape[0]) + 1, (tif_wgt_ar.shape[1]) + 1), dtype=numpy.float32)
    # TODO da sostituire la variabile con Ws_1 ovvero il denom di Ddn
    Wg[1:Wg.shape[0], 1:Wg.shape[1]] = Wg[1:Fd8.shape[0],
                                       1:Wg.shape[1]] + Ws_1  # the weigth to weigth tha flow length
    # adding bottom row and right y axis with zeros
    Wgt = numpy.zeros(shape=((Wg.shape[0]) + 1, (Wg.shape[1]) + 1), dtype=numpy.float32)
    Wgt[:Wgt.shape[0] - 1, :Wgt.shape[1] - 1] = Wg
    #
    start = time.clock()  # for computational time
    # Creating a bigger matrix as large as weight(and all the matrices) to store the weighted flow length values
    W_Fl = numpy.zeros(shape=((Wgt.shape[0]), (Wgt.shape[1])), dtype=numpy.float32)
    W_Fl = W_Fl - 1  # to give -1 to NoData after the while loop calculation
    #
    # Let's go for the search and algo-rhytm for the weighted-Flow-Length
    ND = numpy.where(numpy.isnan(Fdir8) == True)  # fast coordinates all the NoData values, starting from them to go forward and compute flow length
    #
    Y = ND[0]  # rows, NoData indexes
    X = ND[1]  # columns, NoData indexes pay attention not to invert values !!!!!!!!!!!!!!
    #
    # initializing lists for outlet and moving cell coordinates, in function of their position
    YC1 = []
    YC2 = []
    YC3 = []
    YC4 = []
    YC5 = []
    YC6 = []
    YC7 = []
    YC8 = []
    XC1 = []
    XC2 = []
    XC3 = []
    XC4 = []
    XC5 = []
    XC6 = []
    XC7 = []
    XC8 = []
    #
    #   Flow Directions r.watershead
    #   4   3   2
    #   5   -   1
    #   6   7   8
    #
    #   Draining in Direction Matrix
    #   8   7   6
    #   1   -   5
    #   2   3   4
    #
    i1 = Fdir8[Y, X - 1]  # Searching for NoData with cells draining into them, 8 directions
    D1 = numpy.where(i1 == 1)  # l
    YC1.extend(Y[D1])  # coordinates satisfacting the conditions
    XC1.extend(X[D1])
    W_Fl[YC1, XC1] = 0  # initialize flow length at cells draining to NoData
    #
    i2 = Fdir8[Y + 1, X - 1]  # Searching for NoData with cells draining into them, 8 directions
    D2 = numpy.where(i2 == 2)  # lrad2
    YC2.extend(Y[D2])  # coordinates satisfacting the conditions
    XC2.extend(X[D2])
    W_Fl[YC2, XC2] = 0  # initialize flow length at cells draining to NoData
    #
    i3 = Fdir8[Y + 1, X]  # Searching for NoData with cells draining into them, 8 directions
    D3 = numpy.where(i3 == 3)  # l
    YC3.extend(Y[D3])  # coordinates satisfacting the conditions
    XC3.extend(X[D3])
    W_Fl[YC3, XC3] = 0  # initialize flow length at cells draining to NoData
    #
    i4 = Fdir8[Y + 1, X + 1]  # Searching for NoData with cells draining into them, 8 directions
    D4 = numpy.where(i4 == 4)  # lrad2
    YC4.extend(Y[D4])  # coordinates satisfacting the conditions
    XC4.extend(X[D4])
    W_Fl[YC4, XC4] = 0  # initialize flow length at cells draining to NoData
    #
    i5 = Fdir8[Y, X + 1]  # Searching for NoData with cells draining into them, 8 directions
    D5 = numpy.where(i5 == 5)  # l
    YC5.extend(Y[D5])  # coordinates satisfacting the conditions
    XC5.extend(X[D5])
    W_Fl[YC5, XC5] = 0  # initialize flow length at cells draining to NoData
    #
    i6 = Fdir8[Y - 1, X + 1]  # Searching for NoData with cells draining into them, 8 directions
    D6 = numpy.where(i6 == 6)  # lrad2
    YC6.extend(Y[D6])  # coordinates satisfacting the conditions
    XC6.extend(X[D6])
    W_Fl[YC6, XC6] = 0  # initialize flow length at cells draining to NoData
    #
    i7 = Fdir8[Y - 1, X]  # Searching for NoData with cells draining into them, 8 directions
    D7 = numpy.where(i7 == 7)  # l
    YC7.extend(Y[D7])  # coordinates satisfacting the conditions
    XC7.extend(X[D7])
    W_Fl[YC7, XC7] = 0  # initialize flow length at cells draining to NoData
    #
    i8 = Fdir8[Y - 1, X - 1]  # Searching for NoData with cells draining into them, 8 directions
    D8 = numpy.where(i8 == 8)  # lrad2
    YC8.extend(Y[D8])  # coordinates satisfacting the conditions
    XC8.extend(X[D8])
    W_Fl[YC8, XC8] = 0  # initialize flow length at cells draining to NoData
    #
    #start =time.clock()#da cancellare poi.....!!!!!! Solo per check
    count = 1  # "0" passage already done during the previous step
    while len(YC1) or len(YC2) or len(YC3) or len(YC4) or len(YC5) or len(YC6) or len(YC7) or len(YC8) > 0:
        # Converting into array to be able to do operations
        YYC1=numpy.asarray(YC1);XXC1=numpy.asarray(XC1)
        YYC2=numpy.asarray(YC2);XXC2=numpy.asarray(XC2)
        YYC3=numpy.asarray(YC3);XXC3=numpy.asarray(XC3)
        YYC4=numpy.asarray(YC4);XXC4=numpy.asarray(XC4)
        YYC5=numpy.asarray(YC5);XXC5=numpy.asarray(XC5)
        YYC6=numpy.asarray(YC6);XXC6=numpy.asarray(XC6)
        YYC7=numpy.asarray(YC7);XXC7=numpy.asarray(XC7)
        YYC8=numpy.asarray(YC8);XXC8=numpy.asarray(XC8)
        #
        # Now I can do operations and moving towards the right cell!!!!!!!!
        # Weigthing flow length, weights are half sum of pixels weight * travelled length
        # I'm chosing the directions accordingly to Flow_dir step by step going from outlet-nodata to the ridges,
        # each time account for distance (l or l*rad2) multiplied by the half of the weigths of the 2 travelled cells.
        # Then, with variables substitution I'm moving a step further, and adding the prevous pixel value to the new calculated.
        #
        YYC1 = (YYC1);XXC1 = (XXC1 - 1)  # l
        YYC2 = (YYC2 + 1);XXC2 = (XXC2 - 1)  # lrad2
        YYC3 = (YYC3 + 1);XXC3 = (XXC3)  # l
        YYC4 = (YYC4 + 1);XXC4 = (XXC4 + 1)  # lrad2
        YYC5 = (YYC5);XXC5 = (XXC5 + 1)  # l
        YYC6 = (YYC6 - 1);XXC6 = (XXC6 + 1)  # lrad2
        YYC7 = (YYC7 - 1);XXC7 = (XXC7)  # l
        YYC8 = (YYC8 - 1);XXC8 = (XXC8 - 1)  # lrad2
        #
        if count == 1:  # first run zero, like TauDEM, need to check if there is a Nodata pixel receiving flow for all the 8 directions
            if len(YYC1) > 0:
                W_Fl[YYC1, XXC1] = 0
            else:
                pass
            if len(YYC2) > 0:
                W_Fl[YYC2, XXC2] = 0
            else:
                pass
            if len(YYC3) > 0:
                W_Fl[YYC3, XXC3] = 0
            else:
                pass
            if len(YYC4) > 0:
                W_Fl[YYC4, XXC4] = 0
            else:
                pass
            if len(YYC5) > 0:
                W_Fl[YYC5, XXC5] = 0
            else:
                pass
            if len(YYC6) > 0:
                W_Fl[YYC6, XXC6] = 0
            else:
                pass
            if len(YYC7) > 0:
                W_Fl[YYC7, XXC7] = 0
            else:
                pass
            if len(YYC8) > 0:
                W_Fl[YYC8, XXC8] = 0
            else:
                pass
        else:
            W_Fl[YYC1, XXC1] = W_Fl[YC1, XC1] + (cell_s * ((Wgt[YC1, XC1] + Wgt[YYC1, XXC1]) / 2))
            W_Fl[YYC2, XXC2] = W_Fl[YC2, XC2] + (cell_s * math.sqrt(2) * ((Wgt[YC2, XC2] + Wgt[YYC2, XXC2]) / 2))
            W_Fl[YYC3, XXC3] = W_Fl[YC3, XC3] + (cell_s * ((Wgt[YC3, XC3] + Wgt[YYC3, XXC3]) / 2))
            W_Fl[YYC4, XXC4] = W_Fl[YC4, XC4] + (cell_s * math.sqrt(2) * ((Wgt[YC4, XC4] + Wgt[YYC4, XXC4]) / 2))
            W_Fl[YYC5, XXC5] = W_Fl[YC5, XC5] + (cell_s * ((Wgt[YC5, XC5] + Wgt[YYC5, XXC5]) / 2))
            W_Fl[YYC6, XXC6] = W_Fl[YC6, XC6] + (cell_s * math.sqrt(2) * ((Wgt[YC6, XC6] + Wgt[YYC6, XXC6]) / 2))
            W_Fl[YYC7, XXC7] = W_Fl[YC7, XC7] + (cell_s * ((Wgt[YC7, XC7] + Wgt[YYC7, XXC7]) / 2))
            W_Fl[YYC8, XXC8] = W_Fl[YC8, XC8] + (cell_s * math.sqrt(2) * ((Wgt[YC8, XC8] + Wgt[YYC8, XXC8]) / 2))
            #
        #
        # Reconstructing all X and Y of this step and moving on upwards (Downstream if you think in GIS, right?)
        YY = [];XX = []
        YY.extend(YYC1);XX.extend(XXC1)
        YY.extend(YYC2);XX.extend(XXC2)
        YY.extend(YYC3);XX.extend(XXC3)
        YY.extend(YYC4);XX.extend(XXC4)
        YY.extend(YYC5);XX.extend(XXC5)
        YY.extend(YYC6);XX.extend(XXC6)
        YY.extend(YYC7);XX.extend(XXC7)
        YY.extend(YYC8);XX.extend(XXC8)
        #
        YY = numpy.asarray(YY)
        XX = numpy.asarray(XX)
        #
        i1 = Fdir8[YY, XX - 1]  # Searching for cells draining into them, 8 directions
        D1 = numpy.where(i1 == 1)  # l
        YC1 = YY[D1]  # coordinates satisfacting the conditions, HERE i NEED TO ADD ACTUAL LENGTH VALUE + PREVIOUS ONE
        XC1 = XX[D1]
        #
        i2 = Fdir8[YY + 1, XX - 1]  # Searching for cells draining into them, 8 directions
        D2 = numpy.where(i2 == 2)  # lrad2
        YC2 = YY[D2]  # coordinates satisfacting the conditions
        XC2 = XX[D2]
        #
        i3 = Fdir8[YY + 1, XX]  # Searching for cells draining into them, 8 directions
        D3 = numpy.where(i3 == 3)  # l
        YC3 = YY[D3]  # coordinates satisfacting the conditions
        XC3 = XX[D3]
        #
        i4 = Fdir8[YY + 1, XX + 1]  # Searching for cells draining into them, 8 directions
        D4 = numpy.where(i4 == 4)  # lrad2
        YC4 = YY[D4]  # coordinates satisfacting the conditions
        XC4 = XX[D4]
        #
        i5 = Fdir8[YY, XX + 1]  # Searching for cells draining into them, 8 directions
        D5 = numpy.where(i5 == 5)  # l
        YC5 = YY[D5]  # coordinates satisfacting the conditions
        XC5 = XX[D5]
        #
        i6 = Fdir8[YY - 1, XX + 1]  # Searching for cells draining into them, 8 directions
        D6 = numpy.where(i6 == 6)  # lrad2
        YC6 = YY[D6]  # coordinates satisfacting the conditions
        XC6 = XX[D6]
        #
        i7 = Fdir8[YY - 1, XX]  # Searching for cells draining into them, 8 directions
        D7 = numpy.where(i7 == 7)  # l
        YC7 = YY[D7]  # coordinates satisfacting the conditions
        XC7 = XX[D7]
        #
        i8 = Fdir8[YY - 1, XX - 1]  # Searching for cells draining into them, 8 directions
        D8 = numpy.where(i8 == 8)  # lrad2
        YC8 = YY[D8]  # coordinates satisfacting the conditions
        XC8 = XX[D8]
        count = count + 1
    #
    elapsed = (time.clock() - start)  # computational time
    print time.strftime("%d/%m/%Y %H:%M:%S    "), "Process concluded succesfully \n", "%.2f" % elapsed, 'seconds for Weighted-Flow Length calculation with ', int(count), ' iterations'  # truncating the precision

    W_fl = W_Fl[1:W_Fl.shape[0] - 1, 1:W_Fl.shape[1] - 1]#reshaping weigthed flow length, we need this step to homogenize matrices dimensions!!!!!!!!!!
    del W_Fl
    #imposto il valore di zero a 1 per evitare divisioni per zero
    D_down_ar = garray.array()
    W_fl[W_fl == 0] = 1
    D_down_ar[...] = W_fl
    del W_fl
    D_down_ar.write('w_flow_length',null=numpy.nan,overwrite=True)

    """
    --------------------------------------
    WORKING ON D_UP COMPONENT
    --------------------------------------
    """

    grass.run_command('r.watershed',elevation = 'dem',accumulation = 'accMDF',convergence=5, memory=300)
    rasterTemp.append('accMDF')

    tif_dtmsca = garray.array()
    tif_dtmsca.read('acc_watershead_dinf')

    tif_dtmsca = abs(tif_dtmsca)*cell_s

    acc_final_ar = tif_dtmsca / const_ar

    grass.run_command('r.watershed',elevation = 'dem',flow = 'weight',accumulation = "accW",convergence = 5,memory = 300)
    rasterTemp.append('accW')

    acc_W_ar = garray.array()
    acc_W_ar.read('accW')


    grass.run_command('r.watershed', elevation='dem', flow='slope', accumulation="accS", convergence=5, memory=300)
    rasterTemp.append('accS')

    acc_S_ar = garray.array()
    acc_S_ar.read('accS')

    # Computing C_mean as (accW+weigth)/acc_final
    C_mean_ar = (acc_W_ar + tif_wgt_ar) / acc_final_ar
    del (acc_W_ar)  # free memory
    #
    # Computing S mean (accS+s)/acc_final
    S_mean_ar = (acc_S_ar + tif_fdir8_ar) / acc_final_ar
    del (acc_S_ar, tif_fdir8_ar)  # free memory
    #
    # Computing D_up as "%cmean.tif%" * "%smean.tif%" * SquareRoot("%ACCfinal.tif%" * "%resolution.tif%" * "%resolution.tif%")
    cell_area = (const_ar) ** 2  # change of variables, to be sure
    D_up_ar = C_mean_ar * S_mean_ar * numpy.sqrt(acc_final_ar * cell_area)  # to transform from unit values to square units
    #
    # Computing Connectivity index
    ic_ar = numpy.log10(D_up_ar / D_down_ar)

    SDRmax = 0.8;IC0=0.5;k=1
    SDRmap = SDRmax / (1+math.exp((IC0-ic_ar/k)))







'''
    CLEANING TEMP FILE
    '''
    for rast in rasterTemp:
        grass.run_command('g.remove', flags='f', type='raster', name=rast)


if __name__ == "__main__":
    options, flags = grass.parser()
    sys.exit(main())
