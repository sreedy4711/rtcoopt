# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 2017

@author: Reedy
"""
#import modules section
import os  #enables change directory command
import csv #enables reading and writing of csv files
from cvxopt import matrix, solvers #matrix needed to create matrices, solvers needed to actually solve stuff
import time
import sys  #imports reading system arguments
import glob
from datetime import date
from datetime import timedelta
import rtcoopt_step5a


#############################################
############################################
def fasave_file(data, filename):
    with open(filename, 'w',newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        for datum in data:
            spamwriter.writerow([datum[0],datum[1],datum[2]])

def fasave_grd_file(data, filename):
    with open(filename, 'w',newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        for datum in data:
            spamwriter.writerow(datum)

def faget_resource_data(filename):
    startload = time.time()
    badstatuslist = ['OUT','EMR','OFF','OFFNS','Telemetered Resource Status'] #shouldn't be included in optimization
    optdata = {}#initialize dictionary containing needed information for optimiazation
    optdata['TimeStamps'] = set() #initialize set of all timestamps for the dictionary
    optdata['dispatch'] = [] #initialize the dispatch list
    optdata['prices'] = {} #initializes the price list
    optdata['60daygrd'] = [] #initialize output file
    PR = range(75,145,2) # iterable of the 35 different offer segments
    MW = range(74,144,2) # iterable of the 35 different offer segments
    with open(filename, 'r') as csvfile: #open the 60 day SCED report
        offerfilereader = csv.reader(csvfile)
        rowcounter = 0
        for row in offerfilereader: #get a row of data
            rowcounter += 1
            if rowcounter % 10000 == 0:
                print (rowcounter) #just so you have something to watch
            try:
                status = row[151] #look for last blank row
            except IndexError:
                print(rowcounter) #troubleshooting; shouldn't happen  Later stuff will error out
            optdata['60daygrd'].append(row) #build output file
            if status not in badstatuslist: #only look at stuff that should be optimized
                TimeStamp = (row[0],row[1])
                optdata['TimeStamps'].add(TimeStamp)
############ Section that gets HDL, LDL, BP data for the row
                if status == 'ONTEST':
                    HDL = float(row[152])
                    LDL = float(row[152])
                else:
                    HDL = float(row[147])
                    LDL = float(row[150])
                if HDL < 0.05: HDL = 0
                BP = float(row[152])
                try:
                    optdata[TimeStamp]['GTBD'] += BP
                    optdata[TimeStamp]['LDLsum'] += LDL
                except KeyError: #If keyError, optdata['TimeStamp'] hasn't been defined yet, need to initialize
                    counter = 0
                    optdata[TimeStamp] = {}
                    optdata[TimeStamp]['GTBD'] = BP
                    optdata[TimeStamp]['counter'] = 0
                    optdata[TimeStamp]['Resources'] = set()
                    optdata[TimeStamp]['LDLsum'] = LDL
                optdata[TimeStamp][row[2]] = LDL
                optdata[TimeStamp]['Resources'].add(row[2])
                for point in range (1,35,1): #for each segment, grab PROLD, PRNEW, MWOLD, MWNEW (defining points of segment)
                    MWNEW = float(row[MW[point]])
                    MWOLD = float(row[MW[point - 1]])
                    PRNEW = float(row[PR[point]])
                    PROLD = float(row[PR[point - 1]])
                    if (row[2] == 'WOODWRD2_WOODWRD2' and TimeStamp == '7/31/2017 0:00'):
                        print(MWOLD,MWNEW,PROLD,PRNEW)
                    if (MWNEW > MWOLD and MWNEW > LDL and MWOLD < HDL) : #check to se if this is a segment that is within HDL/LDL and should be considered/optimized
                      optdata[TimeStamp][counter] = {} #initialize segment dictionarey
                      if MWOLD < LDL: #check to see if need to adjust MWOLD and PROLD for LDL
                          PROLD = PROLD + (LDL - MWOLD) * (PRNEW - PROLD) / (MWNEW - MWOLD) #linearly interpolated price
                          MWOLD = LDL
                      if MWNEW > HDL: #check to see if need to adjust MWNEW and PRNEW for HDL
                          PRNEW = PRNEW + (HDL - MWNEW) * (PRNEW - PROLD) / (MWNEW - MWOLD) #linearly interpolated price
                          MWNEW = HDL
                      if MWNEW == MWOLD : MWNEW = MWNEW + 0.0001
#### Load up segment dictionary with data
                      optdata[TimeStamp][counter]['MWOLD'] = MWOLD
                      optdata[TimeStamp][counter]['MWNEW'] = MWNEW
                      optdata[TimeStamp][counter]['PROLD'] = PROLD
                      optdata[TimeStamp][counter]['PRNEW'] = PRNEW
                      optdata[TimeStamp][counter]['ResourceName'] = row[2]
                      optdata[TimeStamp]['counter'] += 1

                      counter += 1
    endload = time.time()
    print('load ',endload-startload)
    return optdata

def facreate_matrices(optdata): #### For each TimeStamp, calculate P,q,G,h,A,b
        """ P is a (numvar + numconstraints + undergens + overgens) x (numvar + numconstraints + undergens + overgens) array with diagonal elements
              = price delta/mw delta for the corresponding segment for the first numvar rows/columns (all 0's other)
            q is a (numvar + numconstraints + undergens + overgens) x 1 array/vector =
               for the first numvar rows represent the low level price for each segment,
               the next numconstraint rows represent the violation amounts for each of the modeled network constraints  - the value for each is the Max Shadow Price of the constraint
               the next undergens rows are the steps of the undergeneration power balance penalty curve (up)
               the last overgens rows are the overgeneration powerbalance penalty curve
            G is a (numvar + numconstraints + undergens + overgens) +  (numvar + undergens - 1 + overgens - 1) + numconstraints rows by numvar + numconstraints +  undergens + overgens columns array
              first numvar + numconstraints + undergens + overgens rows are a negative diagonal identity matrix representing low limits for each optimized segment and violation MWs and power balance segments
              second numvar + undergens - 1 + overgens - 1 rows are a positive diagonal identity matrix representing high limits (MW high - MW low for segments, step size for first 8 undergeneration power balance steps) plus a rectangular 0 matrix on right side
              last numconstraints rows are shift factors (segment to constraint) - one row for each modeled constraint
            h is a (numvar + numconstraints + undergens + overgens) +  (numvar + undergens - 1 + overgens - 1) + numconstraints by one column matrix
              first (numvar + numconstraints + undergens + overgens) rows are 0
              second numvar rows are high limits (MW high - MW low) for each segment
              next undergens - 1 + overgens - 1 rows are undergeneration power balance step sizes
              last numconstraints are the mathlimit - LDL flows for each constraint
            A is a 1 row by numvar + numconstraints + undergens + overgens columns matrix.
              First numvar columns are 1.0,
              next numconstraints are 0.0,
              next undergens are 1.0,
              overgens are -1.0
            b is a 1 by 1 matrix with a single entry of GTBD (sum up all the basepoints)"""
        startmatrix = time.time()
        undergens = len(optdata['undergen'])
        overgens = len(optdata['overgen'])
        for TimeStamp in optdata['TimeStamps']:
            numvar = optdata[TimeStamp]['counter']
            try:
                numconstraints = len(optdata[TimeStamp]['networkconstraints'])
            except KeyError:
                numconstraints = 0
            optdata[TimeStamp]['GLists'] = []
            optdata[TimeStamp]['PLists'] = []


            optdata[TimeStamp]['b'] = optdata[TimeStamp]['GTBD'] - optdata[TimeStamp]['LDLsum']

            for counter in range(numvar): #Create P matrix for segment variable columns
                segment = optdata[TimeStamp][counter]
                plist = [0.0] * (numvar + numconstraints + undergens + overgens)
                try:
                    plist[counter] = (segment['PRNEW'] - segment['PROLD']) / (segment['MWNEW'] - segment['MWOLD'])
                except ZeroDivisionError:
                    print(TimeStamp, segment['ResourceName'], segment['MWOLD'], segment['MWNEW'])
                optdata[TimeStamp]['PLists'].append(plist)
            for counter in range(numconstraints + undergens + overgens): #Create P matrix for network constraint, undergen and overgen columns
                plist = [0.0] * (numvar + numconstraints + undergens + overgens)
                plist[numvar + counter] = 0.00001 #kluge-y way of fixing P matrix so that it has appropriate rank - look into making value smaller or going to more generic second order cone program
                optdata[TimeStamp]['PLists'].append(plist)

            optdata[TimeStamp]['qList'] = (numvar + numconstraints + undergens + overgens) * [0.0] #Initialize q vector
            for counter in range(numvar): #linear cost for segements
                segment = optdata[TimeStamp][counter]
                optdata[TimeStamp]['qList'][counter] = segment['PROLD']
            for constraint in range(numconstraints): #linear costs for network constraint violations (max shadow prices)
                optdata[TimeStamp]['qList'][numvar + constraint] = optdata[TimeStamp]['networkconstraints'][constraint][1]
            for undergen in range(undergens): #undergen penalty factors
                optdata[TimeStamp]['qList'][numvar + numconstraints + undergen] = optdata['undergen'][undergen][1]
            for overgen in range(overgens): #overgen penalty factors
                optdata[TimeStamp]['qList'][numvar + numconstraints + undergens + overgen] = optdata['overgen'][overgen][1]

            for counter in range(numvar): #Create G matrix for segment variables
                segment = optdata[TimeStamp][counter]
                glist = [0.0] * ((numvar + numconstraints + undergens + overgens) +  (numvar + undergens - 1 + overgens - 1) + numconstraints)
                glist[counter] = -1.0
                glist[counter + numvar + numconstraints + undergens + overgens] = 1.0
                for constraint in range(numconstraints):
                    constraintid = optdata[TimeStamp]['networkconstraints'][constraint][0]
                    resource = segment['ResourceName']
                    try:
                        glist[constraint + numvar + numconstraints + undergens + overgens + numvar + undergens - 1 + overgens - 1 ] = optdata[TimeStamp]['ShiftFactors'][constraintid][resource]
                        #print(resource,optdata[TimeStamp]['ShiftFactors'][constraintid][resource])
                    except KeyError:
                        optdata['log'].append([TimeStamp,constraintid,resource])
                optdata[TimeStamp]['GLists'].append(glist)
            for constraint in range(numconstraints): #Create G Matrix for network limit violation variables
                glist = [0.0] * ((numvar + numconstraints + undergens + overgens) +  (numvar + undergens - 1 + overgens -1) + numconstraints)
                glist[numvar + constraint] = -1.0
                glist[constraint + numvar + numconstraints + undergens + overgens + numvar + undergens - 1 + overgens - 1] = -1.0
                optdata[TimeStamp]['GLists'].append(glist)
            for undergen in range(undergens): #Create G matrix for undergen variable
                glist = [0.0] * ((numvar + numconstraints + undergens + overgens) +  (numvar + undergens - 1 + overgens -1) + numconstraints)
                glist[numvar + numconstraints + undergen] = -1.0 #undergen >= 0
                if undergen <= undergens - 2:
                    glist[numvar + numconstraints + undergens + overgens + undergen] = 1.0 #first 8 undergen segments <= max
                optdata[TimeStamp]['GLists'].append(glist)
            for overgen in range(overgens): #create G matrix for overgen
                glist = [0.0] * ((numvar + numconstraints + undergens + overgens) +  (numvar + undergens - 1 + overgens -1) + numconstraints)
                glist[numvar + numconstraints + undergens + overgen] = -1.0 #overgen >= 0
                if overgen <= overgens - 2:
                    glist[numvar + numconstraints + undergens + overgens + undergens - 1 + overgen] = 1.0 #first 8 undergen segments <= max
                optdata[TimeStamp]['GLists'].append(glist)

            optdata[TimeStamp]['hList'] = [0.0] * (numvar + numconstraints + undergens + overgens) #set inequality limits for >=0 stuff - segment variables, network limit violations, under/overgen
            for counter in range(numvar):
                segment = optdata[TimeStamp][counter]
                optdata[TimeStamp]['hList'].append(segment['MWNEW'] - segment['MWOLD'])
            for undergen in range(undergens - 1):
                optdata[TimeStamp]['hList'].append(optdata['undergen'][undergen][0])
            for overgen in range(overgens - 1):
                optdata[TimeStamp]['hList'].append(optdata['overgen'][overgen][0])
            for constraint in range(numconstraints):
                optdata[TimeStamp]['hList'].append(optdata[TimeStamp]['networkconstraints'][constraint][2])
                for resource in optdata[TimeStamp]['Resources']:
                    try:
                        constraintid = optdata[TimeStamp]['networkconstraints'][constraint][0]
                        sf = optdata[TimeStamp]['ShiftFactors'][constraintid][resource]
                        optdata[TimeStamp]['hList'][numvar + numconstraints + undergens + overgens + numvar + undergens - 1 + overgens - 1 + constraint] -= sf * optdata[TimeStamp][resource]
                    except KeyError:
                        pass

            optdata[TimeStamp]['AList'] = numvar * [1.0]
            for constraint in range(numconstraints):
                optdata[TimeStamp]['AList'].append(0.0)
            for undergen in range(undergens):
                optdata[TimeStamp]['AList'].append(1.0)
            for overgen in range(overgens):
                optdata[TimeStamp]['AList'].append(-1.0)

            optdata[TimeStamp]['P'] =   matrix(optdata[TimeStamp]['PLists'])
            optdata[TimeStamp]['q'] =   matrix(optdata[TimeStamp]['qList'])
            optdata[TimeStamp]['G'] =   matrix(optdata[TimeStamp]['GLists'])
            optdata[TimeStamp]['h'] =   matrix(optdata[TimeStamp]['hList'])
            optdata[TimeStamp]['A'] =   matrix(optdata[TimeStamp]['AList'],(1,numvar + numconstraints + undergens + overgens))
            optdata[TimeStamp]['b'] =   matrix(optdata[TimeStamp]['b'])
        endmatrix = time.time()
        print('matrix ',endmatrix-startmatrix)
        return optdata

def fasolve(optdata, timestamp):
    print(timestamp)
    startoptimize= time.time()
    test = optdata[timestamp]
    P = test['P']
    #print('P',P.size)
    q = test['q']
    #print('q',q.size)
    G = test['G']
    #print('G',G.size)
    h = test['h']
    #print('h',h.size)
    A = test['A']
    #print('A',A.size)
    b = test['b']
    #print('b',b.size)
    #print(h)
    optdata[timestamp]['Solution']=solvers.qp(P, q, G, h, A, b)
    endoptimize = time.time()
    print('optimize ',endoptimize - startoptimize)
    return optdata

def faconsolidate_results(optdata,timestamp):
    sol = optdata[timestamp]['Solution']
    numvar = optdata[timestamp]['counter']
    try:
        numconstraints = len(optdata[timestamp]['networkconstraints'])
    except KeyError:
        numconstraints = 0
    undergens = len(optdata['undergen'])
    overgens = len(optdata['overgen'])
    #print(sol['status'])
    #print(sol['x']) #solution
    #print(-1*sol['y'][0]) #Shadow Prices for equality constraints
    #print(sol['z'])
    try:
        constraintcounter = 0
        for constraintid in optdata[timestamp]['ShiftFactors']:
            optdata['shadowprices'].append([timestamp,constraintid,sol['z'][numvar + numconstraints + undergens + overgens + numvar + undergens - 1 + overgens - 1 + constraintcounter]])
            constraintcounter += 1
    except KeyError:
        pass
    for resource in optdata[timestamp]['Resources']:
        optdata['prices'][(timestamp,resource)] = -1*sol['y'][0]
        constraintcounter = 0
        try:
            for constraintid in optdata[timestamp]['ShiftFactors']:
                optdata['prices'][(timestamp,resource)] -= optdata[timestamp]['ShiftFactors'][constraintid][resource] * \
                   sol['z'][numvar + numconstraints + undergens + overgens + numvar + undergens - 1 + overgens - 1 + constraintcounter]
                constraintcounter += 1
        except KeyError:
            pass
    #print(sol['z']) #Shadow prices for inequality constraints
    test = optdata[timestamp]
    print(test['counter'])
    klugecost = {}
    for counter in range(numvar):
        resource = optdata[timestamp][counter]['ResourceName']
        x = sol['x'][counter]
        P = optdata[timestamp]['P'][counter,counter]
        q = optdata[timestamp]['q'][counter]
        optdata[timestamp][optdata[timestamp][counter]['ResourceName']] += x
        try:
            klugecost[(timestamp,resource)] += 0.5 * P * x * x + q * x
        except KeyError:
            klugecost[(timestamp,resource)] = 0.5 * P * x * x + q * x

    for row in optdata['60daygrd']:
        grdtimestamp = (row[0],row[1])
        grdresource = row[2]
        if grdtimestamp == timestamp:
            try:
                row.append(klugecost[(grdtimestamp,grdresource)])
                row.append(optdata[grdtimestamp][grdresource])
                row.append(optdata['prices'][(grdtimestamp,grdresource)])
            except KeyError:
                row.append(0)
                row.append(0)
                row.append(0)

    return optdata

def faget_binding_constraint_data(filename,optdata):
    print(filename)
    with open(filename, 'r') as csvfile: #open the binding constraint report
        filereader = csv.reader(csvfile)
        rowcounter = 0
        for row in filereader: #get a row of data
            rowcounter +=1
            #print(row)
            #print((row[0],row[1]),(row[0],row[1]) in optdata['TimeStamps'])
            if (row[0],row[1]) in optdata['TimeStamps'] :
                timestamp = (row[0],row[1])
                constraintid = row[2]
                maxsp = float(row[6])
                limit = float(row[7])
                try:
                    optdata[timestamp]['networkconstraints'].append([constraintid,maxsp,limit]) #list of binding constraints per interval
                except KeyError:
                    optdata[timestamp]['networkconstraints'] = [[constraintid,maxsp,limit]]
                try:
                    optdata['networkconstraints'].append([timestamp,constraintid,maxsp,limit]) #masterlist of binding constraints
                except KeyError:
                    optdata['networkconstraints'] = [[timestamp,constraintid,maxsp,limit]]
    #print(rowcounter)
    return optdata

def faset_up_sf_dicts(optdata):
    for constraint in optdata['networkconstraints']:
        timestamp = constraint[0]
        constraintid = constraint[1]
        try:
            optdata[timestamp]['ShiftFactors'][constraintid] = {}
            #print(timestamp,constraintid)
        except KeyError:
            optdata[timestamp]['ShiftFactors'] = {}
            optdata[timestamp]['ShiftFactors'][constraintid] = {}
            #print(timestamp,constraintid)
    return optdata

def faget_powerbalance(optdata):
    optdata['undergen'] = [(5.0,250.0),(5.0,300.0),(10.0,400.0),(10.0,500.0),(10.0,1000.0),(10.0,2250.0),(50.0,4500.0),(50.0,6000.0),(50.0,7500.0),(100000.0,9001.0)]
    optdata['overgen'] = [(100000.0,250.0)]
    return optdata

def faget_shiftfactor_data(filename,optdata):
    with open(filename, 'r') as csvfile: #open the binding constraint report
        filereader = csv.reader(csvfile)
        rowcounter = 0
        for row in filereader:
            try:
                timestamp = (row[0],row[1])
                constraintid = row[2]
                resourcename = row[5]
                shiftfactor = row[7]
                try:
                    optdata[timestamp]['ShiftFactors'][constraintid][resourcename] = float(shiftfactor)
                except KeyError:
                    pass
            except ValueError:
                    print(shiftfactor)
            except IndexError:
                pass
    return optdata
############################################
############################################

def save_file(data, filename):
    with open(filename, 'w',newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        for datum in data:
            spamwriter.writerow([datum[0],datum[1],datum[2]])
def save_h_file(data,platform):
    filename = 'C:/Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/h_trouble.csv'
    if platform == 'linux': filename = r'/home/reedy' + filename[14:]
    with open(filename, 'w',newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        for datum in data:
            spamwriter.writerow([datum])
def save_resource_file(resourcedict, optdata,timestamp, filename = r'C:/Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/resource_trouble.csv'):
    if optdata['platform'] == 'linux': filename = r'/home/reedy' + filename[14:]
    with open(filename, 'w',newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        for resource in optdata[timestamp]['Resources']:
            spamwriter.writerow([resource,resourcedict[resource]['LSL'],resourcedict[resource]['LDL'],resourcedict[resource]['output'],\
            resourcedict[resource]['HDL'],resourcedict[resource]['HSL'],resourcedict[resource]['status']])

def save_grd_file(data, filename):
    with open(filename, 'w',newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        for datum in data:
            spamwriter.writerow(datum)

def get_resource_data(filename,optdata):
    startload = time.time()
    badstatuslist = ['OUT','EMR','OFF','OFFNS','Telemetered Resource Status'] #shouldn't be included in optimization
    optdata['TimeStamps'] = set() #initialize set of all timestamps for the dictionary
    optdata['dispatch'] = [] #initialize the dispatch list
    optdata['60daygrd'] = [] #initialize output file
    PR = range(75,145,2) # iterable of the 35 different offer segments
    MW = range(74,144,2) # iterable of the 35 different offer segments
    with open(filename, 'r') as csvfile: #open the 60 day SCED report
        offerfilereader = csv.reader(csvfile)
        rowcounter = 0
        for row in offerfilereader: #get a row of data
            rowcounter += 1
            if rowcounter % 10000 == 0:
                print (rowcounter) #just so you have something to watch
            try:
                status = row[151] #look for last blank row
            except IndexError:
                print(rowcounter) #troubleshooting; shouldn't happen  Later stuff will error out
            #optdata['60daygrd'][((row[0],row[1]),row[2])] = row #build output file
            optdata['60daygrd'].append(row) #build output file
            if status not in badstatuslist: #only look at stuff that should be optimized
                TimeStamp = (row[0],row[1])
                Resource = row[2]
                optdata['TimeStamps'].add(TimeStamp)
############ Section that gets HDL, LDL, BP data for the row
                if status == 'ONTEST':
                    HSL = float(row[145])
                    LSL = float(row[148])
                    LDL = float(row[153])
                    HASL = float(row[146])
                    LASL = float(row[149])
                    HDL = float(row[147])
                    output = float(row[153])
                else:
                    #HSL = float(row[145]) - optdata['Resource Parameters'][Resource]['NFRC'] #Round 2 run without NFRC\
                    HSL = float(row[145])
                    LSL = min(float(row[148]),HSL)
                    LDL = min(float(row[150]),HSL)
                    HASL = min(float(row[146]),HSL)
                    LASL = min(float(row[149]),HSL)
                    HDL = min(float(row[147]),HSL)
                    output = float(row[153])

                if optdata['round'] == 3:
                    try:
                        RRSAward = float(row[190])
                        RegAward = float(row[191])
                        try:
                            RegDown = float(row[155])
                        except ValueError:
                            RegDown = 0.0
                    except IndexError:
                        RRSAward = 0.0
                        RegAward = 0.0
                        RegDown = 0.0


############## Correct HDL and LDL for situations that require correction
                if HDL < 0.05: HDL = 0 #accounting for observed SCED behavior
                if (HDL == HASL or status == 'ONREG'): #adjust HDL if HASL limited or ramp sharing
                    HDL = min(max(output + optdata['Resource Parameters'][Resource]['Ramp Rate Up'], HDL), HSL)
                #if ((LDL == LASL and status != 'ONTEST') or status == 'ONREG'): #adjust LDL if LASL limited or ramp sharing
                    #LDL = max(min(output + optdata['Resource Parameters'][Resource]['Ramp Rate Down'], LDL), LSL)
#not adjusting LDL so that can compare costs from one run to the next

                try:
                    optdata[TimeStamp]['LDLsum'] += LDL
                    optdata[TimeStamp]['LDLResources'].add(Resource)
                except KeyError: #If KeyError, optdata['TimeStamp']['LDLsum'] hasn't been defined yet, which means that no grd data has been entered yet; need to initialize
                    counter = 0
                    optdata[TimeStamp]['counter'] = 0
                    optdata[TimeStamp]['Resources'] = set()
                    optdata[TimeStamp]['LDLResources'] = set()
                    optdata[TimeStamp]['LDLResources'].add(Resource)
                    optdata[TimeStamp]['LDLsum'] = LDL
                optdata[TimeStamp][Resource] = {}
                optdata[TimeStamp][Resource]['LDL'] = LDL
                optdata[TimeStamp][Resource]['HDL'] = HDL
                optdata[TimeStamp][Resource]['LSL'] = LSL
                optdata[TimeStamp][Resource]['HSL'] = HSL
                optdata[TimeStamp][Resource]['output'] = output
                optdata[TimeStamp][Resource]['status'] = status
                if optdata['round'] == 3:
                    optdata[TimeStamp][Resource]['Previous Award'] = {}
                    optdata[TimeStamp][Resource]['Previous Award']['RRS'] = RRSAward
                    optdata[TimeStamp][Resource]['Previous Award']['Reg Up'] = RegAward
                    optdata[TimeStamp][Resource]['Previous Award']['Reg Down'] = RegDown

                for point in range (1,35,1): #for each segment, grab PROLD, PRNEW, MWOLD, MWNEW (defining points of segment)
                    MWNEW = float(row[MW[point]])
                    MWOLD = float(row[MW[point - 1]])
                    PRNEW = float(row[PR[point]])
                    PROLD = float(row[PR[point - 1]])
                    if (MWNEW > MWOLD and MWNEW > LDL and MWOLD < HSL) : #check to see if this is a segment that is within HSL/LDL and should be considered/optimized
                      optdata[TimeStamp][counter] = {} #initialize segment dictionary
                      if MWOLD < LDL: #check to see if need to adjust MWOLD and PROLD for LDL
                          PROLD = PROLD + (LDL - MWOLD) * (PRNEW - PROLD) / (MWNEW - MWOLD) #linearly interpolated price
                          MWOLD = LDL
                      if MWNEW > HSL: #check to see if need to adjust MWNEW and PRNEW for HDL
                          PRNEW = PRNEW + (HSL - MWNEW) * (PRNEW - PROLD) / (MWNEW - MWOLD) #linearly interpolated price
                          MWNEW = HSL # I think I can take this out and the constraints will handle this, but leaving language in in case
                      if (MWNEW > MWOLD and MWNEW > LDL and MWOLD < HSL) : #check to see if this is a segment that is within HSL/LDL and should be considered/optimized
#### Load up segment dictionary with data
                          optdata[TimeStamp][counter]['MWOLD'] = MWOLD
                          optdata[TimeStamp][counter]['MWNEW'] = MWNEW
                          optdata[TimeStamp][counter]['PROLD'] = PROLD
                          optdata[TimeStamp][counter]['PRNEW'] = PRNEW
                          optdata[TimeStamp][counter]['ResourceName'] = row[2]
                          optdata[TimeStamp]['Resources'].add(row[2])
                          optdata[TimeStamp]['counter'] += 1
                          counter += 1
    endload = time.time()
    optdata[TimeStamp]['Resources'] = list(optdata[TimeStamp]['Resources'])
    optdata[TimeStamp]['LDLResources'] = list(optdata[TimeStamp]['LDLResources'])
    print('load ',endload-startload)
    return optdata

def make_P_matrix(optdata,TimeStamp):
    """ P is a (numvar + 1 + numvar + 1 + numvar + 1 + numvar + 10 + numconstraints) ^ 2 with diagonal elements
          = price delta/mw delta for the corresponding segment for the first numvar rows/columns (all 0.000001's other)"""
    numvar = optdata[TimeStamp]['counter']
    try:
        numconstraints = len(optdata[TimeStamp]['networkconstraints'])
    except KeyError:
        numconstraints = 0
    optdata[TimeStamp]['PLists'] = []
    print(numvar,numconstraints)
    for counter in range(numvar): #Create P matrix for segment variable columns
        segment = optdata[TimeStamp][counter]
        plist = [0.0] * (numvar + 1 + numvar + 1 + numvar + 1 + numvar + 10 + numconstraints)
        try:
            plist[counter] = (segment['PRNEW'] - segment['PROLD']) / (segment['MWNEW'] - segment['MWOLD'])
        except ZeroDivisionError: #troubleshooting should never happen
            print(TimeStamp, segment['ResourceName'], segment['MWOLD'], segment['MWNEW'])
        optdata[TimeStamp]['PLists'].append(plist)
    for counter in range(1 + numvar + 1 + numvar + 1 + numvar + 10 + numconstraints): #Create P matrix for energy snuff, ns, ns ins, rrs, rrs ins, reg, reg ins (10) and constr viol columns
        plist = [0.0] * (numvar + 1 + numvar + 1 + numvar + 1 + numvar + 10 + numconstraints)
        plist[numvar + counter] = 0.0001 #kluge-y way of fixing P matrix so that it has appropriate rank - look into making value smaller or going to more generic second order cone program
        optdata[TimeStamp]['PLists'].append(plist)
    optdata[TimeStamp]['P'] =   matrix(optdata[TimeStamp]['PLists'])
    optdata[TimeStamp]['PLists'] = []
    return optdata

def make_q_vector(optdata,TimeStamp):
    """q is a (numvar + 1 + numvar + 1 + numvar + 1 + numvar + 10 + numconstraints) x 1 array/vector =
       for the first numvar rows represent the low level price for each segment,
       the next row is the undergeneration penalty
       the next numvar rows are nonspin costs (0)
       the next row is the NonSpin insufficiency penalty ($75)
       the next numvar rows are the rrs costs (0)
       the next row is the rrs insufficiency penalty ($1500)
       the next numvar rows are the reg up costs ($0)
       the next 10 rows are the reg up insufficiency costs
       the next numconstraint rows represent the violation amounts for each of the modeled network constraints  - the value for each is the Max Shadow Price of the constraint
       """
    numvar = optdata[TimeStamp]['counter']
    try:
        numconstraints = len(optdata[TimeStamp]['networkconstraints'])
    except KeyError:
        numconstraints = 0
    optdata[TimeStamp]['qList'] = (numvar + 1 + numvar + 1 + numvar + 1 + numvar + 10 + numconstraints) * [0.0] #Initialize q vector
    #print(len(optdata[TimeStamp]['qList']), 'just initialized')
    for counter in range(numvar): #linear cost for segements
        segment = optdata[TimeStamp][counter]
        optdata[TimeStamp]['qList'][counter] = segment['PROLD']
    #print(len(optdata[TimeStamp]['qList']), 'segments entered')
    optdata[TimeStamp]['qList'][numvar] = 9001 #undergen penalty factors
    optdata[TimeStamp]['qList'][numvar + 1 + numvar] = 75 #nonspin insufficiency penalty factor
    optdata[TimeStamp]['qList'][numvar + 1 + numvar + 1 + numvar] = 1500 #rrs insufficiency penalty factor
    #print(len(optdata[TimeStamp]['qList']), 'undergen, ns and rrs penalty factors entered')
    optdata[TimeStamp]['qList'][(numvar + 1 + numvar + 1 + numvar + 1 + numvar):(numvar + 1 + numvar + 1 + numvar + 1 + numvar + 10)] =\
               [250.0,300.0,400.0,500.0,1000.0,2250.0,4500.0,6000.0,7500.0,9001.0]
    #print(len(optdata[TimeStamp]['qList']), 'reg insuff entered')
    try:
        constraint = 0
        for constraintid in optdata[TimeStamp]['networkconstraints']:
            print('constraint',constraintid,optdata[TimeStamp]['networkconstraints'][constraintid])
            optdata[TimeStamp]['qList'][numvar + 1 + numvar + 1 + numvar + 1 + numvar + 10  + constraint] \
            = optdata[TimeStamp]['networkconstraints'][constraintid][0]
            constraint += 1
    except KeyError:
        pass
    #print(len(optdata[TimeStamp]['qList']), 'max shadow price entered')
    optdata[TimeStamp]['q'] =   matrix(optdata[TimeStamp]['qList'])
    optdata[TimeStamp]['qList'] = []
    return optdata

def make_G_matrix(optdata,TimeStamp):
    """
    G is a (7*numvar + 8*numresources + numconstraints) rows by (numvar + 1 + numvar + 1 + numvar + 1 + numvar + 10 + numvar + 1 + numconstraints) columns matrix
    first numvar rows are bp + regup + rrs + nsrs <= seg mw
    the next numvar rows are bp - reg down >= 0
    the next numvar rows are bp >= 0
    """
    Glists = [] #initialize list of lists
    numvar = optdata[TimeStamp]['counter']
    print(numvar)
    try:
        numconstraints = len(optdata[TimeStamp]['networkconstraints'])
    except KeyError:
        numconstraints = 0
    numresources  = len(optdata[TimeStamp]['Resources'])
    for column in range(numvar): #handle columns for each segment energy awards
        resource = optdata[TimeStamp][column]['ResourceName']
        #print(column, resource)
        glist = (5*numvar + 4*numresources + 2 * numconstraints + 22) * [0.0]
        glist[column] = 1.0 #bp + regup + rrs + nsrs <= seg mw
        glist[numvar + column] = -1.0 #bp  >= 0
        resourcenum = 0
        for rowresource in optdata[TimeStamp]['Resources']:
            if rowresource == resource:
                glist[5 * numvar + resourcenum] = 1.0 # bp + regup + rrs + ns <= eff HSL
                glist[5 * numvar + 3 * numresources + resourcenum] = 1.0 # bp + regup <= eff HDL
            resourcenum += 1
        constraintnum = 0
        try:
            for constraintid in optdata[TimeStamp]['networkconstraints']:
                try:
                    glist[5 * numvar + 4 * numresources + constraintnum] = optdata[TimeStamp]['ShiftFactors'][constraintid][resource]
                    #print(resource,optdata[TimeStamp]['ShiftFactors'][constraintid][resource])
                except KeyError:
                    optdata['log'].append([TimeStamp,constraintid,resource])
                constraintnum += 1
        except KeyError:
            pass
        Glists.append(glist)

    glist = (5*numvar + 4*numresources + 2 * numconstraints + 22) * [0.0] #energy insufficiency column
    glist[5 * numvar + 4 * numresources + 2 * numconstraints] = -1.0
    Glists.append(glist)

    for column in range(numvar): #handle columns for each segment nonspin awards
        resource = optdata[TimeStamp][column]['ResourceName']
        #print(column, resource)
        glist = (5*numvar + 4*numresources + 2 * numconstraints + 22) * [0.0]
        glist[column] = 1.0 #bp + regup + rrs + nsrs <= seg mw
        glist[2 * numvar + column] = -1.0 #ns  >= 0
        resourcenum = 0
        for rowresource in optdata[TimeStamp]['Resources']:
            if rowresource == resource:
                glist[5 * numvar + resourcenum] = 1.0 # bp + regup + rrs + ns <= eff HSL
            resourcenum += 1
        Glists.append(glist)

    glist = (5*numvar + 4*numresources + 2 * numconstraints + 22) * [0.0] #ns insufficiency column
    glist[5 * numvar + 4 * numresources + 2 * numconstraints + 1] = -1.0
    Glists.append(glist)


    for column in range(numvar): #handle columns for each segment rrs awards
        resource = optdata[TimeStamp][column]['ResourceName']
        #print(column, resource)
        glist = (5*numvar + 4*numresources + 2 * numconstraints + 22) * [0.0]
        glist[column] = 1.0 #bp + regup + rrs + nsrs <= seg mw
        glist[3 * numvar + column] = -1.0 #rrs  >= 0
        resourcenum = 0
        for rowresource in optdata[TimeStamp]['Resources']:
            if rowresource == resource:
                glist[5 * numvar + resourcenum] = 1.0 # bp + regup + rrs + ns <= eff HSL
                glist[5 * numvar + numresources + resourcenum] = 1.0 # rrs  <= 20% HSL
            resourcenum += 1
        Glists.append(glist)

    glist = (5*numvar + 4*numresources + 2 * numconstraints + 22) * [0.0] #rrs insufficiency column
    glist[5 * numvar + 4 * numresources + 2 * numconstraints + 2] = -1.0
    Glists.append(glist)

    for column in range(numvar): #handle columns for each segment reg up awards
        resource = optdata[TimeStamp][column]['ResourceName']
        #print(column, resource)
        glist = (5*numvar + 4*numresources + 2 * numconstraints + 22) * [0.0]
        glist[column] = 1.0 #bp + regup + rrs + nsrs <= seg mw
        glist[4 * numvar + column] = -1.0 #reg  >= 0
        resourcenum = 0
        for rowresource in optdata[TimeStamp]['Resources']:
            if rowresource == resource:
                glist[5 * numvar + resourcenum] = 1.0 # bp + regup + rrs + ns <= eff HSL (hsl - nfrc - ldl)
                glist[5 * numvar + 2 * numresources + resourcenum] = 1.0 # reg  <= ramp
                glist[5 * numvar + 3 * numresources + resourcenum] = 1.0 # bp + reg  <= hdl - ldl
            resourcenum += 1
        Glists.append(glist)

    for column in range(10):#regup insufficiency columns greater than zero
        glist = (5*numvar + 4 *numresources + 2 * numconstraints + 22) * [0.0]
        glist[5 * numvar + 4 * numresources + 2 * numconstraints + 3 + column] = -1.0
        if column <= 8:
            glist[5 * numvar + 4 * numresources + 2 * numconstraints + 13 + column] = 1.0
        Glists.append(glist)

    for column in range(numconstraints):#shadow price cap columns
        glist = glist = (5*numvar + 4 *numresources + 2 * numconstraints + 22) * [0.0]
        glist[5 * numvar + 4 * numresources + column] = -1.0
        glist[5 * numvar + 4 * numresources + numconstraints + column] = -1.0
        Glists.append(glist)

    optdata[TimeStamp]['G'] =   matrix(Glists)
    return optdata

def make_h_vector(optdata,TimeStamp):
    optround = optdata['round']
    numvar = optdata[TimeStamp]['counter']
    try:
        numconstraints = len(optdata[TimeStamp]['networkconstraints'])
    except KeyError:
        numconstraints = 0
    numresources  = len(optdata[TimeStamp]['Resources'])
    hlist = (5*numvar + 4 *numresources + 2 * numconstraints + 22) * [0.0]
    for counter in range(numvar): #set h values for bp + nsrs + rrs + reg up <= seg length
        hlist[counter] = optdata[TimeStamp][counter]['MWNEW'] - optdata[TimeStamp][counter]['MWOLD']
    rescounter = 0
    for resource in optdata[TimeStamp]['Resources']:
        try:
            freqaward = (optdata[TimeStamp][resource]['Previous Award']['RRS'] >= 0.1 or \
                         optdata[TimeStamp][resource]['Previous Award']['Reg Up'] >= 0.1 or\
                         optdata[TimeStamp][resource]['Reg Down'])
        except KeyError:
            freqaward = False

        if (optround == 2 or freqaward):
            if optdata[TimeStamp][resource]['status'] == 'ONTEST': hlist[5 * numvar + rescounter] = 0.0
            else:
                hlist[5 * numvar + rescounter] = max(optdata[TimeStamp][resource]['HSL'] - optdata['Resource Parameters'][resource]['NFRC'] - optdata[TimeStamp][resource]['LDL'], 0)
        else:
            hlist[5 * numvar + rescounter] = max(optdata[TimeStamp][resource]['HSL'] - optdata[TimeStamp][resource]['LDL'], 0)

        if ( ( (optround == 2) and (optdata['Resource Parameters'][resource]['RRS Q']) ) or \
             ( (optround == 3) and (optdata['Resource Parameters'][resource]['RRS Q']) and (freqaward or optdata['Resource Parameters'][resource]['NFRC'] <= 0.5)) ):
            hlist[5 * numvar + numresources + rescounter] = 0.2 * (optdata[TimeStamp][resource]['HSL'] )

        if ( ( (optround == 2) and (optdata['Resource Parameters'][resource]['Reg Up Q']) ) or \
             ( (optround == 3) and (optdata['Resource Parameters'][resource]['Reg Up Q']) and (freqaward or optdata['Resource Parameters'][resource]['NFRC'] <= 0.5)) ):
            hlist[5 * numvar + 2 * numresources + rescounter] = \
              max(optdata[TimeStamp][resource]['HDL'] - optdata[TimeStamp][resource]['output'], optdata['Resource Parameters'][resource]['Ramp Rate Up'] )

        hlist[5 * numvar + 3 * numresources + rescounter] =  max(optdata[TimeStamp][resource]['HDL'] - optdata[TimeStamp][resource]['LDL'] , 0.0)

        rescounter += 1

    constraintcounter = 0
    try:
        for constraintid in optdata[TimeStamp]['networkconstraints']:
            hlist[5 * numvar + 4 * numresources + constraintcounter] = optdata[TimeStamp]['networkconstraints'][constraintid][1]
            for resource in optdata[TimeStamp]['LDLResources']:
                try:
                    sf = optdata[TimeStamp]['ShiftFactors'][constraintid][resource]
                    LDL = optdata[TimeStamp][resource]['LDL']
                    hlist[5 * numvar + 4 * numresources + constraintcounter] -= sf * LDL
                except KeyError:
                    pass
            constraintcounter += 1
    except KeyError:
        pass
    for reginsuff in [(5.0,0),(5.0,1),(10.0,2),(10.0,3),(10.0,4),(10.0,5),(50.0,6),(50.0,7),(50.0,8)]:
        hlist[5*numvar + 4 *numresources + 2 * numconstraints + 13 + reginsuff[1]] = reginsuff[0]
    optdata[TimeStamp]['h'] =   matrix(hlist)
    print(numresources,'numresources')
    save_h_file(hlist,optdata['platform'])
    return optdata

def make_A_matrix(optdata,TimeStamp):
    numvar = optdata[TimeStamp]['counter']
    try:
        numconstraints = len(optdata[TimeStamp]['networkconstraints'])
    except KeyError:
        numconstraints = 0
    Alists = []
    for column in range(numvar + 1): #handle energy awards and energy insufficiency
        Alist = [1.0,0.0,0.0,0.0]
        Alists.append(Alist)

    for column in range(numvar + 1): #handle nonspin awards and nonspin insufficiency
        Alist = [0.0,1.0,0.0,0.0]
        Alists.append(Alist)

    for column in range(numvar + 1): #handle rrs awards and insufficiency
        Alist = [0.0,0.0,1.0,0.0]
        Alists.append(Alist)

    for column in range(numvar + 10): #handle reg up awards and insufficiency
        Alist = [0.0,0.0,0.0,1.0]
        Alists.append(Alist)

    for column in range(numconstraints):
        Alist = [0.0,0.0,0.0,0.0]
        Alists.append(Alist)
    optdata[TimeStamp]['A'] =   matrix(Alists)
    return optdata

def make_b_vector(optdata,TimeStamp):
    syscon = optdata[TimeStamp]['System Conditions']
    GTBD = syscon['GTBD']
    LDLSUM = optdata[TimeStamp]['LDLsum']
    NSRS = syscon['NSRS MW']
    RRS = syscon['RRS MW']
    RegUp = syscon['Reg Up MW']
    RegDown = syscon['Reg Down MW']
    optdata[TimeStamp]['b'] = matrix([GTBD - LDLSUM, NSRS,RRS,RegUp])
    print(['GTBD',GTBD,'NSRS MW', NSRS, "RRS MW",RRS, 'Reg Up', RegUp,'LDL Sum',LDLSUM])
    return optdata

def solve(optdata, timestamp):
    print(timestamp)
    AggHSL = 0.0
    AggNFRC = 0.0
    for resource in optdata[timestamp]['Resources']:
        AggHSL += optdata[timestamp][resource]['HSL']
        AggNFRC += optdata['Resource Parameters'][resource]['NFRC']
    print('Aggregate HSL', AggHSL)
    print('Aggregate NFRC', AggNFRC)
    startoptimize= time.time()
    test = optdata[timestamp]
    P = test['P']
    print('P',P.size)
    q = test['q']
    print('q',q.size)
    G = test['G']
    print('G',G.size)
    h = test['h']
    print('h',h.size)
    A = test['A']
    print('A',A.size)
    b = test['b']
    print('b',b.size)

    diag_data = [(q,'q'),(P,'P'),(G,'G'),(h,'h'),(A,'A'),(b,'b')]

    #for (thematrix,thename) in diag_data:
        #filename = r'C:\Users\Admin\Dropbox\PythonCode\RT-Coopt\Matrix{0}_trouble.csv'.format(thename)
        #if optdata['platform'] == 'linux': filename = r'/home/reedy' + filename[14:]
        #print_matrix(thematrix,filename)



    optdata[timestamp]['Solution']=solvers.qp(P, q, G, h, A, b)
    endoptimize = time.time()
    print('optimize ',endoptimize - startoptimize)
    return optdata

def consolidate_results(optdata,timestamp):
    optdata[timestamp]['AS Prices'] = {} #initialize AS Price dict
    sol = optdata[timestamp]['Solution']
    numvar = optdata[timestamp]['counter']
    numresources  = len(optdata[timestamp]['Resources'])
    try:
        numconstraints = len(optdata[timestamp]['networkconstraints'])
    except KeyError:
        numconstraints = 0
    #print(sol['status'])
    #print(sol['x']) #solution
    #print(-1*sol['y'][0]) #Shadow Prices for equality constraints
    #print(sol['z'])
    try:
        constraintcounter = 0  #add shadow price to optdata[timestamp]['networkconstraints'][constraintid] list
        for constraintid in optdata[timestamp]['networkconstraints']:
            optdata[timestamp]['networkconstraints'][constraintid].append(sol['z'][5 * numvar + 4 * numresources + constraintcounter])
            constraintcounter += 1
    except KeyError:
        pass
    for resource in optdata[timestamp]['LDLResources']: #get prices for each resource that was dispatched
        optdata[timestamp][resource]['Energy Price'] = -1*sol['y'][0] #system lambda
        try:
            for constraintid in optdata[timestamp]['ShiftFactors']: #add in the congestion component
                optdata[timestamp][resource]['Energy Price'] -= optdata[timestamp]['ShiftFactors'][constraintid][resource] * \
                  optdata[timestamp]['networkconstraints'][constraintid][2]
        except KeyError:
            pass
    optdata[timestamp]['AS Prices']['NSRS'] = -1*sol['y'][1]  #get the AS prices
    optdata[timestamp]['AS Prices']['RRS'] = -1*sol['y'][2]
    optdata[timestamp]['AS Prices']['Reg Up'] = -1*sol['y'][3]

    test = optdata[timestamp]

    for counter in range(numvar):
        x = sol['x'][counter]
        P = optdata[timestamp]['P'][counter,counter]
        q = optdata[timestamp]['q'][counter]
        try:
            optdata[timestamp] [optdata[timestamp][counter]['ResourceName']]['Energy Award'] += x
            optdata[timestamp][optdata[timestamp][counter]['ResourceName']]['NSRS Award'] += sol['x'][numvar + 1 + counter]
            optdata[timestamp][optdata[timestamp][counter]['ResourceName']]['RRS Award'] += sol['x'][2 * numvar + 2 + counter]
            optdata[timestamp][optdata[timestamp][counter]['ResourceName']]['Reg Up Award'] += sol['x'][3 * numvar + 3 + counter]
            optdata[timestamp] [optdata[timestamp][counter]['ResourceName']]['Production Cost'] += 0.5 * P * x * x + q * x

        except KeyError:
            optdata[timestamp] [optdata[timestamp][counter]['ResourceName']]['Energy Award'] = x + optdata[timestamp] [optdata[timestamp][counter]['ResourceName']]['LDL']
            optdata[timestamp][optdata[timestamp][counter]['ResourceName']]['NSRS Award'] = sol['x'][numvar + 1 + counter]
            optdata[timestamp][optdata[timestamp][counter]['ResourceName']]['RRS Award'] = sol['x'][2 * numvar + 2 + counter]
            optdata[timestamp][optdata[timestamp][counter]['ResourceName']]['Reg Up Award'] = sol['x'][3 * numvar + 3 + counter]
            optdata[timestamp] [optdata[timestamp][counter]['ResourceName']]['Production Cost'] = 0.5 * P * x * x + q * x
    for row in optdata['60daygrd']:
        grdtimestamp = (row[0],row[1])
        grdresource = (row[2])
        if grdtimestamp == timestamp:
            try:
                    row.append(optdata[grdtimestamp][grdresource]['Production Cost'])
                    row.append(optdata[grdtimestamp][grdresource]['Energy Award'])
                    row.append(optdata[grdtimestamp][grdresource]['NSRS Award'])
                    row.append(optdata[grdtimestamp][grdresource]['RRS Award'])
                    row.append(optdata[grdtimestamp][grdresource]['Reg Up Award'])
                    row.append(optdata[grdtimestamp][grdresource]['Energy Price'])
                    row.append(optdata[grdtimestamp]['AS Prices']['NSRS'])
                    row.append(optdata[grdtimestamp]['AS Prices']['RRS'])
                    row.append(optdata[grdtimestamp]['AS Prices']['Reg Up'])
            except KeyError:
                    row.append(0)
                    row.append(0)
                    row.append(0)
                    row.append(0)
                    row.append(0)
                    row.append(0)
                    row.append(0)
                    row.append(0)
                    row.append(0)

    return optdata

def get_binding_constraint_data(filename,optdata):
    print(filename)
    with open(filename, 'r') as csvfile: #open the binding constraint report
        filereader = csv.reader(csvfile)
        rowcounter = 0
        for row in filereader: #get a row of data
            rowcounter +=1
            #print(row)
            #print((row[0],row[1]),(row[0],row[1]) in optdata['TimeStamps'])
            if (row[0],row[1]) in optdata['TimeStamps'] :
                timestamp = (row[0],row[1])
                constraintid = row[2]
                maxsp = float(row[6])
                limit = float(row[7])
                try:
                    optdata[timestamp]['networkconstraints'][constraintid] = [maxsp,limit] #list of binding constraints per interval
                except KeyError:
                    optdata[timestamp]['networkconstraints'] = {}
                    optdata[timestamp]['networkconstraints'][constraintid] = [maxsp,limit]
                try:
                    optdata['networkconstraints'].append([timestamp,constraintid,maxsp,limit]) #masterlist of binding constraints
                except KeyError:
                    optdata['networkconstraints'] = [[timestamp,constraintid,maxsp,limit]]
    #print(rowcounter)
    return optdata

def set_up_sf_dicts(optdata):
    for constraint in optdata['networkconstraints']:
        timestamp = constraint[0]
        constraintid = constraint[1]
        try:
            optdata[timestamp]['ShiftFactors'][constraintid] = {}
            #print(timestamp,constraintid)
        except KeyError:
            optdata[timestamp]['ShiftFactors'] = {}
            optdata[timestamp]['ShiftFactors'][constraintid] = {}
            #print(timestamp,constraintid)
    return optdata

def get_shiftfactor_data(filename,optdata):
    with open(filename, 'r') as csvfile: #open the binding constraint report
        filereader = csv.reader(csvfile)
        rowcounter = 0
        for row in filereader:
            try:
                timestamp = (row[0],row[1])
                constraintid = row[2]
                resourcename = row[5]
                shiftfactor = row[6]
                try:
                    optdata[timestamp]['ShiftFactors'][constraintid][resourcename] = float(shiftfactor)
                except KeyError:
                    pass
            except ValueError:
                    print(shiftfactor)
            except IndexError:
                pass
    return optdata

def get_system_conditions(optdata):
    """
    1)    System Conditions
          0- timestamp
          1- repeated hour flag
          2- Reg up AS obl for that interval
          3- Reg down AS obl for that interval
          4- RRS AS obl for that interval
          5- NSRS AS obl for that interval
          6- Reg up AS price for that interval
          7- Reg down AS price for that interval
          8- RRS AS price for that interval
          9- NSRS AS price for that interval
          10-system lambda for that interval
          11-GTBD for that interval

    """
    filename = 'C:/Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/system_conditions.csv'
    if optdata['platform'] == 'linux': filename = r'/home/reedy' + filename[14:]
    with open(filename, 'r') as csvfile: #open the 60 day SCED report
        filereader = csv.reader(csvfile)
        rowcounter = 0
        for row in filereader:
            if row[0] != 'sced_time_stamp':
                systemconditions = {}
                systemconditions['Reg Down MW'] = float(row[2])
                systemconditions['Reg Up MW'] = float(row[3])
                systemconditions['RRS MW'] = float(row[4])
                systemconditions['NSRS MW'] = float(row[5])
                systemconditions['Reg Down Price'] = float(row[6])
                systemconditions['Reg Up Price'] = float(row[7])
                systemconditions['RRS Price'] = float(row[8])
                systemconditions['NSRS Price'] = float(row[9])
                systemconditions['System Lambda'] = float(row[11])
                systemconditions['GTBD'] = float(row[10])
                optdata[(row[0],row[1])] = {}
                optdata[(row[0],row[1])]['System Conditions'] = systemconditions
    return optdata

def get_resource_parameters(optdata):
    """
    0- Resource
    1- Reg up AS qual status
    2- Reg down AS qual status
    3- RRS AS obl qual status
    4- NSRS AS obl qual status
    5- Ramp Rate Up
    6- Ramp Rate Down
    7- nfrc
    """
    optdata['Resource Parameters'] = {}
    rp = {}
    filename = 'C:/Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/resource_parameters.csv'
    if optdata['platform'] == 'linux': filename = r'/home/reedy' + filename[14:]
    with open(filename, 'r') as csvfile: #open the 60 day SCED report
        filereader = csv.reader(csvfile)
        rowcounter = 0
        for row in filereader:
            if row[0] != 'resource_name':
                resource = row[0]
                rp[resource] = {}
                rp[resource]['Reg Up Q'] = float(row[1])
                rp[resource]['Reg Down Q'] = float(row[2])
                rp[resource]['RRS Q'] = float(row[3])
                rp[resource]['NSRS Q'] = float(row[4])
                rp[resource]['Ramp Rate Up'] = float(row[5])
                rp[resource]['Ramp Rate Down'] = float(row[6])
                rp[resource]['NFRC'] = float(row[7])
    optdata['Resource Parameters'] = rp
    return optdata


def get_grdfilename(wdname,round):
    filename = glob.glob(wdname+'60d_SCED_Gen*.csv')[0]
    return filename

def get_bc_filenames(datefolder):
    filenamelist = glob.glob(datefolder+'*12302.*.csv')
    return filenamelist

def get_shift_factor_filenames(datefolder):
    filenamelist = glob.glob(datefolder+'*12354.*.csv')
    return filenamelist

def get_platform():
    if sys.platform == 'linux':
        platform = 'linux'
    else:
        platform = 'windows'
    return platform

def print_matrix(thematrix,filename):
    rows = thematrix.size[0]
    columns = thematrix.size[1]
    with open(filename, 'w',newline='') as csvfile:
        spamwriter = csv.writer(csvfile, delimiter=',')
        therow = []
        for row in range(rows):
            therow = []
            for column in range(columns):
                therow.append(thematrix[row,column])
            spamwriter.writerow(therow)





def main():

    loaddate = sys.argv[1]
    platform = get_platform()
    start = time.time()
    first_start = time.time()
    wdname = 'C:/Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/'+loaddate+'/'
    if platform == 'linux': wdname = r'/home/reedy' + wdname[14:] #rename for LrwUEb1InunX
    filename = glob.glob(wdname+'60d_SCED_Gen*.csv')[0]
    optdata = faget_resource_data(filename)
    optdata['shadowprices'] = []
    optdata['log'] = []
    optdata = faget_powerbalance(optdata)
    print('got data')
    startbc = time.time()
    for filename in glob.glob(wdname+'cdr.00012302.*.csv'):
        optdata = faget_binding_constraint_data(filename,optdata)
    #for timestamp in optdata['TimeStamps']:print(timestamp)
    optdata = faset_up_sf_dicts(optdata)
    for filename in glob.glob(wdname+'cdr.00012354.*.csv'):
        optdata = faget_shiftfactor_data(filename,optdata)
    endbc = time.time()
    print('binding ',endbc-startbc)
    optdata = facreate_matrices(optdata)
    print('made matrices')
    for timestamp in optdata['TimeStamps']:
        optdata = fasolve(optdata,timestamp)
        print('solved')
        optdata = faconsolidate_results(optdata,timestamp)
        print('consolidated')
    filename = r'C:\Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/{0}/round1grd.csv'.format(loaddate)
    fasave_grd_file(optdata['60daygrd'], wdname + 'round1result.csv')
    print('saved')
    end = time.time()
    print('full run: ',end - start)


############################################end round 1############################################################
















    optdata = {}
    #

    start = time.time()
    optdata['platform'] = get_platform()
    optdata['round'] = 2
    print(optdata['platform'])
    #load system conditions -store in optdata['System Conditions']
    optdata = get_system_conditions(optdata)
    #load resource parameters -store in optdata['Resource Parameters']
    optdata = get_resource_parameters(optdata)
    #begin daily iteration
    #for day in range(31):
        #load grd file -store in optdata['grd file']
    #datefolder = 'C:/Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/'+loaddate.strftime('%Y%m%d')+'/'
    datefolder = r'C:/Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/'+ loaddate + r'/'
    if optdata['platform'] == 'linux': datefolder = r'/home/reedy' + datefolder[14:]
    #grdfilename = get_grdfilename(datefolder)
    grdfilename = r'C:\Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/{0}/round1result.csv'.format(loaddate)
    if optdata['platform'] == 'linux': grdfilename = r'/home/reedy' + grdfilename[14:]
    optdata = get_resource_data(grdfilename,optdata)
    print('got resource data')
    #load binding constraints
    bindingconstraintfilenames = get_bc_filenames(datefolder)
    for bcfilename in bindingconstraintfilenames:
        optdata = get_binding_constraint_data(bcfilename,optdata)
    print('got binding constraint data')
    optdata = set_up_sf_dicts(optdata)
    #load shift factors
    shiftfactorfilenames = get_shift_factor_filenames(datefolder)
    for sffilename in shiftfactorfilenames:
        print(sffilename)
        optdata = get_shiftfactor_data(sffilename,optdata)
    print('got shift factors')

    for timestamp in optdata['TimeStamps']:
    #for garbage in range(1): #just a troubleshooting so that it only does one timestamp
        #timestamp = ("07/20/2016 01:45:28","N")
        #Make P matrix
        optdata = make_P_matrix(optdata,timestamp)
        #Make q vector
        optdata = make_q_vector(optdata,timestamp)
        #Make G matrix
        optdata = make_G_matrix(optdata,timestamp)
        #Make h vector
        optdata = make_h_vector(optdata,timestamp)
        #Make A matrix
        optdata = make_A_matrix(optdata,timestamp)
        #Make b vector
        optdata = make_b_vector(optdata,timestamp)
        print('made matrix for {0}'.format(timestamp[0]))
        #start troubleshooting to be deleted block
        tempdict ={}
        for resource in optdata[timestamp]['Resources']:
            tempdict[resource] = optdata[timestamp][resource]
        save_resource_file(tempdict,optdata,timestamp)

        #optimize
        optdata = solve(optdata,timestamp)
        print('optimized {0}'.format(timestamp[0]))


        #Consolidate Results
        optdata = consolidate_results(optdata,timestamp)
        print('consolidated {0}'.format(timestamp[0]))
    #Publish Results
    filename = r'C:\Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/{0}/newgrd_9.csv'.format(loaddate)
    if optdata['platform'] == 'linux': filename = r'/home/reedy' + filename[14:]
    print(filename)
    save_grd_file(optdata['60daygrd'], filename)
    endofitall=time.time()
    print('total time', endofitall - start)

#############################################################################################################################
#################### end round 2 ############################################################################################
#############################################################################################################################







    optdata = {}
    start = time.time()
    optdata['platform'] = get_platform()
    optdata['round'] = 3
    print(optdata['platform'])
    #load system conditions -store in optdata['System Conditions']
    optdata = get_system_conditions(optdata)
    #load resource parameters -store in optdata['Resource Parameters']
    optdata = get_resource_parameters(optdata)
    #begin daily iteration
    #for day in range(31):
        #load grd file -store in optdata['grd file']
    #datefolder = 'C:/Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/'+loaddate.strftime('%Y%m%d')+'/'
    datefolder = r'C:/Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/'+ loaddate + r'/'
    if optdata['platform'] == 'linux': datefolder = r'/home/reedy' + datefolder[14:]
    #grdfilename = get_grdfilename(datefolder)
    grdfilename = r'C:\Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/{0}/newgrd_9.csv'.format(loaddate)
    if optdata['platform'] == 'linux': grdfilename = r'/home/reedy' + grdfilename[14:]
    optdata = get_resource_data(grdfilename,optdata)
    print('got resource data')
    #load binding constraints
    bindingconstraintfilenames = get_bc_filenames(datefolder)
    for bcfilename in bindingconstraintfilenames:
        optdata = get_binding_constraint_data(bcfilename,optdata)
    print('got binding constraint data')
    optdata = set_up_sf_dicts(optdata)
    #load shift factors
    shiftfactorfilenames = get_shift_factor_filenames(datefolder)
    for sffilename in shiftfactorfilenames:
        print(sffilename)
        optdata = get_shiftfactor_data(sffilename,optdata)
    print('got shift factors')

    for timestamp in optdata['TimeStamps']:
    #for garbage in range(1): #just a troubleshooting so that it only does one timestamp
        #timestamp = ("07/20/2016 01:45:28","N")
        #Make P matrix
        optdata = make_P_matrix(optdata,timestamp)
        #Make q vector
        optdata = make_q_vector(optdata,timestamp)
        #Make G matrix
        optdata = make_G_matrix(optdata,timestamp)
        #Make h vector
        optdata = make_h_vector(optdata,timestamp)
        #Make A matrix
        optdata = make_A_matrix(optdata,timestamp)
        #Make b vector
        optdata = make_b_vector(optdata,timestamp)
        print('made matrix for {0}'.format(timestamp[0]))
        #start troubleshooting to be deleted block
        tempdict ={}
        for resource in optdata[timestamp]['Resources']:
            tempdict[resource] = optdata[timestamp][resource]
        save_resource_file(tempdict,optdata,timestamp)

        #optimize
        optdata = solve(optdata,timestamp)
        print('optimized {0}'.format(timestamp[0]))


        #Consolidate Results
        optdata = consolidate_results(optdata,timestamp)
        print('consolidated {0}'.format(timestamp[0]))
    #Publish Results
    filename = r'C:\Users/Admin/Dropbox/PythonCode/RT-Coopt/FinalData/{0}/newgrd_9.csv'.format(loaddate)

    if optdata['platform'] == 'linux': filename = r'/home/reedy' + filename[14:]
    print(filename)
    save_grd_file(optdata['60daygrd'], filename)
    endofitall=time.time()
    print('total time', (endofitall - first_start)/60.0)




if __name__ == "__main__":
    main()
