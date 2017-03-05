# -*- coding: utf-8 -*-
"""
Created on Sun Aug  7 16:29:51 2016

@author: Daniel
"""

"""import various dependencies"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as dt
import spacepy.time as spyt
import os
"""identifies cdf library path in system, might need to change this on a different system"""
os.environ["CDF_LIB"] = "/Applications/cdf/cdf36_3-dist/lib"
from spacepy import pycdf
import bisect
import datetime
import urllib

def main():
    
    """user interface to identify CDF file to download"""
    pay = raw_input('Enter Payload Here (e.g. 1A): ')
    date = raw_input('Enter Date (e.g. YYYYMMDD): ')
    tr = raw_input('Enter Time Range (e.g. HHMMstart-HHMMend): ')
    
    """creates relevant integers from user input"""
    hstart = int(tr[0:2])
    mstart = int(tr[2:4])
    hend = int(tr[5:7])
    mend = int(tr[7:9])
    y = int(date[0:4])
    m = int(date[4:6])
    d = int(date[6:9])
    
    """downloads ephm cdf file from barreldata server, reads into program"""
    ephmurl = 'http://barreldata.ucsc.edu/data_products/v05/l2/{}/{}/bar_{}_l2_ephm_{}_v05.cdf'.format(pay,date[2:8],pay,date)   
    ephmdata = urllib.urlretrieve(ephmurl)
    cdfephm = pycdf.CDF(ephmdata[0])
    
    print 'EPHM file downloaded.'
    
    """
    ephmpath = '/Users/Daniel/anaconda/workspace/BARREL/{}_{}/bar_{}_l2_ephm_{}_v05.cdf'.format(date[0:4],pay,pay,date)
    cdfephm = pycdf.CDF(ephmpath)
    """
    
    """divides time range into indexed increments"""
    start = datetime.datetime(y,m,d,hstart,mstart,0,0)
    stop = datetime.datetime(y,m,d,hend,mend,0,0)
    start_indephm = bisect.bisect_left(cdfephm['Epoch'], start)
    stop_indephm = bisect.bisect_left(cdfephm['Epoch'], stop)
    
    """reads in time and L values from ephm file"""
    tephm = cdfephm['Epoch'][start_indephm:stop_indephm]
    L_3A = cdfephm['L_Kp2'][start_indephm:stop_indephm]
    """add plot of L value vs MLT (magnetic local time)"""        
    
    """reads in altitude from ephm file"""
    alt = cdfephm['GPS_Alt'][start_indephm:stop_indephm]    
    
    """downloads fspc cdf"""
    fspcurl = 'http://barreldata.ucsc.edu/data_products/v05/l2/{}/{}/bar_{}_l2_fspc_{}_v05.cdf'.format(pay,date[2:8],pay,date)   
    fspcdata = urllib.urlretrieve(fspcurl)
    cdffspc = pycdf.CDF(fspcdata[0])    
    
    print 'FSPC file downloaded.'
    
    """
    fspcpath = '/Users/Daniel/anaconda/workspace/BARREL/{}_{}/bar_{}_l2_fspc_{}_v05.cdf'.format(date[0:4],pay,pay,date)
    cdffspc = pycdf.CDF(fspcpath)
    """
    
    """divides fspc time into indexed increments"""
    start_indfspc = bisect.bisect_left(cdffspc['Epoch'], start)
    stop_indfspc = bisect.bisect_left(cdffspc['Epoch'], stop)    
    
    """reads in time from fpsc file"""
    tfspc = cdffspc['Epoch'][start_indfspc:stop_indfspc]
    
    """look into fspc curve smoothing (smooth 20), maybe ask user if they want to smooth"""    
    
    """checks payload to determine light curve levels, reads in appropriate light curves"""
    if(pay[0] == '1'):    
        light1 = cdffspc['FSPC1'][start_indfspc:stop_indfspc]
    else:
        light1a = cdffspc['FSPC1a'][start_indfspc:stop_indfspc]
        light1b = cdffspc['FSPC1b'][start_indfspc:stop_indfspc] 
        light1c = cdffspc['FSPC1c'][start_indfspc:stop_indfspc]
            
    light2 = cdffspc['FSPC2'][start_indfspc:stop_indfspc]  
    light3 = cdffspc['FSPC3'][start_indfspc:stop_indfspc]  
    light4 = cdffspc['FSPC4'][start_indfspc:stop_indfspc]       
    
    """sets date format"""
    myFmt = dt.DateFormatter('%H:%M')
    
    """plots altitude"""
    plt.subplot(5,1,1)
    plt.plot(tephm, alt)
    plt.title('Payload {}, {}-{}-{}, {}'.format(pay,d,m,y,tr))
    plt.gca().xaxis.set_major_formatter(myFmt)    
    plt.gca().set_ylim([25,40])
    plt.gca().set_xlabel('UT')
    plt.gca().set_ylabel('Altitude')         
    
    """plots light curves"""
    plt.subplot(5,1,2)
    if(pay[0] == '1'):
        plt.plot(tfspc,light1,color = 'red',label = 'level 1')
        plt.legend(loc = 'upper right',frameon = False)        
    else:
        plt.plot(tfspc,light1a,color = 'red',label = '1a')
        plt.plot(tfspc,light1b,color = 'orange',label = '1b')
        plt.plot(tfspc,light1c,color = 'yellow',label = '1c')
        plt.legend(loc = 'upper right',frameon = False)
    plt.gca().xaxis.set_major_formatter(myFmt)
    """plt.yscale('log')"""
    
    plt.subplot(5,1,3)
    plt.plot(tfspc,light2,color = 'green',label = 'level 2')
    plt.legend(loc = 'upper right',frameon = False)
    plt.gca().xaxis.set_major_formatter(myFmt)
    """plt.yscale('log')"""
    
    plt.subplot(5,1,4)
    plt.plot(tfspc,light3,color = 'blue',label = 'level 3')
    plt.legend(loc = 'upper right',frameon = False)
    plt.gca().xaxis.set_major_formatter(myFmt)
    """plt.yscale('log')"""
    
    plt.subplot(5,1,5)
    plt.plot(tfspc,light4,color = 'purple',label = 'level 4')
    plt.legend(loc = 'upper right',frameon = False)
    plt.gca().xaxis.set_major_formatter(myFmt)
    """plt.yscale('log')"""
    plt.show()
    
    """add user interface for time range of event, avg. lat long, avg. alt"""
    
    cdfephm.close()
    cdffspc.close() 

if __name__ == "__main__":
    main()
