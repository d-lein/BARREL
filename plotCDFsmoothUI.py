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
from math import factorial


"""average function"""
def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')


def main():
    
    exitMain = 'n'
    
    """user interface to identify CDF file to download"""
    pay = raw_input('Enter Payload Here (e.g. 1A), must be valid payload: ')
    date = raw_input('Enter Date (e.g. YYYYMMDD), must be valid date for that payload: ')
    
    """iterative loop that allows for user to examine multiple dates without exiting program"""
    while exitMain != 'e':
        readInd = 'x'
        """user interface for file path read vs. online download"""
        readInd = raw_input('Would you like to read in a file from local path or online (p/o)? WARNING:\n'
                            'improper configuration of the file path (in lines 41 and 43 of the source)\n'
                            'will result in a fatal error: ')
        """error loop in case user does not enter appropriate character"""
        while ((readInd != 'o') | (readInd != 'p')):
            if ((readInd == 'o') | (readInd == 'p')):
                break
            readInd = raw_input('ERROR: must enter either "o" or "p".\n'
                                'Would you like to read in a file from local path or online (p/o): ')
            
        if readInd == 'p':
            """reads in files from user determined paths (configure paths in lines 41 and 43)"""
            ephmpath = '/Users/Daniel/anaconda/workspace/BARREL/{}_{}/bar_{}_l2_ephm_{}_v05.cdf'.format(date[0:4],pay,pay,date)
            cdfephm = pycdf.CDF(ephmpath)
            fspcpath = '/Users/Daniel/anaconda/workspace/BARREL/{}_{}/bar_{}_l2_fspc_{}_v05.cdf'.format(date[0:4],pay,pay,date)
            cdffspc = pycdf.CDF(fspcpath)
        
        elif readInd == 'o':
            """downloads ephm cdf file from barreldata server, reads into program"""
            ephmurl = 'http://barreldata.ucsc.edu/data_products/v05/l2/{}/{}/bar_{}_l2_ephm_{}_v05.cdf'.format(pay,date[2:8],pay,date)   
            ephmdata = urllib.urlretrieve(ephmurl)
            cdfephm = pycdf.CDF(ephmdata[0])
            print '\nEPHM file downloaded.'
            """downloads fspc cdf"""
            fspcurl = 'http://barreldata.ucsc.edu/data_products/v05/l2/{}/{}/bar_{}_l2_fspc_{}_v05.cdf'.format(pay,date[2:8],pay,date)   
            fspcdata = urllib.urlretrieve(fspcurl)
            cdffspc = pycdf.CDF(fspcdata[0])    
            print 'FSPC file downloaded.'
            
        """user interface for time range"""
        tr = raw_input('Enter Time Range (e.g. HHMMstart-HHMMend), must be between 0001 and 2359: ')
        
        exitTime = 'n'
        
        """iterative loop that allows user to examine multiple time ranges without exiting the program"""
        while exitTime != 'e':
            """creates relevant integers from user input"""
            hstart = int(tr[0:2])
            mstart = int(tr[2:4])
            hend = int(tr[5:7])
            mend = int(tr[7:9])
            y = int(date[0:4])
            m = int(date[4:6])
            d = int(date[6:9])
    
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
    
            """divides fspc time into indexed increments"""
            start_indfspc = bisect.bisect_left(cdffspc['Epoch'], start)
            stop_indfspc = bisect.bisect_left(cdffspc['Epoch'], stop)    
            
            """reads in time from fpsc file"""
            tfspc = cdffspc['Epoch'][start_indfspc:stop_indfspc]
            
            """checks payload to determine light curve levels, reads in appropriate light curves"""
            if pay[0] == '1':    
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
            plt.figure()
            plt.subplot(5,1,1)
            plt.plot(tephm, alt)
            plt.title('Payload {}, {}-{}-{}, {}'.format(pay,d,m,y,tr))
            plt.gca().xaxis.set_major_formatter(myFmt)    
            plt.gca().set_ylim([25,40])
            plt.gca().set_xlabel('UT')
            plt.gca().set_ylabel('Altitude')         
            
            """user input for fspc curve smoothing"""
            smoothC = raw_input('Would you like to smooth the light curve data (y/n)?'
                                'Note, curve smoothing\nis best used on shorter time'
                                'ranges (<4 hours) to maximize results. A flat\nline '
                                'on a smoothed section indicates no activity: ')
            
            """error loop for improper input"""
            while ((smoothC != 'y') | (smoothC != 'n')):
                if ((smoothC == 'y') | (smoothC == 'n')):
                   break
                readInd = raw_input('ERROR: must enter either "y" or "n"\n'
                                    'Would you like to smooth the light curve data (y/n): ')
                
            """plots light curves according to user input"""
            """without smoothing"""
            if smoothC == 'n':
                plt.subplot(5,1,2)
                """checks to see whether to use FSPC1 or FSPC1a,b,c"""
                if pay[0] == '1':
                    plt.plot(tfspc,light1,color = 'red',label = 'level 1')
                    plt.legend(loc = 'upper right',frameon = False)        
                else:
                    plt.plot(tfspc,light1a,color = 'red',label = '1a')
                    plt.plot(tfspc,light1b,color = 'orange',label = '1b')
                    plt.plot(tfspc,light1c,color = 'yellow',label = '1c')
                    plt.legend(loc = 'upper right',frameon = False)
                plt.gca().xaxis.set_major_formatter(myFmt)
                
                """plots remaining light curves"""
                plt.subplot(5,1,3)
                plt.plot(tfspc,light2,color = 'green',label = 'level 2')
                plt.legend(loc = 'upper right',frameon = False)
                plt.gca().xaxis.set_major_formatter(myFmt)
                
                plt.subplot(5,1,4)
                plt.plot(tfspc,light3,color = 'blue',label = 'level 3')
                plt.legend(loc = 'upper right',frameon = False)
                plt.gca().xaxis.set_major_formatter(myFmt)
                
                plt.subplot(5,1,5)
                plt.plot(tfspc,light4,color = 'purple',label = 'level 4')
                plt.legend(loc = 'upper right',frameon = False)
                plt.gca().xaxis.set_major_formatter(myFmt)
                plt.show()
            
            """with smoothing"""
            if smoothC == 'y':
                plt.subplot(5,1,2)
                """checks to see whether to use FSPC1 or FSPC1a,b,c"""
                if pay[0] == '1':
                    """uses Savitzky Golay filtering (see function below) to smooth curve"""
                    light1sg = savitzky_golay(light1, window_size = 31, order = 4)
                    plt.plot(tfspc,light1sg,color = 'red',label = 'level 1')
                    plt.legend(loc = 'upper right',frameon = False)        
                else:
                    """uses Savitzky Golay filtering (see function below) to smooth curve"""
                    light1asg = savitzky_golay(light1a, window_size = 31, order = 4)
                    light1bsg = savitzky_golay(light1b, window_size = 31, order = 4)
                    light1csg = savitzky_golay(light1c, window_size = 31, order = 4)
                    plt.plot(tfspc,light1asg,color = 'red',label = '1a')
                    plt.plot(tfspc,light1bsg,color = 'orange',label = '1b')
                    plt.plot(tfspc,light1csg,color = 'yellow',label = '1c')
                    plt.legend(loc = 'upper right',frameon = False)
                plt.gca().xaxis.set_major_formatter(myFmt)
                
                """plots remaining light curves"""
                plt.subplot(5,1,3)
                """uses Savitzky Golay filtering (see function below) to smooth curve"""
                light2sg = savitzky_golay(light2, window_size = 31, order = 4)
                plt.plot(tfspc,light2sg,color = 'green',label = 'level 2')
                plt.legend(loc = 'upper right',frameon = False)
                plt.gca().xaxis.set_major_formatter(myFmt)
                
                plt.subplot(5,1,4)
                """uses Savitzky Golay filtering (see function below) to smooth curve"""
                light3sg = savitzky_golay(light3, window_size = 15, order = 4)
                plt.plot(tfspc,light3sg,color = 'blue',label = 'level 3')
                plt.legend(loc = 'upper right',frameon = False)
                plt.gca().xaxis.set_major_formatter(myFmt)
                
                plt.subplot(5,1,5)
                """uses Savitzky Golay filtering (see function above) to smooth curve"""
                light4sg = savitzky_golay(light4, window_size = 31, order = 4)
                plt.plot(tfspc,light4sg,color = 'purple',label = 'level 4')
                plt.legend(loc = 'upper right',frameon = False)
                plt.gca().xaxis.set_major_formatter(myFmt)
                plt.show()
                
            """prints out average lat. and long. over the time range to stdout"""
            lat = cdfephm['GPS_Lat'][start_indephm:stop_indephm]
            lon = cdfephm['GPS_Lon'][start_indephm:stop_indephm]
            """takes means of lat. and long."""
            avlat = mean(lat)
            avlon = mean(lon)
            print 'Average Latitude: {}'.format(avlat)
            print 'Average Longitude: {}'.format(avlon)
            
            """user interface to determine whether user would like to examine another time"""
            print 'Would you like to examine a different time range on this day?'
            exitTime = raw_input('Press return to select a new time range, or enter "e" to select a new date/exit program: ')
            
            if exitTime == 'e':
                break
            tr = raw_input('Enter Time Range (e.g. HHMMstart-HHMMend), must be between 0001 and 2359: ')
        
        """user interface to determine whether user would like to examine another date"""
        print '\nWould you like to examine a different date?'
        exitMain = raw_input('Press return to select a new date, or enter "e" to exit the program: ')
        if exitMain == 'e':
            break
        date = raw_input('Enter Date (e.g. YYYYMMDD), must be valid date for the payload: ')
        """add user interface for time range of event, avg. lat long, avg. alt"""
    
    """closes files"""
    cdfephm.close()
    cdffspc.close() 

if __name__ == "__main__":
    main()
