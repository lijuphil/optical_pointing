# -*- coding: utf-8 -*-
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy.coordinates import AltAz
from datetime import datetime
from operator import itemgetter

import re, sys
import numpy as np
import astropy.units as u
import argParser_sources as argParser
from color import *


def get_telescope_loc(platform):

    """
    Get the astropy earth location object for use in calculating 
    AZ and EL of available sources.
    """

    # currently the same, but may be changed in the future?

    if platform == 'pathfinder':
        
        earth_loc = EarthLocation(lat=37.23333*u.deg, lon=-118.28229*u.deg,
                                  height=1222.0*u.m)

    if platform == 'testbed':

        earth_loc = EarthLocation(lat=37.23333*u.deg, lon=-118.28229*u.deg,
                                  height=1222.0*u.m)

    return earth_loc


def get_sources(source_file):

    """
    Parse the input source file into a dictionary for easier computation

    """

    # RA - 3rd column
    # Dec - 4th column

    source_dict = {}

    with open(source_file, 'r') as sources:
        
        for line in sources.readlines():
            data = line.strip().split()

            source_dict[data[1]] = {'sky_coord': SkyCoord(RA_TO_DEG(data[2]),
                                                            DEC_TO_DEG(data[3]),
                                                            frame='icrs'), 
                                    'magnitude': float(data[7])} 
    
    sources.close()

    return source_dict


def RA_TO_DEG(RA):

    """
    Convert input strings of RA into degrees
    """

    hms = re.split(':', RA)
    hours = float(hms[0])
    minutes = float(hms[1])
    seconds = float(hms[2])

    return ((hours + minutes/60.0 + seconds/3600.0)*15.0)*u.deg


def DEC_TO_DEG(DEC):

    """
    Convert input string of DEC into degrees.
    """

    dms = re.split(':', DEC)
    hours = float(dms[0])
    minutes = float(dms[1])
    seconds = float(dms[2])

    return (hours + minutes/60.0 + seconds/3600.0)*u.deg


def find_sources(earth_loc, source_dict, elevation, azimuth_span):

    """
    Go through source_dict and compute Az/El of each source. Print sources
    that are above specified elevation in degrees. 
    """

    source_list = []

    for name, data in source_dict.items():
        az_el = data['sky_coord'].transform_to(AltAz(obstime=datetime.utcnow(), location=earth_loc))
        if (elevation[1]*u.deg > az_el.alt > elevation[0]*u.deg):
            source_list.append([name, data['magnitude'], az_el])
            
    return source_list


def sort_sources(source_list):

    """
    Sort a list of tuples containing source name, magnitude, and az/el by magnitude.
    """

    source_list.sort(key=itemgetter(1), reverse=True)

    return source_list


def plot_sources(pl_src, pl_mag, pl_alt, pl_az, az_el_range):
    plt.clp()
    col = ['black','violet','iron','lilac','gray','red','basil','tomato','green','yellow','artic','gold','teal','blue','indigo']
    marker = '¤♠♦♣✢✦✧★☆✪✹❊✽❃❋※❄⚙✿❀❁❂●○◉⁌⁍▶▷◀◁▼▽▲△■□◆◇✙✚✖☑☒☀☼➣➤♤♧♡♢'

    if len(pl_src) != 0:
        if len(pl_src) > len(marker):
            marker = marker.join([marker]*int(np.ceil(len(pl_src)/len(marker))))
        if len(pl_src) > len(col):
            col = np.tile(col, int(np.ceil(len(pl_src)/len(col))))

        for i in range(len(pl_src)):
            lab = pl_src[i]
            mag = pl_mag[i]
            plt.scatter([pl_az[i]], [pl_alt[i]], point_marker = marker[i], point_color = col[i], label = str(lab)+','+str(mag))
            plt.legend
        plt.ticks(len(pl_src), len(pl_src))
        plt.xlabel('Azimuth')
        plt.ylabel('Elevation')
        plt.title('Sources in AZ {}--{} deg, EL {}--{}'.format(az_el_range[0],az_el_range[1],az_el_range[2],az_el_range[3]))
        plt.canvas_color('white')
        plt.show()
    

def print_sources(source_list, elevation,azimuth_span, plot):
    print('------------------------------------------------------------------------------')
    print(Color.F_Red,"{} total sources between Elevation {} -- {} degrees.".format(str(len(source_list)), 
                                                                str(elevation[0]), str(elevation[1])),Color.F_Default)

    az_chunks = int(360.0 / azimuth_span)
    
    for i in range(az_chunks):
        print ('------------------------------------------------------------------------------')
        print (Color.F_Green,'Azimuth',azimuth_span*(i),'--',azimuth_span*(i+1),'deg | Elevation', elevation[0],'--', elevation[1],'deg',Color.F_Default)

        pl_src, pl_mag, pl_alt, pl_az, az_el_range = [],[],[],[],[]
        for j in range(len(source_list)):
            alt_az = SkyCoord(source_list[j][2])
            az = alt_az.az*u.deg
            alt = alt_az.alt*u.deg

            if (azimuth_span*(i+1) > az.value > azimuth_span*(i)) and (elevation[1] > alt.value > elevation[0]):
                print("{0} : azimuth is {1.az:.5}, elevation is {2.alt:.4} and has magnitude {3}".format(source_list[j][0], 
                                                                     source_list[j][2], source_list[j][2], 
                                                                     source_list[j][1]))
                pl_src.append(source_list[j][0])
                pl_mag.append(source_list[j][1])
                pl_alt.append(alt.value)
                pl_az.append(az.value)

        az_el_range = [azimuth_span*(i),azimuth_span*(i+1),elevation[0],elevation[1]]
        if (plot == 'True'):
            plot_sources(pl_src, pl_mag, pl_alt, pl_az, az_el_range) # call plotting function

                
        if i+1 == az_chunks and (azimuth_span*(i+1) != 360):
            pl_src, pl_mag, pl_alt, pl_az, az_el_range = [],[],[],[],[]
            print ('------------------------------------------------------------------------------')
            print (Color.F_Green,'Azimuth',azimuth_span*(i+1),'-- 360 deg | Elevation', elevation[0],'--', elevation[1],'deg',Color.F_Default)
            
            for j in range(len(source_list)):
                alt_az = SkyCoord(source_list[j][2])
                az = alt_az.az*u.deg
                alt = alt_az.alt*u.deg
                if (360 >= az.value > azimuth_span*(i+1)) and (elevation[1] > alt.value > elevation[0]):
                    print("{0} : azimuth is {1.az:.5}, elevation is {2.alt:.4} and has magnitude {3}".format(source_list[j][0], 
                                                                         source_list[j][2], source_list[j][2], 
                                                                         source_list[j][1]))
                    pl_src.append(source_list[j][0])
                    pl_mag.append(source_list[j][1])
                    pl_alt.append(alt.value)
                    pl_az.append(az.value)

            az_el_range = [azimuth_span*(i),azimuth_span*(i+1),elevation[0],elevation[1]]
            if (plot == 'True'):
                plot_sources(pl_src, pl_mag, pl_alt, pl_az, az_el_range) # call plotting function


if __name__ == '__main__':

    args = argParser.get_arguments()

    if args.plot == 'True':
        import plotext as plt
    
    for i in range(len(args.elevation)):
        if args.elevation[i] > 90 or args.elevation[i] < 0:
            print ('ERROR: Elevation range should be between 0 and 90 degrees')
            print ('Enter: for e.g. -el 30 65')
            exit()
    
    if args.azimuth_span > 360:
        print ('ERROR: Azimuth span should be between 0 and 360 degrees')
        print ('Enter: for e.g. -az 20')
        exit()

    earth_loc = get_telescope_loc(args.platform)
    
    source_dict = get_sources(args.source_list)

    source_list = find_sources(earth_loc, source_dict, args.elevation, args.azimuth_span)

    source_list = sort_sources(source_list)

    print_sources(source_list, args.elevation, args.azimuth_span, args.plot)