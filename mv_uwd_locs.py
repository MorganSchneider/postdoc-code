# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 13:49:59 2026

@author: mschne28
"""

# ERA5 file times: 00-23 UTC for all files

mv1 = {'name':"Mesomikenda Lake", #1525 EDT
      'date_ymd':[2024,6,13], 'time_utc':'1925', 'lat':47.6513, 'lon':-81.8473,
      'path_length':5140, 'path_width':1110}
mv2 = {'name':"Beauty Lake", #1635 EDT
      'date_ymd':[2024,6,13], 'time_utc':'2035', 'lat':47.5918, 'lon':-80.5881,
      'path_length':7820, 'path_width':800}
mv3 = {'name':"Lac Villebois", #1700 EDT
      'date_ymd':[2024,6,13], 'time_utc':'2100', 'lat':49.0670, 'lon':-78.8273,
      'path_length':2200, 'path_width':500}
mv4 = {'name':"Lac Fricourt", #1942 EDT
      'date_ymd':[2024,6,13], 'time_utc':'2342', 'lat':47.7485, 'lon':-76.4831,
      'path_length':76180, 'path_width':0}

mv5 = {'name':"Lac a Monette", #2055 EDT
      'date_ymd':[2025,6,24], 'time_utc':'0055', 'lat':46.9318, 'lon':-71.2139,
      'path_length':5060, 'path_width':760}

mv6 = {'name':"Thunder Lake", #1800 CDT
      'date_ymd':[2025,7,27], 'time_utc':'2300', 'lat':49.7880, 'lon':-92.5864,
      'path_length':8990, 'path_width':1670}
mv7 = {'name':"Melgund Lake", #1810 CDT
      'date_ymd':[2025,7,27], 'time_utc':'2310', 'lat':49.6493, 'lon':-92.4056,
      'path_length':4620, 'path_width':1020}
mv8 = {'name':"Stormy Lake", #1820 CDT
      'date_ymd':[2025,7,27], 'time_utc':'2320', 'lat':49.3369, 'lon':-92.3048,
      'path_length':14520, 'path_width':1340}
mv9 = {'name':"Kinmoapiku Lake", #1825 CDT
      'date_ymd':[2025,7,27], 'time_utc':'2325', 'lat':49.3184, 'lon':-92.0018,
      'path_length':3230, 'path_width':2500}


locs = {'mv1':mv1, 'mv2':mv2, 'mv3':mv3, 'mv4':mv4, 'mv5':mv5, 'mv6':mv6, 'mv7':mv7, 'mv8':mv8, 'mv9':mv9}


import pickle

dbfile = open("C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/era5/svr_winds/uwd_locs.pkl", 'wb')
pickle.dump(locs, dbfile)
dbfile.close()