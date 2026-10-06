# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 17:00:25 2026

@author: mschne28
"""

from herbie import Herbie
from datetime import datetime


fp = "C:/Users/mschne28/OneDrive - The University of Western Ontario/Documents/qlcs_tornado_outbreaks/"

year = 2026
month = 9
day = 3
model_run = 18
fcst_hour = 3

time = datetime(year, month, day, model_run)
timestr = time.strftime("%Y-%m-%d %H:%M")



H = Herbie(timestr, model='hrrr', product='sfc', fxx=fcst_hour, save_dir=fp)
subset = H.download(r":[U|V]GRD:10 m above|HGT:surface", verbose=True)

H = Herbie(timestr, model='hrrr', product='prs', fxx=fcst_hour, save_dir=fp)
subset = H.download(r":[U|V]GRD:\d+ mb:|:HGT:\d+ mb:", verbose=True)
