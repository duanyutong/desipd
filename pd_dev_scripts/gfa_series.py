import numpy as np
import time
import Pyro4
import sys,os

from DOSlib.application import Application
from DOSlib.discovery import discoverable
from DOSlib.util import dos_parser
from DOSlib.advertise import Seeker

start_num = 1 #for filename
length = 100 #1hr
exp_time = [5,10,15,20,30,40.]
sleep_time = 5.
num_exp = int(length/(exp_time+1))
print("You will take %d GFA exposures with exposure time %f" % (num_exp,exp_time))

role1 = 'GFA0G'
gfa = None

s = Seeker('-dos-','GFAControl')
while gfa==None:
    s.seek()
    if role1 in s.devices:
        gfa=Pyro4.Proxy(s.devices[role1]['pyro_uri'])
        print("GFA connected")
    else:
        print("Not connected to GFA0G")
    time.sleep(1)
gfa.configure()
for i in exptime:
    gfa.prepare_for_exposure(start_num,exptime=i)
    gfa.expose(start_num)
    print("Exposure #%d just completed" % start_num)
    start_num += 1
    gfa.prepare_for_exposure(start_num,exptime=i)
    gfa.expose(start_num)
    start_num +=1
    print("Exposure #%d just completed" % start_num)
"""
for n in range(num_exp):
    num = start_num+n
    print("You are starting exposure %d with an exptime of %f" %(num,exp_time))
    gfa.prepare_for_exposure(num,exptime=exp_time)
    gfa.expose(num)
    print("GFA exposure %d was completed" % num)
    print("You have %f seconds to make your move ...." % sleep_time)
    time.sleep(sleep_time)
    print("Time is up. . . about to take another GFA image")
"""
print("Exposure sequence is complete.")
