import numpy as np
import time,sys,os
import Pyro4

from astropy.io import fits

from DOSlib.application import Application
from DOSlib.discovery import discoverable
from DOSlib.util import dos_parser
from DOSlib.advertise import Seeker


"""Demonstrator script for initializing / moving positioner.
"""

etime = sys.argv[1]
image_id = sys.argv[2]

#etime = 0.3 #sys.argv[1]
data_dir = os.getcwd()
print(data_dir)

# initialization done when starting petal application
role1 = 'FPC'
fpc = None


# Find FPC application
s = Seeker('-dos-','DOStest')
while fpc == None:
    s.seek()
    if role1 in s.devices:
        fpc = Pyro4.Proxy(s.devices[role1]['pyro_uri'])
        print('FPC connected')
    else:
        print('Not connecting to anything using Pyro')
    time.sleep(1)
    
#Configure FPC
fpc.configure()
fpc.set(image_dir=data_dir)

def create_fits(data,num,type):
    image = data['primary']['data']
    Ordered_dict_header = data['primary']['header']
    header_dict = dict(Ordered_dict_header)

    hdr = fits.Header()
    for key, value in header_dict.items():
        hdr[key] = value
        hdr['OBSTYPE'] = type
   
    hdu = fits.PrimaryHDU(image,header=hdr)
    
    hdu.writeto(data_dir+'/'+'FPC_'+str(num)+'.fits')
   
def fpc_images(id):
    """
    Takes a bias, dark, and object image, and spits out a fits file for each
    """
    num = int(id)*100
    fpc.expose(num+1,exptype='bias',exptime=etime,local_copy=True)
    time.sleep(1)
    bias_image = fpc.get_image(num+1)
    time.sleep(1)
    create_fits(bias_image,num+1,'BIAS')
    time.sleep(1)
    fpc.expose(num+2,exptype='dark',exptime=etime,local_copy=True)
    time.sleep(1)
    dark_image = fpc.get_image(num+2)
    time.sleep(1)
    create_fits(dark_image,num+2,'DARK')
    time.sleep(1)
    fpc.expose(num+3,exptype='object',exptime=etime,local_copy=True)
    time.sleep(1)
    object_image = fpc.get_image(num+3)
    time.sleep(1)
    create_fits(object_image,num+3,'OBJECT')

fpc_images(int(image_id))

    


    


