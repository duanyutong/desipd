from petal import posmodel
from petal import posschedule
from petal import posmovetable
from petal import posstate
import posconstants as pc
import numpy as np
import time,sys
import Pyro4

from astropy.io import fits

from DOSlib.application import Application
from DOSlib.discovery import discoverable
from DOSlib.util import dos_parser
from DOSlib.advertise import Seeker


"""Demonstrator script for initializing / moving positioner.
"""
#Script instructions
num_grid = 5 # Number of points in one dimension of grid
command_grid = 'dXdY' # Type of grid want to map out
grid_size = 1
etime = 0.3 #sys.argv[1]
data_dir = '/data/images/fpc/'
local_dir = '/data/images/fpc/'

# initialization done when starting petal application
petal_id = 0
role1 = 'PETAL%d' % petal_id
role2 = 'FPC'
ptl = None
fpc = None


# Find petal & FPC applications
s = Seeker('-dos-','DOStest')
while ptl == None:
    s.seek()
    if role1 in s.devices:
        ptl = Pyro4.Proxy(s.devices[role1]['pyro_uri'])
        print('Petal connected')
    time.sleep(1)

while fpc == None:
    s.seek()
    if role2 in s.devices:
        fpc = Pyro4.Proxy(s.devices[role2]['pyro_uri'])
        print('FPC connected')
    else:
        print('Not connecting to anything using Pyro')
    time.sleep(1)
    
#Configure FPC
fpc.configure()
fpc.set(image_dir=data_dir)

# get configuration info
pos_ids = ptl.get('pos_ids')
print('Using positioners: ', repr(pos_ids))
fid_ids = ptl.get('fid_ids')
print('Using fiducials: ', repr(fid_ids))

print('INITIAL POSITION')
for pos_id in pos_ids:
    print('pos_id: ', ptl.get(posid=pos_id))
    print('expected: ', ptl.expected_current_position_str(pos_id))

def home_positioners(home=True):
    """ Home positioners. If don't want positioners to home, set to False and will fake homing by just saying the positioners are at the homed positions
    """
    if home==True:
        print('MOVE: homing')
        # ptl.set(key='CREEP_TO_LIMITS',value=True) # to force only creeping to hard stops
        ptl.request_homing(pos_ids)
        ptl.schedule_send_and_execute_moves()
    else:
        ptl.set(key=['POS_T','POS_P'],value=[-180,180]) # faking having just homed

def create_fits(data,num):
    image = data['primary']['data']
    Ordered_dict_header = data['primary']['header']
    header_dict = dict(Ordered_dict_header)

    hdr = fits.Header()
    for key, value in header_dict.items():
        hdr[key] = value
   
    hdu = fits.PrimaryHDU(image,header=hdr)
    
    hdu.writeto(local_dir+'FPC_'+str(num)+'.fits')
   

def move_and_fpc_measure(command,targets):
    """ Just a move command that lets you set the type of command. Using only the 'standard syntax'.

        INPUT:    command 'dXdy"  etc.
                  targets: dictionary of dictionaries
                           {pos_id: {'x':0,'y':0}}
    """
    targets = pc.listify2d(targets)
    for i,target in enumerate(targets):
        #print('MOVE: ' + command_grid + ' (' + str(target[0]) + ',' + str(target[1]) + ')')
        log_note = 'demo ' #+ command_grid + ' point ' + str(targets.index(target))

        requests = {}
        for pos_id in pos_ids:                 
            requests[pos_id] = {'command':command_grid, 'target':target, 'log_note':log_note}
        print(requests)

        ptl.request_targets(requests) 
        print("REQUEST targets was complete")
        ptl.schedule_send_and_execute_moves()
        fpc.expose(i,exptype='object',exptime=etime,local_copy=True)
        print("FPC has taken an exposure")
        fpc_image = fpc.get_image(i)
        fpc_image = create_fits(fpc_image,i)


def create_grid(num_point,size):
    """ Makes a target list from a dither grid from one-axis lenght of the dither pattern
        OUTPUT: [[x1,y1],[x2,y2]]
    """
    if isinstance(num_point,int):
        xx = np.linspace(-size,size,num_point)
        yy = np.linspace(-size,size,num_point)
        return [[x0, y0] for x0 in xx for y0 in yy]
    else:
        print('How many points? Must be integer')
    

# This is the script to take a grid of images
targets = create_grid(num_grid,grid_size)
#print(targets)
home_positioners()
move_and_fpc_measure(command_grid,targets)
#print("STARTING EXPOSURE")
#fpc.expose(1,exptype='object',exptime=etime)
#print("FPC image has been taken")
#fpc_image = fpc.get_image(1)
#fpc_image = create_fits(fpc_image,1)

    


    


