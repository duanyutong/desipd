#!/usr/bin/env python3
"""
   FPC.py

   This is the protoDESI FPC control application.

   This application provides the connection between DOS and the fiber photometry camera.
   The application supports both the Michigan sbigcam driver software as well as the INDI software
   For testing and development a simulator can be used to simulate the hardware.
   The configuration variable controller_type selects between these options.

   The following commands are implmented:
       configure       to configure the application and the controller
       get             to read status and configuration information
       set             to update select variables at run time
       expose          to take an image
       get_image       to return image data (in DOS format, ie header dictionary and numpy array)
       write_image     to write image to local disk
       send_image      to send image data to PML interface (in DOS format)

   Several shared variables published by this application:
       STATUS reflects the application's status
       CONNECTED is True when the controller is connected to the hardware and alive (ie responds to version command)

   There is no telemetry information produced by this application.

   Version History
   03/30/2016  v0.2.0     LB     Initial framework
   06/02/2016  v0.5.0     KH     Updated to device framework
   08/28/2016  v0.5.1     DYT    bug fixes and obstype header update
   09/21/2016  v0.6.0     KH     improved exposing status handling to avoid
                                 image shift defects on hardware level
   09/23/2016  v0.6.1     DYT    previous update rolled back earlier updates;
                                 re-apply changes from August
"""
from DOSlib.application import Application
import DOSlib.discovery as discovery
from DOSlib.PML import dos_connection
import DOSlib.mjd
import DOSlib.util
# from DOSlib.util import dos_parser
import DOSlib
import Pyro4
import datetime
import sys
import os
# import math
import time
# import imp
import threading
try:
    import matplotlib.pyplot as plt
    mpl_available = True
except:
    mpl_available = False
import random
import fitsio
from astropy.io import fits
import numpy as np
import glob
from collections import OrderedDict
from sbigcam import SBIGCam
#try:
#    import indiclient
#except:
#    print ('FPC: indi support not available')

tag = 'FPC'

class FPCcam(Application):
    commands = ['configure',
                'get',
                'set',
                'calibrate',
                'expose',
                'get_image',
                'write_image',
                'send_image',
                ]
    defaults = {'fpc_host' : 'decampi21',
                'fpc_port' : 7624,
                'image_dir' : '/data/images',
                'image_name' : 'PROTODESI',
                'controller_type' : 'SBIG',#'SIMULATOR'
                'software' : 'sbigcam',
                'local_copy' : False,
                'auto_purge' : False,
                'subtract_bias' : False,
                'subtract_dark' : False,
                'bias_frame' : None,
                'dark_frame' : None,
                'camera' : 'STi',
                'temp_regulation' : ['on', 'enable_autofreeze'],
                'temp_setpoint' : -10.0,
                'display_image' : False,
                'simulated_count' : 4,
                'status_update_rate' : 5.0,
                }

    # list of keywords NOT to include in FITS header
    fpc_headers = ['local_copy', 'expid', 'image_dir', 'image_name', 'auto_purge', 'cb_name',
                   'cb_function', 'cb_args', 'send_name', 'send_command', 'send_args']

    FITS_ORDER = ['NAXIS', 'NAXIS1', 'NAXIS2', 'EXTNAME', 'DEVICE', 'OBSNUM', 'OBSFRAME',
                  'TILEID', 'OBSTYPE', 'DATE-OBS', 'TIME-OBS',
                  'MJD-OBS', 'OPENSHUT', 'LST', 'RA', 'DEC', 'EXPREQ', 'EXPTIME',
                  'FILENAME', 'CCDTEMP', 'CCDSET', 'CCDTEMP', 'CCDPOWER', 'TECFROZE', 'SHAPE', 'DTYPE']

    def init(self):
        """ Initialization (called directly by the application framework) """
        self.info('Initializing')

        # FPC Status
        self.fpc_status = self.shared_variable('STATUS')
        self.fpc_status.publish()
        self.cam_connected = False
        self.camera_status = self.shared_variable('CAMERASTATUS')
        self.camera_status.publish()
        self.exp_ready = self.shared_variable('FPCREADY',group='EXPOSURE')
        self.exp_ready.publish(allowMultiplePublishers=True)

        # use an event to mark exposure status
        self.exposing = threading.Event()
        self.exposing.clear()
        # Update user information
        self.add_interlock_information(interlock = 'DOS', key ='%s_STATUS' % self.role,
                                       set_condition=['NOTCONNECTED', 'READY'],
                                       enabled=True)

        # configuration information
        self.image_dir = self.config['image_dir']
        if self.image_dir == '.':
            self.image_dir = os.getcwd()
        if not os.access(self.image_dir,os.W_OK):
            self.error('init: image directory not found or insufficient permissions to write image files %s' % str(self.image_dir))
            self.image_dir = None
        self.info('init: using image directory: %s' % self.image_dir)

        self.image_name = self.config['image_name']
        self.local_copy = True if 'T' in str(self.config['local_copy']).upper() else False
        self.auto_purge = True if 'T' in str(self.config['auto_purge']).upper() else False
        self.camera = self.config['camera']
        self.temp_regulation = self.config['temp_regulation'] if isinstance(self.config['temp_regulation'],list) else [self.config['temp_regulation']]
        self.temp_setpoint = self.config['temp_setpoint']
        self.status_update_rate = float(self.config['status_update_rate'])

        self.exposures = {'status' : 'NONE',
                          'expid' : None,
                          'expreq' : None,
                          'local_copy' : self.local_copy,
                          'auto_purge' : self.auto_purge,
                          'obstype' : None,
                          'header' : None,
                          'data' : None}

        try:
            if self.config['controller_type'] == 'SIMULATOR':
                self.client = FakeIndi(display_image = self.config['display_image'] if mpl_available else False,
                                       count = int(self.config['simulated_count']))
            elif self.config['software'] == 'indi':
                try:
                    import indiclient
                except Exception as e:
                    self.error('init: Exception loading indi software: %s' % str(e))
                    sys.exit()
                self.client = indiclient.IndiClient()
            else:
                self.client = SBIG_camera(self.camera)
        except Exception as e:
            self.error('init: Exception accessing camera: %s' % str(e))
            sys.exit()

        # Handle bias and dark frames
        self.bias_frame_type = 'external'
        self.dark_frame_type = 'external'
        self.subtract_bias = True if 'T' in str(self.config['subtract_bias']).upper() else False
        self.subtract_dark = True if 'T' in str(self.config['subtract_dark']).upper() else False
        if str(self.config['dark_frame']).upper() == 'NONE':
            self.dark_frame = None
        else:
            try:
               self.dark_frame = fitsio.read(self.config['dark_frame'])
            except Exception as e:
                rstring = 'init: Exception trying to read dark frame: %s' % str(e)
                self.error(rstring)
                self.dark_frame = None
        if str(self.config['bias_frame']).upper() == 'NONE':
            self.bias_frame = None
        else:
            try:
               self.bias_frame = fitsio.read(self.config['bias_frame'])
            except Exception as e:
                rstring = 'init: Exception trying to read bias frame: %s' % str(e)
                self.error(rstring)
                self.bias_frame = None

        # Connect to the camera (try 10 times)
        for x in range(10):
           self.client.setServer(self.config['fpc_host'], self.config['fpc_port'])
           if (not(self.client.connectServer(model=self.camera))):
               if(self.config["software"] != 'indi'):
                   self.warn("Unable to connect to camera, retrying")
                   self.sleep(0.5)
                   continue
               self.warn("No indiserver running on "+self.client.getHost()+":"+str(self.client.getPort()))
               self.cam_connected = False
           else:
               self.cam_connected = True
               break
        if not self.cam_connected:
            self.error("Cannot connect to the FPC. Will try again in configure().")
            self.fpc_status.write('NOTCONNECTED')

        # Setup discovery, change status to initialized
        if self.connected:
             self._setup_discovery(discovery.discoverable)
        self.info('Initialized')
        self.fpc_status.write('INITIALIZED')
        # call configure() when in device mode
        if self.device_mode == True:
            self.info('calling configure')
            retcode = self.configure()
            if 'FAILED' in retcode:
                raise RuntimeError('init: ' + retcode)


    def configure(self, *args, **kwargs):
        """
        Configure FPC application
        """

        # Parse arguments although currently configure doesn't use any
        try:
            a, kw = DOSlib.util.dos_parser(*args, **kwargs)
        except:
            a = []
            kw = {}
        """ Configure the FPC """
        # Reset app status
        self.fpc_status.write('INITIALIZED')
        self.info('configure: starting')
        self.exposing.clear()

        if not self.cam_connected:
            self.info('Trying to connect to FPC (%s, %s)' % (str(self.config['fpc_host']), str(self.config['fpc_port'])))
            # Connect to the camera (try 10 times)
            for x in range(10):
                self.client.setServer(self.config['fpc_host'], self.config['fpc_port'])
                if (not(self.client.connectServer())):
                    if self.config["software"] != 'indi':
                        self.warn("Unable to connect to camera, retrying")
                        self.sleep(0.5)
                        continue
                    self.warn("No indiserver running on "+self.client.getHost()+":"+str(self.client.getPort()))
                    self.cam_connected = False
                else:
                    self.cam_connected = True
                    break
            if not self.cam_connected:
                rstring = 'configure: Failed to connect to camera.'
                self.error(rstring)
                self.fpc_status.write('ERROR')
                return 'FAILED: ' + rstring
        else:
            # ask server to reset internal state
            try:
                self.client.configure()
            except Exception as e:
                rstring = 'configure: Exception configuring client: %s' % str(e)
                self.error(rstring)
                return 'FAILED: ' + rstring
        self.info("Now verifying FPC connection")
        for x in range(10):
            try:
                time.sleep(0.2)
                if(self.client.getFrame() == None):
                    self.cam_connected = False
                    self.warn("Connection to camera not verified, retry %i"%x)
                    continue
            except Exception as err:
                self.warn("Connection not verified, retry %i"%x)
                self.warn("Error was %s"%repr(err))
                self.cam_connected = False
                continue
            self.cam_connected = True
            break
        if not self.cam_connected:
            self.info('Allowing configuration to complete but all FPC commands are disabled.')
            self.fpc_status.write('NOTCONNECTED')
            return self.SUCCESS

        # get camera status
        self._get_camera_status()

        # Set temperature regulation for SBIG cameras
        if self.config['controller_type'] == 'SBIG':
            try:
                for item in self.temp_regulation:
                    self.info('configure: setting temperature regulation mode: %s, Setpoint %f' % (item, self.temp_setpoint))
                    self.client.set_temperature_regulation(item, self.temp_setpoint)
            except Exception as e:
                self.error('configure: Exception setting temperature regulation: %s' % str(e))

        # Set status, recovery, and, and configuration status and then return
        self.fpc_status.write('READY')
        return self.SUCCESS

    # callback functions for SV and discovery and discovery setup
    def _setup_discovery(self, discoverable):
        # Setup application for discovery and discover other DOS applications
        discoverable(role = self.role, tag = 'FPC', interface = self.role)
        self.info('_setup_discovery: Done')

    def get(self, *args, **kwargs):
        """
        Retrieve information about the FPC application and the FPC at runtime.
        Allowed values for params include:
            status          returns the application status
            exposure        returns information about the exposure
            region, frame   returns ROI information
            size            returns CCD size information
            binning         returns the binning settings
            image_dir       returns the directory used to store images
            image_name      returns the current image base name
            local_copy      returns the state of the local_copy flag
            auto_purge      returns the state of the auto_purge flag
            subtract_bias   Returns True if bias should be subtracted
            subtract_dark   Returns True if a (scaled) dark frame should be subtracted
            bias_frame      Returns filename of bias frame (or "Internal" if self calibrated)
            dark_frame      Returns filename of dark frame (or "Internal" if self calibrated)
            connected       returns true if a camera is connected
            controller      returns the controller type
            software        returns the driver software (indi/sbig)
            display_image   returns the state of the display_image flag
            simulated_count returns the number of simulated spots
            temp_regulation returns regulation modes (SBIG only)
            temp_setpoint   returns temperature setpoint (SBIG only)
            camera_status   return camera status info (SBIG only)
        """
        try:
            a, kw = DOSlib.util.dos_parser(*args, **kwargs)
        except:
            a = []
            kw = {}

        if(len(a) > 2 or len(a) == 0):
            return "FAILED: Invalid arguments %s"%(str(a))
        params = a[0].strip().lower()
        if params in ['state', 'status']:
            return self.fpc_status._value
        elif params.startswith('exposure'):
            copy = dict(self.exposures)
            if 'data' in copy and copy['data'] != None:
                copy['data'] = True
            if 'header' in copy and copy['header'] != None:
                copy['header'] = True
            return copy
        elif params in ['region', 'frame']:
            if(not self.cam_connected):
                return "FAILED: FPC is not connected"
            try:
                return self.client.getFrame()
            except Exception as e:
                self.cam_connected = False
                self.warn("Connection lost: %s", repr(e))
                return "FAILED: Connection lost"
        elif params == 'binning':
            if not self.cam_connected :
                return "FAILED: FPC is not connected"
            try:
                return self.client.getBinning()
            except Exception as e:
                self.cam_connected = False
                self.warn("Connection lost")
                return "FAILED: Connection lost"
        elif params == 'dark_frame':
            if self.dark_frame_type == 'external':
                return self.config['dark_frame']
            elif self.dark_frame is not None:
                return self.dark_frame #'Internal'
            else:
                return None
        elif params == 'bias_frame':
            if self.bias_frame_type == 'external':
                return self.config['bias_frame']
            elif self.bias_frame is not None:
                return self.bias_frame #'Internal'
            else:
                return None
        elif params == 'subtract_bias':
            return self.subtract_bias
        elif params == 'subtract_dark':
            return self.subtract_dark
        elif params == 'temp_regulation':
            return self.temp_regulation
        elif params == 'temp_setpoint':
            return self.temp_setpoint
        elif params == 'fan':
            t = self._get_camera_status()
            r = {}
            for key in ('fan_power', 'fan_speed', 'fan_enabled'):
                if key in t:
                    r[key] = t[key]
            return r
        elif params == 'local_copy':
            return self.local_copy
        elif params.startswith('auto_purge'):
            return self.auto_purge
        elif params == 'image_dir':
            return self.image_dir
        elif params == 'image_name':
            return self.image_name
        elif params == 'connected':
            return self.cam_connected
        elif params.startswith('controller'):
            return self.config['controller_type']
        elif params == 'camera':
            return self.camera
        elif params == 'camera_status':
            return self.camera_status._value
        elif(params.lower() in self.config.keys()):
            return self.config[params]
        else:
            return 'FAILED: Invalid parameter for get command.'

    def set(self, *args, **kwargs):
        """
        Set application configuration parameters including:
            image_dir        sets the directory for local images
            image_name       sets the image base name
            region           sets the ROI. Syntax: region=[<x>,<y>,<width>,<height>]
            binning          sets binning. Syntax: binning=[horizontal,vertical]
            local_copy       if True, a copy of the image is written to disk in FITS format
            auto_purge       if True, image memory and local copies will be purged
            subtract_dark    Sets the subtract dark flag
            subtract_bias    Sets the subtract bias flag
            fan              Turns fan on/off  (SBIG only)
            temp_regulation  Sets temperature control mode (list, e.g. [on, enable_autofreeze])
                             and setpoint
        """
        try:
            a, kw = DOSlib.util.dos_parser(*args, **kwargs)
        except:
            a = []
            kw = {}
        if len(a) == 2:
            key = str(a[0]).lower()
            value = a[1]
        elif len(list(kw.keys())) == 1:
            key = str(list(kwargs.keys())[0]).lower()
            value = kwargs[key]
        else:
            return "FAILED: Invalid parameters for set command."
        if key == 'region':
            if(not self.cam_connected):
                return "FAILED: FPC is not connected"
            try:
                if (type(value) != list and type(value) != tuple):
                    return "FAILED: Region must be a list, [x, y, width, height]"
                if len(value) != 4:
                    return "FAILED: Region must have 4 elements, [x, y, width, height]"
            except Exception as errmsg:
                return "FAILED: %s"%str(errmsg)
            try:
                self.client.setFrame(value[0], value[1], value[2], value[3])
                return self.SUCCESS
            except Exception as e:
                self.cam_connected = False
                self.warn("Connection lost")
                return "FAILED: Connection lost"
        elif key == 'binning':
            if(not self.cam_connected):
                return "FAILED: FPC is not connected"
            try:
                if (type(value) != list and type(value) != tuple):
                    return "FAILED: Parameters must be a list, [horizonal, vertical]"
                if len(value) != 2:
                    return "FAILED: Tuple must have 2 elements, [horizontal, vertical]"
            except Exception as errmsg:
                return "FAILED: %s"%str(errmsg)
            try:
                self.client.setBinning(value[0], value[1])
                return self.SUCCESS
            except Exception as e:
                self.cam_connected = False
                self.warn("Connection lost")
                return "FAILED: Connection lost"
        elif key == 'subtract_bias':
            if value == 'TRUE':
                self.subtract_bias = True
                return self.SUCCESS
            elif value == 'FALSE':
                self.subtract_bias = False
                return self.SUCCESS
            else:
                return 'FAILED: invalid argument for set subtract_bias command'
        elif key == 'subtract_dark':
            if value == 'TRUE':
                self.subtract_dark = True
                return self.SUCCESS
            elif value == 'FALSE':
                self.subtract_dark = False
                return self.SUCCESS
            else:
                return 'FAILED: invalid argument for set subtract_dark command'
        elif 'local_copy' in kw:
            if str(kw['local_copy']).lower() == 'true':
                self.local_copy = True
                return self.SUCCESS
            if str(kw['local_copy']).lower() == 'false':
                self.local_copy = False
                return self.SUCCESS
            else:
                return 'FAILED: invalid argument for set local_copy command'
        elif 'auto_purge' in kw:
            if str(kw['auto_purge']).lower() == 'true':
                self.auto_purge = True
                return self.SUCCESS
            if str(kw['auto_purge']).lower() == 'false':
                self.auto_purge = False
                return self.SUCCESS
            else:
                return 'FAILED: invalid argument for set auto_purge command'
        elif 'image_dir' in kw:
            self.image_dir = kw['image_dir']
            if not os.access(kw['image_dir'],os.W_OK):
                return 'FAILED: image directory not found or insufficient permissions to write image files %s' % str(kw['image_dir'])
            else:
                self.image_dir = kw['image_dir']
                return self.SUCCESS
        elif 'image_name' in kw:
            self.image_name = kw['image_name']
        elif 'fan' in kw:
            if self.config['controller_type'] == 'SBIG':
                if kw['fan'] in ('on', 'ON', True, 1):
                    return self.client.setFan('on')
                elif kw['fan'] in ('off', 'OFF', False, 0):
                    return self.client.setFan('off')
                else:
                    return 'FAILED: Invalid parameter for set fan command'
            else:
                return 'FAILED: set fan command not implemented for this camera'
        elif 'temp_regulation' in kw:
            if self.config['controller_type'] == 'SBIG':
                if 'temp_setpoint' in kw:
                    setpoint = kw['temp_setpoint']
                else:
                    setpoint = None
                l = kw['temp_regulation'] if isinstance(kw['temp_regulation'], list) else [kw['temp_regulation']]
                for item in l:
                    retcode = self.client.set_temperature_regulation(item, setpoint)
                    if retcode == False:
                        return 'FAILED: error setting temperature regulation for item %s' % item
                return self.SUCCESS
            else:
                return 'FAILED: Invalid parameter for set temp_regulation command'
        else:
            return 'FAILED: Invalid arguments for set command'
        return self.FAILED

    def calibrate(self, bias = True, nbias = 5, dark = True, ndark = 5, darktime = 2.0, bias_outfile = None, dark_outfile = None):
        """
        Generate master bias and dark frames
        Images are stacked and median combined
        The master dark frame is scaled to 1 second. darktime is the exposure time used for darks
        """
        if bias == True:
            self.info('calibrate: Preparing master bias frame (count = %s)' % str(nbias))
            bias_data = []
            try:
                for i in range(nbias):
                    print(i)
                    self.expose(1, exptime = 0, auto_purge=False, exptype='bias', local_copy=False)
                    bias_data.append(np.copy(self.exposures['data']))
                bias_stack = np.array(bias_data)
                self.bias_frame = np.median(bias_stack, axis = 0)
                self.bias_frame_type = 'internal'
                if bias_outfile != None:
                    with fitsio.FITS(bias_outfile,'rw',clobber=True) as fits:
                        fits.write(self.bias_frame)
                        fits[-1].write_checksum()                 # add checksum
                        self.info('calibrate: master bias frame generated and written to %s' % bias_outfile)
                else:
                    self.info('calibrate: master bias frame generated')
                # cleanup
                del bias_stack
            except Exception as e:
                rstring = 'calibrate: Exception preparing master bias frame: %s' % str(e)
                self.error(rstring)
                return 'FAILED: ' + rstring
            del bias_data
        if dark == True:
            self.info('calibrate: Preparing master dark frame (count = %s)' % str(ndark))
            dark_data = []
            try:
                for i in range(ndark):
                    self.expose(1, exptype='dark',exptime=darktime, auto_purge=False, local_copy=False)
                    dark_data.append(np.copy(self.exposures['data']))
                dark_stack = np.array(dark_data)
                self.dark_frame = np.median(dark_stack, axis = 0)
                self.dark_frame_type = 'internal'
                # subtract master bias frame before scaling
                if not self.bias_frame is None:
                    self.dark_frame = self.dark_frame - self.bias_frame
                else:
                    self.warn('calibrate: no master bias frame available. Will scale without bias corrections')
                # scale master dark frame
                self.dark_frame = self.dark_frame/float(darktime)
                if dark_outfile != None:
                    with fitsio.FITS(dark_outfile,'rw',clobber=True) as fits:
                        fits.write(self.dark_frame)
                        fits[-1].write_checksum()                 # add checksum
                        self.info('calibrate: master dark frame generated and written to %s' % dark_outfile)
                else:
                    self.info('calibrate: master dark frame generated')
                # cleanup
                del dark_stack
            except Exception as e:
                rstring = 'calibrate: Exception preparing master dark frame: %s' % str(e)
                self.error(rstring)
                return 'FAILED: ' + rstring
            del dark_data
        self.info('calibrate: Done')
        return self.SUCCESS

    def expose(self, expid, *args, **kwargs):
        """
        Take an image  with the FPC.
        Optionally write image in FITS format to file.
        Positional arguments
           expid           exposure number (required)
        Named arguments:
           exptype         object, zero, dark, image, dome_flat
           exptime         exposure time (float seconds)
           local_copy      write a local copy in FITS format if True
           auto_purge      delete old image data in memory and on disk if local_copy is true
           cb_name         Name of PML interface (or proxy) of callback. If None, no callback is performed
           cb_function     Callback function
           cb_args         Arguments passed to callback function
           send_name       Name of PML interface (or proxy) to send image. If None, the image is not sent
           send_command    Command to be used when sending image (send_name and send_command are both required
                           if the image should be sent at the end of digitize and after build (will be forwarded to
                           build()))
           send_args       Arguments passed to send command
           bias_subtract   Subtract bias frame for local_copy or send option
           dark_subtract   Subtract dark frame for local_copy or send option

        The rest of the named arguments passed to the routine but not listed in fpc_headers
        will be added to the header. (expid will be renamed to EXPNUM).
        Returns a dictionary with status information on completion
        """
        try:
            a, kw = DOSlib.util.dos_parser(*args, **kwargs)
        except:
            a = []
            kw = {}
        if(not self.cam_connected):
            rstring = 'expose: FPC is not connected'
            self.error(rstring)
            return "FAILED: " + rstring
        if self.fpc_status._value != 'READY':
            rstring = 'expose: Application is not ready to take an exposure'
            self.error(rstring)
            return "FAILED: " + rstring

        if expid in ['next']:
            last = self.exposures['expid']
            expid = 1 if last == None else last+1

        # setup exposure management
        local_copy = self.config['local_copy'] if 'local_copy' not in kw else kw['local_copy']
        auto_purge = self.config['auto_purge'] if 'auto_purge' not in kw else kw['auto_purge']
        # Callback setup
        cb_name =  None if 'cb_name' not in kw else kw['cb_name']
        cb_function =  None if 'cb_function' not in kw else kw['cb_function']
        cb_args =  None if 'cb_args' not in kw else kw['cb_args']
        send_name = None if 'send_name' not in kw else kw['send_name']
        send_command = None if 'send_command' not in kw else kw['send_command']
        send_args = None if 'send_args' not in kw else kw['send_args']

        if 'exptime' not in kw:
            rstring = 'expose requires an exposure time'
            self.error(rstring)
            return 'FAILED: ' + rstring
        else:
            exptime = float(kw['exptime'])

        if 'exptype' not in kw:
            rstring = 'expose requires a exposure type.'
            self.error(rstring)
            return 'FAILED: ' + rstring
        elif str(kw['exptype']).lower() not in ['object', 'dark', 'image', 'zero', 'bias','dome flat', 'guider', 'focus']:
            rstring = 'expose requires a exposure type. Not %s.' % kw['exptype']
            self.error(rstring)
            return 'FAILED: ' + rstring
        else:
            exptype = str(kw['exptype']).lower()

        self.exposures = {'status' : 'PREPARED',
                          'expid' : expid,
                          'expreq' : exptime,
                          'local_copy' : local_copy,
                          'auto_purge' : auto_purge,
                          'obstype' : exptype,
                          'header' : None,
                          'data' : None}

        # Callback setup
        self.exposures['cb_name'] =  None if 'cb_name' not in kw else kw['cb_name']
        self.exposures['cb_function'] =  None if 'cb_function' not in kw else kw['cb_function']
        self.exposures['cb_args'] =  None if 'cb_args' not in kw else kw['cb_args']
        self.exposures['send_name'] = None if 'send_name' not in kw else kw['send_name']
        self.exposures['send_command'] = None if 'send_command' not in kw else kw['send_command']
        self.exposures['send_args'] = None if 'send_args' not in kw else kw['send_args']

        # setup for header and data
        self.info("expose (%d): exptype %s, exptime %f" % (expid,exptype, exptime))
        self.exposures['status'] = 'EXPOSING'
        if local_copy:
            image_dir = self.image_dir if 'image_dir' not in kw else kw['image_dir']
            image_name = self.image_name if 'image_name' not in kw else kw['image_name']
            image_name.replace('.fits','')
            image_name +='_%s_%08d.fits' % (self.role, expid)
            if not image_name.endswith('.fits'):
                image_name += '.fits'
            filename = os.path.join(image_dir, image_name)

        self.exposures['data'] = None
        self.exposures['header'] = OrderedDict()
        self.exposures['header']['DEVICE'] = self.role
        self.exposures['header']['INSTRUME'] = 'PROTODESI'
        self.exposures['header']['EXPNUM'] = expid
        self.exposures['header']['EXPREQ'] = exptime
        self.exposures['header']['OBSTYPE'] = exptype
        dt = datetime.datetime.utcnow()
        date_obs = dt.isoformat().split('+')[0]
        m = float('%.8f' % DOSlib.mjd.dt2mjd(dt))
        self.exposures['header']['DATE-OBS'] = date_obs
        self.exposures['header']['TIME-OBS'] = dt.time().isoformat()
        self.exposures['header']['MJD-OBS'] =  m
        self.exposures['header']['OPENSHUT'] = date_obs
        self.exposures['header']['LST'] = DOSlib.util.sidereal_time()
        if local_copy:
            self.exposures['header']['FILENAME'] = filename
        # Add telemetry information to header
        s = self._get_camera_status()
        try:
            self.exposures['header']['CCDCOOL']  = s['cooling_enabled']
            self.exposures['header']['CCDTEMP']  = s['imaging_ccd_temperature'] if not isinstance(s['imaging_ccd_temperature'],(float, int)) else round(s['imaging_ccd_temperature'],3)
            self.exposures['header']['CCDPOWER'] = s['imaging_ccd_power'] if not isinstance(s['imaging_ccd_power'],(float, int)) else round(s['imaging_ccd_power'],3)
            self.exposures['header']['CCDSET']   = s['ccd_setpoint']  if not isinstance(s['ccd_setpoint'],(float, int)) else round(s['ccd_setpoint'],3)
            self.exposures['header']['TECFROZE'] = s['tec_frozen']
        except:
            self.debug('expose (%d): some (or all) camera telemetry is not available', expid)

        self.exposing.set()
        try:
            if exptype in ['object', 'image', 'dome flat', 'full', 'guider', 'focus']:
                self.client.setLight()
                self.client.takeExposure(exptime)
            elif exptype == 'dark':
                self.client.setDark()
                self.client.takeExposure(exptime)
            elif exptype in ['zero', 'bias']:
                self.client.setBias()
                self.client.takeExposure(exptime)
            else:
                rstring =  "expose: '%s' is not a valid object type"%exptype
                self.error(rstring)
                return 'FAILED: ' + rstring
            self.exposures['status'] = 'EXPOSED'
        except Exception as e:
            rstring = 'expose: Exception taking exposure %d: %s' % (expid, str(e))
            self.error(rstring)
            self.exposing.clear()
            return 'FAILED: ' + rstring

        # Retrieve data (getData for sbig returns the image, getBlob returns a fits HDUlist)
        blob = None
        try:
            if self.config['software'] == 'indi':
                blob = self.client.getBlob()
            else:
                blob = self.client.getData()
            i = 0
            while(not isinstance(blob, np.ndarray)):
                self.sleep(0.5)
                if self.config['controller_type'] != 'SIMULATOR' and self.config['software'] == 'indi':
                    blob = self.client.getBlob()
                else:
                    blob = self.client.getData()
                i = i+1
                if(i > (exptime*8) + 15):
                    raise TimeoutError("Blob did not return")
        except Exception as e:
            self.cam_connected = False
            self.warn("Connection lost")
            rstring = 'expose: Exception reading data for exposure %d from camera: ' % (expid, str(e))
            self.error(rstring)
            self.exposing.clear()
            return 'FAILED: ' + rstring
        self.exposing.clear()
        self.info("expose (%d): data received" % expid)
        if self.config["controller_type"] == "SIMULATOR":
            self.exposures['data'] = blob
        else:
            if self.config["software"] != 'indi':
                self.exposures['data'] = blob
            else:
                hdulist = fits.HDUList.fromstring(blob)
                self.exposures['data'] = hdulist[0].data
        # Are we supposed to do something with the image:
        if local_copy:
            retcode = self._saving(expid, filename,
                                   bias_subtract = False if 'bias_subtract' not in kw else kw['bias_subtract'],
                                   dark_subtract = False if 'dark_subtract' not in kw else kw['dark_subtract']
                                   )
            if 'FAILED' in retcode:
                self.error('expose: _saving returns: %s' % str(retcode))
        # Send_image?
        if send_name != None and send_command != None:
            try:
                retcode = self.send_image(expid, pml_name = send_name, pml_command = send_command,
                                    bias_subtract = False if 'bias_subtract' not in kw else kw['bias_subtract'],
                                    dark_subtract = False if 'dark_subtract' not in kw else kw['dark_subtract'],
                                    extra_headers = None if 'send_args' not in kw else kw['send_args'])
            except Exception as e:
                self.error('expose: send_image returns %s. Exception: %s' % (repr(retcode), str(e)))
            self.exposures['status'] = 'SENT'

        if self.exposures['cb_name'] != None and self.exposures['cb_function'] != None:
            cb = {'expid' : expid,
                  'source' : self.role,
                  'status' : self.exposures['status'],
                  'cb_args' : cb_args
                  }
            try:
                if str(cb_name).upper().startswith('PYRO'):
                    # Use Pyro4.Proxy directly
                    retcode = Pyro4.Proxy(cb_name).__getattr__(cb_function)(**cb)
                else:
                    retcode = dos_connection(cb_name).execute(cb_function, **cb)
            except Exception as msg:
                self.error('expose: Exception in callback to %s: %s. Arguments: %s, %s' % (str(cb_name), str(msg), repr(args), repr(kw)))


        # Need to clean up memory and disk space?
        if auto_purge:
            # only single copy in memory - nothing to do
            # Remove old image files if local_copy is set
            if local_copy:
                fits_files = sorted(glob.glob(os.path.join(image_dir,'*.fits')), key = self._modified)
                if len(fits_files) > 2:   # set proper threshold (handle sub frames)
                    # newest file is at the end of the list
                    for index in range(len(fits_files) - 2):
                        try:
                            to_remove = os.path.join(image_dir, fits_files[index])
                            os.remove(to_remove)
                            self.debug('expose: removed image file %s' % to_remove)
                        except Exception as e:
                            self.error('expose: Exception removing old fits files: %s' % str(e))

        # We are done
        return self.SUCCESS

    def send_image(self, expid, pml_name = None, pml_command = None, bias_subtract = False, dark_subtract = False, extra_headers = None):
        """
        Send data for exposure = expid to an PML interface (name or proxy)
        Required arguments: expid  (either first positional argument or as a named argument)
                            pml_name (name of PML interface)
                            pml_command (function of PML interface to be called)
        The exposure specified by expid must be in the (internal) exposure status structure.
        """
        if expid != self.exposures['expid']:
            return 'FAILED: send_image: Exposure %s is not available' % str(expid)

        if self.exposures['status'] not in ['EXPOSED', 'SENT']:
#            return 'FAILED: send_image: Exposure %d.%d is not digitized' % (expid, frame)
            return 'FAILED: send_image: Exposure %d is not digitized' % expid
        if extra_headers != None  and isinstance(extra_headers, dict):
            for key in list(extra_headers.keys()):
                self.exposures['header'][str(key).strip().upper()[:8]] = extra_headers[key]    # Careful - overwrites existing values if same keys

        # Detrend?
        image = self._detrend(self.exposures, bias_subtract, dark_subtract)

        # Now send the image(s)
        extras = {'source' : self.role,
                  'status' : self.exposures['status'],
                  'expid' : expid,
                  'primary' : image['primary']
                      }
        try:
            self.info('send_image (%d): sending exposure to %s' % (expid,str(pml_name)))
            # We are setting up a new connection everytime - maybe this can be done more efficiently
            if str(pml_name).upper().startswith('PYRO'):
                # Use Pyro4.Proxy directly
                with Pyro4.Proxy(pml_name) as receiver:
                    retcode = receiver.__getattr__(pml_command)(extras)
            else:
                receiver = dos_connection(pml_name)
                retcode = receiver.execute(pml_command, extras)
                del receiver
            if 'FAILED' in str(retcode):
                self.error('send_image: sending exposure %d returns: %s' % (expid, str(retcode)))
            self.info('send_image (%d): done' % expid)
            return retcode
        except Exception as e:
            rstring = 'send_image: Exception calling %s: %s.' % (str(pml_name), str(e))
            self.error(rstring)
            return 'FAILED: ' + rstring

    def get_image(self, expid, bias_subtract = False, dark_subtract = False, extra_headers = None):
        """
        Returns an image
        Required arguments: expid (int)
        The requested image must be in the internal exposure status structure.
        Any additional headers will be added to the image header.
        (Careful, keys given in extra_headers can potentially overwrite already existing keys)
        Images will be bias subtracted if the flag is true, they will be dark and bias subtracted when
        dark_subtract is true (overwrites bias_subtract setting)
        The call fails if either of the subtract flags is set True but the calibration frames are not available
        """
        if expid in ['last', -1]:
            if self.exposures['status'] not in ['EXPOSED', 'SENT']:
                return 'FAILED: no image available.'
            expid = self.exposures['expid']
        elif expid != self.exposures['expid']:
            raise RuntimeError('Invalid exposure number. Image data not available.')

#        header = {}
        if isinstance(extra_headers, dict):
            for key in list(extra_headers.keys()):
                self.exposures['header'][str(key).strip().upper()[:8]] = extra_headers[key]    # Careful - overwrites existing values if same keys

        # detrend image?
        image = self._detrend(self.exposures, bias_subtract, dark_subtract)

        # Now send the image(s)
        extras = {'source' : self.role,
                  'status' : self.exposures['status'],
                  'expid' : expid,
                  'primary' : image['primary']
                      }
        self.info('Returning exposure %d' % (expid))
        return extras

    def write_image(self, expid, image_dir = None, image_name = None,
                    bias_subtract = False, dark_subtract = False, extra_headers = None, overwrite = True):
        """
        Save image to disk
        Required arguments:    expid (int)
        Optional arguments:    image_name
                               image_dir       Use image_name/image_dir to overwrite current values.
                                               If not given the internal values will be used to construct the filename.
                               overwrite       If True, an existing file with the same name will be overwritten.
                               dark_subtract   a scaled dark frame and a bias frame will be subtracted
                               bias_subtract   a bias frame will be subtracted

        Any additional headers will be added to the image header.
        (Careful, keys given in extra_headers can potentially overwrite already existing keys)
        """
        if expid in ['last', -1]:
            if self.exposures['status'] not in ['EXPOSED', 'SENT']:
                return 'FAILED: no image available.'
            expid = self.exposures['expid']
        elif expid != self.exposures['expid']:
            rstring = 'write_image: invalid exposure number %s. Image data not available.' % str(expid)
            self.error(rstring)
            return 'FAILED: ' + rstring

        if self.exposures['status'] not in ['DIGITIZED', 'EXPOSED', 'SENT']:
#            return 'FAILED: write_image: Exposure %d.%d is not ready' % (expid, frame)
            return 'FAILED: write_image: Exposure %d is not ready' % expid

        dir = image_dir if image_dir != None else self.image_dir
        name = image_name if image_name != None else self.image_name
        complete_name =name + '_%s_%08d.fits' % (self.role, expid)
        filename = os.path.join(dir, complete_name)
        self.info('write_image: using output file %s' % filename)

        # Extra keywords
        if isinstance(extra_headers, dict):
            for key in list(extra_headers.keys()):
                self.expposures['header'][str(key).strip().upper()[:8]] = extra_headers[key]    # Careful - overwrites existing values if same keys
        # call io routine
        try:
            if overwrite and os.path.isfile(filename):   # does the file already exist?
                os.remove(filename)
            self._saving(expid, filename, bias_subtract, dark_subtract)               # output file
        except Exception as e:
            rstring = 'write_image: Exception writing image file: %s' % str(e)
            self.error(rstring)
            return 'FAILED: ' + rstring
        self.info('write_image: Done')
        return self.SUCCESS

    def _saving(self, expid, filename, bias_subtract = False, dark_subtract = False):
        """ Thread to write image to local fits file """

        # Try to order the fits header keywords
        ordered_header = OrderedDict()
        for item in self.FITS_ORDER:
            if item in self.exposures['header']:
                ordered_header[item] = self.exposures['header'][item]
        # get the rest
        ordered_header.update(self.exposures['header'])
        # detrend if necessary
        image = self._detrend(self.exposures, bias_subtract, dark_subtract)
        try:
            if os.path.isfile(filename):                  # does the file already exist?
                os.remove(filename)
            with fitsio.FITS(filename,'rw',clobber=True) as fits:
                fits.write(image['primary']['data'],header=ordered_header)       # write file
                fits[-1].write_checksum()                 # add checksum
            self.info('_saving: Exposure %d written to file %s' % (expid, filename))
            self.exp_ready.write({'expid' : expid, 'source' : self.role, 'image_name' : os.path.basename(filename), 'image_dir' : os.path.dirname(filename)})
        except Exception as e:
            rstring = '_saving: Exception writing FITS file: %s' % str(e)
            self.error(rstring)
            return 'FAILED: ' + rstring
        return self.SUCCESS

    def _detrend(self, exp, bias_subtract, dark_subtract):
        """
        perfrom image detrending
        """
        if (not bias_subtract and not dark_subtract) or self.bias_frame is None or self.dark_frame is None:
            image = {'primary' : {'header' : OrderedDict(exp['header']), 'data' : exp['data']}}
        elif dark_subtract == True:
            image = {'primary' : {'header' : OrderedDict(exp['header'])}}
            exptime = image['primary']['header']['EXPREQ']
            image['primary']['data'] = exp['data'] - exptime * self.dark_frame - self.bias_frame
        elif bias_subtract == True:
            image = {'primary' : {'header' : OrderedDict(exp['header'])}}
            exptime = image['primary']['header']['EXPREQ']
            image['primary']['data'] = exp['data'] - self.bias_frame
        else:
            image = {'primary' : {'header' : OrderedDict(exp['header']), 'data' : exp['data']}}
        return image

# Instance Connection Callbacks
    def about_to_connect_to_instance(self, *args, **kwargs):
        """ Instance management helper function """
        pass

    def did_connect_to_instance(self, *args, **kwargs):
        """ Instance management helper function """
        self.info('connected, setting up discovery stuff')
        discovery.reset()
        discovery.reset_discovered()
        self._setup_discovery(discovery.discoverable)

    def about_to_disconnect_from_instance(self, *args, **kwargs):
        """ Instance management helper function """
        pass

    def did_disconnect_from_instance(self, *args, **kwargs):
        """ Instance management helper function """
        self.info('disconnected, clearing discovery stuff')
        discovery.reset()
        discovery.reset_discovered()

    def _get_camera_status(self):
        if self.exposing.is_set():
            self.warn('get_camera_status: Blocked while exposing')
            return {'STATUS' : 'EXPOSING', 'last_updated' : datetime.datetime.utcnow().isoformat()}
        if self.config['controller_type'] == 'SBIG':
            t = self.client.get_camera_status()
        else:
            t = {}

        # only some keys contain valid information
        tt = {}
        if isinstance(t, dict):
            for key in ('cooling_enabled', 'tec_frozen', 'ccd_setpoint', 'imaging_ccd_temperature', 'imaging_ccd_power'):
                if key in t:
                    tt[key] = t[key]
            tt['status'] = 'SUCCESS'
        else:
            tt['status'] = 'ERROR'
        tt['last_updated'] = datetime.datetime.utcnow().isoformat()
        self.camera_status.write(tt)
        return tt

    def main(self):
        """ Run loop for FPC Application """
        # Wait for application to get configured
        self.info('main: waiting for configure before starting telemetry loop.')
        while not self.shutdown_event.is_set() and self.fpc_status._value != 'READY':
            self.sleep(1.0)

        self.info('main: starting telemetry loop')
        print_once = True
        while not self.shutdown_event.is_set():
            try:
                t = self._get_camera_status()
                print_once = True
            except Exception as e:
                if print_once:
                    self.error('main: Exception getting camera status: %s' % str(e))
                    print_once = False
            if isinstance(t, dict) and 'STATUS' in t and t['STATUS'] != 'SUCCESS':
                self.warn('main: Camera status not available: %r' % t)
                print_once = False
            time.sleep(max(self.status_update_rate, 2.0))

        while not self.shutdown_event.is_set():
            self.sleep(1)
        self.client.disconnectServer()


#######################################################
#                                                     #
# SBIG camera class                                   #
#                                                     #
#######################################################

class SBIG_camera():
    def __init__(self, model = 'STi'):
        self.cam = SBIGCam()
        self.cam.select_camera(model)
        self.model = model
        self.bins = (1,1)
        self.bias = False
        self.image = None
        self.model = None

    def configure(self):
        if self.model ==  'ST8300':
            self.setFrame(0, 0, 3352,2532)
        if self.model == 'STi':
            self.setFrame(0,0,648,484)

    def getFrame(self):
        return (int(self.cam.LEFT.value), int(self.cam.TOP.value), int(self.cam.WIDTH.value), int(self.cam.HEIGHT.value))

    def getBinning(self):
        return (1.0, 1.0)

    def setFrame(self, left, top, width, height):
        self.cam.set_image_size(width, height)
        self.cam.set_window_mode(top, left)

    def setFan(self, state):
        return self.cam.set_fan(state)

    def setBinning(self, a, b):
        self.bins = (a ,b)

    def setLight(self):
        self.cam.set_dark(False)
        self.bias = False

    def setDark(self):
        self.cam.set_dark(True)
        self.bias = False

    def setBias(self):
        self.cam.set_dark(True)
        self.bias = True

    def getBlob(self):
        return fits.PrimaryHDU(self.image)

    def getData(self):
        return self.image

    def get_camera_status(self):
        return self.cam.query_temperature_status()

    def set_temperature_regulation(self, reg, set):
        return self.cam.set_temperature_regulation(reg, set)

    def connectServer(self, *args, model = 'STi', **kwargs):
        self.model = model
        y = self.cam.open_camera()
        x = self.cam.select_camera(model)
        if model == 'STi':
            self.setFrame(0,0,648,484)
        return x and y

    def disconnectServer(self, *args, **kwargs):
        self.cam.close_camera()

    def setServer(self, *args, **kwargs):
        return True

    def getHost(self):
        return "localhost"

    def getPort(self):
        return -1

    def takeExposure(self, t):
        if self.bias:
            t = 0
        #SBIG driver expects exposure time in milliseconds
        self.cam.set_exposure_time(t*1000)
        self.image = self.cam.start_exposure()

#######################################################
#                                                     #
# FPC simulator class                                 #
#                                                     #
#######################################################

class FakeIndi():
    def __init__(self, display_image=False, count = 4):
        self.bins = (1.0, 1.0)
        self.frame = (0, 0, 3352, 2532)
        self.exptype = 'light'
        self.blob = None
        self.display_image = display_image
        self.count = count
        self.model = None

    def configure(self):
        if self.model ==  'ST8300':
            self.setFrame(0, 0, 3352,2532)
        if self.model == 'STi':
            self.setFrame(0,0,648,484)

    def getFrame(self):
        return self.frame

    def getBinning(self):
        return self.bins

    def setFrame(self, a, b, c, d):
        self.frame = (a, b, c, d)

    def setBinning(self, a, b):
        self.bins = (a ,b)

    def setLight(self):
        self.exptype = 'light'

    def setDark(self):
        self.exptype = 'dark'

    def setBias(self):
        self.exptype = 'bias'

    def getBlob(self):
        b = self.blob
        self.blob = None
        return b

    def getData(self):
        if self.blob != None:
            return self.blob[0].data

    def get_camera_status(self):
        return {}

    def set_temperature_regulation(self, reg, set):
        return True

    def connectServer(self, *args, model = 'STi', **kwargs):
        self.model = model
        if self.model ==  'ST8300':
            self.setFrame(0, 0, 3352,2532)
        if self.model == 'STi':
            self.setFrame(0,0,648,484)
        return True

    def disconnectServer(self, *args, **kwargs):
        return True

    def setServer(self, *args, **kwargs):
        return True

    def getHost(self):
        return "SIMULATED"

    def getPort(self):
        return -1

    def takeExposure(self, t, astype='uint16'):
        noise = 500.0
        dark = 5.0
        signal = 10.0
        sigvar = 5.0
        nSpots = self.count
        cov=[[50,0],[0,50]]
        noise_var = 5.0
        if(t != 0):
            time.sleep(t)
        hdulist = fits.HDUList()
        # bias
        n = noise + noise_var * np.random.rand(self.frame[2]-self.frame[0], self.frame[3]-self.frame[1])
        d = n
        # dark
        if self.exptype in ['dark', 'light']:
            i = dark * t + np.random.rand(self.frame[2]-self.frame[0], self.frame[3]-self.frame[1])
            d = d + i
        # objects
        if self.exptype in ['light']:

            for i in range(nSpots):
                mean_x = random.randint(self.frame[0], self.frame[2])
                mean_y = random.randint(self.frame[1], self.frame[3])
                x,y=np.random.multivariate_normal([mean_x, mean_y], cov, 100000).T
                star,xe,ye=np.histogram2d(x,y,bins=[self.frame[2]-self.frame[0],self.frame[3]-self.frame[1]],
                                        range=[[self.frame[0],self.frame[2]],[self.frame[1],self.frame[3]]])
                d += star * (signal + sigvar * random.random())
        if self.display_image:
            try:
                plt.close()
            except:
                pass
            plt.imshow(d, interpolation='nearest')
            plt.show(block=False)

        h = fits.header.Header()
        h.append(("SIMPLE", True, "file does conform to FITS standard "))
        h.append(("SIMPLE", True, "file does conform to FITS standard "))
        if astype == 'uint16':
            h.append(("BITPIX", 16, "number of bits per data pixel"))
        elif astype == 'float64':
            h.append(("BITPIX", -64, "number of bits per data pixel"))
        elif astype == 'float32':
            h.append(("BITPIX", -32, "number of bits per data pixel"))
        h.append(("NAXIS", 2, "number of data axes"))
        h.append(("NAXIS1", int(self.frame[2]-self.frame[0]), "length of data axis 1"))
        h.append(("NAXIS2", int(self.frame[3]-self.frame[1]), "length of data axis 2"))
        h.append(("EXTEND", True, "FITS dataset may contain extensions"))
        h.add_comment("FITS (Flexible Image Transport System) format is defined in 'Astronomy")
        h.add_comment("and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H")
        h.append(("INSTRUME", self.model, "Instrument Name"))
        h.append(("EXPTIME", t, "Total Exposure Time (s)"))
        h.append(("DARKTIME", t, "Total Exposure Time (s)"))
        h.append(("PIXSIZE1", 7.40, "Pixel Size 1 (microns)"))
        h.append(("PIXSIZE2", 7.40, "Pixel Size 2 (microns)"))
        h.append(("XBINNING", int(self.bins[0]), "Binning factor in width"))
        h.append(("YBINNING", int(self.bins[1]), "Binning factor in height"))
        h.append(("FRAME", str(self.exptype), "Frame Type"))
        h.append(("DATAMIN", float(d.min()), "Minimum Value"))
        h.append(("DATAMAX", float(d.max()), "Maximum Value"))
        h.append(("XBAYROFF", int(self.frame[0]), "X offset of Bayer array"))
        h.append(("YBAYROFF", int(self.frame[1]), "Y offset of Bayer array"))
        h.append(("BAYERPAT", "BGGR    ", "Bayer color pattern"))
        h.append(("DATE-OBS", datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), "UTC start date of observation"))
        h.add_comment("Simulated Data")
        phdu = fits.PrimaryHDU(data = d.astype(astype), header = h)
        hdulist.append(phdu)
        self.blob = hdulist

###################################################
if __name__ == '__main__':
    # Are device_mode or service set on the command line?
    if '--device_mode' in sys.argv:
        if '--service' in sys.argv:
            myFPC = FPCcam()
        else:
            myFPC = FPCcam(service='DOStest')
    else:
        if '--service' in sys.argv:
            myFPC = FPCcam(device_mode = True)
        else:
            myFPC = FPCcam(device_mode = True, service='DOStest')
    # Enter run loop
    myFPC.run()
