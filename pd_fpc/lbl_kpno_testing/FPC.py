#!/usr/bin/env python
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
"""
from DOSlib.application import Application
import DOSlib.discovery as discovery
import DOSlib.mjd
import DOSlib.util
import DOSlib
import datetime
import sys
import os
import math
import time
import imp
from astropy.io import fits
import numpy as np
from collections import OrderedDict
from sbigcam import SBIGCam
try:
    import indiclient
except:
    print ('FPC: indi support not available')

class FPCcam(Application):
    commands = ['configure',
                'get', 
                'set', 
                'expose',
                'get_image',
                'write_image',
                'send_image',
                ]
    defaults = {'fpc_host' : 'decampi21',
                'fpc_port' : 7624,
                'image_dir' : '.',
                'image_name' : 'protodesi',
                'controller_type' : 'SBIG',#'SIMULATOR'
                'software' : 'sbigcam',
                'local_copy' : False,
                'auto_purge' : False,
                }

    # list of keywords NOT to include in FITS header
    fpc_headers = ['local_copy', 'expid', 'image_dir', 'image_name', 'auto_purge', 'cb_name', 
                   'cb_function', 'cb_args', 'send_name', 'send_command', 'send_args']
        
    def init(self):
        """ Initialization """
        self.info('Initializing')

        # FPC Status
        self.fpc_status = self.shared_variable('STATUS')  
        self.fpc_status.publish()
        self.cam_connected = False

        # Update user information
        self.add_interlock_information(interlock = 'DOS', key ='FPC_STATUS',
                                       set_condition=['NOTCONNECTED', 'READY'],
                                       enabled=True)

        # configuration information
        self.image_dir = self.config['image_dir']
        if self.image_dir == '.':
            self.image_dir = os.getcwd()
        if not os.access(self.image_dir,os.W_OK):
            self.error('init: image directory not found or insufficient permissions to write image files %s' % str(self.image_dir))
            self.image_dir = None
        self.image_name = self.config['image_name']
        self.local_copy = True if 'T' in str(self.config['local_copy']).upper() else False
        self.auto_purge = True if 'T' in str(self.config['auto_purge']).upper() else False
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
                self.client = FakeIndi()
            elif self.config['software'] == 'indi':
                try:
                    import indiclient
                except Exception as e:
                    self.error('init: Exception loading indi software: %s' % str(e))
                    sys.exit()
                self.client = indiclient.IndiClient()
            else:
                self.client = SBIG_camera()
        except Exception as e:
            self.error('init: Exception accessing camera: %s' % str(e))
            sys.exit()

        # Connect to the camera (try 10 times)
        for x in range(10):
           self.client.setServer(self.config['fpc_host'], self.config['fpc_port'])
           if (not(self.client.connectServer())):
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
        
        if not self.cam_connected:
            self.info('Trying to connect to FPC (%s, %s)' % (str(self.config['fpc_host']), str(self.config['fpc_port'])))
            # Connect to the camera (try 10 times)
            for x in range(10):
                self.client.setServer(self.config['fpc_host'], self.config['fpc_port'])
                if (not(self.client.connectServer())):
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
                rstring = 'configure: Failed to connect to camera.'
                self.error(rstring)
                self.fpc_status.write('ERROR')
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
        status, info, image_dir, image_name, connected,
        region, binning, auto_purge, local_copy, exposure
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
        elif(params.lower() in self.config.keys()):
            return self.config[params]
        else:
            return 'FAILED: Invalid parameter for get command.'

    def set(self, *args, **kwargs):
        """
        Set application configuration parameters including:
        image_dir, image_name, region, binning, local_copy, auto_purge,
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
            return self.SUCCESS
        else:
            return 'FAILED: Invalid arguments for set command'
        return self.FAILED

    def expose(self, expid, *args, **kwargs):
        """
        Take an image  with the FPC.
        Optionally write image in FITS format to file.
        Positional arguments
           expid           exposure number (required)           
        Named arguments:
           exptype         object, zero, dark
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

        # setup exposure management study
        local_copy = self.config['local_copy'] if 'local_copy' not in kw else kw['local_copy']
        auto_purge = self.config['auto_purge'] if 'auto_purge' not in kw else kw['auto_purge']
        # Callback setup        
        cb_name =  None if 'cb_name' not in kw else kw['cb_name']
        cb_function =  None if 'cb_function' not in kw else kw['cb_function']
        cb_args =  None if 'cb_args' not in kw else kw['cb_args']
        send_name = None if 'send_name' not in kw else kw['send_name']
        send_function = None if 'send_function' not in kw else kw['send_function']
        send_args = None if 'send_args' not in kw else kw['send_args']

        if 'exptime' not in kw:
            rstring = 'expose requires an exposure time'
            self.error(rstring)
            return 'FAILED: ' + rstring
        else:
            exptime = float(kw['exptime'])
        if 'exptype' not in kw or kw['exptype'] not in ['object', 'dark', 'zero', 'bias']:
            rstring = 'expose requires an exposure type'
            self.error(rstring)
            return 'FAILED: ' + rstring
        else:
            exptype = str(kw['exptype'])
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
        self.exposures['send_function'] = None if 'send_function' not in kw else kw['send_function']
        self.exposures['send_args'] = None if 'send_args' not in kw else kw['send_args']

        # setup for header and data
        self.info("expose: Exposure %d, exptype %s, exptime %f" % (expid,exptype, exptime))
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
            
        try:
            if exptype == 'object':
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
            return 'FAILED: ' + rstring

        # Retrieve data (getData for sbig returns the image, getBlob returns a fits HDUlist)
        try:
            if self.config['software'] == 'indi':
                blob = self.client.getBlob()
            else:
                blob = self.client.getData()
            i = 0
            while(blob == None):
                self.sleep(0.5)
                if self.config['software'] == 'indi':
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
            return 'FAILED: ' + rstring
        self.info("expose: data received for exposure %d." % expid)
        if(self.config["controller_type"] == "SIMULATOR" or self.config["software"] != 'indi'):
            self.exposures['data'] = blob
        else:
            hdulist = fits.HDUList.fromstring(blob)
            self.exposures['data'] = hdulist[0].data

        # Are we supposed to do something with the image:
        if local_copy:
            retcode = self._saving(auto_purge, image_dir, filename)
            if 'FAILED' in retcode:
                self.error('expose: _savings returns: %s' % str(retcode))
        if self.exposures['send_name'] != None and self.exposures['send_function'] != None:
            extras = {'send_args' : send_args, 'source' : self.role, 'exposure_status' : self.exposures['status'], 'expid' : expid}
            retcode = self.send_image(expid, send_name, send_command, **extras)
            self.exposures['status'] = 'SENT'
        # Is there a callback?
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
        # We are done
        return self.SUCCESS

    def send_image(self, *args, **kwargs):
        """
        Send exposure expid to an PML interface (name or proxy) using pml_command
        """
        try:
            a,kw = dos_parser(*args, **kwargs)
        except:
            return 'FAILED: invalid arguments for send_image command'
        if len(a)>0:
            expid = a[0]
        else:
            if 'expid' in kw:
                expid = kw['expid']
            else:
                return 'FAILED: exposure number missing'
        if expid not in self.exposure_status:
            return 'FAILED: send_image: Exposure %s is not available' % str(expid)
        if len(a)>1:
            pml_name = a[1]
        else:
            if 'pml_name' in kw:
                pml_name = kw['pml_name']
            else:
                return 'FAILED: PML name missing'
        if len(a)>2:
            pml_command = a[2]
        else:
            if 'pml_command' in kw:
                pml_command = kw['pml_command']
            else:
                return 'FAILED: PML command missing'

        header = OrderedDict(self.exposures['header'])
        if 'extra_headers' in kw and isinstance(kw['extra_headers'], dict):
            for key in list(kw['extra_headers'].keys()):
                header[str(key).strip().upper()[:8]] = kw['extra_headers'][key]    # Careful - overwrites existing values if same keys

        exposure = {'header' : header, 'data' : self.exposures['data']}
        try:
            self.info('send_image: sending exposure %d to %s' % (expid, str(pml_name)))
            # We are setting up a new connection everytime - maybe this can be done more efficiently

            if str(pml_name).upper().startswith('PYRO'):
                # Use Pyro4.Proxy directly
                retcode = Pyro4.Proxe(pml_name).__getattr__(pml_command)(**exposure)
            else:
                retcode = dos_connection(pml_name).execute(pml_command, **exposure)
            if 'FAILED' in retcode:
                self.error('send_image: sending exposure %d returns: %s' % (expid, str(retcode)))
            self.info('send_image: Done')
            return retcode
        except Exception as e:
            rstring = 'send_image: Exception calling %s: %s. Arguments: %s, %s' % (str(pml_name), str(msg), repr(args), repr(kw))
            self.error(rstring)
            return 'FAILED: ' + rstring
        
    def get_image(self, expid, extra_headers = None):
        """ 
        Return an image
        """
        if expid in ['last', -1]:
            if self.exposures['status'] not in ['EXPOSED', 'SENT']:
                return 'FAILED: no image available.'
            expid = self.exposures['expid']
        elif expid != self.exposures['expid']:
            raise RuntimeError('Invalid exposure number. Image data not available.')
        
        header = {}
        if isinstance(extra_headers, dict):
            for key in list(extra_headers.keys()):
                header[str(key).strip().upper()[:8]] = extra_headers[key]    # Careful - overwrites existing values if same keys
        # get the rest
        ordered_header = OrderedDict(self.exposures['header'])
        ordered_header.update(header)      
        self.info('Returning exposure %s' % str(expid))
        return {'header' : ordered_header, 'data' : self.exposures['data']}
    
    def write_image(self, expid, image_dir = None, image_name = None, extra_headers = None, overwrite = True):
        """
        Save image to disk
        Use image_name/image_dir to overwrite current values
        add extra headers to fits file
        """
        if expid in ['last', -1]:
            if self.exposures['status'] not in ['EXPOSED', 'SENT']:
                return 'FAILED: no image available.'
            expid = self.exposures['expid']
        elif expid != self.exposures['expid']:
            raise RuntimeError('Invalid exposure number. Image data not available.')
        
        dir = image_dir if image_dir != None else self.image_dir
        name = image_name if image_name != None else self.image_name
        complete_name =name + '_%s_%08d.fits' % (self.role, expid)
        filename = os.path.join(dir, complete_name)
        self.info('write_image: using output file %s' % filename)
        # call io routine
        try:
            if overwrite and os.path.isfile(filename):   # does the file already exist?
                os.remove(filename)
            hdu = fits.PrimaryHDU(self.exposures['data'])
            for key in self.exposures['header']:
                hdu.header.append(key, self.exposures['header'][key])
            if isinstance(extra_headers, dict):
                for key in list(extra_headers.keys()):
                    hdu.header.append(str(key).strip().upper()[:8], extra_headers[key])    # Careful - overwrites existing values if same keys
            hdu.writeto(filename, clobber = overwrite)
        except Exception as e:
            rstring = 'write_image: Exception writing FITS file: %s' % str(e)
            self.error(rstring)
            return "FAILED: " + rstring
        self.info('write_image: Done')
        return self.SUCCESS        

    def _saving(self, auto_purge, image_dir, filename):
        """
        Save image data to local disk
        Clean disk space if auto_purge is set
        """
        
        # Remove old images if auto_purge is set
        if self.auto_purge:
            fits_files = sorted(glob.glob(os.path.join(image_dir,'*.fits')), key = self._modified)
            if len(fits_files) > 1:
                # newest file is at the end of the list
                for index in range(len(fits_files) - 1):
                    try:
                        to_remove = os.path.join(image_dir, fits_files[index])
                        os.remove(to_remove)
                        self.info('_saving: removed image file %s' % to_remove)
                    except Exception as e:
                        self.error('_saving: Exception removing old fits files: %s' % str(e))
        hdu = fits.PrimaryHDU(self.exposures['data'])
        for key in self.exposures['header']:
            hdu.header.append(key, self.exposures['header'][key])
        try:
            hdu.writeto(filename, clobber = True)
        except Exception as e:
            rstring = '_saving: Exception writing FITS file: %s' % str(e)
            self.error(rstring)
            return "FAILED: " + rstring
        return self.SUCCESS

# Instance Connection Callbacks
    def about_to_connect_to_instance(self, *args, **kwargs):
        pass

    def did_connect_to_instance(self, *args, **kwargs):
        self.info('connected, setting up discovery stuff')
        discovery.reset()
        discovery.reset_discovered()
        self._setup_discovery(discovery.discoverable)

    def about_to_disconnect_from_instance(self, *args, **kwargs):
        pass

    def did_disconnect_from_instance(self, *args, **kwargs):
        self.info('disconnected, clearing discovery stuff')
        discovery.reset()
        discovery.reset_discovered()

    def main(self):
        while not self.shutdown_event.is_set():
            self.sleep(1)
        self.client.disconnectServer()

        
#######################################################
#                                                     #
# SBIG camera class                                   #
#                                                     #
#######################################################

class SBIG_camera():
    def __init__(self):
        self.cam = SBIGCam()
        self.bins = (1,1)
        self.bias = False
        self.image = None

    def getFrame(self):
        return (int(self.cam.LEFT.value), int(self.cam.TOP.value), int(self.cam.WIDTH.value), int(self.cam.HEIGHT.value))

    def getBinning(self):
        return (1.0, 1.0)

    def setFrame(self, left, top, width, height):
        self.cam.set_image_size(width, height)
        self.cam.set_window_mode(top, left)

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
    
    def connectServer(self, *args, model = 'ST8300', **kwargs):
        y = self.cam.open_camera()
        x = self.cam.select_camera(model)
        if model == 'ST8300':
            self.setFrame(0,0,3352,2532)
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
        self.cam.set_exposure_time(t*1000)
        self.image = self.cam.start_exposure()

#######################################################
#                                                     #
# FPC simulator class                                 #
#                                                     #
#######################################################

class FakeIndi():
    def __init__(self):
        self.bins = (1.0, 1.0)
        self.frame = (0.0, 0.0, 3362.0, 2504.0)
        self.exptype = 'light'
        self.blob = None

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
        return self.data
    
    def connectServer(self, *args, **kwargs):
        return True

    def disconnectServer(self, *args, **kwargs):
        return True

    def setServer(self, *args, **kwargs):
        return True

    def getHost(self):
        return "SIMULATED"

    def getPort(self):
        return -1

    def takeExposure(self, t):
        if(t != 0):
            time.sleep(t)
        hdulist = fits.HDUList()
        self.d = np.random.rand(self.frame[2], self.frame[3])
        h = fits.header.Header()
        h.append(("SIMPLE", True, "file does conform to FITS standard "))
        h.append(("BITPIX", -32, "number of bits per data pixel"))
        h.append(("NAXIS", 2, "number of data axes"))
        h.append(("NAXIS1", int(self.frame[2]), "length of data axis 1"))
        h.append(("NAXIS2", int(self.frame[3]), "length of data axis 2"))
        h.append(("EXTEND", True, "FITS dataset may contain extensions"))
        h.add_comment("FITS (Flexible Image Transport System) format is defined in 'Astronomy")
        h.add_comment("and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H")
        h.append(("INSTRUME", "'SBIG ST-i Color CCD Camera'", "CCD Name"))
        h.append(("EXPTIME", 1, "Total Exposure Time (s)"))
        h.append(("DARKTIME", 1, "Total Exposure Time (s)"))
        h.append(("PIXSIZE1", 7.40000009536743, "Pixel Size 1 (microns)"))
        h.append(("PIXSIZE2", 7.40000009536743, "Pixel Size 2 (microns)"))
        h.append(("XBINNING", int(self.bins[0]), "Binning factor in width"))
        h.append(("YBINNING", int(self.bins[1]), "Binning factor in height"))
        h.append(("FRAME", str(self.exptype), "Frame Type"))
        h.append(("DATAMIN", float(d.min()), "Minimum Value"))
        h.append(("DATAMAX", float(d.max()), "Maximum Value"))
        h.append(("XBAYROFF", int(self.frame[0]), "X offset of Bayer array"))
        h.append(("YBAYROFF", int(self.frame[1]), "Y offset of Bayer array"))
        h.append(("BAYERPAT", "BGGR    ", "Bayer color pattern"))
        h.append(("DATE-OBS", datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"), "UTC start date of observation"))
        h.add_comment("Generated by INDI (simulated)")
        phdu = fits.PrimaryHDU(data = self.d, header = h)
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
