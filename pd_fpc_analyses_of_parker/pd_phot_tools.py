import os
import numpy as np
from pd_fpc_photometry_paf import pd_fpc_photometry_dict
from photometry_paf import segmentation_photometry
from photutils import properties_table
from astropy.io import fits
import pickle
import matplotlib.pyplot as plt




class PhotTools:
    def __init__(self, sequence_type, file_list, files, match,step_size, fiber_num, directory, name ):
        self.s_type = sequence_type #first_teldither, first_posdither, teldither, posdither, first_stability, stability, one
        self.dither_type = None
        self.dither_num = 5 #size of box
        self.dither_step_size = step_size  #arcsec or mm steps
        self.file_list = np.array(file_list)
        self.files = files #if list are actual file names or just a list of ids for filenames (True: filenames, False: IDs)
        self.img_dir = directory
        self.d = {}
        self.fiber_ids = {3002:(1201, 1373), 3001:(1198, 1059), 3003: (1856, 1059)} #fibre names
        self.DIST_THRESHOLD = 30
        #self.fibers_found = []
        self.fullwell = 60000
        self.name = name
        self.match_type = match
        self.max_fiber_num = fiber_num
    def list_of_files(self):
        """ Creates list of files from the file_list created and starts the dictionary
        """
        if self.files == False:
            print(self.file_list)
            for i,x in enumerate(self.file_list):
                #filename = "PROTODESI_%08d.fits" % x
                filename = 'PROTODESI_FPC_%08d.fits' % x
                self.d[x] = {}
                self.d[x]['path'] = self.img_dir+'/'+filename
                self.d[x]['dither_num'] = i
        elif self.files == True:
            for i, file in enumerate(self.file_list):
                x = os.path.splitext(file)[0][-8:]
                self.d[x] = {}
                self.d[x]['path'] = self.img_dir+'/'+file
                self.d[x]['dither_num'] = i
        else:
            print("Is the file list composed of full filenames (files=True)  or just IDs (files=False)?")

    def simple_photometry(self):
        """
        Takes list of FPC exposure images, makes a folder to hold the segmentation results,    and returns a list of results from the segmentation photometry.

        INPUT: list of exposure IDs that identify the FPC images
        OUTPUT: python dictionary
            KEY: exposure ID
            ITEM: list of information about the photometry
              - path file
              - area:
              - Sum: sum of flux in pixels
              - exptime: from FITS header (requested exposure time)
              - background: mean background of image
              - x: x_centroid pixel value
              - y: x_centroid pixel value
              - flux: sum /exptime
        """
        self.list_of_files()
        for i, values in self.d.items():
            path_file_abs = values['path']

            #Get info from fits header
            exptime = fits.open(path_file_abs)[0].header['EXPREQ']
            date = fits.open(path_file_abs)[0].header['DATE-OBS']
            time = fits.open(path_file_abs)[0].header['MJD-OBS']
            try:
                delta_ra = fits.open(path_file_abs)[0].header['DELTARA']
                delta_dec = fits.open(path_file_abs)[0].header['DELTADEC']
            except:
                delta_ra = None
                delta_dec = None
            self.d[i]['exptime'] = exptime
            self.d[i]['delta_ra'] = delta_ra
            self.d[i]['delta_dec'] = delta_dec
            self.d[i]['date'] = date
            self.d[i]['time'] = time


            #Do photometry
            [segm, props_list_return] = segmentation_photometry(path_file_abs,
                                                bkg_sigma = 3.0,
                                                source_snr = 3.0,
                                                fwhm_kernel = 2.0,
                                                x_size_kernel = 3,
                                                y_size_kernel = 3,
                                                clobber = False)

            #Identify individual fibers
            # there may be more than 3 sources returned, depending on parameters
            # for each of the fibres, look for exactly one match in sources

            fibers_found = []
            if self.match_type == 'location':
                for fiber_num, location in self.fiber_ids.items():
                    for props in props_list_return:

                    # check distance to known fibre position
                    # prop.centroid is in (y, x) format
                        dist = np.linalg.norm(
                                np.fliplr([props.centroid.value])[0]
                                - location)
                        if dist < self.DIST_THRESHOLD:
                            fibers_found.append(fiber_num)
                            print("Fiber %d was identified" % fiber_num)
                            self.d[i][fiber_num] = {}
                            self.d[i][fiber_num]['seg_object'] = props
                            # remove the props of identified source from props list
                            #props_list_return.remove(props)
                            # stop searching for current fibre
                            #break
                for fiber_num in fibers_found:
                    props = self.d[i][fiber_num]['seg_object']
                    data = np.array(properties_table(props,columns=('area','source_sum','background_mean','centroid','max_value')))

                    self.d[i][fiber_num]['area'] = data['area']
                    self.d[i][fiber_num]['sum'] = data['source_sum']
                    self.d[i][fiber_num]['background'] = data['background_mean']
                    self.d[i][fiber_num]['centroid'] = data['centroid']
                    self.d[i][fiber_num]['flux'] = data['source_sum']/exptime
                    self.d[i][fiber_num]['max_value'] = data['max_value']

                    #else:
                    #    # this is not the fibre propwe are looking for
                    #    print('Cannot match fiber %d' % fiber_num)


            elif self.match_type == 'max':
                if self.max_fiber_num == 'all':
                    print("Can only run max for single fibers")
                else:
                    try:
                        max_fiber = int(self.max_fiber_num)
                        self.d[i][max_fiber] = {}
                    except:
                        print("fiber number has to be 3001, 3002, 3003")
                try:
                    all_data = np.array(properties_table(props_list_return,columns=('area')))
                    area = all_data['area']
                    x = np.where(area == max(area))
                    if len(x[0]) > 1:
                        print("Couldn't find a maximum")
                    elif len(x[0]) == 1:
                        prop = props_list_return[x[0]]
                        self.d[i][max_fiber]['seg_object'] = prop
                        #Put info in dictionary
                        props = self.d[i][max_fiber]['seg_object']
                        data = np.array(properties_table(props,columns=('area','source_sum','background_mean','centroid','max_value')))

                        self.d[i][max_fiber]['area'] = data['area']
                        self.d[i][max_fiber]['sum'] = data['source_sum']
                        self.d[i][max_fiber]['background'] = data['background_mean']
                        self.d[i][max_fiber]['centroid'] = data['centroid']
                        self.d[i][max_fiber]['flux'] = data['source_sum']/exptime
                        self.d[i][max_fiber]['max_value'] = data['max_value']
                except:
                    print("no area to speak of")




            #Put all relevant info for each fiber in dictionary




    def save_pickle(self):
        phot_dict = self.d
        for fiber in [3001,3002,3003]:
            for key, values in self.d.items():
                try:
                    del values[fiber]['seg_object']
                except:
                    pass
        #pickle_file = "%s.p" % str(self.name)
        pickle.dump(phot_dict,open("%s.p" % self.name,"wb"))

    def load_pickle(self):
        data_dict = pickle.load(open("%s.p" % self.name, "rb"))
        self.d = data_dict

    def plot_dither(self):
        print(self.dither_type)
        if self.dither_num == 5:
            #dither_array = np.array([[21,20,19,28,17],[22,7,6,5,16],[23,8,1,4,15],[24,9,2,3,14],[24,10,11,12,13]])
            pos_dither_array = np.array([[19,18,17,16,15],[20,6,5,4,14],[21,7,0,3,13],[22,8,1,2,12],[23,24,9,10,11]])
            tel_dither_array = np.array([[11,10,9,24,23],[12,2,1,8,22],[13,3,0,7,21],[14,4,5,6,20],[15,16,17,18,19]])
        new_grids = {}
        name_grids = {}
        for fiber in [3001,3002,3003]: #3002,3003
            new_grids[fiber] = np.zeros([self.dither_num,self.dither_num])
            name_grids[fiber] = np.zeros([self.dither_num,self.dither_num])
            for key, values in self.d.items():
                if self.dither_type == 'telescope':
                    loc = np.where(self.d[key]['dither_num']==tel_dither_array)
                elif self.dither_type == 'positioner':
                    loc = np.where(self.d[key]['dither_num']==pos_dither_array)
                name_grids[fiber][loc] = int(key)
                try:
                    #print(key,self.d[key]['dither_num'],loc)
                    new_grids[fiber][loc] = self.d[key][fiber]['flux']
                    print(self.d[key][fiber]['flux'])
                except:
                    new_grids[fiber][loc] = 0


        new_directory = os.getcwd()+'/'+self.name
        if not os.path.exists(new_directory):
            os.makedirs(new_directory)
        for key, value in new_grids.items():
            names = name_grids[key]
            f,ax = plt.subplots()
            # half = (((self.dither_num-1)/2.)*self.dither_step_size)
            V = np.ma.masked_where(value < 1., value)
            cmap = plt.cm.afmhot
            cmap.set_bad(color='green')
            cls = ax.matshow(V,cmap=cmap)
            for (j,i),label in np.ndenumerate(names):
                ax.text(i,j,int(label),ha='center',va='center')

            ax.set_xlabel('Relative RA offset (arcsec)')
            ax.set_ylabel("Relative DEC offset (arcsec)")
            ax.set_title("Flux for %s pattern: Fiber# %d" % (self.s_type,key))
            ax.set_xticklabels([-3*self.dither_step_size,-2*self.dither_step_size,-self.dither_step_size,0,self.dither_step_size,2*self.dither_step_size])
            ax.set_yticklabels([3*self.dither_step_size,2*self.dither_step_size,self.dither_step_size,0,-self.dither_step_size,-2*self.dither_step_size])
            ax.xaxis.set_ticks_position('bottom')
            f.colorbar(cls)
            plt.savefig("%s/%s_%d.png"% (new_directory,self.name,key))
        plt.show()

    def plot_centroids(self):
        for fiber in self.fibers_found:
            X = []
            Y = []
            for key, values in self.d.items():
                cent = values[fiber]['centroid']
                Y.append(cent[0][0])
                X.append(cent[0][1])
            plt.plot(X,Y,'x')
            plt.xlabel("X pixel value")
            plt.ylabel("Flux")
            plt.title("Centroids of a given fiber" % (fiber,self.name))
            plt.savefig("%s_centroids.png" % self.name)
        #plt.show()

    def plot_stability(self):
        for fiber in self.fibers_found:
            times = []
            fluxes = []
            for key, values in self.d.items():
                times.append(values['time'])
                fluxes.append(values[fiber]['flux'])
            s = np.argsort(np.array(times))
            times = np.array(times)[s]
            fluxes = np.array(fluxes)[s]
            new_directory = os.getcwd()+'/stability/'+self.name
            if not os.path.exists(new_directory):
                os.makedirs(new_directory)
            plt.plot(times,fluxes,'x')
            plt.xlabel("Time image taken (MJD)")
            plt.ylabel("Flux")
            plt.title("Flux stability on fiber %d during guiding - %s" % (fiber,self.name))
            plt.savefig("%s_stability.png" % self.name)
        #plt.show()


    def simple_output(self):
        """ Prints the flux of each fiber identified
        """
        for id,values in self.d.items():
            print("Image %s photometry results:" % id)
            for fiber in [3001,3002,3003]:
                try:
                    fiber_flux = values[fiber]['flux']
                    print("Fiber %d: %f" % (fiber,fiber_flux))
                    max_value = values[fiber]['max_value']
                    if max_value >= self.fullwell:
                        print("Fiber %d is overexposed in this image" % fiber)
                    else:
                        print("This fiber doesn't seem to be overexposed!")
                except:
                    print("No data found for fiber %d" % fiber)

    def run(self):
        if self.s_type == 'first_teldither':
            self.dither_type = 'telescope'
            self.simple_photometry()
            self.save_pickle()
            self.plot_dither()
        if self.s_type == 'first_posdither':
            self.dither_type = 'positioner'
            self.simple_photometry()
            self.save_pickle()
            self.plot_dither()
        if self.s_type == 'teldither':
            self.dither_type = 'telescope'
            self.load_pickle()
            self.plot_dither()
        if self.s_type == 'posdither':
            self.dither_type='positioner'
            self.load_pickle()
            self.plot_dither()
        if self.s_type == 'first_stability':
            self.simple_photometry()
            self.save_pickle()
            self.plot_stability()
        if self.s_type == 'stability':
            self.load_pickle()
            self.plot_stability()
        if self.s_type == 'one':
            self.simple_photometry()
            self.simple_output()
        if self.s_type == 'simple':
            self.load_pickle()
            self.simple_output()






