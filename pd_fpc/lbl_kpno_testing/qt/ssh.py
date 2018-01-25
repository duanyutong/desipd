# -*- coding: utf-8 -*-
"""
SSH utility

Created on Thu Jul 28 09:42:33 2016

@author: givoltage
"""

import paramiko
import os
import stat
#import sys

# define directories
path_remote = r'/home/msdos/SBIG/'
path_local = r'C:\Users\givoltage\Downloads\SBIG'

def paramiko_sftp_progress(transferred, total):
    percentage = transferred/total*100
    total_mb = total/(1024**2)
    print('{0:.2f}%, Total {1:.2f} MB \r'.format(percentage, total_mb))
#    sys.stdout.write('{0:.2f}%, Total {1:.2f} MB \r'.format(percentage, total_mb))
#    sys.stdout.flush()
    
def absolute_full_paths(directory):
    absfullpaths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            fullpath = os.path.abspath(os.path.join(root, file))
            absfullpaths.append(fullpath)
    return absfullpaths  
    
def compare_sets(sftp, dir_remote, dir_local):
    '''
    Compare remote and local files as sets before syncing
    
    '''
    set_remote = set(sftp.listdir(dir_remote))
    list_remote = list(set_remote)
    list_remote_converted = list(list_remote)
    
    for i in range(len(list_remote)):
        list_remote_converted[i] = list_remote_converted[i].replace(':', '_')
    set_remote_converted = set(list_remote_converted)
    
    set_local = set(os.listdir(dir_local))
    set_transfer_converted = \
                    set_remote_converted - (set_remote_converted & set_local)
    set_transfer = set(set_transfer_converted)
    
    for filename in set_transferConverted:
        i = list_remoteConverted.index(filename)
        set_transfer.remove(list_remoteConverted[i])
        set_transfer.add(list_remote[i])
        
    return set_transfer
     
ssh = paramiko.SSHClient()

ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
ssh.connect(
    'desisti.kpno.noao.edu', 
    username='msdos', 
    password='MS-d0s')
    
[stdin, stdout, stderr] = ssh.exec_command("uptime")
stdout.readlines()

# open sftp
# filename = 'sti_test_01_21:17:33.962414.fits'
sftp = ssh.open_sftp()

# list all remote files
#setTransfer = compare_sets(sftp)
#            
# start transfer
def isdir(path):
    try:
        return stat.S_ISDIR(sftp.stat(path).st_mode)
    except IOError:
        #Path does not exist, so by definition not a directory
        return False
        
def sftp_get_all(sftp, path_remote, path_local):
    
    if isdir(path_remote):
        # if remote path is a directory    
        # check that remote path ends with / to ensure compatibility
        # with linux file server
        if path_remote[-1] != '/':
            path_remote = path_remote+'/'
        # create lists of names and full paths
        names = sftp.listdir(path_remote)
        paths = []
        for name in names:
            paths.append(os.path.join(path_remote, name))
        # copy each remote path in paths
        for i in range(len(paths)):
            name = names[i]
            name_local = name.replace(':', '_')
            path = paths[i]
            destination = os.path.join(path_local, name_local)
            print('Processing {}'.format(path))
            if isdir(path):
                print('{} is a directory'.format(name))
                # create local directory
                if not os.path.exists(destination):	# check and create
                    os.makedirs(destination)
            sftp_get_all(sftp, path, destination)
    else:
        # remote path is a file
        print('Transferring file {} ...'.format(path_remote))
        sftp.get(path_remote, path_local,
            callback = paramiko_sftp_progress
            )
            
            


filenameLocal = filename.replace(':', '_')
            
sftp_get_all(sftp, dir_remote, dir_local)
        
#        
#def sftp_transfer_set(set_transfer):
#    
#    for item in setTransfer:
#        
#        # determine whether file or folder
#        
#        print('Transferring {} ...'.format(filename))
#        filenameLocal = filename.replace(':', '_')
#        try:
#            sftp.get(
#                os.path.join(datasetDirRemote, filename), 
#                os.path.join(datasetDirLocal, filenameLocal),
#                callback = paramiko_sftp_progress
#                )
#        except:
#            print('Skipping file: ' + filename)
#            pass
#        print('{} transfer is complete.'.format(filename))

