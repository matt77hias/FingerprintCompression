'''
Configuration
Contains all global parameters and directories
used and accessed in the other .py files
@author     Matthias Moulin & Vincent Peeters
@version    1.0
'''
import os

dir_fingerprints = "data/"  

def get_dir_prefix():
    if (os.environ.get("USERNAME") == "Matthias"):
        return "C:/Users/Matthias/Desktop/Bir/Ma2.1/Wavelets/Project/FingerprintCompression/"
    else:
        return "/Users/vincentpeeters/Wiskundige Ingenieurs/Wavelets/Git/FingerprintCompression/"
        #return "~/Wiskundige Ingenieurs/Wavelets/Git/FingerprintCompression/"

def get_dir_fingerprints():
    return get_dir_prefix() + dir_fingerprints