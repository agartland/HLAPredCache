#! /usr/bin/python

# This script will set appropriate directory settings because some methods such as 'ANN/NetMHC' need this.

import os
import re
from subprocess import call

class Configure(object):
    def __init__(self):
        pass

    def generate_setupinfo(self):
        mainpath = os.path.abspath('.')
        methodpath = os.path.join(mainpath,'mainpath')
        srcpath = os.path.join(mainpath,'src')

        f=open(srcpath+'/setupinfo.template', 'r')
        content = f.read()
        contentb = content % (mainpath,)
        f.close()

        f=open(srcpath + '/setupinfo.py', 'w')
        f.write(contentb)
        f.close()

    def ann(self):
        mainpath      = os.path.abspath('.')
        methodpath    = os.path.join(mainpath,'method')
        netmhcpath = os.path.join(methodpath,'netMHC-3.4')

        fname     = os.path.join(netmhcpath, 'netMHC.template')
        fname_out = os.path.join(netmhcpath, 'netMHC')

        f=open(fname, 'r')
        content = f.read()
        content = re.sub('prediction_directory', netmhcpath, content)
        contentb = content
        f.close()

        print 'netMHC aka ann script:', fname_out
        f=open(fname_out, 'w')
        f.write(contentb)
        f.close()

        os.chmod(fname_out, 0755)


    def netmhcpan(self):
        mainpath      = os.path.abspath('.')
        methodpath    = os.path.join(mainpath,'method')
        netmhcpanpath = os.path.join(methodpath,'netMHCpan-2.8')
        scratchpath   = os.path.join(netmhcpanpath,'scratch')
        
        if (os.path.exists(scratchpath) == False):
            os.mkdir(scratchpath)
        os.chmod(scratchpath, 0777) # So the all can write in this directory.

        fname     = os.path.join(netmhcpanpath, 'netMHCpan.template')
        fname_out = os.path.join(netmhcpanpath, 'netMHCpan')

        f=open(fname, 'r')
        content = f.read()
        contentb = content % (netmhcpanpath, scratchpath)
        f.close()

        print 'netmhcpan script:', fname_out
        f=open(fname_out, 'w')
        f.write(contentb)
        f.close()

        os.chmod(fname_out, 0755)
        
        
    def pickpocket(self):
        mainpath      = os.path.abspath('.')
        methodpath    = os.path.join(mainpath,'method')
        pickpocketpath = os.path.join(methodpath,'pickpocket-1.1')
        scratchpath   = os.path.join(pickpocketpath,'scratch')
        
        if (os.path.exists(scratchpath) == False):
            os.mkdir(scratchpath)
        os.chmod(scratchpath, 0777) # So the all can write in this directory.

        fname     = os.path.join(pickpocketpath, 'PickPocket.template')
        fname_out = os.path.join(pickpocketpath, 'PickPocket')

        f=open(fname, 'r')
        content = f.read()
        contentb = content % (pickpocketpath, scratchpath)
        f.close()

        print 'pickpocket script:', fname_out
        f=open(fname_out, 'w')
        f.write(contentb)
        f.close()

        os.chmod(fname_out, 0755)
    
    def netmhccons(self):
        mainpath      = os.path.abspath('.')
        methodpath    = os.path.join(mainpath,'method')
        netmhcconspath = os.path.join(methodpath,'netMHCcons-1.1')
        netmhcpanpath = os.path.join(methodpath,'netMHCpan-2.8/netMHCpan')
        netmhcpath = os.path.join(methodpath,'netMHC-3.4/netMHC')
        pickpocketpath = os.path.join(methodpath,'pickpocket-1.1/PickPocket')
        scratchpath   = os.path.join(netmhcconspath,'scratch')
        
        if (os.path.exists(scratchpath) == False):
            os.mkdir(scratchpath)
        os.chmod(scratchpath, 0777) # So the all can write in this directory.

        fname     = os.path.join(netmhcconspath, 'netMHCcons.template')
        fname_out = os.path.join(netmhcconspath, 'netMHCcons')

        f=open(fname, 'r')
        content = f.read()
        contentb = content % (netmhcconspath, netmhcpanpath, netmhcpath, pickpocketpath, scratchpath)
        f.close()

        print 'netmhccons script:', fname_out
        f=open(fname_out, 'w')
        f.write(contentb)
        f.close()

        os.chmod(fname_out, 0755)
        
if __name__ == '__main__':
    config = Configure()
    config.generate_setupinfo()
    config.ann()
    config.netmhcpan()
    config.pickpocket()
    config.netmhccons()


