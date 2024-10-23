
import numpy as np

class parameters:

    def __init__(self, dz=1.0, wl=1.0, pw=5e-5, nx=100, qmaxin=1.0):
        self.dz = dz
        self.wl = wl
        self.pw = pw 
        self.nx = nx
        self.qmaxin = qmaxin
        self.calc_qpoints()

    def calc_qpoints(self):
        self.qpoints = np.arange(self.nx)/self.nx
        self.angle = np.arctan( (self.nx)*self.pw/self.dz)
        #self.qmax = 2.0*np.pi*(2.0/self.wl)*np.sin(self.angle/2.0)
        #self.qmax = (2.0/self.wl)*np.sin(self.angle/2.0)
        self.qmax = 2.0*np.pi*self.qmaxin
        self.qpoints *= self.qmax #*1e-10

def extract_required_parameters(pfile):

    sp = parameters()

    with open(pfile,'r') as f:
        for line in f:
            bits = line.split(" ")
            if bits[0]=='z':
                sp.dz = float(bits[2])
            elif bits[0]=='wl':
                sp.wl = float(bits[2])
            elif bits[0]=='pw':
                sp.pw = float(bits[2])          
            elif bits[0]=='nx':
                sp.nx = int(bits[2])          
            elif bits[0]=='saxs_nr':
                sp.saxs_nr = int(bits[2])          
            elif bits[0]=='qmax':
                sp.qmaxin = float(bits[2])          
            elif bits[0]=='elements':
                sp.elements = bits[2:]          
    sp.calc_qpoints()
    return sp

def gaussian(nr, wid):
    r = np.arange(nr) - nr/2
    g = np.exp(- r*r/(2*wid*wid))
    g = np.roll(g, nr//2) 
    g *= 1.0/np.sum(g)
    return g

def convolve_gaussian( data, wid):
    g = gaussian( data.size, wid)
    fg = np.fft.fft(g)
    fdata = np.fft.fft(data)
    fout = np.fft.ifft( fdata*fg.conjugate())
    return fout

def convolve_window( data, wid):
    g = np.zeros(data.size)
    g[:wid//2] = 1.0
    g[-wid//2:] = 1.0
    fg = np.fft.fft(g)
    fdata = np.fft.fft(data)
    fout = np.fft.ifft( fdata*fg.conjugate())
    return fout

