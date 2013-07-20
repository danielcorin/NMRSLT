"""
NMRSLT: Nuclear Magnetic Resonance Spectroscopy Learning Tool
"""

import scipy as sp
import pylab as pl

from matplotlib.widgets import Slider, Button, RadioButtons
from scipy import fftpack as fp

class CurveParams:
    def __init__(self, period, amplitude, env, a, b):
        self.period = period
        self.amplitude = amplitude
        self.env = env
        self.a = a
        self.b = b
    def setAmplitude(self, amp):
        self.amplitude = amp
    def setPeriod(self, per):
        self.period = per
    def setEnvelope(self, env):
        self.env = env
    def getAmplitude(self):
        return self.amplitude
    def getPeriod(self):
        return self.period
    def getEnvelope(self):
        return self.env

class RadioSelection:
    def __init__(self, sel):
        self.selection = sel
    def getSelection(self):
        return self.selection
    def setSelection(self, sel):
        self.selection = sel

def fft(xaxis, real, img, dwell, numPts):
        #check numPts is to a power of 2
        if not float(sp.log2(numPts)).is_integer():
            print 'Number of points is not a power of 2.'
            return
 
        fftReal = fp.fft(real)
        fftImg = fp.fft(img)
        
        #reorganize fourier transform properly
        #1-4 new real
        #2+3 new img
        fftR = sp.real(fftReal[:]) - sp.imag(fftImg[:])
        fftR = fp.fftshift(fftR)
        fftI = sp.imag(fftReal[:]) + sp.real(fftImg[:])
        fftI = fp.fftshift(fftI)
        
        #convert axis time to freq
        freq = timeToFreq(xaxis, dwell, numPts)
        
        return [freq, fftR, fftI]

def timeToFreq(time, dwell, numPts):
        sweepWidth = 1/dwell
        bound = sweepWidth/2
        hzPerPt = (sweepWidth/numPts)
        freq = -bound+(hzPerPt*sweepWidth*time)
        return freq

# will use later if I implement a phase slider
def phase(real, img, degreePhase):
        rad = degreePhase*(sp.pi)/180
        pReal = (real*sp.cos(rad)-img*sp.sin(rad))
        pImg = (img*sp.cos(rad)+real*sp.sin(rad))
        return [pReal, pImg]

def gaussian(x,a,b,c):
    return a * sp.exp(-((x-b)**2)/(2*c**2))
  
# could also use later for more diverse peak types              
def lorentzian(x,a,b,c):
    return a*1/(1+((4/3)*((x-b)/(c))**2))

def multi(x, func, curves):
    out = 0
    for curve in curves:
        out += ( gaussian(x, curve.a, curve.b, curve.getEnvelope()) *
            curve.getAmplitude() * func(2*sp.pi / curve.getPeriod()*x) )
    return out

# params for gaussian noise
mean = 0
std = .025

# data set params
numPts = 16384
dwell = .000001

# create time axis
x = sp.arange(0, numPts * dwell, dwell)

# initial params and equations
a = 1
b = 0
c = .002
period = .002
amp = 1

# build curves
curve1 = CurveParams(period, amp, c, a, b)
curve2 = CurveParams(period / 15, 0.1*amp, c, a, b)
curves = [curve1, curve2]
real = multi(x, sp.cos, curves)
img = multi(x, sp.sin, curves)

# noise to make the data look "real"
noiseX = sp.random.normal(mean, std, x.size)
noiseY = sp.random.normal(mean, std, x.size)
real += noiseX
img += noiseY

# transform data for bottom view port
xfft, realfft, imgfft = fft(x, real, img, x[1]-x[0], x.size)

# create selector object for curve switching with radio buttons
radioSel = RadioSelection(0)

# setup the window
fig = pl.figure(figsize=[9,7])
fig.canvas.set_window_title("NMR transform tool")
ax1 = fig.add_subplot(211)
ax1.set_title("Time domain (seconds)")
pl.ylim([-4,4])
ax2 = fig.add_subplot(212)
ax2.set_title("Frequency domain (Hz)")
pl.xlim([-100000,100000])
pl.ylim([-500,6000])


l, = ax1.plot(x, real, color='blue')
j, = ax2.plot(xfft, realfft, color='red')

pl.subplots_adjust(left=0.1, bottom=0.25, top = 0.95)

# place axes on window
axcolor = 'lightgoldenrodyellow'
axEnvelope = pl.axes([0.1, 0.15, 0.6, 0.03], axisbg=axcolor)
axPeriod = pl.axes([0.1, 0.05, 0.6, 0.03], axisbg=axcolor)
axAmp = pl.axes([0.1, 0.1, 0.6, 0.03], axisbg=axcolor)

# to avoid divide by zero errors   
almostZero = .00000001

# define sliders for each parameter that can vary
envelopeSlider = Slider(axEnvelope, 'Envelope', almostZero, 2*c,
        valinit=c, valfmt = '%1.5f')

periodSlider = Slider(axPeriod, 'Period', almostZero, period * 1.5,
        valinit=period, valfmt = '%1.7f')

amplitudeSlider = Slider(axAmp, 'Amplitude', -2*amp, 2*amp,
        valinit=amp, valfmt = '%1.3f')


# define parameter updating functions for each of the varying parameters
# checks to see which curve is active and changes the corresponding parameter
def updateEnvelope(val):
    curves[radioSel.getSelection()].setEnvelope(envelopeSlider.val)
    updateCurve()
     
def updatePeriod(val):
    curves[radioSel.getSelection()].setPeriod(periodSlider.val)
    updateCurve()

def updateAmp(val):
    curves[radioSel.getSelection()].setAmplitude(amplitudeSlider.val)
    updateCurve()

# general curve update function
def updateCurve():
    # update curves with new parameters
    re = noiseX + multi(x, sp.cos, curves)
    im = noiseY + multi(x, sp.sin, curves)
    
    # re-transform
    xfft, realfft, imgfft = fft(x, re, im, x[1]-x[0], x.size)

    # display new data
    l.set_ydata(re)
    j.set_ydata(realfft)
    pl.draw()

rax = pl.axes([0.85, 0.05, 0.1, 0.1], axisbg=axcolor)
radio = RadioButtons(rax, ('Curve 1', 'Curve 2'), active=0)

def switchCurve(label):
    if label == "Curve 1":
        radioSel.setSelection(0)
    else:
        radioSel.setSelection(1)
    
    # update sliders with selected curve's parameters
    periodSlider.set_val(curves[radioSel.getSelection()].getPeriod())
    amplitudeSlider.set_val(curves[radioSel.getSelection()].getAmplitude())
    envelopeSlider.set_val(curves[radioSel.getSelection()].getEnvelope())

    pl.draw()

# enable interactive controls
periodSlider.on_changed(updatePeriod)
amplitudeSlider.on_changed(updateAmp)
envelopeSlider.on_changed(updateEnvelope)
radio.on_clicked(switchCurve)

# allow window to be closed through the command line
pl.ion()
pl.show()

raw_input("Press enter to close")

