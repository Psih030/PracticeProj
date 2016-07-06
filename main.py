import numpy as np
from numpy.core.umath import sin, cos, pi, exp
import scipy.integrate as integrate
import scipy as sc
import wx

from wx.lib.plot import PlotCanvas, PlotGraphics, PolyLine, PolyMarker

b=1.0
a=2*b
st=5*b
c=a+st
G=1.6668
lbd=8332.2
ym=15
kp=1.0004
list=[1, -2, 3]

def integr(x):
    global a,st
    var=1
    def integrand(bt, x, var):
        global a, st, G, lbd, ym, kp, c
        d=10
        fun=polint(bt)*cos(bt*x)*(lbd+2*G)*(exp(bt*b)*(kp**2-1-2.*bt*b+2.*kp*bt*b)+exp(-bt*b)*(kp**2-1+2*bt*b-2*kp*bt*b))/((kp-1)*(lbd*(kp*exp(2*bt*b)+kp*exp(-2.*bt*b)+2)+G*kp*(exp(2*bt*b)+exp(-2*bt*b)))+2*G*(kp**2+kp+4*bt**2*b**2-2))
        return fun

    return integrate.quad(integrand, 0,10,args=(x,var))

def integr2(x):#int  d bt
    global a,st
    var=1
    def integrand(bt, x, var):
        global a, st, G, lbd, ym, kp, c
        fun=cos(bt*x)*(sin(bt*c)-sin(bt*a))*(lbd+2*G)*(exp(bt*b)*(kp**2-1-2.*bt*b+2.*kp*bt*b)+exp(-bt*b)*(kp**2-1+2*bt*b-2*kp*bt*b))/bt/((kp-1)*(lbd*(kp*exp(2*bt*b)+kp*exp(-2.*bt*b)+2)+G*kp*(exp(2*bt*b)+exp(-2*bt*b)))+2*G*(kp**2+kp+4*bt**2*b**2-2))
        return fun

    return integrate.quad(integrand, 0,100,args=(x,var))


def polint(bt):
    global a,c
    def integrand(ksi,bt):
        global list
        fun=0
        for i in range(0,len(list)):
            fun+=list[i]*ksi**i
        fun=fun*cos(bt*ksi)
        return fun

    return integrate.quad(integrand, a, c, args=(bt))[0]

def drawSinCosWaves():
    d=10
    global a, st, G, lbd, ym, kp, c
    xv=np.arange(0,2*(a+st/2),0.1/2)
    #mass=-ym*(lbd+2*G)*b/G*(kp-1)/(kp+1)
    xv.shape=(xv.size/2, 2)
    nagr=np.arange(xv.shape[0])
    err=np.arange(xv.shape[0])
    for i in range(0,xv.shape[0]):
        x=xv[i,0]
        nagr[i]=-2/pi/G*integr(x)[0]


    xv[:,1]=nagr
    sigma = PolyLine(xv, legend= 'Normal stress', colour='red')


    return PlotGraphics([sigma],"Graph Title", "X Axis", "Y Axis")
 
########################################################################
class MyGraph(wx.Frame):
 
    #----------------------------------------------------------------------
    def __init__(self):
        wx.Frame.__init__(self, None, wx.ID_ANY, 
                          'MyStressPlot')
 
        # Add a panel so it looks the correct on all platforms
        panel = wx.Panel(self, wx.ID_ANY)

        # Add menu
        menuBar = wx.MenuBar()
        menu1=wx.Menu()
        menu1.Append(wx.ID_ANY, "Set load", "")
        menu1.Append(wx.ID_ANY, "Set interval", "")
        menu1.Append(wx.ID_ANY, "Set material", "")
        menu1.Append(wx.ID_ANY, "Set border conditions", "")
        menuBar.Append(menu1, "&Settings")
        menu2=wx.Menu()
        menuBar.Append(menu2, "&Redraw")


        # create some sizers
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        checkSizer = wx.BoxSizer(wx.HORIZONTAL)
 
        # create the widgets
        self.canvas = PlotCanvas(panel)
        self.canvas.Draw(drawSinCosWaves())
        toggleGrid = wx.CheckBox(panel, label="Show Grid")
        toggleGrid.Bind(wx.EVT_CHECKBOX, self.onToggleGrid)
        toggleLegend = wx.CheckBox(panel, label="Show Legend")
        toggleLegend.Bind(wx.EVT_CHECKBOX, self.onToggleLegend)
 
        # layout the widgets
        mainSizer.Add(self.canvas, 1, wx.EXPAND)
        checkSizer.Add(toggleGrid, 0, wx.ALL, 5)
        checkSizer.Add(toggleLegend, 0, wx.ALL, 5)
        mainSizer.Add(checkSizer)
        panel.SetSizer(mainSizer)

        self.SetMenuBar(menuBar)
 
    #----------------------------------------------------------------------
    def onToggleGrid(self, event):
        """"""
        self.canvas.SetEnableGrid(event.IsChecked())
 
    #----------------------------------------------------------------------
    def onToggleLegend(self, event):
        """"""
        self.canvas.SetEnableLegend(event.IsChecked())
 
if __name__ == '__main__':
    app = wx.App(False)
    frame = MyGraph()
    frame.Show()
    app.MainLoop()