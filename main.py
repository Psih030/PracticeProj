import numpy as np
from numpy.core.umath import sin, cos, pi, exp, sinh
import scipy.integrate as integrate
import scipy as sc
import wx

from wx.lib.plot import PlotCanvas, PlotGraphics, PolyLine, PolyMarker

b=1.0
a=4*b
st=1*b
c=a+st
G=1.67
lbd=8332.2
ym=15
kp=1.0004
pol=[1,2,3,4]
var=1
var2=1

def integr(x):
    global a,st
    var=1
    def integrand(bt, x):
        global a, st, G, lbd, ym, kp, c, var
        fun=-3
        if var==1:
            fun=polint(bt)*cos(bt*x)*(lbd+2*G)*(exp(bt*b)*(kp**2-1-2.*bt*b+2.*kp*bt*b)+exp(-bt*b)*(kp**2-1+2*bt*b-2*kp*bt*b))/((kp-1)*(lbd*(kp*exp(2*bt*b)+kp*exp(-2.*bt*b)+2)+G*kp*(exp(2*bt*b)+exp(-2*bt*b)))+2*G*(kp**2+kp+4*bt**2*b**2-2))
        elif var==2:
            fun=polint(bt)*(G*(exp(bt*b)*(kp**2-kp+2*kp*bt*b)-exp(-bt*b)*(kp**2-kp-2*kp*bt*b))+lbd*(kp**2-kp)*(exp(bt*b)-exp(-bt*b)))*cos(bt*x)/((kp**2-kp)*(lbd+G)*sinh(2*bt*b)+4*kp*G*b*bt);
        return fun

    return integrate.quad(integrand, 0,7,args=(x),limit=30)


def polint(bt):
    global a,c
    def integrand(ksi,bt):
        global pol
        fun=0
        for i in range(0,len(pol)):
            fun+=pol[i]*ksi**i
        fun=fun*cos(bt*ksi)
        return fun

    return integrate.quad(integrand, a, c, args=(bt),limit=100)[0]

def polint2(bt):
    global a,c,pol
    summ=0
    for j in range(0,len(pol)):
        sum2=0
        for i in range(0,j):
            sum2+=pohgamer(j,i)/(bt**(1+i))*(-1)*(i/2)*(c**(j-i)*sin(bt*c+(pi/2)**(i%2))-a**(j-i)*sin(bt*a+(pi/2)**(i%2)))
        summ+=pol[j]*sum2
    return summ

def pohgamer(j,i):
    poh=1
    for k in range(0,i):
        poh=poh*(j+1-i+k)
    return poh

def integr2(x):#int  d bt
    global a,st
    var=2
    def integrand(bt, x):
        global a, st, G, lbd, ym, kp, c, var
        fun=-3
        if var == 1:
            fun=cos(bt*x)*(sin(bt*c)-sin(bt*a))*(lbd+2*G)*(exp(bt*b)*(kp**2-1-2.*bt*b+2.*kp*bt*b)+exp(-bt*b)*(kp**2-1+2*bt*b-2*kp*bt*b))/bt/((kp-1)*(lbd*(kp*exp(2*bt*b)+kp*exp(-2.*bt*b)+2)+G*kp*(exp(2*bt*b)+exp(-2*bt*b)))+2*G*(kp**2+kp+4*bt**2*b**2-2))
        elif var ==2:
            fun=(G*(exp(bt*b)*(kp**2-kp+2*kp*bt*b)-exp(-bt*b)*(kp**2-kp-2*kp*bt*b))+lbd*(kp**2-kp)*(exp(bt*b)-exp(-bt*b)))*cos(bt*x)/bt*(sin(bt*c)-sin(bt*a))/((kp**2-kp)*(lbd+G)*sinh(2*bt*b)+4*kp*G*b*bt);
        return fun

    #return integrate.quadrature(integrand, 0,7,args=(x))
    return integrate.quad(integrand, 0,7,args=(x),limit=30)




def drawSinCosWaves():
    d=10000
    global a, st, G, lbd, ym, kp, c
    xv=np.arange(0,2*(a+st/2),0.1/2)
    mass=-ym*(lbd+2*G)*b/G*(kp-1)/(kp+1)
    xv.shape=(xv.size/2, 2)
    nagr=np.arange(xv.shape[0])
    for i in range(0,xv.shape[0]):
        x=xv[i,0]
        #nagr[i]=-2*d/pi/G*integr2(x)[0]
        nagr[i]=-2/pi/G*integr(x)[0]

    if var2==1:
        xv[:,1]=nagr
    if var2==2:
        xv[:,1]=nagr+mass
    sigma = PolyLine(xv, legend= 'Normal stress', colour='red')


    return PlotGraphics([sigma],"Graph Title", "X Axis", "Y Axis")

ID_LOAD=wx.NewId()
ID_INTERVAL=wx.NewId()
ID_MATERIAL=wx.NewId()
ID_COND=wx.NewId()
ID_MAIN_PANEL=wx.NewId()
ID_MAIN_FRAME=wx.NewId()
ID_REFRESH=wx.NewId()
ID_WEIGHT=wx.NewId()
########################################################################
class MyGraph(wx.Frame):
 
    #----------------------------------------------------------------------
    def __init__(self):
        wx.Frame.__init__(self, None, ID_MAIN_FRAME,
                          'MyStressPlot')
 
        # Add a panel so it looks the correct on all platforms
        panel = wx.Panel(self, ID_MAIN_PANEL)

        # Add menu
        menuBar = wx.MenuBar()
        menu1=wx.Menu()
        menu1.Append(ID_LOAD, "Set load", "")
        menu1.Append(ID_INTERVAL, "Set interval", "")
        menu1.Append(ID_MATERIAL, "Set material", "")
        menu1.Append(ID_COND, "Set boundary conditions", "")
        menu1.Append(ID_WEIGHT, "Set or remove dead weight", "")
        menu1.Append(ID_REFRESH, "Refresh plot", "")
        menuBar.Append(menu1, "&Settings")


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

        wx.EVT_MENU(self, ID_REFRESH, self.ReDraw)

        wx.EVT_MENU(self, ID_LOAD, self.ChangeLoad)
        wx.EVT_MENU(self, ID_INTERVAL, self.ChangeInterval)
        wx.EVT_MENU(self, ID_MATERIAL, self.ChangeMaterial)
        wx.EVT_MENU(self, ID_COND, self.ChangeBound)
        wx.EVT_MENU(self, ID_WEIGHT, self.ChangeWeight)
    #----------------------------------------------------------------------
    def onToggleGrid(self, event):
        """"""
        self.canvas.SetEnableGrid(event.IsChecked())
 
    #----------------------------------------------------------------------
    def onToggleLegend(self, event):
        """"""
        self.canvas.SetEnableLegend(event.IsChecked())

    def ReDraw(self,event):
        self.Show(False)
        self.canvas.Draw(drawSinCosWaves())
        self.Show(True)

    def ChangeWeight(self,event):
        dlg = wx.TextEntryDialog(parent=None, message="Write 1 to disable weight and 2 to enable weight", defaultValue='')
        dlg.ShowModal()
        result = dlg.GetValue()
        dlg.Destroy()
        bor=result.split()
        global var2
        var2=int(bor[0])

        self.Show(False)
        self.canvas.Draw(drawSinCosWaves())
        self.Show(True)


    def ChangeBound(self,event):
        dlg = wx.TextEntryDialog(parent=None, message="Write 1 for condition of fixed and 2 for smooth conditions", defaultValue='')
        dlg.ShowModal()
        result = dlg.GetValue()
        dlg.Destroy()
        bor=result.split()
        global var
        var=int(bor[0])

        self.Show(False)
        self.canvas.Draw(drawSinCosWaves())
        self.Show(True)

    def ChangeMaterial(self,event):
        dlg = wx.TextEntryDialog(parent=None, message="Write dead load, Young's modulus and Poisson's ratio, separated with space:", defaultValue='')
        dlg.ShowModal()
        result = dlg.GetValue()
        dlg.Destroy()
        bor=result.split()
        global lbd, G, ym, kp
        ym=float(bor[0])
        E=float(bor[1])
        mu=float(bor[2])
        G=E/(2*(1-mu))
        lbd=mu*E/((1+mu)*(1-2*mu))
        kp=3-4*mu
        self.Show(False)
        self.canvas.Draw(drawSinCosWaves())
        self.Show(True)

    def ChangeLoad(self,event):
        dlg = wx.TextEntryDialog(parent=None, message="Write polinom's coefficients, separated with space (starting with lowest power):", defaultValue='')
        dlg.ShowModal()
        result = dlg.GetValue()
        dlg.Destroy()
        bor=result.split()
        global pol
        l=[]
        for i in range(0,len(bor)):
            l.append(float(bor[i])*b)
        pol=l
        self.Show(False)
        self.canvas.Draw(drawSinCosWaves())
        self.Show(True)

    def ChangeInterval(self,event):
        dlg = wx.TextEntryDialog(parent=None, message='Write interval borders, separated with space:', defaultValue='')
        dlg.ShowModal()
        result = dlg.GetValue()
        dlg.Destroy()
        bor=result.split()
        global a,c,st
        a=float(bor[0])*b
        c=float(bor[1])*b
        st=c-a
        self.Show(False)
        self.canvas.Draw(drawSinCosWaves())
        self.Show(True)



 
if __name__ == '__main__':
    app = wx.App(False)
    frame = MyGraph()
    frame.Show()
    app.SetTopWindow(frame)
    app.MainLoop()