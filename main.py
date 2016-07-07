import numpy as np
from numpy.core.umath import sin, cos, pi, exp, sinh
import scipy.integrate as integrate
import wx
from wx.lib.plot import PlotCanvas, PlotGraphics, PolyLine
b=np.float64(1.0)
a=4*b
st=2*b
c=a+st
G=1.67
lbd=8332.2
ym=15
kp=1.0004
pol=[1,2,3,4,5]
var=1
var2=1

def integr(x):
    global a,st
    var=1
    def integrand(bt, x):
        global a, st, G, lbd, ym, kp, c, var
        fun=-1
        if var==1:
            fun=polint(bt)*cos(bt*x)*(lbd+2*G)*(exp(bt*b)*(kp**2-1-2*bt*b+2*kp*bt*b)+exp(-bt*b)*(kp**2-1+2*bt*b-2*kp*bt*b))/((kp-1)*(lbd*(kp*exp(2*bt*b)+kp*exp(-2*bt*b)+2)+G*kp*(exp(2*bt*b)+exp(-2*bt*b)))+2*G*(kp**2+kp+4*bt**2 * b **2 -2))
        elif var==2:
            fun=polint(bt)*(G*(exp(bt*b)*(kp**2-kp+2*kp*bt*b)-exp(-bt*b)*(kp**2-kp-2*kp*bt*b))+lbd*(kp**2-kp)*(exp(bt*b)-exp(-bt*b)))*cos(bt*x)/((kp**2-kp)*(lbd+G)*sinh(2*bt*b)+4*kp*G*b*bt)
        return fun

    return integrate.quad(integrand, 0,7,args=(x),limit=70,epsabs=1.49e-08,epsrel=1.49e-08,wopts="momcom")


def polint(bt):
    global a,c
    def integrand(ksi,bt):
        global pol
        fun=0
        for i in range(0,len(pol)):
            fun+=pol[i]*ksi**i
        fun=fun*cos(bt*ksi)
        return fun

    return integrate.quad(integrand, a, c, args=(bt),limit=100,epsabs=1.49e-08,epsrel=1.49e-08, wopts="momcom")[0]

def pohgamer(j,i):
    poh=1
    for k in range(0,i):
        poh=poh*(j+1-i+k)
    return poh


def drawSinCosWaves():
    global a, st, G, lbd, ym, kp, c
    xv=np.arange(0,2*(a+st/2),0.1/2)
    mass=-ym*(lbd+2*G)*b/G*(kp-1)/(kp+1)
    xv.shape=(xv.size/2, 2)
    nagr=np.arange(xv.shape[0])
    for i in range(0,xv.shape[0]):
        x=np.float64(xv[i,0])
        nagr[i]=-2/pi/G*integr(x)[0]

    if var2==1:
        xv[:,1]=nagr
    if var2==2:
        xv[:,1]=nagr+mass
    sigma = PolyLine(xv, legend= 'Normal stress', colour='red')


    return PlotGraphics([sigma],"Graph Title", "X Axis", "Y Axis")

class MyChild(wx.Frame):
    def __init__(self,parent):
        wx.Frame.__init__(self, None, title="InfoFrame")
        self.parent=parent
        self.panel=wx.Panel(self)
        self.quotea = wx.StaticText(self.panel, label="a=")
        self.edita = wx.TextCtrl(self.panel, size=(140, -1), value=str(a))

        self.quotec = wx.StaticText(self.panel, label="c=")
        self.editc = wx.TextCtrl(self.panel, size=(140, -1), value=str(c))

        self.quoteE = wx.StaticText(self.panel, label="E=")
        self.editE = wx.TextCtrl(self.panel, size=(140, -1))

        self.quotem = wx.StaticText(self.panel, label="mu=")
        self.editm = wx.TextCtrl(self.panel, size=(140, -1))

        self.quoteg = wx.StaticText(self.panel, label="gamma=")
        self.editg = wx.TextCtrl(self.panel, size=(140, -1), value=str(ym))

        self.quoteG = wx.StaticText(self.panel, label="G="+str(G))

        self.quotel = wx.StaticText(self.panel, label="lambda="+str(lbd))

        strp=""
        for i in range(0,len(pol)):
            strp+=str(pol[i])+" "

        self.quotep = wx.StaticText(self.panel, label="pol=")
        self.editp = wx.TextCtrl(self.panel, size=(140, -1), value=strp)

        self.button = wx.Button(self.panel, label="Save")


        lblList = ['of fixed', 'of smooth']
        self.rbox = wx.RadioBox(self.panel,label = 'Choose boundary condition',  choices = lblList , majorDimension = 1, style = wx.RA_SPECIFY_COLS)
        if var==1:
            self.rbox.SetStringSelection("of fixed")
        elif var==2:
            self.rbox.SetStringSelection("of smooth")

        lblList = [ 'no','yes']
        self.rboxm = wx.RadioBox(self.panel,label = 'Dead load?',  choices = lblList , majorDimension = 1, style = wx.RA_SPECIFY_COLS)
        if var2==1:
            self.rboxm.SetStringSelection("no")
        if var2==2:
            self.rboxm.SetStringSelection("yes")

        self.windowSizer = wx.BoxSizer()
        self.windowSizer.Add(self.panel, 1, wx.ALL | wx.EXPAND)

        self.sizer = wx.GridBagSizer(10, 10)
        self.sizer.Add(self.quotea, (0, 0))
        self.sizer.Add(self.edita, (0, 1))
        self.sizer.Add(self.quotec, (1, 0))
        self.sizer.Add(self.editc, (1, 1))
        self.sizer.Add(self.quotem, (2, 0))
        self.sizer.Add(self.editm, (2, 1))
        self.sizer.Add(self.quoteE, (3, 0))
        self.sizer.Add(self.editE, (3, 1))
        self.sizer.Add(self.quoteG, (4,0))
        self.sizer.Add(self.quotel, (5,0))
        self.sizer.Add(self.quoteg, (6, 0))
        self.sizer.Add(self.editg, (6, 1))
        self.sizer.Add(self.quotep, (7, 0))
        self.sizer.Add(self.editp, (7, 1))
        self.sizer.Add(self.rbox, (8, 0))
        self.sizer.Add(self.rboxm, (8, 1))
        self.sizer.Add(self.button, (9, 0), (1, 2), flag=wx.EXPAND)

        self.border = wx.BoxSizer()
        self.border.Add(self.sizer, 1, wx.ALL | wx.EXPAND, 5)

        self.panel.SetSizerAndFit(self.border)
        self.SetSizerAndFit(self.windowSizer)

        self.button.Bind(wx.EVT_BUTTON, self.OnButton)

    def OnButton(self,event):
        global a, c, G, st, lbd, ym, pol, kp,var,var2
        a=float(self.edita.GetValue())
        c=float(self.editc.GetValue())
        st=c-a
        if self.editE.GetValue()!="":
            E=float(self.editE.GetValue())
        if self.editm.GetValue()!="":
            mu=float(self.editm.GetValue())
        ym=float(self.editg.GetValue())
        if self.editE.GetValue()!="" and self.editm.GetValue()!="":
            G=E/(2*(1-mu))
            lbd=mu*E/((1+mu)*(1-2*mu))
            kp=3-4*mu
        pold=self.editp.GetValue().split()
        for i in range(0,len(pold)):
            print pold[i]
            pol[i]=float(pold[i])
        if self.rbox.GetStringSelection()=="of fixed":
            var=1
        elif self.rbox.GetStringSelection()=="of smooth":
            var=2
        if self.rboxm.GetStringSelection()=="no":
            var2=1
        elif self.rboxm.GetStringSelection()=="yes":
            var2=2
        self.Destroy()

        tp=wx.GetApp().GetTopWindow()
        tp.Show(False)
        tp.canvas.Draw(drawSinCosWaves())
        tp.Show(True)

ID_LOAD=wx.NewId()
ID_INTERVAL=wx.NewId()
ID_MATERIAL=wx.NewId()
ID_COND=wx.NewId()
ID_MAIN_PANEL=wx.NewId()
ID_MAIN_FRAME=wx.NewId()
ID_REFRESH=wx.NewId()
ID_WEIGHT=wx.NewId()
ID_SIDE_PANEL=wx.NewId()
ID_SET_ALL=wx.NewId()

class MyGraph(wx.Frame):

    def __init__(self):
        wx.Frame.__init__(self, None, ID_MAIN_FRAME,
                          'MyStressPlot')

        panel = wx.Panel(self, ID_MAIN_PANEL)

        menuBar = wx.MenuBar()
        menu1=wx.Menu()
        menu1.Append(ID_LOAD, "Set load", "")
        menu1.Append(ID_INTERVAL, "Set interval", "")
        menu1.Append(ID_MATERIAL, "Set material", "")
        menu1.Append(ID_COND, "Set boundary conditions", "")
        menu1.Append(ID_WEIGHT, "Set or remove dead weight", "")
        menu1.Append(ID_REFRESH, "Refresh plot", "")
        menuBar.Append(menu1, "&Settings")
        menu2=wx.Menu()
        menu2.Append(ID_SET_ALL, "Set all", "")
        menuBar.Append(menu2, "&Settings2")


        mainSizer = wx.BoxSizer(wx.VERTICAL)
        checkSizer = wx.BoxSizer(wx.HORIZONTAL)

        self.canvas = PlotCanvas(panel)
        self.canvas.Draw(drawSinCosWaves())
        toggleGrid = wx.CheckBox(panel, label="Show Grid")
        toggleGrid.Bind(wx.EVT_CHECKBOX, self.onToggleGrid)
        toggleLegend = wx.CheckBox(panel, label="Show Legend")
        toggleLegend.Bind(wx.EVT_CHECKBOX, self.onToggleLegend)


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

        wx.EVT_MENU(self, ID_SET_ALL, self.SetAll)

    def SetAll(self,event):
        self.Child= MyChild(self)
        self.Child.Show()


    def onToggleGrid(self, event):
        self.canvas.SetEnableGrid(event.IsChecked())

    def onToggleLegend(self, event):
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