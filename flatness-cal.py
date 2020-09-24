from itertools import combinations 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.colors import Normalize
from matplotlib import cm
import numpy as np
from numpy.linalg import lstsq
import operator
import os
import pandas as pd
from scipy.interpolate import griddata
import sys
import tkinter as tk
from tkinter.filedialog import askopenfilename
################################################################################
PX = 5  # Padx
PY = 1  # Pady
NP = 13 # Max manual input
Pos = 'P'
Neg = 'N'
################################################################################
class CalFunc():
    def SplitGrp(self, Df):
        Grp = []
        for i in range(len(Df.index)):
            G = Pos if Df.Err[i] >= 0 else Neg
            Grp.append(G)
        Df['Grp'] = Grp
        return Df

    def PlaneEq(self, Df):
        P0 = Df.iloc[ 0 , : ]
        P1 = Df.iloc[ 1 , : ]
        P2 = Df.iloc[ 2 , : ]
        ux, uy, uz = [P1.X-P0.X, P1.Y-P0.Y, P1.Z-P0.Z]
        vx, vy, vz = [P2.X-P0.X, P2.Y-P0.Y, P2.Z-P0.Z]
        u_cross_v = [uy*vz-uz*vy, uz*vx-ux*vz, ux*vy-uy*vx]
        Pt  = np.array(P0)
        Norm = np.array(u_cross_v)
        A = Norm[0]
        B = Norm[1]
        C = Norm[2]
        D = 0 if all(flag == 0 for flag in (A, B, C)) else (-Pt.dot(Norm))
        return (A, B, C, D)

    def OrthoD(self, Pt, A, B, C, D, ABS):
        L = (A * Pt.X + B * Pt.Y + C * Pt.Z + D)
        E = (np.sqrt(A * A + B * B + C * C))
        return abs(L/E) if ABS else (L/E)

    def OffsetP(self, A, B, C, Df, Grp):
        if all(flag == 0 for flag in (A, B, C)):
            D = 0
        else:
            Pt = Df.iloc[Df[Grp].argmax()]
            norm = [A, B, C]
            P = np.array([Pt.X, Pt.Y, Pt.Z])
            N = np.array(norm)
            D = -P.dot(N)
        return D

    def IsTruth(self, A, Op, B):
        # operator.eq, operator.gt, operator.ge, operator.lt, operator.le
        return Op(A, B)

    def BestFitPlane(self, Df, Xlim, Ylim):
        Tmp_A = []
        Tmp_B = []
        for i in range(len(Df.X)):
            Tmp_A.append([Df.X[i], Df.Y[i], 1])
            Tmp_B.append(Df.Z[i])
        B = np.matrix(Tmp_B).T
        A = np.matrix(Tmp_A)
        fit, _, _, _ = lstsq(A, B, rcond=None)
        Errors = B - A * fit
        Df['Err'] = Errors  # Insert back into dataframe
        X, Y = np.meshgrid(np.arange(Xlim[0], Xlim[1]), \
                np.arange(Ylim[0], Ylim[1]))
        Z = np.zeros(X.shape)
        for r in range(X.shape[0]):
            for c in range(X.shape[1]):
                Z[r,c] = fit[0] * X[r,c] + fit[1] * Y[r,c] + fit[2]
        return X, Y, Z, Df

    def VPlane(self, Df, Grp):
        Tmp = Df.loc[Df['Grp'] == Grp, ['X', 'Y', 'Z', 'Err']]
        # for Top plane use p and n for bottom plane
        Order = False if (Grp == Pos) else True
        Tmp = Tmp.sort_values(by=['Err'], ascending=Order)
        for Subset in list(combinations(Tmp.index, 3)):
            A, B, C, D = self.PlaneEq(Tmp.loc[Subset, ['X', 'Y', 'Z']])
            # Detect and skip if not valid plane
            if all(flag == 0 for flag in (A, B, C, D)):
                return (Df, A, B, C, D)
            # Start looping
            Exit = False
            Df[Grp] = 0
            Df[Grp + 'C'] = ""
            for i in range(len(Df.index)):
                if i in Subset:
                    Df.loc[i, Grp] = float(0)
                    Df.loc[i, Grp + 'C'] = 'r'
                else:
                    Df.loc[i, Grp + 'C'] = 'b'
                    Dist = self.OrthoD(Df.iloc[ i , : ], A, B, C, D, False)
                    Df.loc[i, Grp] = abs(round(Dist, 3))
                    if Order:
                        Op = operator.lt if C > 0 else operator.gt
                    else:
                        Op = operator.gt if C > 0 else operator.lt
                    if self.IsTruth(round(Dist, 3), Op, 0):
                        Exit = False
                        break
                    else:
                        Exit = True
            if Exit:
                Df.loc[Df[Grp].argmax(), Grp + 'C'] = 'r'
                break
        return (Df, A, B, C, D)
    
    def MeshGrid(self, A, B, C, D, Xlim, Ylim):
        if all(flag == 0 for flag in (A, B, C)):
            return ([], [], [])
        else:
            XX, YY = np.meshgrid(np.arange(Xlim[0], Xlim[1]), np.arange(Ylim[0], Ylim[1]))
            ZZ = ((-A * XX - B * YY -D)* 1.) if C == 0 else ((-A * XX - B * YY -D) * 1. / C)
            return (XX, YY, ZZ)

    def AxLim(self, Df):
        XLow, XHgh = (Df['X'].min(), Df['X'].max())
        YLow, YHgh = (Df['Y'].min(), Df['Y'].max())
        ZLow, ZHgh = (Df['Z'].min(), Df['Z'].max())
        XL, YL, ZL = (XHgh - XLow), (YHgh - YLow), (ZHgh - ZLow)
        if (XL > YL): 
            XLow, XHgh, YLow, YHgh = self.EqualXY(XLow, XHgh, YLow, YHgh)
        else:
            YLow, YHgh, XLow, XHgh = self.EqualXY(YLow, YHgh, XLow, XHgh)           
        Rnd = 0 if ZL > 1 else 3
        ZLow = ZLow - (round(ZL, Rnd) * 0.75)
        ZHgh = ZHgh + (round(ZL, Rnd) * 0.25)
        return ([XLow, XHgh], [YLow, YHgh], [ZLow, ZHgh])

    def EqualXY(self, AL, AH, BL, BH):
        AL = round(AL - round((AH - AL) * 0.15))
        AH = round(AH + round((AH - AL) * 0.15))
        Offset = ((AH - AL) - (BH - BL)) / 2
        BL, BH = (BL - Offset), (BH + Offset )
        return (AL, AH, BL, BH)     

class MenuBar(tk.Menu):
    def __init__(self, parent, controller):
        self.Rd = tk.IntVar()
        self.Rd.set(0)
        self.show_BFP = tk.BooleanVar()
        self.show_BFP.set(False)
        self.ZScale = tk.BooleanVar()
        self.ZScale.set(True)
        self.TP = tk.BooleanVar()
        self.BP = tk.BooleanVar()
        self.TP.set(False)
        self.BP.set(False)

        tk.Menu.__init__(self, parent)
        self.controller = controller
        FileMenu = tk.Menu(self, tearoff=False)
        self.add_cascade(label="File", menu=FileMenu)
        FileMenu.add_command(label="Exit", command=self.quit)

        self.ViewMenu = tk.Menu(self, tearoff=False)
        self.add_cascade(label="View", menu=self.ViewMenu)
        self.ViewMenu.add_checkbutton(label="Top Solution", onvalue=1, \
            offvalue=0, variable=self.TP, selectcolor='lime', \
            command=self.ShowTP)
        self.ViewMenu.add_checkbutton(label="Bot Solution", onvalue=1, \
            offvalue=0, variable=self.BP, selectcolor='lime', \
            command=self.ShowBP)
        self.ViewMenu.add_separator()
        self.ViewMenu.add_radiobutton(label="Virtual Datum", variable=self.Rd, \
            value=0, selectcolor='lime', command=self.RdChg)
        self.ViewMenu.add_radiobutton(label="Surf Gradient", variable=self.Rd, \
            value=1, selectcolor='lime', command=self.RdChg)
        self.ViewMenu.add_separator()
        self.ViewMenu.add_checkbutton(label="Best-Fit Datum", onvalue=1, 
            offvalue=0, variable=self.show_BFP, selectcolor='lime', \
            command=self.ShowBFP)
        self.ViewMenu.add_checkbutton(label="AutoScale", onvalue=1, 
            offvalue=0, variable=self.ZScale, selectcolor='lime', \
            command=self.controller.UpdateFig)

    def ShowBFP(self):
        self.controller.show_BFP(self.show_BFP.get())
    
    def RdChg(self):
        if self.Rd.get():
            if self.TP.get() and self.BP.get():
                self.TP.set(True)
                self.BP.set(False)
        self.controller.UpdateFig()
        
    def ShowTP(self):
        if self.Rd.get():
            self.BP.set(False)
        self.controller.UpdateFig()

    def ShowBP(self):
        if self.Rd.get():
            self.TP.set(False)
        self.controller.UpdateFig()

    def DisableSolu(self, Grp, opt):
        Op = tk.DISABLED if opt else tk.NORMAL
        N = 0 if Grp == Pos else 1
        self.ViewMenu.entryconfig(N, state=Op)
          
    def quit(self):
        sys.exit(0)

class MainApp(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        tk.Tk.wm_title(self, "Flatness Calculator")

        self.menubar = MenuBar(self, controller=self)
        self.config(menu=self.menubar)

        container = tk.Frame(self)
        container.pack(side="left", fill="both")
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        self.frames = {}
        for F in (Manual, Import):
            frame = F(parent=container, controller=self)
            self.frames[F] = frame
            frame.grid(row=0, column=0, sticky="nsew")
        self.show_frame(Import)
        
        self.Fig = Figure(figsize=(6,6), dpi=100)
        self.canvas = FigureCanvasTkAgg(self.Fig, self)
        self.canvas.get_tk_widget().pack(side='left', fill='both', expand=True)
        self.Ax = self.Fig.add_subplot(111, projection="3d", proj_type='ortho')
        self.CbAxes = self.Fig.add_axes([0.15, 0.07, 0.7, 0.01], visible=False)
        self.canvas.draw()

    def show_frame(self, Name):
        frame = self.frames[Name]
        frame.tkraise()

    def ScaleZ(self):
        if self.menubar.ZScale.get():
            self.Ax.set_zlim(self.Zlim[0], self.Zlim[1])
        else:
            XYLim = self.Ax.get_xlim()
            Off = (XYLim[1] -  XYLim[0] - (self.Zlim[1] - self.Zlim[0])) / 2
            self.Ax.set_zlim((self.Zlim[0] - Off), (self.Zlim[1] + Off))

    def UpdateFig(self):
        TP = self.menubar.TP.get()
        BP = self.menubar.BP.get()
        RD = self.menubar.Rd.get()
        T, B, S, G = (False, False, False, True) if RD else \
                        (TP, BP, True, False)
        self.ShowTP(T)
        self.ShowBP(B)
        self.Show_Sup(S)
        self.ShowGd(G)
        self.Ax.set_xlabel('x')
        self.Ax.set_ylabel('y')
        self.Ax.set_zlabel('z')
        self.ScaleZ()
        self.canvas.draw()

    def show_BFP(self, opt):
        if self.BFP:
            self.BFP.set_visible(opt)
            self.canvas.draw()

    def ShowTP(self, opt):
        if not self.menubar.ViewMenu.entrycget(0, 'state') == 'disabled' and \
            hasattr(self, 'P_WF'):
            self.P_WF.set_visible(opt)
            self.PO_WF.set_visible(opt)
            self.P_Sp.set_visible(opt)

    def ShowBP(self, opt):
        if not self.menubar.ViewMenu.entrycget(1, 'state') == 'disabled' and \
            hasattr(self, 'N_WF'):
            self.N_WF.set_visible(opt)
            self.NO_WF.set_visible(opt)
            self.N_Sp.set_visible(opt)

    def ShowGd(self, opt):
        MB = self.menubar
        C = self.Cmap(self.PN(self.PC)) if MB.TP.get() else \
            self.Cmap(self.NN(self.NC))
        CSize = (C.shape[0] - 1) * (C.shape[1] - 1)
        if opt:
            if hasattr(self, 'Surface'):
                self.Surface.set_visible(True)
                self.Surface.set_facecolor(C[:-1, :-1].reshape(CSize, 4))
        else:
            if hasattr(self, 'Surface'):
                self.Surface.set_visible(False)
        self.Show_CBar(opt)
        
    def Show_Sup(self, opt):
        if opt:
            if ((self.menubar.TP.get()) and (self.menubar.BP.get())):
                Str = "Flatness (S1:S2): ({}:{})".format(self.Df[Pos].max(), \
                    self.Df[Neg].max())
            elif self.menubar.TP.get():
                Str = "Flatness S1: {}".format(self.Df[Pos].max())
            elif self.menubar.BP.get():
                Str = "Flatness S2: {}".format(self.Df[Neg].max())
            else:
                Str = ''
            self.SubTitle = self.Fig.suptitle(Str, y=0.05)
        else:
            if hasattr(self, 'SubTitle'):
                self.SubTitle.set_visible(False)

    def Show_Wireframe(self, X, Y, Z, Color):
        if (len(X) and len(Y) and len(Z)):
            return self.Ax.plot_wireframe(X, Y, Z, color=Color, alpha=0.5)

    def Get_Grid(self, Df, X, Y, Grp):
        N = griddata((Df['X'], Df['Y']), Df[Grp], (X, Y), method='cubic')
        if np.isnan(np.sum(N)):
            A = griddata((Df['X'], Df['Y']), Df['Z'], (X, Y), method='nearest')
            N = np.where(np.isnan(N), A, N)
        return N

    def Cal_Surface(self, Df):
        Cmap = cm.get_cmap('jet')
        X1 = np.linspace(Df['X'].min(), Df['X'].max(), \
                len(Df['X'].unique()) * 5)
        Y1 = np.linspace(Df['Y'].min(), Df['Y'].max(), \
                len(Df['Y'].unique()) * 5)
        X2, Y2 = np.meshgrid(X1, Y1)
        Z2 = self.Get_Grid(Df, X2, Y2, 'Z')
        if Pos in Df:
            PC = self.Get_Grid(Df, X2, Y2, Pos)
            PNorm = Normalize(vmin=Df['P'].min(), vmax=Df['P'].max())
        else:
            PC, PNorm = [], []
        if Neg in Df:
            NC = self.Get_Grid(Df, X2, Y2, Neg)
            NNorm = Normalize(vmin=Df['N'].min(), vmax=Df['N'].max())
        else:
            NC, NNorm = [], []
        return (X2, Y2, Z2, Cmap, PC, PNorm, NC, NNorm)

    def Show_CBar(self, opt):
        if opt:
            Norm = self.PN if self.menubar.TP.get() else self.NN
            self.Map = cm.ScalarMappable(cmap=self.Cmap, norm=Norm)
            if not hasattr(self, 'Cbar'):
                self.Cbar = self.Fig.colorbar(self.Map, cax=self.CbAxes, \
                orientation='horizontal')
            else:
                self.Cbar.update_normal(self.Map)
            self.Cbar.ax.tick_params(labelsize='small', rotation=30)
            Tck = self.Cbar.ax.get_xticks()
            Tck = np.append(Tck, Norm.vmax)
            self.Cbar.set_ticks(Tck)
        self.CbAxes.set_visible(opt)

    def Solve(self, pts):
        Cal = CalFunc()
        self.Ax.clear()
        self.Df = pts
        # Get Xlim and Ylim from data
        Xlim, Ylim, self.Zlim = Cal.AxLim(self.Df)
        # Find Best Fit Plane
        BFX, BFY, BFZ, pts = Cal.BestFitPlane(self.Df, Xlim, Ylim)
        # Split all point into 2 group base on vector of the plane
        self.Df = Cal.SplitGrp(self.Df)
        # Calculate 2 possible solution base on group
        self.Df, S1A, S1B, S1C, S1D = Cal.VPlane(self.Df, Pos)
        Cal.OffsetP(S1A, S1B, S1C, pts, Pos)
        # Create mesh grid for S1
        S1X, S1Y, S1Z = Cal.MeshGrid(S1A, S1B, S1C, S1D, Xlim, Ylim)
        S1OX, S1OY, S1OZ = Cal.MeshGrid(S1A, S1B, S1C, \
                            Cal.OffsetP(S1A, S1B, S1C, self.Df, Pos), \
                            Xlim, Ylim)
        self.Df, S2A, S2B, S2C, S2D = Cal.VPlane(self.Df, Neg)
        S2X, S2Y, S2Z = Cal.MeshGrid(S2A, S2B, S2C, S2D, Xlim, Ylim)
        S2OX, S2OY, S2OZ = Cal.MeshGrid(S2A, S2B, S2C, \
                            Cal.OffsetP(S2A, S2B, S2C, self.Df, Neg), \
                            Xlim, Ylim) 
        # Plotting section
        if 'PC' in self.Df:
            self.P_Sp = self.Ax.scatter(self.Df.X, self.Df.Y, self.Df.Z, \
                c=self.Df.PC)
        if 'NC' in self.Df:
            self.N_Sp = self.Ax.scatter(self.Df.X, self.Df.Y, self.Df.Z, \
                c=self.Df.NC)
        # Label scatter point with number
        for i in range(len(self.Df)):
            self.Ax.text(self.Df.X[i], self.Df.Y[i], self.Df.Z[i], \
                '%s' % (str(i+1)), size=8, zorder=1, color='k')
        # Plot Surface
        self.SX, self.SY, self.SZ, self.Cmap, self.PC, self.PN, self.NC, \
            self.NN = self.Cal_Surface(self.Df)
        self.Surface = self.Ax.plot_surface(self.SX, self.SY, self.SZ, \
                        rstride=1, cstride=1, linewidth=0, antialiased=False, \
                        alpha=0.7, visible=False)
        # Plot Wireframe
        self.BFP = self.Show_Wireframe(BFX, BFY, BFZ, 'yellow')
        self.BFP.set_visible(self.menubar.show_BFP.get())
        self.P_WF = self.Show_Wireframe(S1X, S1Y, S1Z, 'green')
        self.PO_WF = self.Show_Wireframe(S1OX, S1OY, S1OZ, 'green')
        self.N_WF = self.Show_Wireframe(S2X, S2Y, S2Z, 'cyan')
        self.NO_WF = self.Show_Wireframe(S2OX, S2OY, S2OZ, 'cyan')
        # Calculate best result base on Max dist
        if (Pos in self.Df) and (Neg in self.Df):
            self.menubar.DisableSolu(Pos, False)
            self.menubar.DisableSolu(Neg, False)
            T, B = (True, False) if self.Df[Pos].max() < self.Df[Neg].max() \
                    else (False, True)
            self.menubar.TP.set(T)
            self.menubar.BP.set(B)
        else:
            T, B, DS = (True, False, Neg) if Pos in self.Df else (False, True, Pos)
            self.menubar.TP.set(T)
            self.menubar.BP.set(B)
            self.menubar.DisableSolu(DS, True)
        # Update canvas
        self.UpdateFig()
           
class Import(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        self.controller = controller

        BtnFrame = tk.Frame(self)
        ManualBtn = tk.Button(BtnFrame, text="Manual", relief='raised', \
                                command=lambda: controller.show_frame(Manual))
        ManualBtn.pack(side='left', fill='both', expand=True)
        ImportBtn = tk.Button(BtnFrame, text="Import", relief='sunken', \
                                command=lambda: controller.show_frame(Import))
        ImportBtn.pack(side='left', fill='both', expand=True)
        BtnFrame.pack(side='top', fill='x', padx=PX)

        tk.Label(self, text="Import file:").pack(side="top", anchor='w', \
            padx=PX)
        self.AddEntry = tk.Entry(self)
        self.AddEntry.pack(side="top", fill="x", padx=PX, pady=5)
        ImportBtn = tk.Button(self, text='Import', command=self.OpenFile)
        ImportBtn.pack(side="top", anchor="w", padx=PX)

    def OpenFile(self):
        InitDir = 'shell:MyDocumentsFolder' if (os.name == "nt") else \
                  (os.path.expanduser('~'))
        FileName = askopenfilename(initialdir = InitDir, \
            filetypes = (("csv files","*.csv"),("all files","*.*")))
        try:
            if FileName is not None:
                self.AddEntry.delete(0, "end")
                self.AddEntry.insert(0, FileName)
                self.pts = self.ReadFile(FileName)
                self.controller.Solve(self.pts)
        except:
            pass
    
    def ReadFile(self, FileName):
        pts = pd.read_table(FileName, sep=',', dtype='float64', \
            float_precision='high')
        return pts

class Manual(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)

        Name = [("Pts", 4), ("X", 8), ("Y", 8), ("Z", 8)]
        self.entries = []
        self.controller = controller

        BtnFrame = tk.Frame(self)
        ManualBtn = tk.Button(BtnFrame, text="Manual", relief='sunken', \
                                command=lambda: controller.show_frame(Manual))
        ManualBtn.pack(side='left', fill='both', expand=True)
        ImportBtn = tk.Button(BtnFrame, text="Import", relief='raised', \
                                command=lambda: controller.show_frame(Import))
        ImportBtn.pack(side='left', fill='both', expand=True)
        BtnFrame.pack(fill='x', padx=PX)

        tk.Label(self, text="Point Input:").pack(side="top", anchor='w', \
            padx=PX)
        self.PtsRow = tk.Frame(self)
        for name in Name:
            tk.Label(self.PtsRow, text=name[0], width=name[1]).pack(fill='x', \
                side='left', expand=True)
        self.PtsRow.pack(fill='both', padx=PX)
        for i in range(NP):
            self.PtsRow = tk.Frame(self)
            tk.Label(self.PtsRow, text=f"P{i + 1}", width=4)\
                .pack(fill='x', side='left', expand=True)
            self.entries.append([])
            for j in range(3):
                Entry = tk.Entry(self.PtsRow, width=8, validate='all', \
                    vcmd=(self.register(self.validate_float), '%S', '%d'))
                self.entries[i].append(Entry)
                self.entries[i][j].pack(fill='x', side='left', expand=True)
                self.entries[i][j].bind('<Return>', self.EntryEnter)
                self.entries[i][j].bind('<FocusOut>', self.LooseFocus)
            self.PtsRow.pack(fill='both', padx=PX)
            
        self.SolveFrame = tk.Frame(self, padx=PX, pady=PY)
        self.ClearBtn = tk.Button(self.SolveFrame, text='Clear', \
            command=self.ClearEntry)
        self.ClearBtn.pack(fill='both', side='left', expand=True)
        self.CalculateBtn = tk.Button(self.SolveFrame, text='Calculate', \
            command=self.GenPts)
        self.CalculateBtn.pack(fill='both', side='left', expand=True)
        self.SolveFrame.pack(fill='both', padx=PX, pady=(PY*10, PY))
   
    def ClearEntry(self):
        for i in range(13):
            for j in range(3):
                self.entries[i][j].config(bg = 'white')
                self.entries[i][j].delete(0, tk.END)
        self.controller.ClearAX(self)

    def validate_float(self, Char, Act):
        return True if ((Char in '0123456789.-+') or (Act == '0')) else False

    def EntryEnter(self, event):
        print('Value: %s' % (event.widget.get()))
        event.widget.tk_focusNext().focus()

    def IsFloat(self, Value):
        try:
            float(Value)
            return True
        except ValueError:
            return False

    def LooseFocus(self, event):
        widget=event.widget
        Value = widget.get()
        if Value != "":
            try:
                float(Value)
                widget['bg'] = 'green'
            except ValueError:
                print('Error')
                widget['bg'] = 'red'
        else:
            widget['bg'] = 'white'
    
    def GenPts(self):
        pts = []
        for i in range(NP):
            PX = self.entries[i][0].get()
            PY = self.entries[i][1].get()
            PZ = self.entries[i][2].get()
            if (self.IsFloat(PX) and self.IsFloat(PY) and self.IsFloat(PZ)):
                pts.append([PX, PY, PZ])
        if len(pts) >= 4:
            df = pd.DataFrame(data=np.array(pts), columns=['X', 'Y', 'Z'], \
                dtype='float64')
            self.controller.Solve(df)

def GetIcon():
    datafile = "app.ico"
    datafile = os.path.join(os.path.dirname(__file__), datafile) \
               if not hasattr(sys, "frozen") else \
               os.path.join(sys.prefix, datafile)
    return datafile

if __name__ == "__main__":
    app = MainApp()
    app.attributes('-alpha', 0.0)
    app.eval('tk::PlaceWindow . center')
    app.attributes('-alpha', 1.0)
    app.iconbitmap(default=GetIcon())
    app.mainloop()
