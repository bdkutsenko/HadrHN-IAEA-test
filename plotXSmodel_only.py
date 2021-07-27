#!/usr/bin/env python
#
# Plot XS pictures 
# A. Bagulya
# 23.03.2019
#

import os, sys
from ROOT import TROOT, gROOT, gStyle, gPad, TCanvas, TStyle, TPad
from ROOT import TFile, TH1, TLegend, TGraphAsymmErrors
from array import array
from math import log10

def plotXS():
    gROOT.Reset()
    gROOT.SetStyle('Plain')
    gROOT.SetBatch(True)

    #HADSRC = readEnvVar('HADSRC')
    testpath = os.path.join('./','')
    plotpath = os.path.join('./','plotXS')
    datapath = os.path.join('./','data')

    if len(sys.argv) < 3:
        print('''Usage: python plotXS.py file.root 0/1''')
        sys.exit()
    else:
        fname = sys.argv[1]
        flag  = int(sys.argv[2])

    project = os.path.splitext(os.path.basename(fname))[0].split('_')[0]
    target = os.path.splitext(os.path.basename(fname))[0].split('_')[1]
    model_name = os.path.splitext(os.path.basename(fname))[0].split('_')[2]
    Z = os.path.splitext(os.path.basename(fname))[0].split('_')[3]
    A = os.path.splitext(os.path.basename(fname))[0].split('_')[4]
    if model_name == 'kNS':
        model_name = 'NS'
        
    if int(flag) == 1 and project != 'gamma':
        print('Flag 1 (Inelastic cross-section from 0 to 200 MeV) implemented only for gamma')
        exit()
        
    part = {'proton':'p', 'neutron':'n', 'pi+':'#pi^{+}', 'pi-':'#pi^{-}', 'kaon+':'K^{+}', 'kaon-':'K^{-}', 'gamma':'#g' , 'deuteron':'d', 'triton':'t', 'He3':'he3', 'alpha':'alpha'}

    global PLAB, SIG, errPLABl, errPLABh, errSIGl, errSIGh

    PLAB, SIG = array ( 'd' , [] ), array ( 'd' , [] )
    errPLABl, errPLABh, errSIGl, errSIGh = array ( 'd' , [] ), array ( 'd' , [] ), array ( 'd' , [] ), array ( 'd' , [] )

    c1 = TCanvas('c1', 'c1',200,100,900,600)
    gStyle.SetOptStat(0)
    c1.SetFillColor(0)
    c1.SetBorderMode(0)
    c1.SetBorderSize(0)
    c1.SetFrameBorderMode(0)

    hasDataEl = 0;
    hasDataIn = 0;
    if project == 'gamma':
        hasDataEl = 0;
    
    if project == 'alpha':
        hasDataEl = 0
        hasDataIn = 0;

    Name = 'Total, elastic and inelastic cross sections ' + project + ' + ' + target + ' in mb, model = ' + model_name
    if project == 'gamma':
        Name = 'Total cross sections gamma' + ' + '+target + A+ ' in mb, model = ' + model_name
        
    crossampl = [2,   1,   6,   10,  6,  20,  24,  26,  29,  30,  30,  30,  45,  45,  50,  60,  70,  80,  65,  90, 100,  100,  110,  110,  120,  120,  130,  130,  140,  140, 150,  150,  160,   170,   170,  200,   210,   220,  230,  240,   250, 260,   270,  280,  290, 300, 310, 320, 330, 340, 350,   360, 370,   380, 390, 400,   410,   420,   430,   440,  450,   460,   470, 480,   490,   500,   510,   520,   530, 540,  550,   560, 600, 600,   610,   620,   630, 650, 660,   680, 700, 720,  740,   760,   780,   800,   820,   840,   860, 880, 900]
    

    lbound = 1.;
    rbound = 100000;
    if project == 'neutron':
        hhh = gPad.DrawFrame(1, 1., 7, 10000. ,'cross sections in mb')
    elif project == 'kaon+':
        hhh = gPad.DrawFrame(1, 1., 7, 100. ,'cross sections in mb')
    elif project == 'gamma' and int(flag) == 0:
        hhh = gPad.DrawFrame(1., 1., 7., crossampl[int(Z)] ,'cross sections in mb')
    elif project == 'gamma' and int(flag) == 1:
        hhh = gPad.DrawFrame(lbound, 0.01, rbound, 4*crossampl[int(Z)] ,'cross sections in mb')
    elif project == 'proton':
        hhh = gPad.DrawFrame(1, 1., 12, 1000. ,'cross sections in mb')
    else:
        hhh = gPad.DrawFrame(1, 1., 7, 1000. ,'cross sections in mb')

    hhh.SetTitle(Name)
    if int(flag) == 0:
        hhh.GetXaxis().SetTitle('log10(P/(MeV/c))')
    else:
        hhh.GetXaxis().SetTitle('P, MeV/c')
    
    hhh.GetYaxis().SetTitle('#sigma, mb')
    hhh.GetXaxis().SetTitleOffset(1.2)
    hhh.GetYaxis().SetTitleOffset(1.2)
    hhh.Draw('AXIS SAME')
    if int(flag) == 0:
        gPad.SetLogy(1)
    else :
        gPad.SetLogx(1)
        gPad.SetLogy(1)
    gPad.Update()

    leg = TLegend(.77, .75, 0.97, .95)
    leg.SetTextSize(.03)
    leg.SetHeader('Cross sections') 

    rootflname = os.path.join(testpath, fname)
    print('File to open ', rootflname)
    f1 = TFile(rootflname)
    if int(flag) == 1:
        AddH1(f1, leg, 2, 'h4', 'Photonuclear Low Energy region ' + model_name)
    else:
        AddH1(f1, leg, 2, 'h0', 'Inelastic ' + model_name)
    c1.Update()
    if project != 'gamma':
        AddH1(f1, leg, 3, 'h1', 'Elastic ' + model_name)
        c1.Update()
        AddH1(f1, leg, 4, 'h2', 'Total ' + model_name)
        c1.Update()

    if int(flag) == 0:
       fout = os.path.join(plotpath, 'A_' + os.path.splitext(os.path.basename(fname))[0]) + '.png'
    else:
       fout = os.path.join(plotpath, 'A_' + os.path.splitext(os.path.basename(fname))[0])+'_'+str(lbound)+'_'+str(rbound) + 'MeV.png'
        
    c1.Print(fout)
    c1.Close()

def readEnvVar(envar):
    if os.environ.get(envar) is None:
        print('''No environment variable %s''' % envar)
        sys.exit()
    else:
        return os.environ.get(envar)

def openFile(fileName):
    global infile
    try:
        infile = open(fileName, 'r')
    except IOError:
        print('Input file <',fileName, '> does not exist! Exit')
        sys.exit(2)

def AddH1(ff, leg, idx, hh, title):
    h1 = gROOT.FindObject(hh)
    h1.SetLineColor(idx)
    h1.SetLineWidth(2)
    h1.Draw('HISTO SAME C')
    leg.AddEntry(h1, title, 'l')

def readFile(ff):
    openFile(ff)
    lines = infile.readlines()
    num = int(lines[10])
    for line in lines[11:11+num]:
        plab = 1000.*float(line.split()[1])
        pmin = 1000.*float(line.split()[2])
        pmax = 1000.*float(line.split()[3])
        sig = float(line.split()[4])
        errsigh = float(line.split()[5])
        errsigl = float(line.split()[6])

        PLAB.append(log10(plab))
        errPLABl.append(0.)
        errPLABh.append(0.)
        SIG.append(sig)
        errSIGl.append(errsigh)
        errSIGh.append(errsigl)
    infile.close()

def AddGraph(ff, leg, title, idx):
     openFile(ff)
     readFile(ff)
     globals()['gr' + str(idx)] = TGraphAsymmErrors(len(PLAB),PLAB,SIG,errPLABl,errPLABh,errSIGl,errSIGh)
     globals()['gr' + str(idx)].SetMarkerStyle(19+idx)
     globals()['gr' + str(idx)].SetMarkerColor(1)
     globals()['gr' + str(idx)].SetMarkerSize(0.4)
     globals()['gr' + str(idx)].Draw('SAME P')
     leg.AddEntry(globals()['gr' + str(idx)],title,'p')

     del PLAB[:]
     del errPLABl[:]
     del errPLABh[:]
     del SIG[:]
     del errSIGl[:]
     del errSIGh[:]

###______________________________
if __name__ == '__main__':
    plotXS()
