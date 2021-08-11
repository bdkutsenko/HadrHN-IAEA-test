#!/usr/bin/env python
#
# Plot XS pictures 
# A. Bagulya
# 23.03.2019
#
# Updated for phtotonuclear
# Kutsenko Bogdan
# 20.07.2021




import os, sys
from ROOT import TROOT, gROOT, gStyle, gPad, TCanvas, TStyle, TPad
from ROOT import TFile, TH1, TLegend, TGraphAsymmErrors
from array import array
from math import log10
import errno

def plotXS():
    gROOT.Reset()
    gROOT.SetStyle('Plain')
    gROOT.SetBatch(True)

    if len(sys.argv) < 3:
        print('''Usage: python plotXS.py file.root 0/1/2''')
        sys.exit()
    else:
        fname = sys.argv[1]
        fname2 = sys.argv[2]
        

        flag  = int(sys.argv[3])

    project = os.path.splitext(os.path.basename(fname))[0].split('_')[0]
    target = os.path.splitext(os.path.basename(fname))[0].split('_')[1]
    model_name = os.path.splitext(os.path.basename(fname))[0].split('_')[2]
    Z = os.path.splitext(os.path.basename(fname))[0].split('_')[3]
    A = os.path.splitext(os.path.basename(fname))[0].split('_')[4]

    project2 = os.path.splitext(os.path.basename(fname2))[0].split('_')[0]
    target2 = os.path.splitext(os.path.basename(fname2))[0].split('_')[1]
    model_name2 = os.path.splitext(os.path.basename(fname2))[0].split('_')[2]
    Z2 = os.path.splitext(os.path.basename(fname2))[0].split('_')[3]

    if model_name == 'kNS':
        model_name = 'NS'

    #HADSRC = readEnvVar('HADSRC')
    testpath = os.path.join('./','')
    dirName = os.path.join('plotXS/',str(target))

    
    try:
    # Create target Directory
        os.mkdir(dirName)
        print("Directory " , dirName ,  " Created ") 
    except OSError as e:
        if e.errno == errno.EEXIST:
            print("Directory " , dirName ,  " already exists")
        else:
            raise

        
    plotpath = os.path.join('./','plotXS/'+str(target))
    datapath = os.path.join('./','data')
    expDatapath = os.path.join('./','expData')

    
    if (int(flag) == 0 or int(flag == 1)) and project != 'gamma':
        print('Flag 0/1 (Inelastic cross-section from 0 to 200 MeV) implemented only for gamma')
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
        Name = 'Total cross sections gamma' + ' + '+target + A+ ' in mb, model = ' + model_name + ' + ' + model_name2
         
    crossampl = [0,   4,   6,   10,  6,  20,  24,  26,  29,  30,  30,  30,  70,  90,  50,  60,  70,  80,  65,  90, 100,  100,  110,  110,  120,  120,  130,  130,  140,  140, 150,  150,  160,   170,   170,  200,   210,   220,  230,  240,   250, 260,   270,  280,  290, 300, 310, 320, 330, 340, 350,   360, 370,   380, 390, 400,   410,   420,   430,   440,  450,   460,   470, 480,   490,   500,   510,   520,   530, 540,  550,   560, 600, 600,   610,   620,   630, 650, 660,   680, 700, 720,  740,   760,   780,   800,   820,   840,   860, 880, 900]
    lbound0=1.
    rbound0=100000.
    lbound1=1.
    rbound1=160.
    if project == 'neutron':
        hhh = gPad.DrawFrame(1, 1., 7, 10000. ,'cross sections in mb')
    elif project == 'kaon+':
        hhh = gPad.DrawFrame(1, 1., 7, 100. ,'cross sections in mb')
    elif project == 'gamma' and int(flag) == 2:
        hhh = gPad.DrawFrame(0.1, 0.01, 13., crossampl[int(Z)] ,'cross sections in mb')
    elif project == 'gamma' and int(flag) == 1 :
        hhh = gPad.DrawFrame(lbound0, 0.01, rbound0, 4*crossampl[int(Z)] ,'cross sections in mb')
    elif project == 'gamma' and int(flag) == 0:
        hhh = gPad.DrawFrame(lbound1, 0., rbound1, crossampl[int(Z)] ,'cross sections in mb')
    elif project == 'proton':
        hhh = gPad.DrawFrame(1, 1., 12, 1000. ,'cross sections in mb')
    else:
        hhh = gPad.DrawFrame(1, 1., 7, 1000. ,'cross sections in mb')

    hhh.SetTitle(Name)
    if int(flag) == 2:
        hhh.GetXaxis().SetTitle('log10(P/(MeV/c))')
    else:
        hhh.GetXaxis().SetTitle('P, MeV/c')
    
    hhh.GetYaxis().SetTitle('#sigma, mb')
    hhh.GetXaxis().SetTitleOffset(1.2)
    hhh.GetYaxis().SetTitleOffset(1.2)
    hhh.Draw('AXIS SAME')
    if int(flag) == 2:
        gPad.SetLogy(1)
    elif int(flag) == 1 :
        gPad.SetLogx(1)
        gPad.SetLogy(1)
    gPad.SetGrid(1, 1)
    gPad.Update()

    leg = TLegend(.45, .75, 0.97, .9)
    leg.SetTextSize(.03)
    leg.SetHeader('Cross sections') 

    rootflname = os.path.join(testpath, fname)
    print('File to open ', rootflname)
    f1 = TFile(rootflname)
    if int(flag) == 0:
       AddH1(f1, leg, 2, 3, 'h5', 'Photonuclear Lowest Energy region ' + model_name) 
    elif int(flag) == 1:
        AddH1(f1, leg, 2, 3, 'h4', 'Photonuclear Low Energy region ' + model_name)
    else:
        AddH1(f1, leg, 2, 3, 'h0', 'Inelastic ' + model_name)
    c1.Update()
    if project != 'gamma':
        AddH1(f1, leg, 3, 2, 'h1', 'Elastic ' + model_name)
        c1.Update()
        AddH1(f1, leg, 4, 2, 'h2', 'Total ' + model_name)
        c1.Update()
    rootflname = os.path.join(testpath, fname2)
    print('File to open ', rootflname)
    f2 = TFile(rootflname)
    if int(flag) == 0:
        AddH1(f2, leg, 4, 2, 'h5', 'Photonuclear Lowest Energy region ' + model_name2)
    elif int(flag) == 1:
        AddH1(f2, leg, 4, 2, 'h4', 'Photonuclear Low Energy region ' + model_name2)
    else:
        AddH1(f2, leg, 4, 2, 'h0', 'Inelastic ' + model_name2 )


    datafile = os.path.join(expDatapath, 'rpp2016-' + str(project) + str(target) +'_'+str(Z)+'_'+str(A)+ '_total.dat')
    if os.path.isfile(datafile) and flag == 2:
        AddGraph(datafile, leg, 'data Tot', 2, 0)
    datafile = os.path.join(expDatapath, str(project) + str(target) +'_'+str(Z)+'_'+str(A)+ '_total.dat')
    if os.path.isfile(datafile) and flag != 2:
        AddGraph(datafile, leg, 'data Tot', 2, 1)
        
    c1.Update()
    leg.Draw()
    if int(flag) == 2:
       fout = os.path.join(plotpath, 'A_' + model_name + '_' + model_name2 +'_'+str(Z)+'_'+str(A)+ '.png')
    elif int(flag) == 0:
       fout = os.path.join(plotpath, 'A_' + model_name + '_' + model_name2 + '_' +str(lbound1)+'_'+str(rbound1) + 'MeV_'+str(Z)+'_'+str(A)+'.png')
    elif int(flag) == 1:
        fout = os.path.join(plotpath, 'A_' + model_name + '_' + model_name2 +'_'+str(lbound0)+'_'+str(rbound0) +'MeV_'+str(Z)+'_'+str(A) +'.png')
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

def AddH1(ff, leg, idx, wdth, hh, title):
    h1 = ff.Get(hh)
    h1.SetLineColor(idx)
    h1.SetLineWidth(wdth)
    h1.Draw('HISTO SAME C')
    leg.AddEntry(h1, title, 'l')
    
def readFileEXFOR(ff):    
    openFile(ff)
    lines=infile.readlines()
    for x in lines:
        li=x.strip()
        if not li.startswith("#"):
            PLAB.append(float(x.split()[0]))
            errPLABl.append(0.)
            errPLABh.append(0.)
            SIG.append(float(x.split()[1]))
            errSIGl.append(float(x.split()[2]))
            errSIGh.append(float(x.split()[2]))

    infile.close()

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

def AddGraph(ff, leg, title, idx, fileType):
     openFile(ff)
     if fileType == 0:
         readFile(ff)
     if fileType == 1:
         readFileEXFOR(ff)
         
     globals()['gr' + str(idx)] = TGraphAsymmErrors(len(PLAB),PLAB,SIG,errPLABl,errPLABh,errSIGl,errSIGh)
     globals()['gr' + str(idx)].SetMarkerStyle(19+idx)
     globals()['gr' + str(idx)].SetMarkerColor(1)
     globals()['gr' + str(idx)].SetMarkerSize(0.5)
     globals()['gr' + str(idx)].Draw('SAME P')
     leg.AddEntry(globals()['gr' + str(idx)],title,'ep')

     del PLAB[:]
     del errPLABl[:]
     del errPLABh[:]
     del SIG[:]
     del errSIGl[:]
     del errSIGh[:]

###______________________________
if __name__ == '__main__':
    plotXS()
