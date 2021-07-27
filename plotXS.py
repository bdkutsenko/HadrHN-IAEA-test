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

    HADSRC = readEnvVar('HADSRC')
    testpath = os.path.join(HADSRC,'HadrHN_class_test','test')
    plotpath = os.path.join(HADSRC,'HadrHN_class_test','plotXS')
    datapath = os.path.join(HADSRC,'HadrHN_class_test','data')

    if len(sys.argv) < 3:
        print '''Usage: python plotXS.py file.root 0/1'''
        sys.exit()
    else:
        fname = sys.argv[1]
        flag  = sys.argv[2]

    project = os.path.splitext(os.path.basename(fname))[0].split('_')[0]
    target = os.path.splitext(os.path.basename(fname))[0].split('_')[1]
    model_name = os.path.splitext(os.path.basename(fname))[0].split('_')[2]
    if model_name == 'kNS':
        model_name = 'NS'

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

    hasDataEl = 1;
    hasDataIn = 1;
    if project == 'gamma':
        hasDataEl = 0;
    
    if project == 'alpha':
        hasDataEl = 0
        hasDataIn = 0;

    Name = 'Total, elastic and inelastic cross sections ' + part[project] + ' + ' + part[target] + ' in mb, model = ' + model_name
    if project == 'gamma':
        Name = 'Total cross sections gamma' + part[target] + ' in mb, model = ' + model_name

    if project == 'neutron':
        hhh = gPad.DrawFrame(1, 1., 7, 10000. ,'cross sections in mb')
    elif project == 'kaon+':
        hhh = gPad.DrawFrame(1, 1., 7, 100. ,'cross sections in mb')
    elif project == 'gamma':
        hhh = gPad.DrawFrame(1, 0.01, 7, 1. ,'cross sections in mb')
    elif project == 'proton' and flag == '1':
        hhh = gPad.DrawFrame(1, 1., 12, 1000. ,'cross sections in mb')
    else:
        hhh = gPad.DrawFrame(1, 1., 7, 1000. ,'cross sections in mb')

    hhh.SetTitle(Name)
    hhh.GetXaxis().SetTitle('log10(P/(MeV/c))')
    hhh.GetYaxis().SetTitle('#sigma, mb')
    hhh.GetXaxis().SetTitleOffset(1.2)
    hhh.GetYaxis().SetTitleOffset(1.2)
    hhh.Draw('AXIS SAME')

    gPad.SetLogy(1)
    gPad.Update()

    leg = TLegend(.77, .75, 0.97, .95)
    leg.SetTextSize(.03)
    leg.SetHeader('Cross sections') 

    rootflname = os.path.join(testpath, fname)
    print 'File to open ', rootflname
    f1 = TFile(rootflname)

    AddH1(f1, leg, 2, 'h0', 'Inelastic ' + model_name)
    c1.Update()
    if project != 'gamma':
        AddH1(f1, leg, 3, 'h1', 'Elastic ' + model_name)
        c1.Update()
        AddH1(f1, leg, 4, 'h2', 'Total ' + model_name)
        c1.Update()
#    AddH1(f1, leg, 5, 'h3', 'EM Elastic')
#    c1.Update()

    if project == 'pi+':
        datafile = os.path.join(datapath, 'rpp2016-pip' + part[target] + '_elastic.dat')
    elif project == 'pi-':
        datafile = os.path.join(datapath, 'rpp2016-pim' + part[target] + '_elastic.dat')
    elif project == 'kaon+':
        datafile = os.path.join(datapath, 'rpp2016-kp' + part[target] + '_elastic.dat')
    elif project == 'kaon-':
        datafile = os.path.join(datapath, 'rpp2016-km' + part[target] + '_elastic.dat')
    else:
        datafile = os.path.join(datapath, 'rpp2016-' + part[project] + part[target] + '_elastic.dat')

##    if project != 'gamma':
    if hasDataEl == 1:
        AddGraph(datafile, leg, 'data El PDG16', 1)

    if project == 'pi+':
        datafile = os.path.join(datapath, 'rpp2016-pip' + part[target] + '_total.dat')
    elif project == 'pi-':
        datafile = os.path.join(datapath, 'rpp2016-pim' + part[target] + '_total.dat')
    elif project == 'kaon+':
        datafile = os.path.join(datapath, 'rpp2016-kp' + part[target] + '_total.dat')
    elif project == 'kaon-':
        datafile = os.path.join(datapath, 'rpp2016-km' + part[target] + '_total.dat')
    elif project == 'gamma':
        datafile = os.path.join(datapath, 'rpp2016-gamma' + part[target] + '_total.dat')
    else:
        datafile = os.path.join(datapath, 'rpp2016-' + part[project] + part[target] + '_total.dat')

    if hasDataIn == 1:
        AddGraph(datafile, leg, 'data Tot PDG16', 2)

    leg.Draw('SAME')
    c1.Update() 

    if flag == '1':
        fout = os.path.join(plotpath, 'A_' + os.path.splitext(os.path.basename(fname))[0]) + '_cosmic.png' 
    else:
        fout = os.path.join(plotpath, 'A_' + os.path.splitext(os.path.basename(fname))[0]) + '.png'
    c1.Print(fout)
    c1.Close()

def readEnvVar(envar):
    if os.environ.get(envar) is None:
        print '''No environment variable %s''' % envar
        sys.exit()
    else:
        return os.environ.get(envar)

def openFile(fileName):
    global infile
    try:
        infile = open(fileName, 'r')
    except IOError:
        print 'Input file <',fileName, '> does not exist! Exit'
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
     globals()['gr' + str(idx)].SetMarkerSize(0.8)
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
