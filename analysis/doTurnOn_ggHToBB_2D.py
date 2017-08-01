import ROOT
import copy
from math import *
from array import array

#import tdrstyle
import CMS_lumi, tdrstyle

#set the tdr style
#tdrstyle.setTDRStyle()

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.lumi_7TeV = "4.8 fb^{-1}"
CMS_lumi.lumi_13TeV = ""
#CMS_lumi.lumi_13TeV = "7.1 fb^{-1} (RunG)"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
#CMS_lumi.lumi_sqrtS = "RunG 13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPeriod =0
iPos = 11

runName      = "2017"
fitName      = ""

ffName = "fittedFunctions3_.h"
f = open( ffName , 'w')

ROOT.gROOT.LoadMacro("tdrstyleTrigger.C")
ROOT.gROOT.LoadMacro("test.h")
ROOT.setTDRStyle()

minTurnOn_funct   = 0
maxTurnOn_funct   = 1.05

minRatio    = 0.8
maxRatio    = 1.2

functionMin = 80
functionMax = 500

ped = "function"

title = "aaa"

var           ="Jet_pt[3]"
trigger       ="ntrgObjects_hltQuadPFCentralJetLooseID45>=4"
binning       =(40,0,120)
preselection  ="1"

maxev=100000000

def getTitle(fileName):
    file_ = ROOT.TFile.Open(fileName)
    file_.cd()
    name = file_.GetListOfKeys().First().GetName()
    canvas =  file_.Get(name)
    pad =  canvas.GetPrimitive("unten")
    title =  pad.GetPrimitive("Ratio").GetXaxis().GetTitle()
    return title
    
def DivideTGraph(num,den):
    Ns_den   = den.GetN()
    Xs_den   = den.GetX()
    Ys_den   = den.GetY()
    EXLs_den = den.GetEXlow()
    EXHs_den = den.GetEXhigh()
    EYLs_den = den.GetEYlow()
    EYHs_den = den.GetEYhigh()

    print "den.GetN()",den.GetN()
    print "num.GetN()",num.GetN()

    Ys_num   = num.GetY()
    EYLs_num = num.GetEYlow()
    EYHs_num = num.GetEYhigh()

    print "DivideTGraph: new"
    
    bins   = [ i for i in range( Ns_den) if Ys_den[i]>0]

    print "Xs_den",Xs_den
    Xs_new   = [ Xs_den[i] for i in bins]
    print "Xs_new",Xs_new
    Ys_new   = [ Ys_num[i]/(Ys_den[i]) for i in bins]
    EXLs_new = [ EXLs_den[i] for i in bins]
    EXHs_new = [ EXHs_den[i] for i in bins]
    [ EYLs_num[i] for i in bins]
    [ ((EYLs_num[i]/(Ys_num[i]+1E-3))**2) for i in bins]
    [ Ys_new[i]*sqrt((EYLs_num[i]/(Ys_num[i]+1E-3))**2) for i in bins]
    EYLs_new = [ Ys_new[i]*sqrt((EYLs_num[i]/(Ys_num[i]+1E-3))**2+(EYHs_den[i]/(Ys_den[i]+1E-3))**2) for i in bins]
    EYHs_new = [ Ys_new[i]*sqrt((EYHs_num[i]/(Ys_num[i]+1E-3))**2+(EYLs_den[i]/(Ys_den[i]+1E-3))**2) for i in bins]

    print "DivideTGraph: len"

    n = len(Xs_new)

    print "DivideTGraph: array"

    Xs_new = array.array('f',Xs_new)
    Ys_new = array.array('f',Ys_new)
    EXLs_new = array.array('f',EXLs_new)
    EXHs_new = array.array('f',EXHs_new)
    EYLs_new = array.array('f',EYLs_new)
    EYHs_new = array.array('f',EYHs_new)

    print "DivideTGraph: ratio"

    ratio = ROOT.TGraphAsymmErrors(n, Xs_new, Ys_new, EXLs_new, EXHs_new, EYLs_new, EYHs_new)
    print "DivideTGraph: done"

    return ratio

def makeHistos():
    tempNum = ROOT.TH1F("num", "", binning[0], binning[1], binning[2])
    tempDen = ROOT.TH1F("den", "", binning[0], binning[1], binning[2])
    tempNum.SetDirectory(0)
    tempDen.SetDirectory(0)
    for (i,weight) in fileNames:
	file_ = ROOT.TFile(i)
	NEvents = file_.Get("Count").GetBinContent(1)
	tree = ROOT.TChain("tree")
	tree.Add(i)
	print "fileName=",i
	print "var=",var
	print "trigger=",trigger
	print "preselection=",preselection
	print "binning=",str(binning)
	print "weight=",weight
	tree.Draw(var+">>num"+str(binning),str("("+preselection+"&&"+trigger+")*"+weight+"/"+NEvents),"",maxev)
	print "Draw:\t",var+">>num"+str(binning),str("("+preselection+"&&"+trigger+")*"+weight+"/"+NEvents)
	num = ROOT.gDirectory.Get("num")
	num = copy.copy(num)
	tree.Draw(var+">>den"+str(binning),str("("+preselection+")*"+weight+"/"+NEvents),"",maxev)
	print "Draw:\t",var+">>den"+str(binning),str("("+preselection+")*"+weight+"/"+NEvents)
	den = ROOT.gDirectory.Get("den")
	den = copy.copy(den)
	print "num,den = ",num.Integral(),den.Integral()
	tempNum.Add(num)
	tempDen.Add(den)
    print "TOTAL num,den = ",tempNum.Integral(),tempDen.Integral()
    return tempNum,tempDen

def make2DHistos():
    tempNum = ROOT.TH2F("num", "", binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
    tempDen = ROOT.TH2F("den", "", binning[0], binning[1], binning[2], binning[3], binning[4], binning[5])
    tempNum.SetDirectory(0)
    tempDen.SetDirectory(0)
    for (i,weight) in fileNames:
        file_ = ROOT.TFile(i)
        NEvents = str(file_.Get("Count").GetBinContent(1))
        tree = ROOT.TChain("tree")
        tree.Add(i)
        print "fileName=",i
        print "var=",var
	print "var2=",var2
        print "trigger=",trigger
        print "preselection=",preselection
        print "binning=",str(binning)
        print "weight=",weight
        tree.Draw(var+":"+var2+">>num"+str(binning),str("("+preselection+"&&"+trigger+")"),"",maxev)
        print "Draw:\t",var+":"+var2+">>num"+str(binning),str("("+preselection+"&&"+trigger+")")
        num = ROOT.gDirectory.Get("num")
        num = copy.copy(num)
        tree.Draw(var+":"+var2+">>den"+str(binning),str("("+preselection+")"),"",maxev)
        print "Draw:\t",var+":"+var2+">>den"+str(binning),str("("+preselection+")")
        den = ROOT.gDirectory.Get("den")
        den = copy.copy(den)
        print "num,den = ",num.Integral(),den.Integral()
        tempNum.Add(num)
        tempDen.Add(den)
    print "TOTAL num,den = ",tempNum.Integral(),tempDen.Integral()
    return tempNum,tempDen
    
    
def getMCAndData(fileName):
    file_ = ROOT.TFile.Open(fileName)
    file_.cd()
    name = file_.GetListOfKeys().First().GetName()
    canvas =  file_.Get(name)
    pad =  canvas.GetPrimitive("oben")
    mystack =  pad.GetPrimitive(name)
    MC_tmp =  mystack.GetStack().Last()
    data_tmp =  pad.GetPrimitive("noData")
    MC =  MC_tmp.Clone("MC")
    data =  data_tmp.Clone("data")
    MC.SetTitle("MC")
    data.SetTitle("data")
    MC.GetXaxis().SetTitle(title)
    data.GetXaxis().SetTitle(title)
    MC.SetMarkerStyle(20)
    data.SetMarkerStyle(20)
    MC.SetMarkerColor(ROOT.kBlack)
    data.SetMarkerColor(ROOT.kBlack)
    MC = copy.copy(MC)
    data = copy.copy(data)
    file_.Close()
    return MC,data

def doRatio(num, den, option=""):
    mratio =  den.Clone("mratio")
    mratio.SetTitle("Ratio")
    mratio.Reset()
    if option is "b":
        mratio.Divide(num,den,1,1,"b")
    else:
        mratio.Divide(num,den)
    for i in range(mratio.GetNbinsX()):
	for j in range(mratio.GetNbinsY()):
		if mratio.GetBinContent(i,j) > 1:
			print "Bin: ", i,j,mratio.GetBinContent(i,j)
			print "Num Bin: ", i,j,num.GetBinContent(i,j)
			print "Den Bin: ", i,j,den.GetBinContent(i,j)
    return mratio
#        mratio = ROOT.TGraphAsymmErrors()
#        mratio.SetMarkerColor(ROOT.kBlack)
#        mratio.SetLineColor(ROOT.kBlack)
#        mratio.SetName("ratio")
#        mratio.GetXaxis().SetTitle(title)
#        mratio = histo.Clone(triggerName+"_eff")
#        mratio.Divide(histo,inclusive,1.,1.,"B")
#        mratio.Divide(histo,inclusive,1.,1.,"cl=0.683 b(1,1) mode")
#        print num.GetNbinsX(),num.GetXaxis().GetXmin(),num.GetXaxis().GetXmax()
#        print den.GetNbinsX(),den.GetXaxis().GetXmin(),den.GetXaxis().GetXmax()
#        for i in range(num.GetNbinsX()):print num.GetBinLowEdge(i),
#        print ""
#        for i in range(den.GetNbinsX()): print den.GetBinLowEdge(i),
#        print ""

#        for i in range(den.GetNbinsX()+2):
#            if den.GetBinContent(i)<=0:
#                den.SetBinError(i,1.)
#                den.SetBinContent(i,0)
#                num.SetBinContent(i,0)
#        for i in range(num.GetNbinsX()+2):
#            if num.GetBinContent(i)<=0:
#                num.SetBinError(i,10.)
#                num.SetBinContent(i,1.E-7)
#        for i in range(den.GetNbinsX()+2):
#            if num.GetBinContent(i)>den.GetBinContent(i):
#                print "WARNING!"
#                print num.GetBinContent(i),den.GetBinContent(i)
#                num.SetBinContent(i,den.GetBinContent(i))
#            num.SetBinContent(i,num.GetBinContent(i))
#            den.SetBinContent(i,den.GetBinContent(i))

#        mratio.Divide(num,den,"cl=0.683 b(1,1) mode")
#        print "End ratio. bins:",mratio.GetN()," num:",num.GetNbinsX()," den:",den.GetNbinsX()
#        return mratio

def confidenceInterval(graph, function):
    fit = function.Clone("fit")
    fitUp = function.Clone("fitUp")
    fitUp.SetLineColor(ROOT.kRed)
    fitUp.SetLineStyle(2)
    fitDown = function.Clone("fitDown")
    fitDown.SetLineStyle(2) 
    print "Fit1"
    fit.FixParameter(4,0)
    fit.FixParameter(5,0)
    graph.Fit(fit,"","",fit.GetXmin(),fit.GetXmax())
#    fit.ReleaseParameter(4)
#    fit.ReleaseParameter(5)
    print "Fit2"
    graph.Fit(fit,"","",fit.GetXmin(),fit.GetXmax()) #was WW
    print "Fit3"
    graph.Fit(fit,"","",fit.GetXmin(),fit.GetXmax())
    parameters = [0]*function.GetNpar()
    for i in range(len(parameters)):
        parameters[i] = fit.GetParameter(i)

    parametersUp = parameters[:]
    parametersDown = parameters[:]

    looseRange=10.
    tightRange=10.

    print "Up/down fit"

    ## FixParameters
    for i in range(len(parameters)):
#        fit.ReleaseParameter(i)
#        function.SetParLimits(4,0,10)

#        if i in [0]: #  x0 can go down
#            pass
#        elif i in [3]: # global efficiency can go up
#            pass
        if i in [0,2]: #  x0 can go down
            pass
        else: # fix the other parameters
            fit.FixParameter( i, parameters[i] )
    fitResult = graph.Fit(fit,"SEV0","",fit.GetXmin(),fit.GetXmax())
    for i in range(len(parameters)):

        parameters[i] = fit.GetParameter(i)
        if i in [0]: #  x0 can go down
            parametersUp[i] = fit.GetParameter(i)   + fitResult.LowerError(i)
            parametersDown[i] = fit.GetParameter(i) + fitResult.UpperError(i)
    #    elif i in [1]: # check-me!
    #        parametersUp[i] = fit.GetParameter(i)   + fitResult.UpperError(i)
    #        parametersDown[i] = fit.GetParameter(i) + fitResult.LowerError(i)
    #    elif i in [2]: # check-me!
    #        parametersUp[i] = fit.GetParameter(i)   + min(fitResult.LowerError(i),-0.02)
    #        parametersDown[i] = fit.GetParameter(i) + max(fitResult.UpperError(i),+0.02)
        elif i in [3]: # global efficiency can go up
            parametersUp[i] = fit.GetParameter(i)   + min(fitResult.LowerError(i),-0.02)
            parametersDown[i] = fit.GetParameter(i) + max(fitResult.UpperError(i),+0.02)
        else: # fix the other parameters
            parametersUp[i] = fit.GetParameter(i)
            parametersDown[i] = fit.GetParameter(i)

    if (fit.GetParameter(3)<0): ##if [3]<0, I have to exchange [0],[1] min/max
        for i in [0,1]:
            print "CHANGE"
            print parametersUp[i],parametersDown[i]
            tmp = parametersUp[i]
            parametersUp[i] = parametersDown[i]
            parametersDown[i] = tmp
            print parametersUp[i],parametersDown[i]

    for i in range(len(parameters)):
        print "i=,",i,"\t",parameters[i],"\t",parametersUp[i],"\t",parametersDown[i]
        fit.SetParameter(i,parameters[i])
        fitUp.SetParameter(i,parametersUp[i])
        fitDown.SetParameter(i,parametersDown[i])

    return fit,fitUp,fitDown

def getPavetext():
#  addInfo = ROOT.TPaveText(0.17,0.835,0.25,0.935,"NDC")
  addInfo = ROOT.TPaveText(0.17,0.735,0.25,0.835,"NDC")
  addInfo.SetFillColor(0)
  addInfo.SetLineColor(0)
  addInfo.SetFillStyle(0)
  addInfo.SetBorderSize(0)
  addInfo.SetTextFont(42)
  addInfo.SetTextSize(0.030)
  addInfo.SetTextAlign(12)
  return addInfo


def doPlots():

    num,den =  make2DHistos()
    turnOn = doRatio(num,den,"b")

    #function = ROOT.TF1("turnonPt","1-(0.5-0.5*erf( (x-[0])/[1]))*([3])-[2] ",functionMin,functionMax)
#    function = ROOT.TF1("turnonPt","1-(0.5-0.5*TMath::Erf( (x-[0])/[1]))*[3]-[2] ",functionMin,functionMax)
#    function = ROOT.TF1("turnonPt","1-(0.5-0.5*TMath::Erf( (x-[0])/[1]))*[3]-[2] ",functionMin,functionMax)
#    function = ROOT.TF1("turnonPt","1-(0.5-0.5*TMath::Erf( (x-[0])/([1])*(x-[0]>[5]) + ((x-[0])/([1]+[4]) + [5]*(1/[1]+[4]-1/[1]))*(x-[0]<=[5]) ))*[3]-[2] ",functionMin,functionMax)
#    function = ROOT.TF1("turnonPt","1-(0.5-0.5*TMath::Erf( (x-[0])*(1+[4]*x**2)/([1]+[5]*x**2)))*([3])-[2] ",functionMin,functionMax)
#    function = ROOT.TF1("turnonPt","(0.5+0.5*TMath::Erf( (x-[0])*(x-[0]>[5])/[1] + (x-[0])*(x-[0]<[5])/[2] + [5]*(1/[1]-1/[2])*(x-[0]<[5]) ))*[4]+[3] ",functionMin,functionMax)
#    function = ROOT.TF1("turnonPt","1-(0.5-0.5*TMath::Erf( (x-[0])/[1]))*([3])-(0.5-0.5*TMath::Erf( (x-[4])/[5]))*([6])-[2]  ",functionMin,functionMax)

#    function = ROOT.TF1("turnonPt","TMath::Erf( (x-[0])/[1]*(([4]+[5]*x)) )*[2]-[3]",functionMin,functionMax)

#    print "Using:",parametersTurnOn_funct
#    function.SetParameters(*parametersTurnOn_funct)
#    function.SetParLimits(1,0,100)
#    function.SetParLimits(2,0,1)
#    function.SetParLimits(3,-1,2)
#    function.SetParLimits(4,0,1)
#    function.SetParLimits(5,0,1)
#    function.SetLineWidth(2)


#    TurnOn_funct = function.Clone("TurnOn_funct")

    c1 = ROOT.TCanvas("c1","",800,600)
    c1.cd()
#    TurnOn_funct,TurnOn_functUp,TurnOn_functDown = confidenceInterval(turnOn,TurnOn_funct)
#
#    turnOn.SetMaximum(maxTurnOn_funct)
#    turnOn.SetMinimum(minTurnOn_funct)
    turnOn.GetXaxis().SetTitle(xtitle)
    turnOn.GetYaxis().SetTitle(ytitle)
    turnOn.GetZaxis().SetTitle("Efficiency")
#    turnOn.GetZaxis().SetRangeUser(0.9,1)

    turnOn.Draw("COLZ")
#    TurnOn_funct.Draw("same")
#    TurnOn_functUp.Draw("same")
#    TurnOn_functDown.Draw("same")

    addInfo = getPavetext()
#    addInfo.AddText("p_{T} > 400 GeV,")
    addInfo.AddText("p_{T} > 400 GeV, |#eta| < 2.5")
    addInfo.AddText("M_{SD} > 40 GeV")
    addInfo.Draw("same")

    #c1.SaveAs("turnOn_"+ped+".C")
    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
    c1.SaveAs("turnOn2D_"+ped+"_"+runName+fitName+".pdf")
    #c1.SaveAs("turnOn_"+ped+".root")
#    f.write(('TF1* %s = new TF1("%s","'%(ped,ped)           + str(TurnOn_funct.GetExpFormula("P"))+'");\n'))
#    f.write(('TF1* %sUp = new TF1("%sUp","'%(ped,ped)       + str(TurnOn_functUp.GetExpFormula("P"))+'");\n'))
#    f.write(('TF1* %sDown = new TF1("%sDown","'%(ped,ped)   + str(TurnOn_functDown.GetExpFormula("P"))+'");\n'))
    
	

ROOT.gROOT.SetBatch()
#ROOT.gROOT.Reset()
#ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptFit(0)

#stops = [ 0.0, 1.0]
#red =   [ 1.0, 0.3]
#green = [ 1.0, 0.3]
#blue =  [ 1.0, 1.0]

#s = array('d', stops)
#r = array('d', red)
#g = array('d', green)
#b = array('d', blue)

#npoints = len(s)
#ROOT.TColor.CreateGradientColorTable(npoints, s, r, g, b, 999)

minRatio    = 0.5
maxRatio    = 1.5

#fileName    ="../GluGluHToBB_partial.root"
fileNames    = [("../GluGluHToBB.root","1.")]
#fileName     = "root://cmseos.fnal.gov//store/user/lpchbb/HeppyNtuples/V23/SingleMuon" + runName + ".root"
#fileName    = "ZvvHighPt_V20_SingleMuon.root"
#fileName    = "ZvvHighPt_V20_TT_TuneCUETP8M1_13TeV-powheg-pythia8.root"
#fileData    = "/gpfs/ddn/srm/cms/store/user/arizzi/VHBBHeppyV20/SingleMuon/VHBB_HEPPY_V20_SingleMuon__Run2015D-16Dec2015-v1/160210_081323/0000/tree*.root"

#preselection = "HLT_BIT_HLT_IsoMu24_v  && Vtype==2  && CSVsorted[2]>0.6"
preselection = "AK8PuppiJets_pt[0]>400&&abs(AK8PuppiJets_eta[0])<2.5&&AK8SoftDropMass[0]>40"
#preselection = "AK8PFJetsCHS_pt[0]>350&&abs(AK8PFJetsCHS_eta[0])<2.5&&AK8SoftDropMass[0]>40"

#&& diHiggs"#s && LPt_mass > 400"
parametersTurnOn_funct = ()
#################### HLT_AK8PFJet550_PFAK8BTagCSV_v2 #########################
parametersTurnOn_funct = (200,100,0.01,1,1E-3,1E-3)

Nbins       = 8#50 
functionMin = 0#100
functionMax = 1#200
NbinsY      = 6
yMin        = 400
yMax        = 1000
var             = "AK8PFJetsCHS_DoubleSV[0]"
var2            = "AK8PuppiJets_pt[0]"
#var             = "Jet_pt[1]+Jet_pt[2]+Jet_pt[3]"
#trigger         = "HLT_AK8PFJet360_TrimMass30_v9==1"
trigger         = "HLT_AK8PFJet360_TrimMass30_v9==1||((l1Jets_pt[0]>170 && AK8CaloJetsIDPassed_pt[0]>290 && abs(AK8CaloJetsIDPassed_eta[0])<5.0 && HighestAK8PFJetForBTag_csv>0.3 && AK8PFJetsMatchedToCaloJets280_pt[0]>340 && abs(AK8PFJetsMatchedToCaloJets280_eta[0])<5.0 && AK8TrimModJets_mass[0] > 30))"
#trigger         = "(l1Jets_pt[0]>170 && AK8CaloJetsIDPassed_pt[0]>300 && abs(AK8CaloJetsIDPassed_eta[0])<5.0 && HighestAK8PFJetForBTag_csv>0.24 && AK8PFJetsMatchedToCaloJets300_pt[0]>350 && abs(AK8PFJetsMatchedToCaloJets300_eta[0])<5.0)"
#trigger         = "(l1Jets_pt[0]>170 && AK8CaloJetsIDPassed_pt[0]>350 && abs(AK8CaloJetsIDPassed_eta[0])<5.0 && HighestAK8PFJetForBTag_csv>0.4 && AK8PFJetsMatchedToCaloJets350_pt[0]>400 && abs(AK8PFJetsMatchedToCaloJets350_eta[0])<5.0)"
#trigger         = "l1Jets_pt>170 && AK8CaloJetsIDPassed_pt>500 && abs(AK8CaloJetsIDPassed_eta)<5.0 && AK8PFJetsMatchedToCaloJets500_pt>550 && abs(AK8PFJetsMatchedToCaloJets500_eta)<5.0 && AK8PFJetForBTag_csv>0.63"
#trigger         = "HLT_AK8PFJet550_PFAK8BTagCSV_v2==1"
binning         = (NbinsY,yMin,yMax,Nbins,functionMin,functionMax)
ped             = "Signal_PFAK8BTagCSV_Trig"
ytitle          = "Leading AK8 double-b"
xtitle		= "Leading AK8 p_{T}"
doPlots()

#################### L1 high#########################
'''
parametersTurnOn_funct = (200,100,0.01,1,1E-3,1E-3)
Nbins       = 15#50
functionMin = 100#200
functionMax = 400#350
var             = "Jet_pt[0]+Jet_pt[1]+Jet_pt[2]+Jet_pt[3]"
#var             = "Jet_pt[1]+Jet_pt[2]+Jet_pt[3]"
trigger         = "ntrgObjects_hltQuadCentralJet30>=1"
binning         = (Nbins,functionMin,functionMax)
ped             = "QaD_DoubleJet_L1h"
title           = "p^{T}_{2}+p^{T}_{3}+p^{T}_{4}"
#doPlots()
'''
##################### CaloPt4 low ###et######################
parametersTurnOn_funct = (100,20,0.01,1,1E-3,1E-3)

Nbins       = 30#60
functionMin = 15#25
functionMax = 145#85
var             = "Jet_pt[3]"
#var             = "Sum$(Pt4(Jet_pt,Jet_eta,Jet_puId,3,Iteration$,Length$))"
preselection    = preselection + "&&"+ trigger
trigger         = "ntrgObjects_hltQuadCentralJet45>=4"
binning         = (Nbins,functionMin,functionMax)
ped             = "QaD_DoubleJet_CaloPt4"
title           = "p^{T}_{4}"
#doPlots()

##################### CaloPt4 high ###et######################
'''
parametersTurnOn_funct = (100,20,0.01,1,1E-3,1E-3)
Nbins       = 20
functionMin = 20
functionMax = 150
var             = "Jet_pt[3]"
#var             = "Sum$(Pt4(Jet_pt,Jet_eta,Jet_puId,3,Iteration$,Length$))"
#preselection    = preselection + "&&"+ trigger
trigger         = "ntrgObjects_hltQuadCentralJet30>=4"
binning         = (Nbins,functionMin,functionMax)
ped             = "QaD_DoubleJet_CaloPt4h"
title           = "p^{T}_{4}"
#doPlots()
'''


##################### CSV3 #########################
parametersTurnOn_funct = (100,20,0.01,1,1E-3,1E-3)
Nbins       = 24
functionMin =  0.4 #CSVL =  0.460 
functionMax = 1
#var             = "-log(1-Jet_btagCSV[aJCidx[0]])"
#var             = "Jet_btagCSV[aJCidx[0]]"
var             = "MaxIf$(Jet_btagCSV,Jet_btagCSV!=Max$(Jet_btagCSV)&&Jet_btagCSV!=MaxIf$(Jet_btagCSV,Jet_btagCSV!=Max$(Jet_btagCSV)))"
#var             = "Sum$(CSV(Jet_btagCSV,Jet_eta,Jet_puId,2,Iteration$,Length$))"
preselection    = preselection + "&&"+ trigger
#trigger         = "ntrgObjects_hltTripleCSV0p67>=3"
trigger = "ntrgObjects_hltBTagCaloCSVp087Triple>=3"
binning         = (Nbins,functionMin,functionMax)
ped             = "QaD_DoubleJet_CSV3"
title           = "CSV_{3}"
#doPlots()

###################### PFPt4 ########################
parametersTurnOn_funct = (0,20,0.01,1,1E-3,1E-3)

Nbins       = 50
functionMin = 20#25
functionMax = 140#125
var             = "Jet_pt[3]"
#var             = "Sum$(Pt4(Jet_pt,Jet_eta,Jet_puId,3,Iteration$,Length$))"
preselection    = preselection + "&&"+ trigger
trigger         = "ntrgObjects_hltQuadPFCentralJetLooseID45>=4"
binning         = (Nbins,functionMin,functionMax)
ped             = "QaD_DoubleJet_PFPt4"
title           = "p^{T}_{4}"
#doPlots()



################### L1 low ###################
parametersTurnOn_funct = (200,100,0.01,1)

Nbins           =  20#50
functionMin     =  250#100
functionMax     =  600#200
var             = "Jet_pt[0]+Jet_pt[1]+Jet_pt[2]+Jet_pt[3]"
#var             = "Jet_pt[1]+Jet_pt[2]+Jet_pt[3]"
trigger         = "ntrgObjects_hltQuadCentralJet30>=1"
binning         = (Nbins,functionMin,functionMax)
ped             = "QaD_QuaJet_L1"
title           = "p^{T}_{1}+p^{T}_{2}+p^{T}_{3}+p^{T}_{4}"
#doPlots()

#################### L1 high #########################
'''
parametersTurnOn_funct = (200,100,0.01,1)
Nbins           = 20#50
functionMin     =  300#120
functionMax     =  600#270
var             = "Jet_pt[0]+Jet_pt[1]+Jet_pt[2]+Jet_pt[3]"
#var             = "Jet_pt[1]+Jet_pt[2]+Jet_pt[3]"
trigger         = "ntrgObjects_hltQuadCentralJet45>=1"
binning         = (Nbins,functionMin,functionMax)
ped             = "QaD_QuaJet_L1h"
title           = "p^{T}_{1}+p^{T}_{2}+p^{T}_{3}+p^{T}_{4}"
doPlots()
'''
##################### CaloPt4 low #########################
parametersTurnOn_funct = (100,20,0.01,1)

Nbins       = 30#140#70
functionMin = 40#30#35
functionMax = 150#105
var             = "Jet_pt[3]"
#var             = "Sum$(Pt4(Jet_pt,Jet_eta,Jet_puId,3,Iteration$,Length$))"
preselection    = preselection + "&&"+ trigger
trigger         = "ntrgObjects_hltQuadCentralJet30>=4"
binning         = (Nbins,functionMin,functionMax)
ped             = "QaD_QuaJet_CaloPt4"
title           = "p^{T}_{4}"
#doPlots()


##################### CaloPt2 #########################
parametersTurnOn_funct = (100,20,0.01,1,1E-3,1E-3)

Nbins       = 50
functionMin = 60
functionMax = 250
var             = "Jet_pt[1]"
#var             = "Sum$(Pt4(Jet_pt,Jet_eta,Jet_puId,1,Iteration$,Length$))"
preselection    = preselection + "&&"+ trigger
trigger         = "ntrgObjects_hltDoubleCentralJet90>=2"
binning         = (Nbins,functionMin,functionMax)
ped             = "QaD_DoubleJet_CaloPt2"
title           = "p^{T}_{2}"
#doPlots()


##################### CaloPt4 high #########################
'''
parametersTurnOn_funct = (100,20,0.01,1)
Nbins       = 20#140#70
functionMin = 40#30#35
functionMax = 100#100#105
var             = "Jet_pt[3]"
#var             = "Sum$(Pt4(Jet_pt,Jet_eta,Jet_puId,3,Iteration$,Length$))"
#preselection    = preselection + "&&"+ trigger
trigger         = "ntrgObjects_hltQuadCentralJet45>=4"
binning         = (Nbins,functionMin,functionMax)
ped             = "QaD_QuaJet_CaloPt4h"
title           = "p^{T}_{4}"
doPlots()
'''

##################### CSV3 #########################
parametersTurnOn_funct = (100,20,0.01,1)
Nbins       = 10#80#24
functionMin = 0.6#0.2#0.4 #CSVL =  0.460 
functionMax = 1.0#1.0#1.0
#var             = "-log(1-Jet_btagCSV[aJCidx[0]])"
#var             = "Jet_btagCSV[aJCidx[0]]"
var             = "MaxIf$(Jet_btagCSV,Jet_btagCSV!=Max$(Jet_btagCSV)&&Jet_btagCSV!=MaxIf$(Jet_btagCSV,Jet_btagCSV!=Max$(Jet_btagCSV)))"
#var              = "CSVsorted[2]"
#var              = "Jet_btagCSV[2]"
#var             = "Sum$(CSV(Jet_btagCSV,Jet_eta,Jet_puId,2,Iteration$,Length$))"
preselection    = preselection + "&&"+ trigger
# old ! trigger         = "ntrgObjects_hltTripleCSV0p67>=3"
trigger = "ntrgObjects_hltBTagCaloCSVp087Triple>=3"
binning         = (Nbins,functionMin,functionMax)
ped             = "QaD_QuaJet_CSV3"
title           = "CSV_{3}"
#doPlots()

###################### PFPt4 ########################
parametersTurnOn_funct = (0,20,0.01,1)

Nbins       = 15#100#50
functionMin = 45#30#40
functionMax = 205#130#140
var             = "Jet_pt[3]"
#var             = "Sum$(Pt4(Jet_pt,Jet_eta,Jet_puId,3,Iteration$,Length$))"
preselection    = preselection + "&&"+ trigger
trigger         = "ntrgObjects_hltQuadPFCentralJetLooseID30>=4"
binning         = (Nbins,functionMin,functionMax)
ped             = "QaD_QuaJet_PFPt4"
title           = "p^{T}_{4}"
#doPlots()

##############################################

###################### PFPt2 ########################
parametersTurnOn_funct = (0,20,0.01,1,1E-3,1E-3)

Nbins       = 30
functionMin = 80
functionMax = 400
var             = "Jet_pt[1]"
#var             = "Sum$(Pt4(Jet_pt,Jet_eta,Jet_puId,1,Iteration$,Length$))"
preselection    = preselection + "&&"+ trigger
trigger         = "ntrgObjects_hltDoublePFCentralJetLooseID90>=2"
binning         = (Nbins,functionMin,functionMax)
ped             = "QaD_DoubleJet_PFPt2"
title           = "p^{T}_{2}"
#doPlots()


##f.Close()
