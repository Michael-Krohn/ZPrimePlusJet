import ROOT
from VBFutils import Sort,GetVariablesToFill

PU_max = 60.
PU_min = 50.
lumi = 2.E34

class Jet:
    def __init__(self, pt, eta, phi, mass, csv):
        self.pt = pt
        self.eta = eta
        self.phi = phi
        self.mass = mass
        self.csv = csv

def getNum(tree):
    num = 0.
    den = 0.
    num2 = 0.
    num3 = 0.
    tree.SetBranchStatus("*",0)
    tree.SetBranchStatus("maxPUptHat",1)
    tree.SetBranchStatus("l1Jets_num",1)
    tree.SetBranchStatus("l1Jets_pt[l1Jets_num]",1)
    tree.SetBranchStatus("l1HT",1)
    tree.SetBranchStatus("caloJets_num",1)
    tree.SetBranchStatus("caloJets_pt[caloJets_num]",1)
    tree.SetBranchStatus("caloJets_eta[caloJets_num]",1)
    tree.SetBranchStatus("caloJets_phi[caloJets_num]",1)
    tree.SetBranchStatus("caloJets_mass[caloJets_num]",1)
    tree.SetBranchStatus("caloJets_csv[caloJets_num]",1)
    tree.SetBranchStatus("pfJets_pt[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_eta[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_phi[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_mass[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_csv[pfJets_num]",1)
    tree.SetBranchStatus("pfJets_num",1)
    tree.SetBranchStatus("run",1)
    tree.SetBranchStatus("HLT_HT300_Beamspot_v2",1)
    tree.SetBranchStatus("HLT_PFJet500_v12",1)
    for ev in tree:
#        if not ev.run==297292: continue
        den+=1

#	if not (ev.HLT_HT300_Beamspot_v2==1): continue
	if not (ev.HLT_PFJet500_v12==1): continue

        num += 1

    print "Numerator: ", num
    print "Numerator2: ", num2
    print "Numerator3: ", num3
    print "Denominator: ", den    
    return num,den

def getFraction(fileName):
    file_ = ROOT.TFile(fileName)
    tree = file_.Get("tree")
#    den = file_.Get("Count").GetBinContent(1)
#    den *= (PU_max-PU_min)/(63-28)
    num,den = getNum(tree)
    if num<=0: return (0.,0.)
    file_.Close()
    fraction = num/den
    return (fraction , (fraction*(1-fraction)/den)**0.5)

def getRate(fileName, xsection, lumi):
    sampleRate = 1.E-36*xsection*lumi
    (fraction, errFraction) = getFraction(fileName)
    return (sampleRate * fraction , sampleRate * errFraction)

## https://github.com/cms-steam/RateEstimate/blob/master/datasetCrossSections/datasetCrossSectionsSummer16.py#L66
bkgCrossSection = [
#    ("QCD15to30.root"   , 1837410000.),
#    ("QCD30to50.root"   , 140932000.),
#    ("QCD50to80.root"   , 19204300.),
    ("../QCD80to120.root"  , 2762530.),
    ("../QCD120to170.root" , 471100.),
    ("../QCD170to300.root" , 117276.),
    ("../QCD300to470.root" , 7823.),
    ("../QCD470to600.root" , 648.2),
]

totalRate = 0
totalRateError = 0

#for (fileName,xsection) in bkgCrossSection:
#    rate , rateError = getRate(fileName, xsection, lumi)
#    print fileName,"\t",rate," +/- ",rateError
#    totalRate += rate
#    totalRateError += rateError

#print "totalRate\t",totalRate," +/- ",totalRateError

for signal in ["ParkedData_Sklims/HLTPhysics_2017C.root"]:
#for signal in ["ParkedData_Sklims/HLTPhysics_2017C.root","ParkedData_Sklims/HLTPhysics_2017B_Run297292to297296.root","ParkedData_Sklims/HLTPhysics_2017B_Run297424to297435.root"]:
#for signal in ["ParkedData_Sklims/HLTPhysics_2017B_Run297292to297296.root"]:
    (fraction, errFraction) = getFraction(signal)
    print signal,"\t",fraction," +/- ",errFraction
    print signal, "Rate(54kHz L1 Rate)\t",fraction*54000.,"+/-",errFraction*54000.
