import ROOT 

def GetVariablesToFill(b1,b2,q1,q2):
    return ( abs(q1.Eta() - q2.Eta()), abs(b1.DeltaPhi(b2)), (q1+q2).M(), (b1+b2).M())

def SetVariable(tree,name,option='F',lenght=1):
    if option is 'F': arraytype='f'
    elif option is 'O': arraytype='i'
    elif option is 'I': arraytype='i'
    else:
        print 'option ',option,' not recognized.'
        return

    variable = array(arraytype,[0]*lenght)
    if lenght>1: name = name + '['+str(lenght)+']'
    tree.Branch(name,variable,name+'/'+option)
    return variable

def SortByEta(vect):
    vect.sort(key=lambda x: x.Eta(), reverse=True)

def SortByPt(vect):
    vect.sort(key=lambda x: x.Pt(), reverse=True)

def SortByCSV(vect):
    vect.sort(key=lambda x: x.csv, reverse=True)

def Sort(lorentzvectorswithcsv,method=''):
    (b1,b2,q1,q2) = [ROOT.TLorentzVector()]*4
    if len(lorentzvectorswithcsv)>=4:
        if method is '1BTagAndEta':

            SortByPt(lorentzvectorswithcsv)
            maxElementsPt = 4
            lorentzvectorswithcsv = lorentzvectorswithcsv[:maxElementsPt]
            SortByCSV(lorentzvectorswithcsv)
            b1 = lorentzvectorswithcsv[0]
            del lorentzvectorswithcsv[0]

            SortByEta(lorentzvectorswithcsv)
            q1 = lorentzvectorswithcsv[0]
            b2 = lorentzvectorswithcsv[1]
            q2 = lorentzvectorswithcsv[2]

        elif method is '2BTagAndPt':

            SortByPt(lorentzvectorswithcsv)
            maxElementsPt = 6
            lorentzvectorswithcsv = lorentzvectorswithcsv[:maxElementsPt]
            SortByCSV(lorentzvectorswithcsv)
            b1 = lorentzvectorswithcsv[0]
            b2 = lorentzvectorswithcsv[1]
            del lorentzvectorswithcsv[1]
            del lorentzvectorswithcsv[0]

            SortByPt(lorentzvectorswithcsv)
            q1 = lorentzvectorswithcsv[0]
            q2 = lorentzvectorswithcsv[1]

        elif method is 'Eta':

            SortByPt(lorentzvectorswithcsv)
            maxElementsPt = 4
            lorentzvectorswithcsv = lorentzvectorswithcsv[:maxElementsPt]
            SortByEta(lorentzvectorswithcsv)
            q1 = lorentzvectorswithcsv[0]
            b1 = lorentzvectorswithcsv[1]
            b2 = lorentzvectorswithcsv[2]
            q2 = lorentzvectorswithcsv[3]

        elif method is 'Gen':

            SortByPt(lorentzvectorswithcsv)
            for i,jet in enumerate(lorentzvectorswithcsv):
                if jet.mcMatchId==25 and jet.mcFlavour==5:
                    break

            b1_ = lorentzvectorswithcsv[i]
            del lorentzvectorswithcsv[i]

            for i,jet in enumerate(lorentzvectorswithcsv):
                if jet.mcMatchId==25 and jet.mcFlavour==-5:
                    break

            b2_ = lorentzvectorswithcsv[i]
            del lorentzvectorswithcsv[i]

            mjj = 0
            for i,jet in enumerate(lorentzvectorswithcsv):
                if jet.mcFlavour!=0 and jet.mcFlavour!=21 and (abs(jet.mcFlavour)<=2):
                    for j,jet2 in enumerate(lorentzvectorswithcsv):
                        if j>i and jet2.mcFlavour!=0 and jet2.mcFlavour!=21 and (abs(jet2.mcFlavour)<=2):
                            if (jet+jet2).M()>mjj:
                                mjj = (jet+jet2).M()
                                q1_ = jet
                                q2_ = jet2

            try:
                (b1,b2,q1,q2) = (b1_,b2_,q1_,q2_)
            except:
#                print "excption"
                pass

        else:
            print 'method:',method,' not found'
    else:
#        print 'less than 4 jets!'
        pass

    if b2.Pt()>b1.Pt(): (b1,b2) = (b2,b1)
    if q2.Pt()>q1.Pt(): (q1,q2) = (q2,q1)
    return (b1,b2,q1,q2)
