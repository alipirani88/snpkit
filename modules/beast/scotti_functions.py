'''Functions to generate a SCOTTI xml file. Modified from:
https://github.com/Taming-the-BEAST/SCOTTI-Tutorial/blob/master/scripts/SCOTTI_generate_xml.py
'''

# Import modules
import sys
import os
import math
import csv
from collections import defaultdict
from datetime import datetime
import time
import re

#does it represent an integer?
def IsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

#does it represent a float?
def IsFloat(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False

#make tree string into format for mugration
def removeStates(treeString):
        l=len(treeString)
        newString=""
        i=0
        while i < l :
                if treeString[i]=="[":
                        while treeString[i]!="]":
                                i+=1
                        i+=1
                newString+=treeString[i]
                i+=1
        return newString                

#scale branch lengths 
def scaleTreeString(treeString,scale):
        l=len(treeString)
        newString=""
        i=0
        while i < l :
                if treeString[i]==":":
                        newString+=":"
                        numb=""
                        i+=1
                        while treeString[i]!="," and treeString[i]!=")" and treeString[i]!=";":
                                numb+=treeString[i]
                                i+=1
                        newString+=str(float(numb)*scale)
                newString+=treeString[i]
                i+=1
        return newString

def generate_scotti_xml(fasta, dates, hosts, hostTimes, output=None, overwrite=False, maxHosts=None,  unlimLife=False, penalizeMigration=False, numIter=100000000, mutationModel='HKY', fixedAs=0, fixedCs=0, fixedGs=0, fixedTs=0):
    """Generates an xml file that can be run in BEAST2.

    Args:
        fasta:          Input fasta alignment. If not specified, looks in the
                        working directory for SCOTTI_aligment.fas
        output:		    output file that will contain the newly created SCOTTI
                        xml to be run with BEAST2. (default: fasta file base name)
        overwrite:      Overwrite output files.
        dates:          Input csv file with sampling dates associated to each
                        sample name (same name as fasta file). If not
                        specified, looks for SCOTTI_dates.csv
        hosts:		Input csv file with sampling hosts associated to each
                        sample name (same name as fasta file). If not
                        specified, looks for SCOTTI_hosts.csv
        hostTime:       Input csv file with the earliest times in which hosts
                        are infectable, and latest time when hosts are
                        infective. Use same host names as those used for the
                        sampling hosts file. If not specified, looks for
                        SCOTTI_hostTimes.csv
        maxHosts:       Maximum number of hosts (between sampled, non-sampled,
                        and non-observed) allowed.
        unlimLife:      span (unlimited contribution in time to the outbreak).
                        By default non-sampled non-observed hosts have limited
                        duration. Usually asymptomaic patients have limited
                        contributions, while to model environmental
                        contamination it's better to use unlimited hosts.
                        Unlimited life span can also decrease computational
                        demand.
        penalizeMigration: penalize lineages that are still in a deme after deme
                        closure. By default, don't (original version).
        numIter:        Number of iterations of the MCMC. Default=100,000,000.
        mutationModel:  String represnting the mutation model to be used.
                        Possibilities are "HKY" or "JC". Further alternaties
                        can be specified by modifying the output xml manually.
        fixedAs:        Number of sites fixed for A in the genome and not
                        included in the alignment.
        fixedCs:        Number of sites fixed for C in the genome and not
                        included in the alignment.
        fixedGs:        Number of sites fixed for G in the genome and not
                        included in the alignment.
        fixedTs:        Number of sites fixed for T in the genome and not
                        included in the alignment.

    Output: 
        XML file that can be used to run SCOTTI in BEAST2.
    """
    # Set output file name if none specified
    if output is None:
        output = fasta.split('.')[0]
    
    # Check if input arguments make sense
    if (not (os.path.isfile(fasta))):
        print("Error: input fasta file not found, please specify it correctly in input using the option -f .\n")
        return()
    if (os.path.isfile(output + ".xml")) and (not overwrite):
        print("Output xml file already exists, to allow overwriting of output file use option -ov, otherwise specify a different output file.\n")
        return()
    if not (os.path.isfile(dates)):
        print("Error: input dates csv file not found, please specify it correctly in input using the option -d .\n")
        return()
    if not (os.path.isfile(hosts)):
        print("Error: input hosts csv file not found, please specify it correctly in input using the option -s .\n")
        return()
    if not (os.path.isfile(hostTimes)):
        print("Error: input hosts csv file not found, please specify it correctly in input using the option -i .\n")
        return()
    if (numIter<=0):
        print("Error: the number of MCMC iterations specified is negative or null.\n")
        return()
    elif (numIter<=1000):
        print("Very small number of MCMC iterations, might want to specifiy something bigger.\n")
    if (mutationModel!="HKY" and mutationModel!="JC"):
        print("Mutation model not supported. Please specify something among \"HKY\" or \"JC\".\n")
        return()
    if (fixedAs<0 or fixedCs<0 or fixedGs<0 or fixedTs<0):
        print("Error: the number of fixed sites for each base type must be positive.\n")
        return()

    #read fasta file
    sqFile=open(fasta)
    seqs={}
    line=sqFile.readline()
    linesplit=line.split()
    while line!="" and (len(linesplit)<1 or linesplit[0][0]!=">"):
        line=sqFile.readline()
        linesplit=line.split()
    if line=="":
        print("Wrong format in fasta file\n")
        return()
    while line!="":
        if len(linesplit[0])>1:
            name=linesplit[0].replace(">","")
        else:
            name=linesplit[1]
        seq=""
        line=sqFile.readline()
        linesplit=line.split()
        while line!="" and (len(linesplit)>0  and linesplit[0][0]!=">") and line!="\n":
            seq+=linesplit[0]
            line=sqFile.readline()
            linesplit=line.split()
        if name in list(seqs.keys()):
            print("Error, multiple entries with same name "+name+" in fasta file. Only first part of the name line is used for the name.")
            return()
        seqs[name]=seq
        while line!="" and  (len(linesplit)<1 or linesplit[0][0]!=">"):
            line=sqFile.readline()
            linesplit=line.split()
    sqFile.close()

    #read in file with host-sample information
    hFile=open(hosts)
    hosts={}
    line=hFile.readline()
    linesplit=("".join(line.split())).split(",")  
    while line!="" and len(linesplit)>1: 
        sam=linesplit[0]
        hos=linesplit[1]
        if not (sam in list(seqs.keys())):
            print("Error: sample "+sam+" found in sample-host file but not found in the fasta file.")
            return()
        if sam in list(hosts.keys()):
            print("Error, multiple entries with same name "+sam+" in hosts file. ")
            return()
        hosts[sam]=hos
        line=hFile.readline()
        linesplit=("".join(line.split())).split(",")
        while line=="\n":
            line=hFile.readline()
            linesplit=("".join(line.split())).split(",")
    hFile.close()

    #read in sampling dates
    dFile=open(dates)
    dates={}
    line=dFile.readline()
    linesplit=("".join(line.split())).split(",")  
    while line!="" and len(linesplit)>1:  
        sam=linesplit[0]
        dat=linesplit[1]
        if not (IsFloat(dat)):
            print("Error: date "+str(dat)+" for sample "+sam+" cannot be converted to a number")
            return()
        if not (sam in list(seqs.keys())):
            print("Error: sample "+sam+" found in sampling dates file but not found in the fasta file.")
            return()
        if sam in list(dates.keys()):
            print("Error, multiple entries with same name "+sam+" in dates file. ")
            return()
        dates[sam]=float(dat)
        line=dFile.readline()
        linesplit=("".join(line.split())).split(",")
        while line=="\n":
            line=dFile.readline()
            linesplit=("".join(line.split())).split(",")
    dFile.close()
        
    #read in file with host times information
    hFile=open(hostTimes)
    hostT={}
    line=hFile.readline()
    linesplit=("".join(line.split())).split(",")  
    while line!="" and len(linesplit)>2: 
        hos=linesplit[0]
        dat=float(linesplit[1])
        dat2=float(linesplit[2])
        if (not (IsFloat(dat))) or (not (IsFloat(dat2))):
            print("Error: date "+str(dat)+" or "+str(dat2)+" for sample "+sam+" cannot be converted to a number")
            return()
        #if not (hos in seqs.keys()):
        #    print "Error: sample "+hos+" found in host times file but not found in the fasta file."
        #    return()
        if hos in list(hostT.keys()):
            print("Error, multiple entries with same name "+hos+" in hosts-times file. ")
            return()
        hostT[hos]=[dat, dat2]
        line=hFile.readline()
        linesplit=("".join(line.split())).split(",")
        while line=="\n":
            line=hFile.readline()
            linesplit=("".join(line.split())).split(",")
    hFile.close()

    #check that the details are correct
    for s in list(seqs.keys()):
        if not (s in list(dates.keys())):
            print("Error, no date specified for sample "+s)
            return()
        if not (s in list(hosts.keys())):
            print("Error, no host specified for sample "+s)
            return()
    for s in list(hosts.keys()):
        if not (s in list(seqs.keys())):
            print("Error, no sequence specified for sample "+s)
            return()
        if not (s in list(dates.keys())):
            print("Error, no date specified for sample "+s)
            return()
    for s in list(dates.keys()):
        if not (s in list(seqs.keys())):
            print("Error, no sequence specified for sample "+s)
            return()
        if not (s in list(hosts.keys())):
            print("Error, no host specified for sample "+s)
            return()
    for s in list(hosts.keys()):
        if not (hosts[s] in list(hostT.keys())):
            print("Warning: no dates specified for host "+hosts[s]+". I will assume this host is available throghtout")
            hostT[hosts[s]]=[-1000000000000,1000000000000]


    for h in list(hostT.keys()):
        found=0
        for s in list(hosts.keys()):
            if hosts[s]==h:
                found=1
        if found==0:
            print("Warning: no sample specified for host "+h+". I will assume this host could have been infected or not. If you are sure this host was infected, then please include a non-informative sample (all \"N\"s) from it.")

    xml=open(output+'.xml',"w")
    xml.write("<beast version=\'2.0\' namespace=\'beast.evolution.alignment:beast.core:beast.core.parameter:beast.evolution.tree:beast.evolution.tree.coalescent:beast.core.util:beast.evolution.operators:beast.evolution.sitemodel:beast.evolution.substitutionmodel:beast.evolution.likelihood:beast.evolution.tree:beast.math.distributions:multitypetreeVolz.distributions:multitypetreeVolz.operators:multitypetreeVolz.util\'>\n\n"+"  <data id=\"alignmentVar\" dataType=\"nucleotide\">\n")
    for s in list(seqs.keys()):
        xml.write("    <sequence taxon=\'"+s+"\' value=\'"+seqs[s]+"\'/>\n")
    xml.write("  </data>\n")
    st=""
    if fixedAs>0 and fixedCs>0 and fixedGs>0 and fixedTs>0:
        st="constantSiteWeights=\'"+str(fixedAs)+" "+str(fixedCs)+" "+str(fixedGs)+" "+str(fixedTs)+"\'"
    xml.write("  <data id=\'alignment\' spec=\'FilteredAlignment\' filter=\'-\' data=\'@alignmentVar\' "+st+"/>\n\n")

    #hosts
    xml.write("   <typeTraitSet spec=\'TraitSetWithLimits\' id=\'typeTraitSet\' traitname=\"type\" value=\"")
    hL=list(hosts.keys())
    for s in hL:
        xml.write(s+"="+hosts[s])
        if s!=hL[len(hL)-1]:
            xml.write(",")
    xml.write("\">\n     <taxa spec=\'TaxonSet\' alignment=\'@alignment\'/>\n   </typeTraitSet>\n\n")

    #dates
    xml.write("   <timeTraitSet spec=\'TraitSetWithLimits\' id=\'timeTraitSet\' traitname=\"date-forward\" value=\"")
    dL=list(dates.keys())
    for s in range(len(dL)):
        xml.write(dL[s]+"="+str(dates[dL[s]]))
        if (s!=(len(dL)-1)):
            xml.write(",")
    xml.write("\">\n     <taxa spec=\'TaxonSet\' alignment=\'@alignment\'/>\n   </timeTraitSet>\n\n")

    #end infectability
    xml.write("   <endInfectionsTraitSet spec=\'InfectionTraitSet\' id=\'endInfectionsTraitSet\' type=\'end\' traitname=\"date-forward\" value=\"")
    dL=list(hostT.keys())
    for s in dL:
        xml.write(s+"="+str(hostT[s][1]))
        if s!=dL[len(dL)-1]:
            xml.write(",")
    xml.write("\">\n          <taxaDates idref=\'timeTraitSet\'/>\n           <taxaTypes idref=\'typeTraitSet\'/>\n   </endInfectionsTraitSet>\n\n")

    #start infectability
    xml.write("   <startInfectionsTraitSet spec=\'InfectionTraitSet\' id=\'startInfectionsTraitSet\' type=\'start\' traitname=\"date-forward\" value=\"")
    dL=list(hostT.keys())
    for s in dL:
        xml.write(s+"="+str(hostT[s][0]))
        if s!=dL[len(dL)-1]:
            xml.write(",")
    xml.write("\">\n          <taxaDates idref=\'timeTraitSet\'/>\n           <taxaTypes idref=\'typeTraitSet\'/>\n   </startInfectionsTraitSet>\n\n")

    #check that times fit
    for s in list(seqs.keys()):
        h=hosts[s]
        if dates[s]<=hostT[hosts[s]][0]:
            print("Error, sampling time "+str(dates[s])+" of sample "+s+" must be after start of infectability time "+str(hostT[hosts[s]][0])+" of host "+h)
            return()
        if dates[s]>=hostT[hosts[s]][1]:
            print("Error, sampling time "+str(dates[s])+" of sample "+s+" must be before end of infectiousness time "+str(hostT[hosts[s]][0])+" of host "+h)
            return()

    #include mutation model
    xml.write("   <siteModel spec=\"SiteModel\" id=\"siteModel\">\n     <mutationRate spec=\'RealParameter\' id=\"mutationRate\" value=\"1.0\"/>\n     <substModel spec=\"")
    if mutationModel=="HKY":
       xml.write("HKY\">\n       <kappa spec=\'RealParameter\' id=\"hky.kappa\" value=\"1.0\"/>\n")
    elif mutationModel=="JC":
       xml.write("jc69\">\n")
    elif mutationModel=="GTR":
       xml.write("GTR\" rateAC=\"@rateAC\" rateAG=\"@rateAG\" rateAT=\"@rateAT\" rateCG=\"@rateCG\" rateGT=\"@rateGT\" >\n       <parameter estimate=\'false\' name=\"rateCT\" lower=\"0.0\" value=\"1.0\" id=\"rateCTfixed\"/>\n")
    else:
        print("Error, mutation model "+mutationModel+" not recognized.")
        return()
    xml.write("       <frequencies estimate=\"false\" spec=\'Frequencies\'>\n        <frequencies spec=\'RealParameter\' id=\""+mutationModel+".freq\" value=\"0.25 0.25 0.25 0.25\"/>\n       </frequencies>\n     </substModel>\n   </siteModel>\n\n")
                  
    #migration model
    nHos=len(list(hostT.keys()))
    meanL=0.0
    for h in list(hostT.keys()):
        meanL+=hostT[h][1]-hostT[h][0]
    meanL=meanL/len(list(hostT.keys()))
    ma=maxHosts
    if ma<nHos:
        ma=nHos+5
        print("Maximum number of hosts specified is less than the minimum number of hosts. Setting it to "+str(int(ma)))
    if ma<(nHos+2):
        print("Warning: the maximum number of hosts is less than "+str(int(nHos+2))+". If there are problems in the execution of the BEAST2 analysis, try a higher values of maxHosts.")
    start=(nHos+ma)/2
    if start<(nHos+2) and ma>=(nHos+2):
        start=nHos+2
    xml.write("   <migrationModelUniform spec=\'MigrationModelUniform\' id=\'migModel\' minDemes=\'"+str(int(nHos))+"\'>\n     <rate spec=\'RealParameter\' value=\""+str(1.0/(meanL*nHos))+"\" dimension=\"1\" id=\"rate\"/>\n     <popSize spec=\'RealParameter\' value=\"1.0\" dimension=\"1\" id=\"popSize\"/>\n     <numDemes spec=\'IntegerParameter\' value=\""+str(int(start))+"\" dimension=\"1\" id=\"numDemes\" lower=\'"+str(int(nHos))+"\' upper=\'"+str(int(ma))+"\'/>\n     <trait idref=\'startInfectionsTraitSet\' type=\'start\'/>\n     <trait idref=\'endInfectionsTraitSet\' type=\'end\'/>\n   </migrationModelUniform>\n\n")            

    #priors
    xml.write("   <input spec=\'CompoundDistribution\' id=\'parameterPriors\'>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@mutationRate\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"-3.0\" S=\"6.0\"/>\n     </distribution>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@popSize\">\n       <distr spec=\"LogNormalDistributionModel\"  M=\"0.0\" S=\"6.0\"/>\n     </distribution>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@rate\">\n       <distr spec=\'LogNormalDistributionModel\' M=\""+str(math.log(1.0/(meanL*nHos)))+"\" S=\"0.2\"/>\n     </distribution>\n")
    if mutationModel=="HKY":
        xml.write("     <distribution spec=\'beast.math.distributions.Prior\' x=\"@hky.kappa\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"4.0\"/>\n     </distribution>\n   </input>\n\n")
    elif mutationModel=="JC":
       xml.write("   </input>\n\n")
    elif mutationModel=="GTR":
        xml.write("     <distribution spec=\'beast.math.distributions.Prior\' x=\"@rateAC\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"4.0\"/>\n     </distribution>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@rateAG\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"4.0\"/>\n     </distribution>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@rateAT\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"4.0\"/>\n     </distribution>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@rateCG\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"4.0\"/>\n     </distribution>\n     <distribution spec=\'beast.math.distributions.Prior\' x=\"@rateGT\">\n       <distr spec=\'LogNormalDistributionModel\' M=\"0.0\" S=\"4.0\"/>\n     </distribution>\n   </input>\n\n")      
                  
    #Likelihoods
    xml.write("<input spec=\'TreeLikelihood\' id=\"treeLikelihood\">\n     <data idref=\"alignment\"/>\n     <tree idref=\"tree\"/>\n     <siteModel idref=\'siteModel\'/>\n   </input>\n\n   <input spec=\'StructuredCoalescentTreeDensityNew\' id=\'treePrior\' limitedLifespan=\'"+str(unlimLife)+"\' penalizeMigration=\'"+str(penalizeMigration)+"\'>\n     <multiTypeTreeConcise idref=\"tree\"/>\n     <migrationModelUniform idref=\"migModel\"/>\n   </input>\n\n   <run spec=\"MCMC\" id=\"mcmc\" chainLength=\""+str(numIter)+"\" storeEvery=\"10000\">\n\n")

    #migration model
    xml.write("     <init spec=\'StructuredCoalescentMultiTypeTreeConcise\' id=\'tree\' nTypes=\""+str(int(start))+"\">\n         <migrationModelUniform spec=\'MigrationModelUniform\' minDemes=\'"+str(int(nHos))+"\'>\n             <rate spec=\'RealParameter\' value=\""+str(1.0/(meanL*nHos))+"\" dimension=\"1\"/>\n             <popSize spec=\'RealParameter\' value=\"1.0\" dimension=\"1\"/>\n             <numDemes spec=\'IntegerParameter\' value=\""+str(int(start))+"\" dimension=\"1\"/>\n         </migrationModelUniform>\n         <trait idref=\'typeTraitSet\'/>\n         <trait idref=\'timeTraitSet\'/>\n     </init>\n\n")    

    #state nodes              
    xml.write("     <state>\n       <stateNode idref=\"tree\"/>\n       <stateNode idref=\"rate\"/>\n       <stateNode idref=\"popSize\"/>\n       <stateNode idref=\"numDemes\"/>\n       <stateNode idref=\"mutationRate\"/>\n")
    if mutationModel=="HKY":
        xml.write("       <stateNode idref=\"hky.kappa\"/>\n")
    elif mutationModel=="GTR":              
        xml.write("       <stateNode idref=\"rateAC\"/>\n       <stateNode idref=\"rateAG\"/>\n       <stateNode idref=\"rateAT\"/>\n       <stateNode idref=\"rateCG\"/>\n       <stateNode idref=\"rateGT\"/>\n")
    xml.write("       <stateNode idref=\""+mutationModel+".freq\"/>\n     </state>\n\n")     
                  
    #compound distribution
    xml.write("     <distribution spec=\'CompoundDistribution\' id=\'posterior\'>\n       <distribution idref=\"treeLikelihood\"/>\n       <distribution idref=\'treePrior\'/>\n       <distribution idref=\"parameterPriors\"/>\n     </distribution>\n\n")
                  
    #operators
    xml.write("     <operator spec=\'ScaleOperator\' id=\'RateScaler\' parameter=\"@rate\" scaleFactor=\"0.8\" weight=\"1\"/>\n     <operator spec=\"ScaleOperator\" id=\"PopSizeScaler\" parameter=\"@popSize\" scaleFactor=\"0.8\" weight=\"1\"/>\n     <operator spec=\"ShortRangeUniformOperator\" id=\"NumDemesScaler\" parameter=\"@numDemes\" weight=\"1\"/>\n     <operator spec=\"ScaleOperator\" id=\"muRateScaler\" parameter=\"@mutationRate\" scaleFactor=\"0.8\" weight=\"1\"/>\n     <operator spec=\"DeltaExchangeOperator\" id=\"freqExchanger\" parameter=\"@"+mutationModel+".freq\" delta=\"0.01\" weight=\"0.1\"/>\n     <operator id=\'treeScaler.t\' spec=\'ScaleOperator\' scaleFactor=\"0.5\" weight=\"3\" tree=\"@tree\"/>\n     <operator id=\'treeRootScaler.t\' spec=\'ScaleOperator\' scaleFactor=\"0.5\" weight=\"3\" tree=\"@tree\" rootOnly=\'true\'/>\n     <operator id=\'UniformOperator.t\' spec=\'Uniform\' weight=\"30\" tree=\"@tree\"/>\n     <operator id=\'SubtreeSlide.t\' spec=\'SubtreeSlide\' weight=\"15\" gaussian=\"true\" size=\"1.0\" tree=\"@tree\"/>\n     <operator id=\'narrow.t\' spec=\'Exchange\' isNarrow=\'true\' weight=\"15\" tree=\"@tree\"/>\n     <operator id=\'wide.t\' spec=\'Exchange\' isNarrow=\'false\' weight=\"3\" tree=\"@tree\"/>\n     <operator id=\'WilsonBalding.t\' spec=\'WilsonBalding\' weight=\"3\" tree=\"@tree\"/>\n")
    if mutationModel=="HKY":
        xml.write("     <operator spec=\'ScaleOperator\' id=\'kappaScaler\' parameter=\"@hky.kappa\" scaleFactor=\"0.8\" weight=\"0.1\"/>\n")
    elif mutationModel=="GTR": 
        xml.write("     <operator spec=\'ScaleOperator\' id=\'ACScaler\' parameter=\"@rateAC\" scaleFactor=\"0.8\" weight=\"0.1\"/>\n     <operator spec=\'ScaleOperator\' id=\'AGScaler\' parameter=\"@rateAG\" scaleFactor=\"0.8\" weight=\"0.1\"/>\n     <operator spec=\'ScaleOperator\' id=\'ATScaler\' parameter=\"@rateAT\" scaleFactor=\"0.8\" weight=\"0.1\"/>\n     <operator spec=\'ScaleOperator\' id=\'CGScaler\' parameter=\"@rateCG\" scaleFactor=\"0.8\" weight=\"0.1\"/>\n     <operator spec=\'ScaleOperator\' id=\'GTScaler\' parameter=\"@rateGT\" scaleFactor=\"0.8\" weight=\"0.1\"/>\n")
    xml.write("\n")

    #loggers
    xml.write("     <logger logEvery=\""+str(int(numIter/10000))+"\" fileName=\""+output+".log\">\n       <model idref=\'posterior\'/>\n       <log idref=\"posterior\"/>       <log idref=\"treeLikelihood\"/>\n       <log idref=\"migModel\"/>\n       <log idref=\"mutationRate\"/>\n       <log spec=\'TreeHeightLogger\' tree=\'@tree\'/>\n       <log idref=\""+mutationModel+".freq\"/>\n")
    if mutationModel=="HKY":
        xml.write("       <log idref=\"hky.kappa\"/>\n")
    elif mutationModel=="GTR":               
        xml.write("       <log idref=\"rateAC\"/>\n       <log idref=\"rateAG\"/>\n       <log idref=\"rateAT\"/>\n       <log idref=\"rateCG\"/>\n       <log idref=\"rateGT\"/>\n")
    xml.write("     </logger>\n")
         
    xml.write("     <logger logEvery=\""+str(int(numIter/5000))+"\" fileName=\""+output+".trees\" mode=\"tree\">\n       <log idref=\'treePrior\'/>\n     </logger>\n")

    xml.write("     <logger logEvery=\""+str(int(numIter/1000))+"\">\n       <model idref=\'posterior\'/>\n       <log idref=\"posterior\"/>\n       <log idref=\"treeLikelihood\"/>\n       <log spec=\'TreeHeightLogger\' tree=\'@tree\'/>\n       <log idref=\"migModel\"/>\n       <log idref=\"mutationRate\"/>\n       <ESS spec=\'ESS\' name=\'log\' arg=\"@treePrior\"/>\n       <ESS spec=\'ESS\' name=\'log\' arg=\"@posterior\"/>\n     </logger>\n\n   </run>\n </beast>\n\n")

    return(output+'.xml')


def check_list_int(li):
    """Checks a list to see if any of the values can be converted to int.

    If any value in the input list can be converted to an int, True is 
    returned. Otherwise, False  is returned.

    Args:
        li: list to check

    Returns:
        Boolean value. True if at least one element in the list can be 
        converted to int, false otherwise.    
    """
    no_numbers = True
    for i in li:
        try:
            int(i)
            no_numbers = False
            break
        except:
            continue
    return not no_numbers

def add_prefixes(hosts, dates, hostTimes, fasta, id_pref='S', host_pref='H'):
    """Changes sequence id and host name prefixes in scotti input.

    Checks to see if the sequence ids and host names in the scotti input
    files are numeric. If they are, adds a character to the beginning of each
    sample id ('S') and host ('H') and writes output files with renamed.
    This fixes an error that occurs when numeric sequences are used to 
    generate the MCC tree.

    Args:
        hosts: csv where the first column is the sample id and the 
            second column is the host
        dates: csv where the first column is the sample id and the
            second column is the sampling date 
        hostTimes: csv where the first column is the host, the second
            column is the earliest time the host is infectable, and
            the third column is the latest time the host is infective
        fasta: list of fasta file(s) to make scotti xml(s)
        id_pref: the prefix to prepend to the sample names, if needed 
            (default: 'S')
        host_pref: the prefix to prepend to the host names, if needed
            (default: 'H') 

    Output:
        Nothing if all sample ids and host names already have a letter in them.
        Up to 4 renamed files if some sample ids and/or host names do not have 
        a letter in them. These are the original files renamed as .renamed.csv
        and .renamed.fa (for the fasta file)

    Returns:
        Nothing if all sample ids and host names already have a letter in them.
        A list with the names of the output files if they are generated.
        The order of the files is: hosts, dates, hostTimes, fasta
    """

    hosts_dict = defaultdict(list)

    with open(hosts) as f:
        reader = csv.reader(f)
        for row in reader:
            for (i,v) in enumerate(row):
                hosts_dict[i].append(v)    

    i_hosts = hosts_dict[0] 
    h_hosts = hosts_dict[1]

    num_in_i = check_list_int(i_hosts)
    num_in_h = check_list_int(h_hosts)

    if all([not num_in_i, not num_in_h]):
        return

    dates_dict = defaultdict(list)

    with open(dates) as f:
        reader = csv.reader(f)
        for row in reader:
            for (i,v) in enumerate(row):
                dates_dict[i].append(v)

    i_dates = dates_dict[0] 
    d_dates = dates_dict[1]

    times_dict = defaultdict(list)

    with open(hostTimes) as f:
        reader = csv.reader(f)
        for row in reader:
            for (i,v) in enumerate(row):
                times_dict[i].append(v)        

    h_times = times_dict[0] 
    t1_times = times_dict[1]
    t2_times = times_dict[2]

    hosts_renamed = hosts.split('.')[0] + '.renamed.csv'
    dates_renamed = dates.split('.')[0] + '.renamed.csv'
    times_renamed = hostTimes.split('.')[0] + '.renamed.csv'
    fasta_renamed = [fa.split('.')[0] + '.renamed.fa' for fa in fasta]
  
    if num_in_i:
        i_hosts_renamed = [id_pref + x for x in i_hosts]
        i_dates_renamed = [id_pref + x for x in i_dates]
        for i,fa in enumerate(fasta):
            os.system('sed \'s/>/>S/g\' ' + fa + ' > ' + fasta_renamed[i])
    else:
        i_hosts_renamed = i_hosts
        i_dates_renamed = i_dates
    
    if num_in_h:
        h_hosts_renamed = [host_pref + x for x in h_hosts]
        h_times_renamed = [host_pref + x for x in h_times]
    else:
        h_hosts_renamed = h_hosts
        h_times_renamed = h_times

    with open(hosts_renamed, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(i_hosts_renamed, h_hosts_renamed))

    with open(dates_renamed, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(i_dates_renamed, d_dates))

    with open(times_renamed, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(h_times_renamed, t1_times, t2_times))
    
    print('Changed prefixes for ids and/or hosts. New files:\n' + 
        hosts_renamed + '\n' + dates_renamed + '\n' + times_renamed + '\n' + ' '.join(fasta_renamed))

    return [hosts_renamed, dates_renamed, times_renamed, fasta_renamed]

def get_seconds_since_epoch(date):
    """Gets econds since epoch of a datetime date.

    Args: 
        date: datetime date

    Returns:
        seconds since eopch

    """
    return time.mktime(date.timetuple())


def convert_ymd_to_decimal(dates):
    """Converts a list of dates in yyyy-mm-dd or yyyy/mm/dd format to decimal format.

    Modified from: https://stackoverflow.com/questions/6451655/python-how-to-convert-datetime-dates-to-decimal-years

    Args:
        dates: list of dates to convert
    
    Returns:
        list of dates in decimal format

    """

    decimal_dates = []

    for d in dates:
        if re.search('-',d) is not None:
            date = datetime.strptime(d, '%Y-%m-%d')
        elif re.search('/',d) is not None:
            date = datetime.strptime(d, '%Y/%m/%d')
        else:
            print('Dates not in correct format. Format should be yyyy-mm-dd or yyyy/mm/dd')
            return

        start_this_year = datetime(year=date.year, month=1, day=1)
        start_next_year = datetime(year=date.year+1, month=1, day=1)

        year_elapsed = get_seconds_since_epoch(date) - get_seconds_since_epoch(start_this_year)
        year_duration = get_seconds_since_epoch(start_next_year) - get_seconds_since_epoch(start_this_year)
        fraction = year_elapsed/year_duration

        decimal_dates.append(date.year + fraction)

    return decimal_dates
        

def generate_csvs_for_scotti(csv_file, id_col, host_col, culture_date_col,  
                             first_date_col, last_date_col, 
                             alt_first_date_col=None, outdir='.',outpref=None):
    """Takes as input a csv file with sample ids, host ids, and dates, and 
    makes csvs to be used in scotti.

    The input csv should containt a column for each of the following: sample 
    id, host id, first date, last date. More columns are okay. The output csvs
    will be (id,host), (id,culture_date), (host,first_date,last_date)

    Args:
        csv_file: csv file with sample id, host id, first date, last date
        id_col: column name for sample ids
        host_col: column name for host ids
        culture_date_col: column name for culture dates
        first_date_col: column name for first date
        last_date_col: column name for last date
        alt_first_date_col: alternative first date (if first_date_col value is NA)
        outdir: output directory (default: current working directory)
        outpref: output file prefix (defaut: input file base name)

    Output:
        host csv: (id,host)
        culture date csv: (id,culture_date)
        first last date csv: (host,first_date,last_date)
    """ 

    # Read in data
    columns = defaultdict(list)

    with open(csv_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            for (k,v) in row.items():
                columns[k].append(v)
    
    # Get columns of interest
    ids = columns[id_col]
    hosts = columns[host_col]
    culture_dates = columns[culture_date_col]
    first_date = columns[first_date_col]
    last_date = columns[last_date_col]
    # If some first dates are NA (ex. no first negative culture), 
    # use alternative first date (ex. first positive culture)
    if alt_first_date_col is not None:
        alt_first_date = columns[alt_first_date_col]
        for i,d in enumerate(first_date):
            if d == 'NA':
                first_date[i] = alt_first_date[i]

    # Convert dates to decimal format
    culture_dates_dec = convert_ymd_to_decimal(culture_dates)
    first_date_dec = convert_ymd_to_decimal(first_date)
    last_date_dec = convert_ymd_to_decimal(last_date)

    # Set output file names
    if outpref is None:
        outpref = csv.split('.')[0]
    hosts_csv = outdir + '/' + outpref + '_hosts.csv'
    dates_csv = outdir + '/' + outpref + '_dates.csv'
    hostTimes_csv = outdir + '/' + outpref + '_hostTimes.csv'

    # Write hosts output file
    with open(hosts_csv, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(ids, hosts))

    # Write dates output file
    with open(dates_csv, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(ids, culture_dates_dec))

    # Get unique hosts and first/last date
    idxs = [idx for idx, item in enumerate(hosts) if item not in hosts[:idx]]
    uniq_hosts = [hosts[i] for i in idxs]
    uniq_first_date  = [first_date_dec[i] for i in idxs]
    uniq_last_date  = [last_date_dec[i] for i in idxs]

    # Write hostTimes output file
    with open(hostTimes_csv, 'w') as f:
        writer = csv.writer(f)
        writer.writerows(zip(uniq_hosts, uniq_first_date, uniq_last_date))

    return([hosts_csv, dates_csv, hostTimes_csv])


def subset_scotti_csvs(hosts, dates, hostTimes, fasta):
    """Subsets scotti csvs to only include samples and hosts in the fasta file.

    csv files for each input fasta will be generated.

    Args:
        hosts: csv where the first column is the sample id and the 
            second column is the host
        dates: csv where the first column is the sample id and the
            second column is the sampling date 
        hostTimes: csv where the first column is the host, the second
            column is the earliest time the host is infectable, and
            the third column is the latest time the host is infective
        fasta: list of fasta file(s) to make scotti xml(s)

    Output:

    Returns:

    """

    hosts_dict = defaultdict(list)

    with open(hosts) as f:
        reader = csv.reader(f)
        for row in reader:
            for (i,v) in enumerate(row):
                hosts_dict[i].append(v)

    i_hosts = hosts_dict[0]
    h_hosts = hosts_dict[1]

    dates_dict = defaultdict(list)

    with open(dates) as f:
        reader = csv.reader(f)
        for row in reader:
            for (i,v) in enumerate(row):
                dates_dict[i].append(v)

    i_dates = dates_dict[0]
    d_dates = dates_dict[1]

    times_dict = defaultdict(list)

    with open(hostTimes) as f:
        reader = csv.reader(f)
        for row in reader:
            for (i,v) in enumerate(row):
                times_dict[i].append(v)

    h_times = times_dict[0]
    t1_times = times_dict[1]
    t2_times = times_dict[2]

    hosts_renamed = hosts.split('.')[0] + '.renamed.csv'
    dates_renamed = dates.split('.')[0] + '.renamed.csv'
    times_renamed = hostTimes.split('.')[0] + '.renamed.csv'
    fasta_renamed = [fa.split('.')[0] + '.renamed.fa' for fa in fasta]

    # For each fasta file
    for fa in fasta:
        # Get ids from fasta file
        fa_ids = []
        with open(fa) as f:
            for line in f:
                if re.search('^>', line) is not None:
                    fa_ids.append(re.sub('^>', '', line))
        
         
     
