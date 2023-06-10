# read genes that have a ko assignment from A.kegg50.txt file in the annotations directory
# process genes/A.contig2gene.gff and contigs/A.contig2read.gff for getting the counts for each gene

#
from __future__ import division
import re
import os
import sys
import time
import glob


KEGG_folder = "./kegg_files/" 
output_folder = "./metabolic_nets/"
data_folder = "./DATA/"

excluded_pathways = ["01100", "01110", "01120"]

exclude_set = []



print "reading reaction file \n"
reactionFile = open(KEGG_folder+"reaction", "r")
RNEnzyme = {}
uniqueEnzyme = {}
reactionRPAIR = {}
RPAIRLINE = False
ENZYMELINE = False
for line in reactionFile:
    #line = line.strip()
    ENTRYLINE = re.match('ENTRY\s+', line)
    if ENTRYLINE:
        RPAIRLINE = False
        RN_NUMBER = re.search('ENTRY\s+(\S+)\s+.*', line)
        RN_NUMBER = RN_NUMBER.group(1)

    if not ENZYMELINE:
        ENZYMELINE = re.match('ENZYME\s+', line)
    else:
        ENZYMELINE = re.match('\s+', line)

    if ENZYMELINE:
        RPAIRLINE = False
        #EZ_NUMBER = re.search('ENZYME\s+(\S+)\s*.*', line)
        #EZ_NUMBER = EZ_NUMBER.group(1)
        if re.match('ENZYME\s+', line):
            EZ_NUMBERS = re.search('ENZYME\s+([\S+\s*]+)', line)
        else:
            EZ_NUMBERS = re.search('\s+([\S+\s*]+)', line)
        EZ_NUMBERS = EZ_NUMBERS.group(1)
        EZ_NUMBERS = EZ_NUMBERS.split()
        for EZ_NUMBER in EZ_NUMBERS:
            if EZ_NUMBER not in uniqueEnzyme:
                #print 'New EC number: '+ EZ_NUMBER+ '\n'
                uniqueEnzyme[ EZ_NUMBER ] = 0
            uniqueEnzyme[ EZ_NUMBER ] += 1
            if EZ_NUMBER not in RNEnzyme:
                RNEnzyme[ EZ_NUMBER ] = set()
        
            if (RN_NUMBER != ''):
                RNEnzyme[ EZ_NUMBER ].add( RN_NUMBER )

        #print "%s %s\n" % (EZ_NUMBER , RNEnzyme[ EZ_NUMBER ])
        RN_NUMBER = ''
        continue

    if not RPAIRLINE:
        RPAIRLINE = re.match('RPAIR\s+', line)
    else:
        RPAIRLINE = re.match('\s+', line)

    if RPAIRLINE:
        #print line + '\n'
        if re.match('RPAIR\s+', line):
            RPAIR_SEARCH = re.search('RPAIR\s+\S+\s+(\S+)\s+(\S+).*', line)
        else:
            RPAIR_SEARCH = re.search('\s*\S+\s+(\S+)\s+(\S+).*', line)
        RPAIR = RPAIR_SEARCH.group(1)
        RPAIR_TYPE = RPAIR_SEARCH.group(2)
        #if RPAIR_TYPE == 'main':
        if RN_NUMBER not in reactionRPAIR:
            reactionRPAIR[RN_NUMBER] = set()
        #[RPAIR_Product, RPAIR_Substrate] = RPAIR.split('_')
        reactionRPAIR[RN_NUMBER].add(RPAIR)
#       print 'Reaction ' + RN_NUMBER + ' has a RPAIR: ' + RPAIR +  ' of type '+  RPAIR_TYPE + '\n' 
        
reactionFile.close()

print "reaction2EQ file\n"
#reaction2EQFile = open(KEGG_folder+"reaction_list", "r")
reaction2EQFile = open(KEGG_folder+"/reaction_mapformula.lst", "r")
reaction2EQ = {}
for line in reaction2EQFile:
    reactionNumber = line.strip().split(":")[0].strip()
    pathwayID = line.strip().split(":")[1].strip()
    reactionFormula = line.strip().split(":")[2].strip()
    #for currency_metabolite in currency_compounds:
    #    if currency_metabolite in reactionFormula:
            
    if reactionNumber not in reaction2EQ: reaction2EQ[reactionNumber] = {}
    reaction2EQ[reactionNumber][pathwayID] = reactionFormula


all_EC_Counts_files = glob.glob( "%s/*.EC.abundance" % data_folder)

for sample_file_name in all_EC_Counts_files:
    sample_name = os.path.splitext(os.path.split(str(sample_file_name))[1])[0]
    #sample_name = sample_name.split(".")[0]
    sample_name = sample_name[:-3] # removing .EC fron the end
    if sample_name in exclude_set: continue
    print 'processing '+ sample_name + '\n'
    print "reading gene abundance\n"
    EC_abundance = {} 
    all_EC_abundance_info = open(sample_file_name)
    for line in all_EC_abundance_info:
        One_EC_abundance_info = line.split(":")
        EC = One_EC_abundance_info[0].strip()
        #EC = ''.join(EC.split(' '))
        EC = EC.split(' ')[1].strip()
        if EC not in EC_abundance: EC_abundance[EC] = 0
        EC_abundance[EC] = float(One_EC_abundance_info[1].strip())

    print "writing metabolic_net_file\n"
    Reaction_output = open(output_folder + sample_name + ".metabolic.net", "w")

    AllReactions = {}
                
    #wrirte reactions to file
    for EC_NUMBER in EC_abundance:
        if EC_NUMBER in RNEnzyme:
            Reactions = RNEnzyme[ EC_NUMBER]
        else:
            print 'no reactions for EC number: '+ EC_NUMBER+ '\n'
            continue
                    
        #count_adjustment = len(KO_EC[KO]) * len(Reactions) # split the weight between all posible reactions and KOs
        count_adjustment = 1
        for RN_NUMBER in Reactions:
            if RN_NUMBER not in AllReactions:
                AllReactions[RN_NUMBER] = 0;
            AllReactions[RN_NUMBER] += EC_abundance[EC_NUMBER] / count_adjustment

    for RN_NUMBER in AllReactions:
        if RN_NUMBER in reactionRPAIR: # if ( RN_NUMBER in reaction2EQ) and (RN_NUMBER in reactionRPAIR): #

            for i in range(int(round(AllReactions[RN_NUMBER] ))):
                for RPAIR in reactionRPAIR[RN_NUMBER]:
                    [RPAIR_Substrate, RPAIR_Product] = RPAIR.split('_')
                    Reaction_output.write(RPAIR_Substrate+','+RN_NUMBER+','+RPAIR_Product+'\n')
                Reaction_output.write('\n')

        
    Reaction_output.close()

