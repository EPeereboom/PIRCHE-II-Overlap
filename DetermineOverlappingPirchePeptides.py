import argparse
from os import makedirs
from os.path import join, isdir


def extractPeptidesFromColumn(peptideColumn=None, sampleID=None):
    peptides = []
    peptideLookup = {}

    #print('Peptide Column:' + str(peptideColumn))
    if (len(peptideColumn.strip()) > 0):
        # Lots of peptides in these columns, separated by |
        presentedPeptides = peptideColumn.split('|')
        for presentedPeptide in presentedPeptides:
            # print('presentedPeptide:' + str(presentedPeptide))
            # print('it has ' + str(len(presentedPeptide.strip().split(' '))) + ' tokens.')

            # Peptides have 3 parts.
            ninemer, fifteenmer, score, weight = presentedPeptide.strip().split(' ')
            score = float(score.replace('(','').replace(')',''))
            weight = float(weight.replace('[', '').replace(']', ''))

            # detect duplicate peptides
            if(ninemer in peptideLookup.keys() and str(peptideLookup[ninemer]) != str((ninemer,fifteenmer,score,weight))):
                pass
                #print('Warning!: duplicate 9mer found for Sample=' + str(sampleID) + ' and ninemer=' + str(ninemer))
                (beforeNinemer, beforeFifteenmer, beforeScore, beforeWeight) = peptideLookup[ninemer]
                #print('Ninemer1: ' + str((beforeNinemer, beforeFifteenmer, beforeScore, beforeWeight)))
                #print('Ninemer2: ' + str((ninemer, fifteenmer, score, weight)))

                # Use the fifteenmer and ic50 from the peptide w bigger weight.
                # Sum the weights
                if(float(beforeWeight) > float(weight)):
                    combinedFifteenmer=beforeFifteenmer
                    combinedScore= beforeScore
                    combinedWeight = beforeWeight + weight
                else:
                    combinedFifteenmer=fifteenmer
                    combinedScore= score
                    combinedWeight = beforeWeight + weight

                #print('CombinedNinemer: ' + str((ninemer, combinedFifteenmer, combinedScore, combinedWeight)))
                if(combinedWeight > 1):
                    raise Exception ('Weight should not be > 1)')
                peptideLookup[ninemer] = (ninemer, combinedFifteenmer, combinedScore, combinedWeight)

            else:
                peptideLookup[ninemer] = (ninemer, fifteenmer, score, weight)

    for ninemer in peptideLookup.keys():
        peptides.append(peptideLookup[ninemer])

    else:
        pass
    return peptides

def getSharedPeptides(immunizerPeptides=None, recallPeptides=None):
    allOverlap9mers=set()
    immunizer15mers=set()
    recall15mers=set()

    for immunizerPeptide in immunizerPeptides:
        (immunizerNinemer, immunizerFifteenmer, immunizerScore, immunizerWeight) = immunizerPeptide
        for recallPeptide in recallPeptides:
            (recallNinemer, recallFifteenmer, recallScore, recallWeight) = recallPeptide

            if(immunizerNinemer == recallNinemer):

                if(immunizerWeight == recallWeight):
                    # Great!
                    pass
                else:
                    # This is not a problem:
                    # but that is perfectly fine. The weight in the immunizer file is always 1,
                    # as we provided the high res patient DRB1 typing and the donor (bead) is already high res.
                    # For the recall epitope, it's also based on the uncertainty of the donor's typing. That's why we
                    # also take the weight of the recall epitope and not the weight of the immunizer
                    #print('Warning, the immunizer has a different weight than the recall epitope for peptide (' + immunizerNinemer + ':' + immunizerFifteenmer + '/' + recallFifteenmer + ')')
                    #print('immunizerWeight:' + str(immunizerWeight))
                    #print('recallWeight:' + str(recallWeight))
                    pass

                allOverlap9mers.add(immunizerNinemer + ':' + str(recallWeight))
                immunizer15mers.add(str(immunizerFifteenmer) + ':' + str(immunizerWeight))
                recall15mers.add(str(recallFifteenmer) + ':' + str(recallWeight))

    # For the excel sheet we need
    # 9mers;immunizer15mers;recall15mers
    sharedPeptideString = '"'
    for overlap9mer in sorted(list(allOverlap9mers)):
        sharedPeptideString = sharedPeptideString + overlap9mer + '\n'
    sharedPeptideString = sharedPeptideString.strip('\n') + '"' + separator + '"'

    for immunizer15mer in sorted(list(immunizer15mers)):
        sharedPeptideString = sharedPeptideString + immunizer15mer + '\n'
    sharedPeptideString = sharedPeptideString.strip('\n') + '"' + separator + '"'

    for recall15mer in sorted(list(recall15mers)):
        sharedPeptideString = sharedPeptideString + recall15mer + '\n'
    sharedPeptideString = sharedPeptideString.strip('\n') + '"'

    return sharedPeptideString, allOverlap9mers, immunizer15mers,recall15mers

def writePeptideList(outputFileName=None, recallEpitopes_drb1_1=None, recallEpitopes_drb1_2=None,immunizers_drb1_1=None,immunizers_drb1_2=None):
    print('Writing list of peptides:' + str(outputFileName))
    with open(outputFileName, 'w') as output:
        output.write('sample_id' + separator + 'file' + separator + 'drb1_?' + separator +  '9mer'
                     + separator + '15mer' + separator + 'ic50' + separator + 'weight' + newline)

        uniqueRows = set()

        for sampleID in recallEpitopes_drb1_1.keys():
            for recallEpitope in recallEpitopes_drb1_1[sampleID]:
                uniqueRows.add((sampleID,'recall','drb1_1',str(recallEpitope[0]),str(recallEpitope[1]),str(recallEpitope[2]),str(recallEpitope[3])))

        for sampleID in recallEpitopes_drb1_2.keys():
            for recallEpitope in recallEpitopes_drb1_2[sampleID]:
                uniqueRows.add((sampleID, 'recall', 'drb1_2', str(recallEpitope[0]), str(recallEpitope[1]),str(recallEpitope[2]), str(recallEpitope[3])))

        for sampleID in immunizers_drb1_1.keys():
            for hlaAllele in immunizers_drb1_1[sampleID]:
                for immunizerEpitope in immunizers_drb1_1[sampleID][hlaAllele]:
                    uniqueRows.add((sampleID, 'immunizer', 'drb1_1', str(immunizerEpitope[0]), str(immunizerEpitope[1]), str(immunizerEpitope[2]), str(immunizerEpitope[3])))

        for sampleID in immunizers_drb1_2.keys():
            for hlaAllele in immunizers_drb1_2[sampleID]:
                for immunizerEpitope in immunizers_drb1_2[sampleID][hlaAllele]:
                    uniqueRows.add((sampleID, 'immunizer', 'drb1_2', str(immunizerEpitope[0]), str(immunizerEpitope[1]),str(immunizerEpitope[2]), str(immunizerEpitope[3])))

        for uniqueRow in sorted(list(uniqueRows)):
            #print('uniqueRow:' + str(uniqueRow))
            output.write(separator.join(uniqueRow) + newline)

def loadImmunizersCSV(inputFile=None):
    print('Loading Immunizer File' + str(inputFile))
    # Load Immunizers CSV (Patient, Bead, Peptides)
    # CSV with Patient, HLA type, Peptides(could be spread on lots of columns)
    # Store as a set of peptides in a nested dictionary. immunizers[sampleID][allele]=set()

    immunizers_drb1_1 = {}
    immunizers_drb1_2 = {}
    immunizersBeadLookup = {}
    with open(inputFile, "r") as input:
        rowCount = 0
        for row in input:
            rowCount += 1
            if rowCount > 1:
                # print('processing row ' + str(row))
                columns = row.split(separator)
                sampleID = columns[0].strip()
                beadID = columns[1].strip()
                hlaAllele = columns[2].strip()
                hlaLocus = columns[3].strip()
                valueLum = columns[4].strip()
                lra = columns[5].strip()
                lraMfi = columns[6].strip()

                # MFIs and whatnot...

                if (sampleID is not None and len(sampleID) > 0):
                    # The pirch II peptides are in the remaining columns, starting with 2
                    drb1_1_peptides = extractPeptidesFromColumn(peptideColumn=columns[7])
                    drb1_2_peptides = extractPeptidesFromColumn(peptideColumn=columns[8])

                    # Have we seen this sample ID before?
                    if (sampleID not in immunizers_drb1_1):
                        immunizers_drb1_1[sampleID] = {}
                        immunizers_drb1_2[sampleID] = {}
                        immunizersBeadLookup[sampleID] = {}

                    # Store this set of peptides for later.
                    immunizers_drb1_1[sampleID][hlaAllele] = drb1_1_peptides
                    immunizers_drb1_2[sampleID][hlaAllele] = drb1_2_peptides
                    # immunizersBeadLookup[sampleID][hlaAllele]=(beadID,pircheII,pircheIISd)
                    immunizersBeadLookup[sampleID][hlaAllele] = (beadID, hlaLocus, valueLum, lra, lraMfi)
                    if (verbose):
                        print('Stored Immunizer for sample ' + sampleID + ' and hla allele ' + hlaAllele + ' :' + str(
                            drb1_1_peptides))
                        print('Stored Immunizer for sample ' + sampleID + ' and hla allele ' + hlaAllele + ' :' + str(
                            drb1_2_peptides))
                else:
                    print('Row without a sample ID:' + str(row))
    return immunizers_drb1_1, immunizers_drb1_2, immunizersBeadLookup

def loadRecallEpitopes(inputFile=None):
    print('Loading Recall File' + str(inputFile))
    # Load Recall Epitopes
    # It's a CSV of Sample ID, and a list of Peptides.
    recallEpitopes_drb1_1 = {}
    recallEpitopes_drb1_2 = {}
    with open(inputFile, "r") as input:
        rowCount = 0
        for row in input:
            rowCount += 1
            if rowCount > 1:
                columns = row.split(separator)
                sampleID = columns[0].strip()

                if(sampleID is not None and len(sampleID) > 0):
                    # The recall peptides are in the remaining columns, starting with 1
                    # Store this set of peptides for later.
                    recallEpitopes_drb1_1[sampleID] = extractPeptidesFromColumn(peptideColumn=columns[1])
                    recallEpitopes_drb1_2[sampleID] = extractPeptidesFromColumn(peptideColumn=columns[2])
                    #print('recallEpitopes_drb1_1 for sample ' + sampleID + ' :' + str(recallEpitopes_drb1_1[sampleID]))
                    #print('recallEpitopes_drb1_2 for sample ' + sampleID + ' :' + str(recallEpitopes_drb1_2[sampleID]))

    return recallEpitopes_drb1_1,recallEpitopes_drb1_2

def loadHomozygosityFile(inputFile=None):
    print('Loading Homozygosity File' + str(inputFile))
    homozygosityLookup = {}
    with open(inputFile, "r") as input:
        rowCount = 0
        for row in input:
            rowCount += 1
            if rowCount > 1:
                columns = row.split(separator)
                sampleID = columns[0].strip()
                # TODO: Do something with they typings? It's here if i need it.
                DRB1_1 = columns[1].strip()
                DRB1_2 = columns[2].strip()
                if(columns[3].strip() == '1'):
                    homozygosityLookup[sampleID] = True
                elif(columns[3].strip() == '0'):
                    homozygosityLookup[sampleID] = False
                else:
                    print('Strange Row!:' + str(row))
                    raise Exception ('Expected to see either a 1 or 0 for homozygosity in the input file.')


    return homozygosityLookup

def writeOverlappingPeptides(outputFileName=None, immunizerSampleIDs=None, immunizers_drb1_1=None, immunizers_drb1_2=None, recallEpitopes_drb1_1=None, recallEpitopes_drb1_2=None, immunizersBeadLookup=None):
    print('Finding the overlapping peptides between immunizers and recall:' + str(outputFileName))
    # Output, a new CSV
    # So the output columns would be:
    # Sample_ID (same as input file)
    # Bead_ID (same as input file)
    # Mol_LUM (same as input file)
    # PIRCHE_II (same as input file)
    # PIRCHE_II_SD (same as input file)
    # Peptides (so the overlapping peptides, can be in 1 cell. Each cell in this column will look something like this: GRCVDGLRR VYLEGRCVDGLRRYL [0.9797] | FRXMYGLRX VYLEFRXMYGLRXYL [0.8397] | GRCVDGXRR VYLEGRCVDGXRRYL [0.9292])
    # Sum_weights
    with open(outputFileName, 'w') as output:
        # header line
        #output.write('Sample_ID' + separator + 'Bead_ID' + separator + 'HLA' + separator + 'PIRCHE_II' + separator + 'PIRCHE_II_SD' + separator + 'Peptides' + separator + 'Sum_Weights' + newline)
        output.write('sampleID' + separator + 'BeadID' + separator + 'Mol_LUM dimers'
            + separator + 'Locus_LUM' + separator + 'Value_LUM'
            + separator + 'LRA'	+ separator + 'MFI/LRA_MFI' + separator + 'DRB1_1_Overlap_9mers'
            + separator + 'DRB1_1_Immunizer_15mers'	+ separator + 'DRB1_1_Recall_15mers'
            + separator + 'DRB1_2_Overlap_9mers'
            + separator + 'DRB1_2_Immunizer_15mers' + separator + 'DRB1_2_Recall_15mers'
            + newline)

        drb1_1_sharedPeptides = {}
        drb1_2_sharedPeptides = {}

        for sampleID in sorted(list(immunizerSampleIDs)):
            drb1_1_sharedPeptides[sampleID] = {}
            drb1_2_sharedPeptides[sampleID] = {}
            for hlaAllele in sorted(list(immunizers_drb1_1[sampleID])):
                sharedDrb1_1_String, Drb1_1_Overlap9mers, Drb1_1_immunizer15mers, Drb1_1_recall15mers = getSharedPeptides(immunizerPeptides=immunizers_drb1_1[sampleID][hlaAllele], recallPeptides=recallEpitopes_drb1_1[sampleID])
                sharedDrb1_2_String, Drb1_2_Overlap9mers, Drb1_2_immunizer15mers, Drb1_2_recall15mers = getSharedPeptides(immunizerPeptides=immunizers_drb1_2[sampleID][hlaAllele], recallPeptides=recallEpitopes_drb1_2[sampleID])

                drb1_1_sharedPeptides[sampleID][hlaAllele] = (Drb1_1_Overlap9mers, Drb1_1_immunizer15mers, Drb1_1_recall15mers)
                drb1_2_sharedPeptides[sampleID][hlaAllele] = (Drb1_2_Overlap9mers, Drb1_2_immunizer15mers, Drb1_2_recall15mers)

                #(beadID, pircheII, pircheIISd) = immunizersBeadLookup[sampleID][hlaAllele]
                (beadID,hlaLocus,valueLum,lra,lraMfi) = immunizersBeadLookup[sampleID][hlaAllele]

                outputLine = (sampleID + separator + beadID + separator + hlaAllele + separator + hlaLocus + separator
                    + valueLum + separator + lra + separator + lraMfi + separator)

                outputLine = outputLine + sharedDrb1_1_String + separator
                outputLine = outputLine + sharedDrb1_2_String + newline

                output.write(outputLine)

        return drb1_1_sharedPeptides, drb1_2_sharedPeptides

def writeCombinedOverlappingPeptides(outputFileName=None,drb1_1_sharedPeptides=None, drb1_2_sharedPeptides=None,immunizerSampleIDs=None):
    print('Combining the overlapping peptides per-bead for sample:' + str(outputFileName))
    drb1_1_combinedPeptides = {}
    drb1_2_combinedPeptides = {}
    with open(outputFileName, 'w') as output:
        # header line
        # output.write('Sample_ID' + separator + 'Bead_ID' + separator + 'HLA' + separator + 'PIRCHE_II' + separator + 'PIRCHE_II_SD' + separator + 'Peptides' + separator + 'Sum_Weights' + newline)
        output.write('sampleID' + separator + 'DRB1_1_Overlap_9mers'
                     + separator + 'DRB1_1_Immunizer_15mers' + separator + 'DRB1_1_Recall_15mers'
                     + separator + 'DRB1_2_Overlap_9mers'
                     + separator + 'DRB1_2_Immunizer_15mers' + separator + 'DRB1_2_Recall_15mers'
                     + newline)



        for sampleID in sorted(list(immunizerSampleIDs)):

            drb1_1_combined9mers=set()
            drb1_1_combined_immunizer_15mers=set()
            drb1_1_combined_recall_15mers=set()
            drb1_2_combined9mers=set()
            drb1_2_combined_immunizer_15mers=set()
            drb1_2_combined_recall_15mers=set()
            for hlaAllele in sorted(list(drb1_1_sharedPeptides[sampleID])):
                #print('tuple:' + str(drb1_1_sharedPeptides[sampleID][hlaAllele]))
                #print('9mers' + str(drb1_1_sharedPeptides[sampleID][hlaAllele][0]))
                drb1_1_combined9mers = drb1_1_combined9mers.union(drb1_1_sharedPeptides[sampleID][hlaAllele][0])
                drb1_1_combined_immunizer_15mers = drb1_1_combined_immunizer_15mers.union(drb1_1_sharedPeptides[sampleID][hlaAllele][1])
                drb1_1_combined_recall_15mers = drb1_1_combined_recall_15mers.union(drb1_1_sharedPeptides[sampleID][hlaAllele][2])
                drb1_2_combined9mers = drb1_2_combined9mers.union(drb1_2_sharedPeptides[sampleID][hlaAllele][0])
                drb1_2_combined_immunizer_15mers = drb1_2_combined_immunizer_15mers.union(drb1_2_sharedPeptides[sampleID][hlaAllele][1])
                drb1_2_combined_recall_15mers = drb1_2_combined_recall_15mers.union(drb1_2_sharedPeptides[sampleID][hlaAllele][2])

            drb1_1_combinedPeptides[sampleID] = (drb1_1_combined9mers,drb1_1_combined_immunizer_15mers,drb1_1_combined_recall_15mers)
            drb1_2_combinedPeptides[sampleID] = (drb1_2_combined9mers,drb1_2_combined_immunizer_15mers,drb1_2_combined_recall_15mers)


            # Write a line
            outputLine = (sampleID + separator + '"')

            for overlap9mer in sorted(list(drb1_1_combined9mers)):
                outputLine = outputLine + overlap9mer + ','
            outputLine = outputLine.strip(',') + '"' + separator + '"'

            for immunizer15mer in sorted(list(drb1_1_combined_immunizer_15mers)):
                outputLine = outputLine + immunizer15mer + ','
            outputLine = outputLine.strip(',') + '"' + separator + '"'

            for recall15mer in sorted(list(drb1_1_combined_recall_15mers)):
                outputLine = outputLine + recall15mer + ','
            outputLine = outputLine.strip(',') + '"' + separator + '"'

            for overlap9mer in sorted(list(drb1_2_combined9mers)):
                outputLine = outputLine + overlap9mer + ','
            outputLine = outputLine.strip(',') + '"' + separator + '"'

            for immunizer15mer in sorted(list(drb1_2_combined_immunizer_15mers)):
                outputLine = outputLine + immunizer15mer + ','
            outputLine = outputLine.strip(',') + '"' + separator + '"'

            for recall15mer in sorted(list(drb1_2_combined_recall_15mers)):
                outputLine = outputLine + recall15mer + ','
            outputLine = outputLine.strip(',') + '"' + newline

            output.write(outputLine)


    return drb1_1_combinedPeptides, drb1_2_combinedPeptides

def writeScores(outputFileName=None, immunizerSampleIDs=None, drb1_1_combinedPeptides=None, drb1_2_combinedPeptides=None):
    print('Writing Scores:' + str(outputFileName))

    with open(outputFileName, 'w') as output:
        # header line
        # output.write('Sample_ID' + separator + 'Bead_ID' + separator + 'HLA' + separator + 'PIRCHE_II' + separator + 'PIRCHE_II_SD' + separator + 'Peptides' + separator + 'Sum_Weights' + newline)
        output.write('sampleID' + separator + 'DRB1_1_Overlap_9mers_Sum_Weights'
                     + separator + 'DRB1_1_Immunizer_Sum_Weights' + separator + 'DRB1_1_Recall_Sum_Weights'
                     + separator + 'DRB1_2_Overlap_9mers_Sum_Weights'
                     + separator + 'DRB1_2_Immunizer_Sum_Weights' + separator + 'DRB1_2_Recall_Sum_Weights'
                     + newline)

        for sampleID in sorted(list(immunizerSampleIDs)):
            # Write a line
            drb1_1_9mer_sum = 0#len(list(drb1_1_combinedPeptides[sampleID][0]))
            drb1_1_immunizer_sum = 0
            drb1_1_recall_sum = 0
            drb1_2_9mer_sum = 0#len(list(drb1_2_combinedPeptides[sampleID][0]))
            drb1_2_immunizer_sum = 0
            drb1_2_recall_sum = 0

            for peptide9mer in list(drb1_1_combinedPeptides[sampleID][0]):
                peptide9merSeq, peptide9merWeight = peptide9mer.split(':')
                drb1_1_9mer_sum += float(peptide9merWeight)

            for immunizer15mer in list(drb1_1_combinedPeptides[sampleID][1]):
                immSeq, immWeight = immunizer15mer.split(':')
                drb1_1_immunizer_sum += float(immWeight)

            for recall15mer in list(drb1_1_combinedPeptides[sampleID][2]):
                recSeq, recWeight = recall15mer.split(':')
                drb1_1_recall_sum += float(recWeight)

            for peptide9mer in list(drb1_2_combinedPeptides[sampleID][0]):
                peptide9merSeq, peptide9merWeight = peptide9mer.split(':')
                drb1_2_9mer_sum += float(peptide9merWeight)

            for immunizer15mer in list(drb1_2_combinedPeptides[sampleID][1]):
                immSeq, immWeight = immunizer15mer.split(':')
                drb1_2_immunizer_sum += float(immWeight)

            for recall15mer in list(drb1_2_combinedPeptides[sampleID][2]):
                recSeq, recWeight = recall15mer.split(':')
                drb1_2_recall_sum += float(recWeight)

            outputLine = (sampleID + separator
                + str(drb1_1_9mer_sum) + separator + str(drb1_1_immunizer_sum) + separator + str(drb1_1_recall_sum) + separator
                + str(drb1_2_9mer_sum) + separator + str(drb1_2_immunizer_sum) + separator + str(drb1_2_recall_sum) + newline
                )

            output.write(outputLine)

def writeFinalScore(outputFileName=None, immunizerSampleIDs=None, drb1_1_combinedPeptides=None, drb1_2_combinedPeptides=None, homozygosityLookup=None):
    print('Writing Final Score:' + str(outputFileName))

    with open(outputFileName, 'w') as output:
        output.write('sampleID' + separator + 'homozygous' + separator + 'Combined_recall_Score'
            + separator + 'DRB1_1_9mers' + separator + 'DRB1_2_Recall_9mers'
            + newline)

        for sampleID in sorted(list(immunizerSampleIDs)):
            # Write a line
            drb1_1_9mer_sum = 0
            drb1_2_9mer_sum = 0


            isHomozygous = homozygosityLookup[sampleID]

            # TODO: If homozygous, use one column. If not calculate 2 scores and sum them.
            for recall15mer in list(drb1_1_combinedPeptides[sampleID][0]):
                recSeq, recWeight = recall15mer.split(':')
                drb1_1_9mer_sum += float(recWeight)

            for recall15mer in list(drb1_2_combinedPeptides[sampleID][0]):
                recSeq, recWeight = recall15mer.split(':')
                drb1_2_9mer_sum += float(recWeight)

            if(isHomozygous):
                recallScore = drb1_1_9mer_sum
            else:
                recallScore = drb1_1_9mer_sum + drb1_2_9mer_sum



            outputLine = (sampleID + separator + str(isHomozygous) + separator +
                str(recallScore) + separator +
                ','.join(sorted(list(drb1_1_combinedPeptides[sampleID][0]))) + separator +
                ','.join(sorted(list(drb1_2_combinedPeptides[sampleID][0]))) + newline)

            output.write(outputLine)

def filterImmunizers(immunizers_drb1=None, immunizersBeadLookup=None, filterMethod=None):
    print('Filtering Immunizers by method: ' + str(filterMethod))
    filterMethod = str(filterMethod).upper()
    #("NONE", "5PERCENT", "10PERCENT", "CUTOFFRATIO", "HIGHESTMFI"):
    filteredImmunizers= {}

    if(filterMethod=='NONE'):
        filteredImmunizers=immunizers_drb1

    elif(filterMethod=='HIGHESTRATIO'):
        #print('Input Keys:' + str())
        for sampleID in immunizers_drb1:
            maxAlleles=[]
            maxRatio=0
            for hlaAllele in (immunizers_drb1[sampleID]):
                currentRatio = float(immunizersBeadLookup[sampleID][hlaAllele][4])
                if(currentRatio > maxRatio):
                    maxAlleles=[hlaAllele]
                    maxRatio=currentRatio
                elif(currentRatio == maxRatio):
                    maxAlleles.append(hlaAllele)
                else:
                    pass
            if(len(maxAlleles) > 1):
                pass
                print('for sample ' + str(sampleID) + ' two beads were found with the same max mfi ratio:' + str(maxAlleles))
            if(maxRatio==0):
                raise Exception('No maximum MFI bead was found for sample ' + str(sampleID))
            filteredImmunizers[sampleID] = {}
            for maxAllele in maxAlleles:
                filteredImmunizers[sampleID][maxAllele]=immunizers_drb1[sampleID][maxAllele]

    elif (filterMethod in['5PERCENTMFI','1PERCENTMFI','10PERCENTMFI']):
        if(filterMethod=='5PERCENTMFI'):
            percentageCutoff=.05
        elif (filterMethod == '1PERCENTMFI'):
            percentageCutoff = .01
        elif (filterMethod == '10PERCENTMFI'):
            percentageCutoff = .10
        else:
            raise Exception ('Unexpected cutoff Value!' + str(filterMethod))

        for sampleID in immunizers_drb1:
            # Calculate a cutoff, per sample
            maxMFI = 0
            minMFI = 9999999999 #Arbitrarily large, this is bad practice. Oh well.
            for hlaAllele in (immunizers_drb1[sampleID]):
                currentMFI = float(immunizersBeadLookup[sampleID][hlaAllele][2])
                if (currentMFI > maxMFI):
                    maxMFI = currentMFI
                if (currentMFI < minMFI):
                    minMFI = currentMFI

            cutoffMFI = minMFI + (maxMFI - minMFI)*(1-percentageCutoff)
            print('in sample ' + str(sampleID) + ' For a min MFI ' + str(minMFI) + ' and max MFI ' + str(maxMFI) + ' i calculated cuttoff ' + str(cutoffMFI))

            filteredImmunizers[sampleID] = {}
            # Add beads if they are >= cutoff
            for hlaAllele in (immunizers_drb1[sampleID]):
                currentMFI = float(immunizersBeadLookup[sampleID][hlaAllele][2])

                if currentMFI >= cutoffMFI:
                    filteredImmunizers[sampleID][hlaAllele] = immunizers_drb1[sampleID][hlaAllele]
                    #print('in sample ' + str(sampleID) + ' one was bigger than cutoff')

    elif (filterMethod in['5PERCENTRATIO','1PERCENTRATIO','10PERCENTRATIO']):
        if(filterMethod=='5PERCENTRATIO'):
            percentageCutoff=.05
        elif (filterMethod == '1PERCENTRATIO'):
            percentageCutoff = .01
        elif (filterMethod == '10PERCENTRATIO'):
            percentageCutoff = .10
        else:
            raise Exception ('Unexpected cutoff Value!' + str(filterMethod))

        for sampleID in immunizers_drb1:
            # Calculate a cutoff, per sample
            maxRatio = 0
            minRatio = 9999999999 #Arbitrarily large, this is bad practice. Oh well.
            for hlaAllele in (immunizers_drb1[sampleID]):
                currentRatio = float(immunizersBeadLookup[sampleID][hlaAllele][4])
                if (currentRatio > maxRatio):
                    maxRatio = currentRatio
                if (currentRatio < minRatio):
                    minRatio = currentRatio

            cutoffRatio = minRatio + (maxRatio - minRatio)*(1-percentageCutoff)
            print('in sample ' + str(sampleID) + ' For a minRatio ' + str(minRatio) + ' and maxRatio ' + str(maxRatio) + ' i calculated cuttoff ' + str(cutoffRatio))

            filteredImmunizers[sampleID] = {}
            # Add beads if they are >= cutoff
            for hlaAllele in (immunizers_drb1[sampleID]):
                currentRatio = float(immunizersBeadLookup[sampleID][hlaAllele][4])

                if currentRatio >= cutoffRatio:
                    filteredImmunizers[sampleID][hlaAllele] = immunizers_drb1[sampleID][hlaAllele]
                    #print('in sample ' + str(sampleID) + ' one was bigger than cutoff')

    elif (filterMethod in['5PERCENTOVERALLRATIO','1PERCENTOVERALLRATIO','10PERCENTOVERALLRATIO']):
        if(filterMethod=='5PERCENTOVERALLRATIO'):
            percentageCutoff=.05
        elif (filterMethod == '1PERCENTOVERALLRATIO'):
            percentageCutoff = .01
        elif (filterMethod == '10PERCENTOVERALLRATIO'):
            percentageCutoff = .10
        else:
            raise Exception ('Unexpected cutoff Value!' + str(filterMethod))

        allRatios = []

        for sampleID in immunizers_drb1:
            for hlaAllele in (immunizers_drb1[sampleID]):
                currentRatio = float(immunizersBeadLookup[sampleID][hlaAllele][4])
                allRatios.append(currentRatio)

        print('Found a total of ' + str(len(allRatios)) + ' ratios.')
        sampleCount = int(len(allRatios) * percentageCutoff)
        topRatios = sorted(allRatios, reverse=True)[0:sampleCount]
        minRatio = min(topRatios)
        print('The top ' + str(percentageCutoff) + ' samples is ' + str(sampleCount))
        print('topsamples:' + str(topRatios))
        print('minRatio=' + str(minRatio))

        for sampleID in immunizers_drb1:
            filteredImmunizers[sampleID] = {}
            for hlaAllele in (immunizers_drb1[sampleID]):
                currentRatio = float(immunizersBeadLookup[sampleID][hlaAllele][4])
                if(currentRatio >= minRatio):
                    filteredImmunizers[sampleID][hlaAllele] = immunizers_drb1[sampleID][hlaAllele]
                    #print('Just added immunizer(ratio=' + str(currentRatio) + ') ' + str(immunizers_drb1[sampleID][hlaAllele]))


    elif (filterMethod in['5RATIOCUTOFF', '2RATIOCUTOFF','3RATIOCUTOFF']):
        if(filterMethod=='5RATIOCUTOFF'):
            ratioCutoff=5
        elif (filterMethod == '2RATIOCUTOFF'):
            ratioCutoff = 2
        elif (filterMethod == '3RATIOCUTOFF'):
            ratioCutoff = 3
        else:
            raise Exception ('Unexpected ratio cutoff Value!' + str(filterMethod))

        for sampleID in immunizers_drb1:
            filteredImmunizers[sampleID] = {}
            for hlaAllele in (immunizers_drb1[sampleID]):
                #print('lookup:' + str(immunizersBeadLookup[sampleID][hlaAllele]))
                currentRatio = float(immunizersBeadLookup[sampleID][hlaAllele][4])
                if(currentRatio >= ratioCutoff):
                    filteredImmunizers[sampleID][hlaAllele] = immunizers_drb1[sampleID][hlaAllele]
                    #print('in sample ' + str(sampleID) + ' one was bigger than cutoff')

    else:
        raise Exception('Unknown filter method:' + str(filterMethod))

    return filteredImmunizers


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", help="verbose operation", action="store_true")
    parser.add_argument("-i", "--immunizers", help="immunizers input file", required=True)
    parser.add_argument("-r", "--recall", help="recall epitopes input file", required=True)
    parser.add_argument("-o", "--output", help="output directory", required=True)
    parser.add_argument("-f", "--filter", help="Method to filter the beads", required=True)
    parser.add_argument("-H", "--homozygosity", help="homozygosity input file", required=True)


    args = parser.parse_args()

    # newline and separator characters, change if necessary.
    newline='\r\n'
    separator=';'

    verbose = args.verbose
    if verbose:
        print("verbose mode active")

    print('We are merging data now.')

    if not isdir(args.output):
        makedirs(args.output)

    # Load the data.
    immunizers_drb1_1, immunizers_drb1_2, immunizersBeadLookup= loadImmunizersCSV(inputFile=args.immunizers)
    recallEpitopes_drb1_1, recallEpitopes_drb1_2 = loadRecallEpitopes(inputFile=args.recall)


    # Double check, our sample IDs should be exactly the same.
    immunizerSampleIDs = set(immunizers_drb1_1.keys()).union(set(immunizers_drb1_2.keys()))
    recallSampleIDs = set(recallEpitopes_drb1_1.keys()).union(set(recallEpitopes_drb1_2.keys()))
    #print('immunizerSamples:' + str(sorted(list(immunizerSampleIDs))))
    #print('recallSampleIDs:' + str(sorted(list(recallSampleIDs))))
    if(len(immunizerSampleIDs - recallSampleIDs) > 0):
        print('Warning! These sample ids are Missing in the recall epitopes file' + str((immunizerSampleIDs - recallSampleIDs)))
    if(len(recallSampleIDs - immunizerSampleIDs) > 0):
        print('Warning! These sample ids are Missing in the immunizer file' + str((recallSampleIDs - immunizerSampleIDs)))

    homozygosityLookup = loadHomozygosityFile(inputFile=args.homozygosity)

    # Filter the Immunizers, select beads based on MFI
    filtered_immunizers_drb1_1 = filterImmunizers(immunizers_drb1=immunizers_drb1_1, immunizersBeadLookup=immunizersBeadLookup, filterMethod=args.filter)
    filtered_immunizers_drb1_2 = filterImmunizers(immunizers_drb1=immunizers_drb1_2, immunizersBeadLookup=immunizersBeadLookup, filterMethod=args.filter)

    # output a simple list of peptides
    writePeptideList(outputFileName=join(args.output, '4A.peptideList.csv'),recallEpitopes_drb1_1=recallEpitopes_drb1_1,recallEpitopes_drb1_2=recallEpitopes_drb1_2,immunizers_drb1_1=filtered_immunizers_drb1_1,immunizers_drb1_2=filtered_immunizers_drb1_2)
    drb1_1_sharedPeptides, drb1_2_sharedPeptides = writeOverlappingPeptides(outputFileName=join(args.output, '4B.overlapping_pirche_peptides.csv'), immunizerSampleIDs=immunizerSampleIDs, immunizers_drb1_1=filtered_immunizers_drb1_1, immunizers_drb1_2=filtered_immunizers_drb1_2, recallEpitopes_drb1_1=recallEpitopes_drb1_1, recallEpitopes_drb1_2=recallEpitopes_drb1_2, immunizersBeadLookup=immunizersBeadLookup)
    drb1_1_combinedPeptides, drb1_2_combinedPeptides = writeCombinedOverlappingPeptides(outputFileName=join(args.output, '4C.combined_overlapping_pirche_peptides.csv'),immunizerSampleIDs=immunizerSampleIDs,drb1_1_sharedPeptides=drb1_1_sharedPeptides, drb1_2_sharedPeptides=drb1_2_sharedPeptides)
    writeScores(outputFileName=join(args.output, '4D.scores_per_patient.csv'),immunizerSampleIDs=immunizerSampleIDs,drb1_1_combinedPeptides=drb1_1_combinedPeptides, drb1_2_combinedPeptides=drb1_2_combinedPeptides)
    writeFinalScore(outputFileName=join(args.output, '4E.final_score_per_patient.csv'), immunizerSampleIDs=immunizerSampleIDs, drb1_1_combinedPeptides=drb1_1_combinedPeptides, drb1_2_combinedPeptides=drb1_2_combinedPeptides, homozygosityLookup=homozygosityLookup)

    print('Done, output written to ' + args.output)



