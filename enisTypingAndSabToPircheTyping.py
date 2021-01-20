import argparse
#import sys
import re
#from os.path import isfile, join
#from os import listdir

newline = '\r\n'

def hlaListFromString(input):
    ags = input.split(" ")
    asList = []
    for ag in ags:
        asList.append(ag.strip())
    return asList

def removeBroads(bsMap, input):
    output = []
    for ag in input:
        output.append(ag)
    for ag in input:
        if ag in bsMap.keys():
            if bsMap[ag] in output:
                output.remove(bsMap[ag])
    return output

def readMFICsv(inputCsv=None,sampleIDCol=0,alleleCol=2):
    sampleAlleles={}
    with open(inputCsv, "r") as input:
        rowCount = 0
        for row in input:
            rowCount += 1
            if rowCount > 1:
                cols = row.split(";")
                sampleID = cols[sampleIDCol]
                alleles = cols[alleleCol]
                #alleleList = alleles.split(',')
                alleleList = [alleles]
                # heterodimers (DAB1/DQA1) are split by a comma. Add both. Not sure if this is correct.
                for allele in alleleList:
                    if(sampleID in sampleAlleles.keys()):
                        sampleAlleles[sampleID].append(allele)
                    else:
                        sampleAlleles[sampleID] = [allele]
    return sampleAlleles

def pircheAgString(input, separator):
    if(input is None):
        return ''

    output = ""
    for ag in input:
        match = re.match(r"([A-Za-z]+)([0-9]+)", ag, re.I)
        if match:
            items = match.groups()
            locus = items[0]
            allele = items[1]
            print('found locus ' + str(locus))
            if(locus == 'DRB'):
                output = output + ag + separator
            elif not locus == "Bw" and not locus == "BW" and not (locus == "DR" and (allele == "51" or allele == "52" or allele == "53")):
                if locus == "Cw":
                    locus = "C"
                if locus == "CW":
                    locus = "C"
                if locus == "DQ":
                    locus = "DQB1"
                if locus == "DR":
                    locus = "DRB1"
                output = output + locus + "*" + allele + separator

    print('returning output ' + str(output))

    return output

def pirchAgStringWithSABAlleles(input, separator, patient):
    if(verbose):
        print('Getting Pirch String from a list of length' + str(len(input)))
    output = ""
    for alleleIndex, alleleString in enumerate(sorted(list(set(input)))):
        alleleString=alleleString.replace('\r','').replace('\n','')
        output = output + patient + '-Bead' + str(alleleIndex+1) + separator + alleleString + newline

    return output

def removeMolecular(input):
    output = []
    for ag in input:
        if not "*" in ag and not ":" in ag and not ag == "DQB1":
            output.append(ag)
    return output


def extrapolateDrb1Typings(id, patientTypings, extrapolatedDrb1Typings):
    #print('patient ' + str(id) + ' typing before:' + str(patientTypings))
    typingsWithExtrapolatedDrb1 = []

    # If there is an extrapolated typing for this patient
    if(id in extrapolatedDrb1Typings.keys()):
        #print('patient ' + str(id) + ' is in the list of patients with extrapolated typings.')
        for ag in patientTypings:
            #print('ag:' + str(ag))
            if not str(ag).startswith('DR'):
                typingsWithExtrapolatedDrb1.append(ag)
        extrapolatedDr1, extrapolatedDr2 = extrapolatedDrb1Typings[id].split(';')
        typingsWithExtrapolatedDrb1.append(extrapolatedDr1)
        typingsWithExtrapolatedDrb1.append(extrapolatedDr2)

        #print('patient ' + str(id) + ' typing after:' + str(typingsWithExtrapolatedDrb1))
        return typingsWithExtrapolatedDrb1

    # otherwise this patient probably didn't reach the cutoff for confidence in the extrapolated typing
    else:
        #print('patient ' + str(id) + ' is NOT in the list of patients with extrapolated typings.')
        return None


    #return typingsWithExtrapolatedDrb1
    return patientTypings


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", help="verbose operation", action="store_true")
    parser.add_argument("-i", "--input", help="input file")
    parser.add_argument("-b", "--broad", help="broad-split file")
    parser.add_argument("-d", "--donorcsv", help="donor csv bead file")
    parser.add_argument("-p", "--patient", help="patient column")
    parser.add_argument("-I", "--id", help="id column")
    parser.add_argument("-e", "--extrapolated_drb1", help="extrapolated drb1 typings input file")
    parser.add_argument("-o", "--output", help="output file")
    parser.add_argument("-O", "--outputPirche", help="output file for PIRCHE format")

    args = parser.parse_args()

    verbose = args.verbose
    if verbose:
        print("verbose mode active")

    if (verbose):
        print("Input " + args.input)
        print("Broad-Split-File " + args.broad)
        print("Output " + args.output)

    extrapolatedDrb1Typings={}
    with open(args.extrapolated_drb1, "r") as extFile:
        rowCount = 0
        for row in extFile:
            rowCount += 1
            if rowCount > 1:
                sampleId, drb1_1, drb1_2, frequency = row.split(";")
                sampleId = sampleId.replace('patient','')
                #print('read sample:' + str(sampleId))
                extrapolatedDrb1Typings[sampleId] = (drb1_1 + ';' + drb1_2)

    splitToBroad = {}
    with open(args.broad, "r") as bsfile:
        for row in bsfile:
            if not row.startswith("#"):
                cols = row.split(";")
                locus = cols[0]
                broad = cols[1]
                splits = cols[2].split("/")

                for split in splits:
                    splitToBroad[locus + split] = locus + broad

    if (verbose):
        print("Split to broad: " + str(splitToBroad))

    donorTypings = readMFICsv(inputCsv=args.donorcsv)

    with open(args.outputPirche, "w") as outputPirche:
        with open(args.output, "w") as output:

            with open(args.input, "r", encoding="utf-8", errors="ignore") as input:
                rowCount = 0
                for row in input:
                    rowCount += 1
                    if rowCount > 1:
                        cols = row.split(";")
                        id = cols[int(args.id)]
                        patientTyping = hlaListFromString(cols[int(args.patient)])
                        #donorTyping = hlaListFromString(cols[int(args.donor)])
                        donorTyping = donorTypings[id]

                        patientTypingWithoutBroad = removeBroads(splitToBroad, patientTyping)
                        donorTypingWithoutBroad = removeBroads(splitToBroad, donorTyping)

                        patientTypingWithoutMolecular = removeMolecular(patientTypingWithoutBroad)
                        patientTypingExtrapolatedDrb1 = extrapolateDrb1Typings(id, patientTypingWithoutMolecular, extrapolatedDrb1Typings)
                        #donorTypingWithoutMolecular = removeMolecular(donorTypingWithoutBroad)
                        donorTypingWithoutMolecular=donorTypingWithoutBroad
                        if verbose:
                            # TODO: I think there's an issue here, double check how the broad and split being removed.
                            print("Patient typing: " + str(patientTyping) + newline + " without broads: " + str(patientTypingExtrapolatedDrb1))
                            print("Donor typing: " + str(donorTyping) + newline + " without broads: " + str(donorTypingWithoutMolecular))

                        patientTypingPirche = pircheAgString(patientTypingExtrapolatedDrb1, ",")
                        #donorTypingPirche = pircheAgString(donorTypingWithoutMolecular, ";")
                        donorTypingPirche = pirchAgStringWithSABAlleles(donorTypingWithoutMolecular, ",", "patient" + id)

                        if verbose:
                            print("Patient Typing String: " + str(patientTypingPirche))
                            print("Donor Typing String: " + str(donorTypingPirche))

                        if donorTypingPirche != "" and patientTypingPirche != "":
                            combinedTypings = ("patient" + id + "," + patientTypingPirche
                                + newline + donorTypingPirche + "," + newline)
                            outputPirche.write(combinedTypings)
                            outputPirche.flush()

                        output.write(row.strip() + ";\"" + "patient" + id + ";" + patientTypingPirche + "\";\"" + donorTypingPirche + "\";")
                    else:
                        output.write(row.strip())
                    output.write(newline);
                    output.flush()