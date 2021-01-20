import argparse
import sys
import re
from os.path import isfile, join
from os import listdir

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
    parser.add_argument("-d", "--donor", help="donor column")
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
                        donorTyping = hlaListFromString(cols[int(args.donor)])


                        patientTypingWithoutBroad = removeBroads(splitToBroad, patientTyping)
                        donorTypingWithoutBroad = removeBroads(splitToBroad, donorTyping)

                        patientTypingWithoutMolecular = removeMolecular(patientTypingWithoutBroad)
                        donorTypingWithoutMolecular = removeMolecular(donorTypingWithoutBroad)

                        # Only extrapolating drb1 to high res for the patients, not the donors.
                        patientTypingExtrapolatedDrb1 = extrapolateDrb1Typings(id, patientTypingWithoutMolecular, extrapolatedDrb1Typings)
                        print('Patient ' + str(id) + ' with extrapolated typings:' + str(patientTypingExtrapolatedDrb1))

                        if verbose:
                            print("Patient typing: " + str(patientTyping) + "; without broads: " + str(patientTypingExtrapolatedDrb1))
                            print("Donor typing: " + str(donorTyping) + "; without broads: " + str(donorTypingWithoutMolecular))

                        patientTypingPirche = pircheAgString(patientTypingExtrapolatedDrb1, ",")
                        donorTypingPirche = pircheAgString(donorTypingWithoutMolecular, ",")

                        if donorTypingPirche != "" and patientTypingPirche != "":
                            combinedTypings = "patient" + id + "," + patientTypingPirche + "\n" + "donor" + id + "," + donorTypingPirche
                            outputPirche.write(combinedTypings + "\n,\n")
                            outputPirche.flush()

                        output.write(row.strip() + ";\"" + "patient" + id + ";" + patientTypingPirche + "\";\"" + "donor" + id + ";" + donorTypingPirche + "\";")
                    else:
                        output.write(row.strip())
                    output.write("\n");
                    output.flush()
