import argparse
import sys
import re
from os.path import isfile, join
from os import listdir


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--verbose", help="verbose operation", action="store_true")
    parser.add_argument("-i", "--input", help="typings input file", required=True)
    parser.add_argument("-o", "--output", help="output file for extrapolated typing", required=True)
    parser.add_argument("-4", "--reject", help="rejected sample list output file", required=True)
    parser.add_argument("-c", "--cutoff", help="cutoff value for confidence of typing, default .65", default='.65')


    args = parser.parse_args()

    # newline and separator characters, change if necessary.
    newline='\r\n'
    separator=';'

    verbose = args.verbose
    if verbose:
        print("verbose mode active")

    patientIDColumn = 0
    drb1_1Column = 7
    drb1_2Column = 8
    frequencyColumn = 19
    cutoffFrequency = float(args.cutoff)

    potentialTypings={}

    print('We are extrapolating typing data now.')
    print('opening input file:' + str(args.input))
    with open(args.input, "r") as input:
        for row in input:
            #print('Row:' + str(row))
            columns = row.split(separator)
            if('#' in columns[patientIDColumn]):
                patientID, otheridentifier = columns[patientIDColumn].split('#')
                drTypings = separator.join(sorted([columns[drb1_1Column],columns[drb1_2Column]]))
                # In the input file, a value of 10,000 indicates a single value. Not sure why it is 10,000 but that means there is only one option.
                if(str(columns[frequencyColumn]) == '10,000'):
                    frequency = 1
                else:
                    frequency = float(columns[frequencyColumn])

                # print('This is an interesting data row:' + str(row))
                # print('ID:' + str(columns[patientIDColumn]))
                # print('drb1_1:' + str(columns[drb1_1Column]))
                # print('drb1_2:' + str(columns[drb1_2Column]))
                # print('Frequency:' + str(columns[frequencyColumn]))
                # print('sorted typing:' + drTypings)

                if(patientID not in potentialTypings.keys()):
                    potentialTypings[patientID] = {}

                if(drTypings not in potentialTypings[patientID].keys()):
                    potentialTypings[patientID][drTypings] = []

                # Store this freq
                potentialTypings[patientID][drTypings].append(frequency)

    extrapolatedTypings = []
    rejectedSampleList = []

    # Loop through potential typings to get the extrapolated typing
    for patientID in sorted(potentialTypings.keys()):
        totalPatientFrequencies = 0
        extrapolatedTypingFound = False
        for drTyping in sorted(potentialTypings[patientID].keys()):
            totalTypingFrequencies = 0
            #print('checking patient ' + str(patientID) + ' and typings ' + str(drTyping))
            for frequency in potentialTypings[patientID][drTyping]:
                totalPatientFrequencies += float(frequency)
                totalTypingFrequencies += float(frequency)

            if(totalTypingFrequencies >= cutoffFrequency):
                #print('Patient ' + str(patientID) + ' and typings ' + str(drTyping) + ' frequency ' + str(totalTypingFrequencies) + ' exceeded the cutoff value(' + str(cutoffFrequency) + ')')
                extrapolatedTypings.append((patientID,drTyping,str(totalTypingFrequencies)))
                extrapolatedTypingFound = True

        #print('Total patient frequencies:' + str(totalPatientFrequencies))
        # sanity check
        if(totalPatientFrequencies > 1):
            print('Warning! Total frequencies patient ' + str(patientID) + ' was greater than 1:' + str(totalPatientFrequencies))
        # record rejected sample IDs
        if(not extrapolatedTypingFound):
            rejectedSampleList.append(patientID)


    # write output
    print('opening output file')
    with open(args.output, "w") as output:
        output.write('Patient_ID' + separator + 'DRB1_1' + separator + 'DRB1_2' + separator + 'Sum_Frequency(>=' + str(cutoffFrequency) + ')' + newline)
        for extrapolatedTyping in extrapolatedTypings:
            output.write(separator.join(extrapolatedTyping) + newline)

    with open(args.reject, "w") as output:
        output.write('Patient_ID(No frequency >=' + str(cutoffFrequency) + ')' + newline)
        for rejectedSample in rejectedSampleList:
            output.write(rejectedSample + newline)

    print('Done.')
