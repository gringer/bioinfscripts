#!/usr/bin/python

# room_label.py -- generate room labels

import cgi # for cgi forms
import cgitb # for cgi trace-back
import re # for regular expression parsing
import os # file existence, urandom [sessionIDs]
import base64 # for encoding sessionIDs
import subprocess # for running external programs
import csv # for parsing csv files (e.g. blast output)
import time # for results file cleanup
import tempfile # for blast results
import sys # for command-line argument parsing
import shutil # for copying files

def printFile(fileName, parameters, printContent):
    if(printContent):
        # HTTP header
        print('Content-type: text/html\n')
    if(not(os.path.exists(fileName))):
        print('File does not exist: %s' % fileName);
        return
    readFile = open(fileName, 'r')
    for line in readFile:
        paramMatches = re.findall("%\((.*?)\)", line)
        for param in paramMatches:
            # doing a 'manual' replacement because the usual %(param)
            # notation doesn't seem to work
            if(param in parameters):
                line = line.replace('%(' + param + ')',
                                    parameters[param])
            else:
                line = line.replace('%(' + param + ')', '')
        seenFields = re.findall('<(?:input|textarea|select).*?name="(.*?)"', line)
        for field in seenFields:
            parameters['seenFields'].append(field)
        print(line.rstrip())

def printHiddenValues(lastForm, parameters):
    # make sure runProgram state isn't preserved across multiple submits
    # [don't want it to try running more than once]
    parameters['seenFields'].append('runProgram')
    # add fields from 'addFields' to hidden values
    for field in parameters['addFields']:
        if(not(field in parameters['seenFields']) and
           field in parameters):
            print('<input type="hidden" name="%s" value="%s">' %
                  (field, parameters[field]))
    # add fields from previous form to hidden values
    for field in lastForm:
        if(not(field in parameters['seenFields'])):
            print('<input type="hidden" name="%s" value="%s">' %
                  (field, lastForm.getfirst(field)))

def loadForm(lastForm, parameters):
    # retrieves data values from the previous form
    for field in lastForm:
        parameters[field] = lastForm.getfirst(field, '')

def loadDefaults(fileName, parameters):
    # loads default form values from a file
    reader = csv.reader(open(fileName))
    for row in reader:
        if(len(row) == 2):
            parameters[row[0]] = row[1]

def writeError(errorMessage, parameters):
    printFile('../room_label/header.html', myparams, True)
    print("<h1>Error</h1>")
    print("<p>%s</p>" % errorMessage)
    printHiddenValues(form, myparams)
    printFile('../room_label/footer.html', myparams, False)
    exit(0)

def runGenerator(lastForm, parameters, docType):
    roomID = parameters['ID']
    if(len(roomID) == 0):
        roomID = "X.XX"
    outputFileName = 'room_label/tmp/room_label_%s.pdf' % roomID
    svgFileName = 'room_label/tmp/room_label_%s.svg' % roomID
    resultSVG = open(svgFileName, mode='w+b')
    shutil.copy('../room_label/bg_default.png',
                'room_label/tmp/bg_%s.png' % roomID)
    for line in open('../room_label/MIMR_Room_template.svg', 'r'):
        if('bg_default.png' in line):
            if(('inputFile' in parameters) and
               (lastForm['inputFile'].filename)):
                fileItem = lastForm['inputFile']
                fName, fExt = os.path.splitext(fileItem.filename)
                open('room_label/tmp/bg_%s.%s' % (roomID, fExt),
                     'wb').write(fileItem.file.read())
                line = line.replace('bg_default.png',
                                    'bg_%s.%s' % (roomID,fExt))
            else:
                line = '<!--' + line.strip() + '-->\n'
        if('@' in line):
            paramMatches = re.findall("@(.*?)@", line)
            for param in paramMatches:
                if(param in parameters):
                    line = line.replace('@' + param + '@', parameters[param])
                else:
                    line = line.replace('@' + param + '@', '['+param+']')
        resultSVG.write(line)
    resultSVG.close()
    if(docType == "pdf"):
        exportLine = '--export-pdf=%s' % outputFileName
        commandLine = list(('inkscape',
                            exportLine,
                            resultSVG.name))
        runProcess = subprocess.call(commandLine)
        jamFileName = outputFileName.replace(".pdf","-pdfjam.pdf")
        commandLine = list(('pdfjam','--landscape','--a4paper',
                            '--scale','0.71',
                            '--outfile', jamFileName,
                            outputFileName))
        runProcess = subprocess.call(commandLine)
        print('Content-type: application/pdf')
        print('Content-Disposition: attachment; ' +
              'filename=room_label_%s.pdf\n' % roomID)
        sys.stdout.write(open(jamFileName,'rb').read())
    else:
        print('Content-type: image/svg+xml')
        print('Content-Disposition: attachment; ' +
              'filename=room_label_%s.svg\n' % roomID)
        sys.stdout.write(open(svgFileName,'rb').read())
    exit(0)

### Begin Actual Program ###

cgitb.enable() # make errors visible on web pages

form = cgi.FieldStorage()   # FieldStorage object to
                            # hold the form data
myparams = {
    "class_query"  : "taboff",
    "class_params" : "taboff",
    "class_results": "taboff",
    "seenFields"   : list(),
    "addFields"    : ('resultsExist','sessionID','blastCommand'),
    }

#blastn -db db/3alln_smed -query /tmp/tmpGFk3hs -outfmt 5 -task blastn -evalue 10 -max_target_seqs 100 -word_size 11

currentProgram = form.getfirst("selectProgram","blastn")
currentTab = form.getfirst("selectTab","query")

# overwrite default values with previous form values
loadForm(form, myparams)

# add sessionID if it doesn't already exist
if(not('sessionID' in myparams)):
    myparams['sessionID'] = base64.b64encode(os.urandom(16))

myparams['program'] = currentProgram
# activate current tab
myparams['class_' + currentTab] = "tabon"

if('REQUEST_URI' in os.environ):
    myparams['request_uri'] = os.environ['REQUEST_URI']
else:
    myparams['request_uri'] = 'shell'
    for arg in sys.argv:
        if('=' in arg):
            (opt,val) = arg.split('=',2)
            myparams[opt] = val

if(not ('ID' in myparams)):
    myparams['ID'] = "X.XX"

# run BLAST (if requested)
docType = form.getfirst("runProgram","")
if(len(docType) > 0):
    runGenerator(form, myparams, docType)

if((not('resultsExist' in myparams)) or (myparams['resultsExist'] != 'True')):
    myparams['class_results'] += " tabdisabled"
else:
    # retrieve result file, and display on tab (if tab is visible)
    if(currentTab == 'results'):
        myparams['results'] = getResults(myparams)

printFile('../room_label/labeler.html', myparams, True)
if('errors' in myparams):
    print('<h3>Errors:</h3><pre>%s</pre>' % myparams['errors'])
printHiddenValues(form, myparams)
printFile('../room_label/footer.html', myparams, False)
