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
import math # for label placement

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
            parameters['seenFields'].add(field)
        print(line.rstrip())

def printHiddenValues(lastForm, parameters):
    # make sure runProgram state isn't preserved across multiple submits
    # [don't want it to try running more than once]
    parameters['seenFields'].add('runProgram')
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
    inputFile = '../room_label/MIMR_Room_template.svg'
    hbase = 'unece/'
    if('ANGELS' in parameters):
        parameters['ANGELS'] = parameters['ANGELS'].replace('\n','</flowPara><flowPara>')
    if (('selectClass' in parameters) and (parameters['selectClass'] == 'lab')):
        inputFile = '../room_label/MIMR_Lab_template.svg'
    for line in open(inputFile, 'r'):
        if('bg_default.png' in line):
            if(('inputFile' in parameters) and
               (lastForm['inputFile'].filename)):
                fileItem = lastForm['inputFile']
                fName, fExt = os.path.splitext(fileItem.filename)
                fHead, fTail = os.path.split(fileItem.filename)
                open('room_label/tmp/bg_%s%s' % (roomID, fExt),
                     'wb').write(fileItem.file.read())
                if(fExt == ".pdf"):
                    commandLine = list(('convert',
                                        '-density','150',
                                        'room_label/tmp/bg_%s%s' %
                                        (roomID, fExt),
                                        'room_label/tmp/bg_%s.jpg' %
                                        roomID))
                    runProcess = subprocess.call(commandLine)
                    fExt = ".jpg"
                if(docType == "svg"):
                    fHead, fTail = os.path.split(fileItem.filename)
                    line = line.replace('bg_default.png', fTail)
                else:
                    line = line.replace('bg_default.png',
                                        'bg_%s%s' % (roomID,fExt))
            else:
                line = '<!--' + line.strip() + '-->\n'
        if('@' in line):
            paramMatches = re.findall("@(.*?)@", line)
            for param in paramMatches:
                if(param in parameters):
                    line = line.replace('@' + param + '@', parameters[param])
                else:
                    line = line.replace('@' + param + '@', '['+param+']')
        if('${HAZARDS}' in line):
            hazardStr = "    "
            if('warnBoxes' in lastForm):
                boxHeight = 200
                boxWidth = 440
                boxSY = 823
                hazardItems = list(lastForm.getvalue('warnBoxes'))
                numHazards = len(hazardItems)
                hLines = 2 if (numHazards > 3) else 1
                hWidth = boxWidth / (int(math.ceil(numHazards / hLines)))
                if(hWidth > (boxHeight/hLines)):
                    hWidth = (boxHeight/hLines)
                hHeight = hWidth
                nextHazard = 0
                lPos = boxWidth - (hWidth * numHazards) / (2 * hLines)
                for item in hazardItems:
                    xPos = (int(nextHazard / hLines) * hWidth +
                            hWidth*0.05 + lPos + 7.5 +
                            (numHazards % 2) * (nextHazard % hLines) * (hWidth / 2))
                    yPos = int(nextHazard % hLines) * hHeight + hWidth*0.05 + boxSY + 7.5
                    hazardStr = ((hazardStr + '<image xlink:href="%s%s.svg" ' % (hbase,item)) +
                                 'height="%f" width="%f" x="%f" y="%f" />' %
                                 (hWidth * 0.9, hHeight * 0.9, xPos, yPos))
                    nextHazard += 1
            line = hazardStr + '\n'
        line = line.replace('&','&amp;')
        resultSVG.write(line)
    resultSVG.close()
    if(docType == "pdf"):
        exportLine = '--export-pdf=%s' % outputFileName
        commandLine = list(('inkscape',
                            exportLine,
                            resultSVG.name))
        runProcess = subprocess.call(commandLine)
        jamFileName = outputFileName.replace(".pdf","-pdfjam.pdf")
        commandLine = list(('pdfjam','--landscape',
                            '--a4paper',
#                            '--a5paper',
#                            '--preamble',
#                            '\usepackage[cam,a4,center]{crop}'
                            '--scale', '0.71',
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
    "seenFields"   : set(),
    "addFields"    : ('resultsExist','sessionID','blastCommand','selectClass'),
    }

#blastn -db db/3alln_smed -query /tmp/tmpGFk3hs -outfmt 5 -task blastn -evalue 10 -max_target_seqs 100 -word_size 11

currentType = form.getfirst("selectClass","office")
currentTab = form.getfirst("selectTab","query")

# overwrite default values with previous form values
loadForm(form, myparams)

# add sessionID if it doesn't already exist
if(not('sessionID' in myparams)):
    myparams['sessionID'] = base64.b64encode(os.urandom(16))

myparams['type'] = currentType
# activate current tab
myparams['class_' + currentTab] = "tabon"

if('REQUEST_URI' in os.environ):
    myparams['request_uri'] = os.environ['REQUEST_URI'].split("?")[0]
else:
    myparams['request_uri'] = 'shell'
    for arg in sys.argv:
        if('=' in arg):
            (opt,val) = arg.split('=',2)
            myparams[opt] = val

if(not ('ID' in myparams)):
    myparams['ID'] = "X.XX"

# run generator program (if requested)
docType = form.getfirst("runProgram","")
if(len(docType) > 0):
    runGenerator(form, myparams, docType)

printFile('../room_label/%s-labeler.html' % currentType, myparams, True)
if('errors' in myparams):
    print('<h3>Errors:</h3><pre>%s</pre>' % myparams['errors'])
printHiddenValues(form, myparams)
printFile('../room_label/footer.html', myparams, False)
