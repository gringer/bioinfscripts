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
        print('File does not exist: %s' % fileName)
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
    pc2FileName = '@@fail@@'
    svgFileName = 'room_label/tmp/room_label_%s.svg' % roomID
    resultSVG = open(svgFileName, mode='w+b')
    pc2SVG = None
    inputFile = '@@fail@@'
    hbase = 'unece/'
    doPC2 = False
    if('ANGELS' in parameters):
        parameters['ANGELS'] = parameters['ANGELS'].replace('\n','</flowPara><flowPara>')
    if (not ('type' in parameters)):
        parameters['type'] = 'office'
    # this file name should really be extracted from the HTML file
    if (parameters['type'] == 'lab'):
        inputFile = '../room_label/MIMR_Lab_template.svg'
        if(lastForm.getvalue('pc2cb') == "pc2cb"):
            doPC2 = True
    if (parameters['type'] == 'office'):
        inputFile = '../room_label/MIMR_Room_template.svg'
    if (parameters['type'] == 'office2'):
        inputFile = '../room_label/MIMR_Office2_template.svg'
    if (parameters['type'] == 'plain'):
        inputFile = '../room_label/MIMR_Plain_template.svg'
    if (parameters['type'] == 'BRU'):
        inputFile = '../room_label/MIMR_BRU_template.svg'
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
                    if((param == "ID") and (parameters[param] == "")):
                        line = line.replace('Room ','')
                    line = line.replace('@' + param + '@', parameters[param])
                else:
                    line = line.replace('@' + param + '@', '['+param+']')
        if('${HAZARDS}' in line):
            hazardStr = "    "
            if('warnBoxes' in lastForm):
                boxHeight = 200
                boxWidth = 440
                boxCY = 823 + 215/2
                boxCX = 221 + 460/2
                formHazards = lastForm.getvalue('warnBoxes')
                hazardItems = list(formHazards)
                if(isinstance(formHazards, str)):
                    hazardItems = list([formHazards])
                numHazards = len(hazardItems)
                hLines = 2 if (numHazards > 3) else 1
                hWidth = boxWidth / (int(math.ceil(numHazards / hLines)))
                if(hWidth > (boxHeight/hLines)):
                    hWidth = (boxHeight/hLines)
                hHeight = hWidth
                nextHazard = 0
                lPos = boxCX - (hWidth * numHazards) / (2 * hLines)
                tPos = boxCY - (hHeight * hLines / 2)
                for item in hazardItems:
                    xPos = (int(nextHazard / hLines) * hWidth +
                            hWidth*0.05 + lPos +
                            (numHazards % 2) * (nextHazard % hLines) * (hWidth / 2))
                    yPos = int(nextHazard % hLines) * hHeight + hWidth*0.05 + tPos
                    hazardStr = ((hazardStr + '<image xlink:href="%s%s.svg" ' % (hbase,item)) +
                                 'height="%f" width="%f" x="%f" y="%f" />' %
                                 (hWidth * 0.9, hHeight * 0.9, xPos, yPos))
                    nextHazard += 1
            line = hazardStr + '\n'
        line = line.replace('&','&amp;')
        resultSVG.write(line)
    resultSVG.close()
    if(doPC2):
        pc2FileName = 'room_label/tmp/pc2_label_%s.svg' % roomID
        pc2SVG = open(pc2FileName, mode='w+b')
        for line in open('../room_label/PC2_instructions_template.svg', 'r'):
            if('@' in line):
                paramMatches = re.findall("@(.*?)@", line)
                for param in paramMatches:
                    if(param in parameters):
                        line = line.replace('@' + param + '@',
                                            parameters[param])
                    else:
                        line = line.replace('@' + param + '@',
                                            '['+param+']')
            line = line.replace('&','&amp;')
            pc2SVG.write(line)
        pc2SVG.close()
    if(docType == "pdf"):
        exportLine = '--export-pdf=%s' % outputFileName
        commandLine = list(('inkscape',
                            exportLine,
                            resultSVG.name))
        runProcess = subprocess.call(commandLine)
        if(doPC2):
            exportLine = '--export-pdf=%s' % pc2FileName
            commandLine = list(('inkscape',
                                exportLine,
                                pc2SVG.name))
            runProcess = subprocess.call(commandLine)
        jamFileName = outputFileName.replace(".pdf","-pdfjam.pdf")
        commandLine = list(('pdfjam','--landscape',
                            '--preamble',
                            '\usepackage[cross,axes,a4,center,noinfo]{crop}',
                            '--papersize', '{15.8cm,21.03cm}',
                            '--outfile', jamFileName,
                            outputFileName))
        #sys.stderr.write("Running '%s'...\n" % (" ".join(commandLine)))
        runProcess = subprocess.call(commandLine)
        if(doPC2):
            catFileName = jamFileName.replace("-pdfjam.pdf","-pdfjam-PC2.pdf")
            commandLine = list(('pdftk','A=%s' % jamFileName,
                                'B=%s' % pc2FileName, 'cat',
                                'Awest','B',
                                'output',
                                catFileName))
            runProcess = subprocess.call(commandLine)
            print('Content-type: application/pdf')
            print('Content-Disposition: attachment; ' +
                  'filename=room_label_%s.pdf\n' % roomID)
            sys.stdout.write(open(catFileName,'rb').read())
        else:
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

# overwrite default values with previous form values
loadForm(form, myparams)

currentType = form.getfirst("selectClass","office")

# add sessionID if it doesn't already exist
if(not('sessionID' in myparams)):
    myparams['sessionID'] = base64.b64encode(os.urandom(16))

myparams['type'] = currentType

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
