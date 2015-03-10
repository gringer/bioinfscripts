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
    printFile('room_label/header.html', myparams, True)
    print("<h1>Error</h1>")
    print("<p>%s</p>" % errorMessage)
    printHiddenValues(form, myparams)
    printFile('room_label/footer.html', myparams, False)
    exit(0)

def runGenerator(lastForm, parameters):
    roomID = parameters['ID']
    if(len(roomID) == 0):
        roomID = "X.XX"
    outputFileName = 'room_label/tmp/room_label_%s.pdf' % roomID
    resultSVG = open('room_label/tmp/room_label_%s.svg' % roomID, mode='w+b')
    shutil.copy('room_label/bg_default.png',
                'room_label/tmp/bg_%s.png' % roomID)
    for line in open('room_label/MIMR_Room_template.svg', 'r'):
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
    print('Content-Disposition: attachment; ' +
          'filename=room_label_%s.pdf\n' % roomID)
    sys.stdout.write(open(jamFileName,'rb').read())
    exit(0)

def getResults(parameters):
    resultStorageName = 'room_label/results/resultsFiles.csv'
    mostRecentFileName = None
    mostRecentTime = 0
    mostRecentErrorFileName = None
    mostRecentErrorTime = 0
    formattedPreResult = ''
    formattedSummaryResult = ''
    formattedFullResult = ''
    # sort out GBrowse pattern replacements
    gbrowsePatterns = dict()
    if('gbrowse_patterns' in parameters):
        patternList = parameters['gbrowse_patterns'].split(';')
        for pattern in patternList:
            components = pattern.split(':',1)
            gbrowsePatterns[components[0]] = components[1]
    # find location of results files and error output
    reader = csv.reader(open(resultStorageName))
    for row in reader:
        sessionID = row[0]
        resultFileName = row[1]
        creationTime = float(row[2])
        if((sessionID == parameters['sessionID']) and
           (creationTime >= mostRecentTime)):
            mostRecentFileName = resultFileName
            mostRecentTime = creationTime
        if((sessionID == parameters['sessionID']+".err") and
           (creationTime >= mostRecentErrorTime)):
            mostRecentErrorFileName = resultFileName
            mostRecentErrorTime = creationTime
    resultsFound = False
    if(mostRecentFileName != None):
        if((mostRecentErrorFileName != None) and (os.path.getsize(mostRecentErrorFileName) > 0)):
            f = open(mostRecentErrorFileName, 'r')
            errorStr = "<pre>"
            for line in f:
                errorStr += line
            errorStr += "</pre>"
            writeError(errorStr, parameters)
            return
        if(os.path.getsize(mostRecentFileName) == 0):
            return('The results file is empty. Have a sip of your favourite beverage then click "Results" again.')
        resultFile = open(mostRecentFileName, 'r')
        blast_records = NCBIXML.parse(resultFile)
        formattedPreResult += ('<p>BLAST Run started: %s</p>\n'
                               % time.strftime('%Y-%b-%d %H:%M:%S',
                                               time.localtime(mostRecentTime)))
        formattedFullResult += '<h2>Match Details</h2>'
        formattedFullResult += '<pre>\n'
        numAlignments = 0
        queries = set()
        summaryTable = list()
        for blast_record in blast_records:
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    resultsFound = True
                    query = blast_record.query
                    subject = alignment.hit_def
                    if(" " in query):
                        query = query[0:query.find(" ")]
                    if(" " in subject):
                        subject = subject[0:subject.find(" ")]
                    identity = float(hsp.identities) / (hsp.align_length) * 100
                    coverage = float(abs(hsp.query_end - hsp.query_start)) / blast_record.query_length * 100
                    subjCoverage = float(abs(hsp.sbjct_end - hsp.sbjct_start)) / alignment.length * 100
                    # place appropriate hyperlinks into subject names
                    for subPattern in gbrowsePatterns:
                        if(subPattern in subject):
                            subject = (
                                gbrowsePatterns[subPattern] %
                                (subject,
                                 min(hsp.sbjct_start, hsp.sbjct_end),
                                 max(hsp.sbjct_start, hsp.sbjct_end),
                                 subject,
                                 min(hsp.sbjct_start, hsp.sbjct_end),
                                 max(hsp.sbjct_start, hsp.sbjct_end),
                                 subject))
                    alignmentText = '<a name="%d" href="#summary">**** Alignment %d ****</a>\n' % (
                        numAlignments, numAlignments)
                    alignmentText += 'query: %s\n' % query
#                    alignmentText += str(blast_record.__dict__)
                    alignmentText += 'query length: %s\n' % blast_record.query_length
                    alignmentText += 'subject: %s\n' % subject
                    alignmentText += 'subject length: %s\n' % alignment.length
                    alignmentText += 'align length: %s\n' % hsp.align_length
                    alignmentText += 'score: %s\n' % hsp.score
                    alignmentText += 'bits: %s\n' % hsp.bits
                    alignmentText += 'identity: %0.2f%%\n' % identity
                    alignmentText += 'query coverage: %0.2f%%\n' % coverage
                    alignmentText += 'subject coverage: %0.2f%%\n' % subjCoverage
                    alignmentText += 'e value: %g\n' % hsp.expect
                    querySpos = oldQPos = queryPos = hsp.query_start
                    sbjctSpos = oldSPos = sbjctPos = hsp.sbjct_start
                    alignSpos = alignPos = 0
                    queryDir =  1 if (hsp.query_start < hsp.query_end) else -1
                    sbjctDir = 1 if (hsp.sbjct_start < hsp.sbjct_end) else -1
                    incQuery = False
                    incSbjct = False
                    for hsp.char in hsp.match:
                        if(hsp.query[alignPos] != '-'):
                            oldQPos = queryPos
                            if(not incQuery):
                                querySpos = queryPos
                            incQuery = True
                            queryPos += queryDir
                        if(hsp.sbjct[alignPos] != '-'):
                            oldSPos = sbjctPos
                            if(not incSbjct):
                                sbjctSpos = sbjctPos
                            incSbjct = True
                            sbjctPos += sbjctDir
                        alignPos += 1
                        if((alignPos % 100 == 0) or (alignPos >= len(hsp.match))):
                            alignmentText += '\n'
#                           alignmentText += '      %8s %s\n' % ('', ''.join('         *' * 10))
                            alignmentText += 'Query %8d %s %-8d\n' % (
                                querySpos, hsp.query[alignSpos:alignPos], oldQPos)
                            alignmentText += '      %8s %s\n' % ('', hsp.match[alignSpos:alignPos])
                            alignmentText += 'Sbjct %8s %s %-8d\n' % (
                                sbjctSpos, hsp.sbjct[alignSpos:alignPos], oldSPos)
                            incQuery = False
                            incSbjct = False
                            alignSpos = alignPos
                    # alignmentText += '\n'
                    # alignmentText += 'Query %8d %s %-8d\n' % (hsp.query_start, hsp.query, hsp.query_end)
                    # alignmentText += '      %8s %s\n' % ('', hsp.match)
                    # alignmentText += 'Sbjct %8s %s %-8d\n' % (hsp.sbjct_start, hsp.sbjct, hsp.sbjct_end)
                    formattedFullResult += alignmentText + '\n'
                    queries.add(blast_record.query)
                    summaryTable.append((query, subject, hsp.score,
                        coverage, identity, hsp.expect))
                    numAlignments += 1
        formattedFullResult += '</pre>\n'
        formattedPreResult += ('<p>Number of alignments: %d</p>\n'
                               % numAlignments)
        formattedSummaryResult += '<h2><a name="summary"></a>Summary</h2>\n'
        formattedSummaryResult += '<table class="sortable">\n'
        formattedSummaryResult += ('<thead>\n' +
                                   '  <tr>' +
                                   '<th>Alignment</th>' +
                                   '<th>Query</th>' +
                                   '<th>Subject</th>' +
                                   '<th>Bitscore</th>' +
                                   '<th>Coverage %</th>' +
                                   '<th>Identity %</th>' +
                                   '<th>E value</th>' +
                                   '</tr>\n' +
                                   '</thead>\n')
        formattedSummaryResult += '<tbody>\n'
        for alignment in range(numAlignments):
            formattedSummaryResult += (('  <tr><td><a href="#%d">%d</a></td><td>%s</td><td>%s</td>' +
                                        '<td>%0.2f</td><td>%0.2f</td><td>%0.2f</td><td>%5g</td></tr>\n') % (
                    alignment,
                    alignment,
                    (summaryTable[alignment])[0],
                    (summaryTable[alignment])[1],
                    (summaryTable[alignment])[2],
                    (summaryTable[alignment])[3],
                    (summaryTable[alignment])[4],
                    (summaryTable[alignment])[5]))
        formattedSummaryResult += '</tbody>\n'
        formattedSummaryResult += '</table>\n'
    if(resultsFound):
        return(formattedPreResult + formattedSummaryResult + formattedFullResult)
    else:
        return('No hits were found')

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
if(form.getfirst("runProgram","") == "TRUE"):
    runGenerator(form, myparams)

if((not('resultsExist' in myparams)) or (myparams['resultsExist'] != 'True')):
    myparams['class_results'] += " tabdisabled"
else:
    # retrieve result file, and display on tab (if tab is visible)
    if(currentTab == 'results'):
        myparams['results'] = getResults(myparams)

printFile('room_label/labeler.html', myparams, True)
if('errors' in myparams):
    print('<h3>Errors:</h3><pre>%s</pre>' % myparams['errors'])
printHiddenValues(form, myparams)
printFile('room_label/footer.html', myparams, False)
