#!/usr/bin/python

# fastafetch.py -- extract sequences from a fasta file, web UI for samtools faidx

import cgi # for cgi forms
import cgitb # for cgi trace-back
import re # for regular expression parsing
import os # file existence, urandom [sessionIDs]
import base64 # for encoding sessionIDs
import subprocess # for running external programs
import csv # for parsing csv files (e.g. blast output)
import time # for results file cleanup
import tempfile # for blast results
from Bio.Blast import NCBIXML # for XML parsing
from decimal import Decimal # for scientific notation

# set-up variables
fastaDBdir = 'db/fasta/'
runFetch = False

def printHiddenValues(lastForm, parameters):
    # make sure runBlast state isn't preserved across multiple submits
    # [don't want it to try running more than once]
    parameters['seenFields'].append('runBlast')
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

def printOptions(lastForm, parameters):
    print('<p><label accesskey=f>FASTA file: <select name="queryDB">')
    fastaNames = list()
    for filename in sorted(os.listdir(fastaDBdir)):
        if(filename.endswith(".fai")):
            baseName = filename[:-4]
            if(('queryDB' in parameters) and
               (baseName == parameters['queryDB'])):
                fastaNames[0:0] = [baseName]
            else:
                fastaNames.append(baseName)
    for baseName in fastaNames:
        print('<option value="%s">%s</option>' % (baseName, baseName))
    print('</select></p>')
    print('')
    print('<p class="textOption"><span style="color: white">or </span>' +
          '<input type="radio" value="getSeq" name="seqOpt" id="optGS" checked />')
    print('<label>Sequence ID(s): <textarea cols=40 rows=4 id="fastaResult" name="seqID"></textarea>' +
          '</label></p>')
    print('<p>or <input type="radio" value="getList" name="seqOpt" id="optGL" />')
    print('<label for="optGL">List sequences</label></p>')
    print('<button type="submit" class="blastbutton" name="fetch" value="fetchFASTA">' +
          'Fetch</button>')

def loadForm(lastForm, parameters):
    # retrieves data values from the previous form
    for field in lastForm:
        parameters[field] = lastForm.getfirst(field, '')

### Begin Actual Program ###

cgitb.enable() # make errors visible on web pages

form = cgi.FieldStorage()   # FieldStorage object to
                            # hold the form data
myparams = {
    "class_query"  : "taboff",
    "class_params" : "taboff",
    "class_results": "taboff",
    "seenFields"   : list(),
    "addFields"    : ('resultsExist','sessionID' ,'blastCommand'),
    }

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

if(myparams['program'] in ('blastn', 'blastx', 'tblastx')):
    myparams['inputType'] = 'nucleotide';
if(myparams['program'] in ('blastp', 'tblastn')):
    myparams['inputType'] = 'protein';

if('REQUEST_URI' in os.environ):
    myparams['request_uri'] = os.environ['REQUEST_URI']
else:
    myparams['request_uri'] = '[THE FILE YOU RAN]'

# run BLAST (if requested)
if(form.getfirst("fetch","") == "fetchFASTA"):
    runFetch = True

if((not('resultsExist' in myparams)) or (myparams['resultsExist'] != 'True')):
    myparams['class_results'] += " tabdisabled"
else:
    # retrieve result file, and display on tab (if tab is visible)
    if(currentTab == 'results'):
        myparams['results'] = getResults(myparams)

print('Content-type: text/html\n')
print '''<html><head>
<title>FASTA Sequence Fetcher</title>
<style id="pageStyle" type="text/css">
  .textOption *{
    vertical-align: top;
  }
</style>
</head><body>
<h1>FASTA Sequence Fetcher</h1>
<em>[Interface to <tt>SAMtools faidx</tt>]</em>
<form method="post" action="%s" enctype="multipart/form-data">
''' % (myparams['request_uri'])

printOptions(form, myparams)

printHiddenValues(form, myparams)

print('</form>')

if(runFetch):
    queryPath = fastaDBdir + myparams['queryDB']
    if(myparams['seqOpt'] == 'getList'):
        seqs = list()
        with open(queryPath + '.fai') as f:
            print('<h2>Sequence list</h2>')
            print('<p><em>[From "%s"]</em></p>' % myparams['queryDB'])
            print('<textarea cols=80 rows=10 id="fastaResult">')
            for line in f:
                line = line.split()[0]
                seqs.append(line)
        for seq in sorted(seqs):
            print seq
        print('</textarea>')
    else:
        resultFile = tempfile.TemporaryFile()
        errorFile = tempfile.TemporaryFile()
        cmdArgs = ['samtools','faidx', queryPath]
        seqIDs = myparams['seqID'].split()
        for seqID in seqIDs:
            if(not ':' in seqID):
                # prevent too-large sequences from being accidentally included
                cmdArgs.append(seqID + ":1-19999")
            else:
                cmdArgs.append(seqID)
        resData = subprocess.call(cmdArgs, stdout=resultFile, stderr=errorFile)
        errorFile.seek(0)
        resultFile.seek(0)
        print('<h2>Sequence</h2>')
        print('<p><em>[From "%s"]</em></p>' % myparams['queryDB'])
        print('<textarea cols=80 rows=10 id="fastaResult">')
        errorAdded = False
        for line in resultFile:
            line = line.rstrip()
            if(line.endswith(":1-19999")): # remove too-large indicator
                line = line[:-8]
            print line
        for line in errorFile:
            if(not errorAdded):
                print '** Error **'
                errorAdded = True
            print line.rstrip()
        resultFile.close()
        errorFile.close()
        print('</textarea>')

print '''</body>
</html>'''
