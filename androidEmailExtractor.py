#!/usr/bin/python3

## androidEmailExtractor.py -- extract emails from Android mail databases

## Copyright (C) 2018 David Eccles (gringer) <bioinformatics@gringene.org>

#############################################################################
##                                                                         ##
##  This program is free software: you can redistribute it and/or modify   ##
##  it under the terms of the GNU General Public License as published by   ##
##  the Free Software Foundation, either version 3 of the License, or      ##
##  (at your option) any later version.                                    ##
##                                                                         ##
##  This program is distributed in the hope that it will be useful,        ##
##  but WITHOUT ANY WARRANTY; without even the implied warranty of         ##
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          ##
##  GNU General Public License for more details.                           ##
##                                                                         ##
##  You should have received a copy of the GNU General Public License      ##
##  along with this program.  If not, see <http://www.gnu.org/licenses/>.  ##
##                                                                         ##
#############################################################################

## This expects the files 'EmailProvider.db' and 'EmailProviderBody.db' to
## be in the current directory. It will write to standard output the emails
## from those SQLite databases (as SMTP-ish format), assuming a specific
## mail storage format that matches at least one version of Android

## Note: only the text version of emails will be displayed

## Note: any emails that have a mailBoxKey of 6 will be ignored;
##       this hopefully corresponds with the "TRASH" folder

## The inspiration for this code was from here:
## http://code.activestate.com/recipes/578106-android-gmail-export-messages-from-sqlite-database/

from datetime import datetime
from textwrap import wrap, fill
import sqlite3

## Load message text bodies into memory

conn = sqlite3.connect('EmailProviderBody.db')
cursor = conn.cursor()
cursor.execute("select messageKey, textContent, textReply from Body")
rows = cursor.fetchall()
textReplies = dict()
textContents = dict()

for row in rows:
    textReplies[row[0]] = row[1]
    textContents[row[0]] = row[2]

cursor.close()
conn.close()

## Load message metadata into memory

conn = sqlite3.connect('EmailProvider.db')
cursor = conn.cursor()
cursor.execute("select _rowid_, fromList, toList, ccList, bccList, replyToList, subject, snippet, timeStamp, mailboxKey from Message")
rows = cursor.fetchall()

fnum = 0

for row in rows:
    fname = str(fnum)
    fnum += 1
    if(int(row[9]) == 6): ## ignore things in TRASH folder
        continue
    #print(fname)
    print("-- Message #%s --" % row[0])
    fromList = row[1].split("")
    print("From: \"%s\" <%s>" % (
        (fromList[1] if (len(fromList) > 1) else ""),
        fromList[0]))
    if((row[2] is not None)):
        toList = row[2].split("");
        toStrList = list()
        for toAddr in toList:
            toAddrList = toAddr.split("")
            toStrList.append("\"%s\" <%s>" % (
                (toAddrList[1] if (len(toAddrList) > 1) else ""),
                toAddrList[0]))
        print("To: %s" % "; ".join(toStrList))
    if(row[3] is not None):
        ccList = row[3].split("");
        ccStrList = list()
        for ccAddr in ccList:
            ccAddrList = ccAddr.split("")
            ccStrList.append("\"%s\" <%s>" % (
                (ccAddrList[1] if (len(ccAddrList) > 1) else ""),
                ccAddrList[0]))
        print("Cc: %s" % "; ".join(toStrList))
    if(row[4] is not None):
        bccList = row[4].split("");
        bccStrList = list()
        for bccAddr in bccList:
            bccAddrList = bccAddr.split("")
            bccStrList.append("\"%s\" <%s>" % (
                (bccAddrList[1] if (len(bccAddrList) > 1) else ""),
                bccAddrList[0]))
        print("Bcc: %s" % "; ".join(toStrList))
    if(row[5] is not None):
        rtList = row[5].split("");
        rtStrList = list()
        for rtAddr in rtList:
            rtAddrList = rtAddr.split("")
            rtStrList.append("\"%s\" <%s>" % (
                (rtAddrList[1] if (len(rtAddrList) > 1) else ""),
                rtAddrList[0]))
        print("Reply to: %s" % "; ".join(toStrList))
    print("Date: %s" % datetime.fromtimestamp(row[8] / 1000).strftime("%a, %d %B %Y %H:%M:%S"))
    print("Subject: %s" % str(row[6]))
    if(row[0] in textReplies and (textReplies[row[0]] is not None)
       and (textContents[row[0]] is not None)):
        repList = textReplies[row[0]].splitlines()
        for repLine in repList:
            repLineChunks = wrap(repLine, 90)
            for chunk in repLineChunks:
                print("> %s" % chunk)
        print()
    if(row[0] in textContents and (textContents[row[0]] is not None)):
        strList = textContents[row[0]].splitlines()
        for strLine in strList:
            print(fill(strLine, 90))
    if(row[0] in textReplies and (textContents[row[0]] is None)):
        strList = textReplies[row[0]].splitlines()
        for strLine in strList:
            print(fill(strLine, 90))
    print()
    
cursor.close()
conn.close()

print("-- End of messages --")
