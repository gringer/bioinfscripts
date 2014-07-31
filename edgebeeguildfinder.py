#!/usr/bin/python

# Edgebee guild finder -- given a list of user names,
# work out what edgebee guild those users come from. Output is csv,
# containing User, Game, and Guild

# also found at http://pastebin.com/W6sYanuQ

# Author: David Eccles (gringer) <coding@gringer.org>

# usage: echo -e "user1\nuser2" | ./edgebeeguildfinder.py

import cookielib, urllib2
import fileinput
import re
import sys
import json
import csv

userNames = []
addedUsers = set()

for line in fileinput.input():
    if(line.startswith("{")):
        jsonObj = json.loads(line)
        if(('result' in jsonObj) and ('players' in jsonObj['result'])):
            for userName in map(lambda x: x['name'],jsonObj['result']['players']):
                if(not userName in addedUsers):
                    userNames.append(userName)
                    addedUsers.add(userName)
    else:
        userName = line.rstrip()
        if(not userName in addedUsers):
            userNames.append(userName)
            addedUsers.add(userName)

cj = cookielib.CookieJar()
opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cj))
f = opener.open('http://www.edgebee.com/signin?' +
                'username=gringerscripts&password=gringerscripts&remember=1')

cr = csv.writer(sys.stdout)
cr.writerow(['User','Game','Guild','LastOnline','Registered'])
for user in userNames:
    foundSandP = False
    try:
       f = opener.open('http://www.edgebee.com/user?name=%s' % (user))
    except:
       continue
    inItem = False
    gameName = None
    extraStat = None
    extraStats = {}
    extraStats['lastOn'] = None
    extraStats['regOn'] = None
    extraStats['guildName'] = None
    for line in f.read().splitlines(True):
        line = line.rstrip()
        if('playerListItem' in line):
            inItem = True
            gameName = None
            extraStats['guildName'] = None
        if(inItem):
            if('h2' in line):
                gameName = re.compile(r'<[^>]+>').sub('', line).lstrip()
            if('h4' in line):
                inItem = False
        elif('h4' in line):
            if('Registered on' in line):
                extraStat = 'regOn'
            elif('Last online' in line):
                extraStat = 'lastOn'
            else:
                extraStat = None
        elif((extraStat is not None) and ('<br>' in line)):
            extraStats[extraStat] = re.compile(r'<[^>]+>').sub('', line).lstrip()
        elif('Guild:' in line):
            guildName = re.compile(r'<[^>]+>').sub('', line).lstrip()
            guildName = guildName.replace('Guild: ','')
            extraStats['guildName'] = guildName
        if(('<br style="clear:both"/>' in line)  and (gameName is not None)):
            cr.writerow([user, gameName, extraStats['guildName'],
                                extraStats['lastOn'],extraStats['regOn']])
