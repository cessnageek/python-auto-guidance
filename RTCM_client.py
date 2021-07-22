#!/usr/bin/env python3 -i

import time
import serial
import base64
import socket
from errno import ENETUNREACH, EPIPE
from userPassVariable import userPass

url = '165.206.203.10'
mountpoint = 'RTCM3_IMAX'

userPass64 = base64.b64encode(userPass).decode('ascii')
header =\
"GET /RTCM3_IMAX HTTP/1.1\r\n" +\
"Host http://165.206.203.10:10000\r\n" +\
"Ntrip-Version: Ntrip/2.0\r\n" +\
"User-Agent: NTRIP Client 1.0\r\n\r\n" +\
"Authorization: Basic {}\r\n".format(userPass64)

connectToServer = False


try:
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.settimeout(10.0)
    s.connect((url,10000))
    print("connected")
    s.send(header.encode('ascii'))
except socket.timeout:
    print("Initial connection timed out.. trying again")
    connectToServer = True
except IOError as e:
    if e.errno == ENETUNREACH:
        print("Internet unreachable... trying again")
        connectToServer = True
    elif e.errno == EPIPE:
        print("Broken pipe... trying again")
        connectToServer = True
    else:
        print("Unknown issue... trying again")
        connectToServer = True

RTCMDelay = 2
RTCMTimer = time.time()
RTCMToServer = ''

gps = serial.Serial(port = "/dev/ttyUSB0", baudrate=38400,\
    bytesize=8, timeout=10, stopbits=serial.STOPBITS_ONE)

print("Before loop")

def getNMEA(gps):
    gps.flushInput()
    gpggaFound = False
    endOfGPGGA = False
    outStrChars = []
    while(not endOfGPGGA):
        end = len(outStrChars)-1
        outStrChars.append(gps.read(1).decode('utf-8'))
        if(len(outStrChars) > 4 and (outStrChars[end] == 'A')):
               if(outStrChars[end-1] == 'G'):
                   if(outStrChars[end-2] == 'G'):
                       if(outStrChars[end-3] == 'P' or outStrChars[end-3] == 'N'):
                           if(outStrChars[end-4] == 'G'):
                               gpggaFound = True; 
        if(gpggaFound and outStrChars[end] == '$'):
            endOfGPGGA = True
    return ''.join(outStrChars)

while(True):
    if(time.time() - RTCMTimer > RTCMDelay):
        if(connectToServer):
            print("Connecting...")
            connectToServer = False

            try:
                s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
                s.settimeout(10.0)
                s.connect((url,10000))
                print("connected")
                s.send(header.encode('ascii'))
            except socket.timeout:
                print("Initial connection timed out.. trying again")
                connectToServer = True
            except IOError as e:
                if e.errno == ENETUNREACH:
                    print("Internet unreachable... trying again")
                    connectToServer = True
                elif e.errno == EPIPE:
                    print("Broken pipe... trying again")
                    connectToServer = True
                else:
                    print("Unknown issue... trying again")
                    connectToServer = True

        if(not connectToServer):
            dataFromServer = b''
            #print("Delta: " + str(time.time() - RTCMTimer))
            RTCMToServer = getNMEA(gps)
            print(RTCMToServer)
            try:
                s.send(RTCMToServer.encode('ascii'))
                dataFromServer = s.recv(1024)
            except socket.timeout:
                print("Socket timed out... trying again.")
                connectToServer = True
            except IOError as e:
                if e.errno == ENETUNREACH:
                    print("Internet unreachable... trying again")
                    connectToServer = True
                elif e.errno == EPIPE:
                    print("Broken pipe... trying again")
                    connectToServer = True
                else:
                    print("Unknown issue... trying again")
                    connectToServer = True

            gps.write(dataFromServer)

        RTCMTimer = time.time()

input("Press enter to end: ")