#!/usr/bin/env python3 -i

# Ryley Hindman, 2021-2022
# This is the main autotrac file.

## Import Statements
from PyQt5 import QtGui
from PyQt5.QtCore import QObject, pyqtSignal, QThread

from PyQt5.QtWidgets import QApplication, QMainWindow, QGraphicsScene,\
 QGraphicsView, QGraphicsItem, QGraphicsRectItem, QLabel, QVBoxLayout, QWidget,\
 QPushButton, QMenuBar, QMenu, QAction, QFileDialog
from PyQt5.QtGui import QPen, QBrush
from PyQt5.Qt import Qt
 
import numpy

from mpmath import *
mp.dps = 20
print(mp)

import sys

import time
import serial
import numpy
import math
import copy

from threading import Lock
from queue import Queue

## Global Variable Initialization

data_lock = Lock()
quit_lock = Lock()
point_lock = Lock()

newPost = False
pl1ToCalc = 0
pl2ToCalc = 0
pl3ToCalc = 0
newplToCalc = False
quit = False
keepRunning = True
pointButtonPushed = False
endRowTurnaround = False
fileQueue = Queue()

## Function: sign(x)
# This function returns:
# 1 if the input x is larger than 0
# 0 if the input x equals 0
# -1 if the input x is smaller than 0
def sign(x):
    if(x > 0):
        return 1
    elif(x == 0):
        return 0
    elif(x < 0):
        return -1

## Function: addToQueue(x,y,z,lat,longi)
#Add x,y,z,lat,longi to the fileQueue to put
#the values into the file
def addToQueue(x, y, z, lat, longi):
    global fileQueue
    #print("Adding to queue")
    fileQueue.put((x,y,z,lat,longi))

## Class: fileWorker(QObject)
# This class encapsulates the file reading and writing for saving the track
# log.
class fileWorker(QObject):
    global fileQueue
    finished = pyqtSignal()

    ## Function: runFileWorker(self)
    # This function is called once on the file save thread intialization.  It then
    # enters the infinite loop until the main thread requests a shutdown.
    def runFileWorker(self):
        # Need to clean up this run variable.  There is a run variable used in
        # the main loop, and I think it is for a different purpose.
        run = True
        self.openFileWrite('saveFile.txt')
        print("Opened file")
        while(run):
            # If the interuption is requested from the main thread, this will
            # attempt a clean shutdown of the file worker.
            if(QThread.currentThread().isInterruptionRequested()):
                run = False
                print("Interrupted")

            if(not fileQueue.empty()):
                data = fileQueue.get(False)
                #print("Committing1!")

                if(not data == None):
                #    print("Committing2!")
                    self.commitToFile(data)
            
        self.closeFile()

    ## Function: openFileRead(self, named)
    # This function opens a save file to read.  It ensures that the filename
    # ends with the .txt extension.
    def openFileRead(self, name):
        if(not name[-4:] == '.txt'):
            raise ReferenceError("Filename must end in \'.txt\'")
        self.saveFile = open(name, "r")
        print("Opened file read")

    ## Function: openFileWrite(self, name)
    # This function opens a save file to write.  It ensures that the filename
    # ends with the .txt extension.
    def openFileWrite(self, name):
        if(not name[-4:] == '.txt'):
            raise ReferenceError("Filename must end in \'.txt\'")
        self.saveFile = open(name, "w")
        print("Opened file write")

    ## Function commitToFile(self, data)
    # This function commits the array 'data' to the currently open file.  There
    # is no protection for a non-opened save file, so there is a possibility for
    # an exception to be thrown here.
    def commitToFile(self,data):
        if(not len(data) == 5):
            raise ReferenceError("Data length is not 5")

        self.saveFile.write(str(data[0]))
        self.saveFile.write(',')
        self.saveFile.write(str(data[1]))
        self.saveFile.write(',')
        self.saveFile.write(str(data[2]))
        self.saveFile.write(',')
        self.saveFile.write(str(data[3]))
        self.saveFile.write(',')
        self.saveFile.write(str(data[4]))
        self.saveFile.write('\n')

    ## Function closeFile(self)
    # This function closes the currently open text file.  There is no protection
    # for a non-opened save file, so there is a possibility for an exception to
    # be thrown here.
    def closeFile(self):
        self.saveFile.close()
        print("Closed file")


## Function: doesIntersect(x1,y1,x2,y2,p1x,p1y,p2x,p2y)
# This function calculates if the line defined by
# (p1x,p1y) -> (p2x,p2y) intersects the line defined
# by (x1,y1) -> (x2,y2).
# Currently, this function is unused, but could potentially be used to improve
# the graphics rendering efficiency.  It has not been tested.
def doesIntersect(x1,y1,x2,y2,p1x,p1y,p2x,p2y):
    if((p1x == p2x and p1y == p2y) or (x1 == x2 and y1 == y2)):
        return False
    v1x = x2 - x1
    v1y = y2 - y1

    v2x = p2x - p1x
    v2y = p2y - p1y

    denom = (v1x*v2y) - (v1y*v2x)

    #If both are parallel
    if(denom == 0):
        return False

    t = ((y1*v2x) - (p1y*v2x) - (x1*v2y) + (p1x*v2y)) / (denom)
    #print(t)
    if(not v2y == 0):
        r = ((t*v1y)+(y1 - p1y))/v2y
    else:
        r = ((t*v1x)+(x1-p1x))/v2x

    #print(r)

    if(t > 1 or t < 0 or r > 1):
        return False
    else:
        return True

## Function: doesIntersectRect(x1,y1,x2,y2,x3,y3,x4,y4,p1x,p1y,p2x,p2y)
# This function determines if the line defined by (p1x,p1y) -> (p2x,p2y)
# intersects the rectangle defined by (x1,y1)->(x2,y2)->(x3,y3)->(x4,y4).
# Returns True if the input line intersects the rectangle.
# Returns False if the input line does not intersect the rectangle.
# Requres doesIntersect() to work.
def doesIntersectRect(x1,y1,x2,y2,x3,y3,x4,y4,p1x,p1y,p2x,p2y):
    if(doesIntersect(x1,y1,x2,y2,p1x,p1y,p2x,p2y)):
        return True
    if(doesIntersect(x2,y2,x3,y3,p1x,p1y,p2x,p2y)):
        return True
    if(doesIntersect(x3,y3,x4,y4,p1x,p1y,p2x,p2y)):
        return True
    if(doesIntersect(x4,y4,x1,y1,p1x,p1y,p2x,p2y)):
        return True
    return False

## Function: getNMEA(gps)
# This function takes a gps input as a Serial object.  It then flushes the 
# input buffer, and waits for the whole GNGGA or GPGGA string to arrive.
# The function returns the string ending with GNGGA or GPGGA.
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

## Function: calcDist(x,y,z,v1,v2,v3)
# This function takes in an (x,y,z) position vector and a (v1,v2,v3) direction 
# vector.  It then calculates the length of (x,y,z) projected along(v1,v2,v3)
# from the origiin of (v1,v2,v3).  The use case is in the main loop, where we
# are calculating the distance from the plane.
def calcDist(x,y,z,v1,v2,v3):
    if(not(v1 == 0 and v2 == 0 and v3 == 0)):
        dist = mpf(-(x * v1 + y * v2 + z * v3)/mpf(mp.sqrt(mp.power(v1,2) +\
         mp.power(v2,2) + mp.power(v3,2))))
    else:
        dist = 0
    return dist

## Function: getECEF(lat,longi,height)
# This function uses the ECEF conversion functions to get earth-centered
# earth-fixed coordinates from the latitude, longitude, and height data.
def getECEF(lat,longi,height):
    latd = mp.radians(lat)
    longd = mp.radians(longi)
    a = 6378137
    b = float(6356752.3142)
    eSq = float(0.0066943799901377997230156324803829)
    N = a/(mp.sqrt(1-(eSq * mp.sin(latd)**2)))
    xOut = mpf((N + height) * mp.cos(latd) * mp.cos(longd))
    yOut = mpf((N + height) * mp.cos(latd) * mp.sin(longd))
    zOut = mpf(((mp.power(b,2) * N / mp.power(a,2)) + height) * mp.sin(latd))
    return xOut, yOut, zOut

## Function: parseNMEA(inputStr)
# This function takes an inputStr argument and parses it into the latitude,
# longitude, and height data.  It also provides the precision of the input
# data in the qualOut variable.
def parseNMEA(inputStr):
    GPGGAref = inputStr.rfind('$GNGGA')
    latOut = 0
    longOut = 0
    qualOut = 0
    geoidSepOut = 0
    if(GPGGAref != -1):
        inputStr = inputStr[GPGGAref+1:len(inputStr)-1]
        endDelim = inputStr.find('$')
        if(endDelim != -1):
            inputStr = inputStr[0:endDelim]
            data = inputStr.split(',')
            if(len(data) == 15):
                latStr = data[2]
                latDir = data[3]
                longStr = data[4]
                longDir = data[5]
                qualStr = data[6]
                altStr = data[9]
                geoidSepStr = data[11]
                if(latDir == 'S'):
                    latDMS = mpf(latStr)
                    latDeg = mp.floor(latDMS) * -0.01
                    latMinutes = latDMS - (latDeg * 100)
                    latOut = mpf(latDeg - (latMinutes / 60))
                else:
                    latDMS = mpf(latStr)
                    latDeg = mp.floor(latDMS) * 0.01
                    latMinutes = latDMS - (latDeg * 100)
                    latOut = mpf(latDeg + (latMinutes / 60))

                if(longDir == 'W'):
                    longDMS = mpf(longStr)
                    longDeg = mp.floor(longDMS) * -0.01
                    longMinutes = longDMS - (longDeg * 100)
                    longOut = mpf(longDeg - (longMinutes / 60))
                else:
                    longDMS = mpf(longStr)
                    longDeg = mp.floor(longDMS) * 0.01
                    longMinutes = longDMS - (longDeg * 100)
                    longOut = mpf(longDeg + (longMinutes / 60))

                qualOut = float(qualStr)
                altOut = float(altStr)
                geoidSepOut = float(geoidSepStr)
            else:
                qualOut = -2 #malformed input

        else:
            qualOut = -2 #malformed input

    else:
        qualOut = -1 #GPGGA not available

    #print("Lat: " + str(latOut) + " Long: " + str(longOut))
    return latOut, longOut, qualOut, altOut, geoidSepOut

## Function updatePosition(gps)
# This function takes a gps Serial object as an input.  It then gets the NMEA
# data from this gps input, parses it, and then returns the x,y,z values
# from the ECEF converting function, as well as the GPS precision
def updatePosition(gps):
    inputNMEA = getNMEA(gps)
    lat,longi,qual,alt,geoidSep = parseNMEA(inputNMEA)
    #print("Latitude: ")
    #print(lat)
    #print("Longitude: ")
    #print(longi)
    if (qual >= 0):
        height = alt + geoidSep
        x,y,z = getECEF(lat, longi, height)
    return mpf(x),mpf(y),mpf(z),qual

## Function calcPlaneThroughOrigin(x1,y1,z1,x2,y2,z2)
# This function calculates a plane through point (x1,y1,z1) and (x2,y2,z2) and
# the origin of the world.  It then returns a normal vector to this plane.
def calcPlaneThroughOrigin(x1,y1,z1,x2,y2,z2):
    #initial cross final position.  Normal points to
    # the right of the initial direction of travel
    v1 = mpf(y1 * z2 - z1 * y2)
    v2 = mpf(z1 * x2 - x1 * z2)
    v3 = mpf(x1 * y2 - y1 * x2)
    return v1, v2, v3

def calcVelocity(a,b,x1,y1,z1,x2,y2,z2,delT):
    t = mpf(math.sqrt((a**2 * b**2)/(x2**2 * b**2 + y2**2 * b**2 + z2**2 * a*2)))
    xEllipse = x2 * t
    yEllipse = y2 * t
    zEllipse = z2 * t
    n1 = mpf(2 * xEllipse / (a**2))
    n2 = mpf(2 * yEllipse / (a**2))
    n3 = mpf(2 * zEllipse / (b**2))
    normal = [n1, n2, n3]
    velocity = [x2 - x1, y2 - y1, z2 - z1]
   
    v1 = velocity[0] - normal[0] * (normal[0] * velocity[0] + normal[1] * velocity[1] + normal[2]\
     * velocity[2])/mpf(math.sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2))
    v2 = velocity[1] - normal[1] * (normal[0] * velocity[0] + normal[1] * velocity[1] + normal[2]\
     * velocity[2])/mpf(math.sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2))
    v3 = velocity[1] - normal[1] * (normal[0] * velocity[0] + normal[1] * velocity[1] + normal[2]\
     * velocity[2])/mpf(math.sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2))
    v1 = v1/delT
    v2 = v1/delT
    v3 = v1/delT

    return v1, v2, v3

# This function calculates the direction of travel projected onto the 
# tangent plane to the WGS84 ellipsoid at the current position, as well as
# the direction of the line projected onto the tangent plane.
# Inputs:
# a, b from WGS84.  x1, y1, z1, x2, y2, z2 from previous and current location
# n1, n2, n3 from line plane
# e1, e2 last and current error
def calcDirAndNormalProjection(a,b,x1,y1,z1,x2,y2,z2,n1,n2,n3,e1,e2):
    direction = numpy.array([x2-x1, y2-y1, z2-z1])
    norm = numpy.array([n1,n2,n3])

    t = mpf(math.sqrt((a**2 * b**2)/(x2**2 * b**2 + y2**2 * b**2 + z2**2 * a*2)))
    xEllipse = x2 * t
    yEllipse = y2 * t
    zEllipse = z2 * t
    t1 = mpf(2 * xEllipse / (a**2))
    t2 = mpf(2 * yEllipse / (a**2))
    t3 = mpf(2 * zEllipse / (b**2))

    tangent = numpy.array([t1, t2, t3])

    lineProj = numpy.cross(tangent, norm)
    dirProj = numpy.cross(numpy.cross(tangent, direction), tangent)

    lineProjMag = numpy.linalg.norm(lineProj)
    directionMag = numpy.linalg.norm(direction)
    if(not (lineProjMag == 0 or directionMag == 0)):
        sign = 0
        crossProduct = numpy.cross(lineProj, direction)
        dotProduct = numpy.dot(lineProj, direction)

        if(e1 > e2):
            sign = -1
        elif(e1 < e2):
            sign = 1

        phi = numpy.arcsin(float(sign * numpy.linalg.norm(crossProduct)/ \
            (lineProjMag * directionMag)))

        return phi
    else:
        return 0


# This object defines a world point, which is an x, y, z coordinate in the ECEF
# coordinate system.  This is a linked list data structure, so this object has
# iterate methods.  The iterate method just returns self.next, even if self.next
# does not exist.  The iterateSafe method returns self.next and True if self.next
# does exist.  Else, it returns self and False if self.next does not exist.
# The get function returns the current x, y, z coordinates for this object.
# clearNext clears self.next.
class WorldPoint():
    def __init__(self, xWorld, yWorld, zWorld, empty=False):
        if (not empty):
            self.x = xWorld
            self.y = yWorld
            self.z = zWorld
            self.next = None
        elif(empty):
            self.x = None
            self.y = None
            self.z = None
            self.next = None

    def add(self, point):
        self.next = copy.deepcopy(point)
        return self.next

    def iterate(self):
        return self.next

    def iterateSafe(self):
        if(not self.next is None):
            return self.next, True
        else:
            return self, False

    def get(self):
        return self.x, self.y, self.z

    def clearNext(self):
        self.next = None
        return self

# This object holds the conversion from world ECEF to scene coordinates.
class SceneCoordinates():
    def __init__(self,x0, y0, z0, x1, y1, z1, n1, n2, n3):
        self.basis = numpy.array([[0,0,0],
                                [0,0,0],
                                [0,0,0]])

        self.inverseBasis = numpy.array([[0,0,0],
                                        [0,0,0],
                                        [0,0,0]])

        self.readyToTransform = False

        self.origin = numpy.array([[0],[0],[0]])

        self.calcBasis(x0, y0, z0, x1, y1, z1, n1, n2, n3)


    #This function calculates the basis from (x0, y0, z0)
    #(x1, y1, z1) and <n1, n2, n3>
    def calcBasis(self, x0, y0, z0, x1, y1, z1, n1, n2, n3):
        v3 = numpy.array([n1,n2,n3])/numpy.linalg.norm([n1, n2, n3])
        v2 = numpy.array([x1-x0,y1-y0,z1-z0])/numpy.linalg.norm([x1- x0, y1-y0, z1-z0])

        v1 = numpy.cross(v2, v3)

        self.origin = numpy.array([[x0],[y0],[z0]])

        self.setBasis(v1, v2, v3)


    #This function sets the basis matrix to the 
    #basis vectors <v1> <v2> <v3>
    def setBasis(self, v1, v2, v3):
        self.basis = numpy.float64(numpy.array([[v1[0],v2[0],v3[0]],
                                [v1[1],v2[1],v3[1]],
                                [v1[2],v2[2],v3[2]]]))
        self.calcInverseBasis()


    #This function calculates the inverse matrix of the basis vectors
    #<v1> <v2> <v3>
    def calcInverseBasis(self):
        self.readyToTransform = True
        try:
            self.inverseBasis = numpy.linalg.inv(self.basis)
        except numpy.linalg.LinAlgError as e:
            #print(self.basis)
            #print(e)
            #print("Unable to calculate inverse basis")
            self.readyToTransform = False


    #This function takes care of the coordinate transformation from
    #world coordinates to local coordinates, with the origin at 
    #(Xb, Yb, Zb) and basis vectors <Xv, Yv, Zv>^(-1) (input inverse matrix)
    def calcSceneLocalCoordinates(self, x, y, z):
        localCoord = numpy.array([[0],[0],[0]])
        success = False
        if(self.readyToTransform):
            relativeCoord = numpy.array([[x],[y],[z]]) - self.origin
            #print("localcoord: ")
            #print(localCoord)
            localCoord = numpy.matmul(self.inverseBasis, numpy.float64(relativeCoord))
            success = True

        return localCoord, success

    


class Worker(QObject):
    finished = pyqtSignal()
    errorOut = pyqtSignal(float)
    steerOut = pyqtSignal(int)
    distOut = pyqtSignal(float)
    normVectOut = pyqtSignal(tuple)
    qualOut = pyqtSignal(int)
    localPointOut = pyqtSignal(WorldPoint)

    def runMainLoop(self):
        a = 6378137
        b = mpf(6356752.3142)
        run = 1
        firstPos = False
        secondPos = False
        thirdPos = False
        fileDialog = False
        x1 = mpf(0)
        y1 = mpf(0)
        z1 = mpf(0)
        x2 = mpf(0)
        y2 = mpf(0)
        z2 = mpf(0)
        x3 = mpf(0)
        y3 = mpf(0)
        z3 = mpf(0)
        x = mpf(0)
        y = mpf(0)
        z = mpf(0)
        lastX = mpf(0)
        lastY = mpf(0)
        lastZ = mpf(0)
        prevXPost = mpf(0)
        prevYPost = mpf(0)
        prevZPost = mpf(0)
        v1 = 0
        v2 = 0
        v3 = 0
        pl1 = mpf(0.0)
        pl2 = mpf(0.0)
        pl3 = mpf(0.0)
        offset = 1.803 # 5ft in meters
        rowNum = 0 # start on row 0
        rate = 10 # computation rate in hz
        sideOfLine = 1
        initDirV1 = 0
        initDirV2 = 0
        initDirV3 = 0
        sideOfLine = 0 # 1 if right, -1 if left
        saveFile = ""
        velTimer = time.time()
        first = 1
        quality = 0
        listBeginning = WorldPoint(0,0,0,True)
        lastError = 0
        kX = 0.25
        kPhi = 0.65
        kD = 0.03#0.13/4*1.3
        kDPhi = 0.1
        thirdPosCount = 0

        lastPoints = [[0], [0], [0]]

        mode = 1 #1 = straight track
        #2 = headland turn
        #3 = headland
        #4 = error

        travelDirection = 1 #1 = go to the right of initial track
        #2 = go tot he lef to initial track

        lineStartY = 0
        lastRowNum = 0
        wheelAngle = 0

        lastPhiNormalized = 0

        every4Phi = 0
        every4X = 0
        every4Y = 0
        every4Z = 0

        wheelAngleCounter = 0

        steeringOffset = 0

        steeringAngle = 0

        global newPost
        global pl1ToCalc
        global pl2ToCalc
        global pl3ToCalc
        global newplToCalc
        global keepRunning
        global pointButtonPushed
		global endRowTurnaround

        serialNumber = 0
        print('trying /dev/ttyACM0')
        retry = True
        while(retry):
            try:
                gps = serial.Serial(port = ('/dev/ttyACM'+str(serialNumber)), baudrate=9600,\
                 bytesize=8, timeout=10, stopbits=serial.STOPBITS_ONE)
                retry = False
            except serial.SerialException:
                serialNumber = serialNumber + 1
                if(serialNumber > 0):
                    raise serial.SerialException('Unable to connect to GPS.  Check the serial connection.')
                print('trying /dev/ttyACM'+str(serialNumber))

        arduino = serial.Serial(port = '/dev/ttyACM1', baudrate=9600,\
            bytesize=8, timeout=10, stopbits=serial.STOPBITS_ONE, write_timeout=0)

        f = open("fencePosts.txt", "a")

        while (run == 1):
            mainTimer = time.time()

            if(QThread.currentThread().isInterruptionRequested()):
                run = 0
                print("Interrupted")

            if(not thirdPos):
                point_lock.acquire()
                if(pointButtonPushed):
                    if(not firstPos):
                        x1, y1, z1, quality = updatePosition(gps)
                        print("First push")
                        print(str(x1) + " " + str(y1) + " " + str(z1) + " " + str(quality))
                        addToQueue(x1,y1,z1,0,0)

                        firstPos = True
                        pointButtonPushed = False
                    elif(not secondPos):
                        x2, y2, z2, quality = updatePosition(gps)
                        print("Second push")
                        pl1, pl2, pl3 = calcPlaneThroughOrigin(x1,y1,z1,x2,y2,z2) # Normal vector to plane through origin
                        addToQueue(x2,y2,z2,0,0)

						#Calculate normal to ellipse at point 2
                        t = mpf(math.sqrt((a**2 * b**2)/(x2**2 * b**2 + y2**2 * b**2 + z2**2 * a*2)))
                        xEllipse = x2 * t
                        yEllipse = y2 * t
                        zEllipse = z2 * t
                        n1 = mpf(2 * xEllipse / (a**2))
                        n2 = mpf(2 * yEllipse / (a**2))
                        n3 = mpf(2 * zEllipse / (b**2))
                        
                        localScene = SceneCoordinates(x1, y1, z1, x2, y2, z2, n1, n2, n3)

                        #Add 0,0,0 and the local coordinates for x2, y2, z2 to our linked list
                        pointsLinkedList = listBeginning.add(WorldPoint(0,0,0,True))
                        points, success = localScene.calcSceneLocalCoordinates(x2,y2,z2)
						lineStartY = points[1][0]
                        if(not success):
                            raise Exception('Local world coordinate calculation unsuccessful')
                        pointsLinkedList = pointsLinkedList.add(WorldPoint(points[0][0], points[1][0], points[2][0]))
                        #self.localPointOut.emit(pointsLinkedList)
						self.localPointOut.emit(WorldPoint(points[0][0], points[1][0], points[2][0]))

                        secondPos = True
                        pointButtonPushed = False
                    elif(not thirdPos):
                        initDirV1, initDirV2, initDirV3 = calcVelocity(a,b,x1,y1,z1,x2,y2,z2,1)
                        x3, y3, z3, quality = updatePosition(gps)
                        print("Third push")
                        distOff = calcDist(x3,y3,z3,pl1,pl2,pl3)
                        sideOfLine = int(distOff/abs(distOff))
                        thirdPos = True
                        pointButtonPushed = False
                else:
                    x, y, z, quality = updatePosition(gps)

                point_lock.release()

            if(thirdPos):
                row = 0
                error = 0
                effort = 0 #- for turning left + for turning right
                kp = 0 #unused for now
                x, y, z, quality = updatePosition(gps)
                #thirdPosCount = thirdPosCount + 1
                #print(thirdPosCount)
                points, success = localScene.calcSceneLocalCoordinates(x,y,z)
                print("x")
                print(points[0][0])
                print("y")
                print(points[1][0])
                print("z")
                print(points[2][0])
                if(not success):
                    raise Exception('Local world coordinate calculation unsuccessful')

                pointsLinkedList = pointsLinkedList.add(WorldPoint(points[0][0], points[1][0], points[2][0]))
                #self.localPointOut.emit(pointsLinkedList)

                distFromLast = sqrt((x-lastX)**2 + (y-lastY)**2 + (z-lastZ)**2)
                if(distFromLast > 0.5):
					self.localPointOut.emit(WorldPoint(points[0][0], points[1][0], points[2][0]))
                    addToQueue(x,y,z,0,0)
                    lastX = x
                    lastY = y
                    lastZ = z
                data_lock.acquire()
                if(newPost == True):
                    prevXPost = x
                    prevYPost = y
                    prevZPost = z
                    print("NEW POST")
                    newPost = False
                if(newplToCalc == True):
                    pl1 = pl1ToCalc
                    pl2 = pl2ToCalc
                    pl3 = pl3ToCalc
                    newplToCalc = False
                if(endRowTurnaround):
                    mode = 2
                    endRowTurnaround = False
                data_lock.release()

                if(first == 1):
                    xlast = x
                    ylast = y
                    zlast = z
                    first = 0
                    prevXPost = x
                    prevYPost = y
                    prevZPost = z

                # this is our actual distance from the original plane along
                # a normal to the plane
                actualValue = calcDist(x,y,z,pl1,pl2,pl3)

                rowNum = round(actualValue / offset)
                #print("rowNum")
                #print(rowNum)

                # this is a simple calculation for our desired distance
                # from the original plane
                rowNum = 0
                desiredValue = sideOfLine * offset * rowNum
                #print("desiredValue")
                #print(desiredValue)


                # this is our 3D earth-referenced velocity vector projected onto
                # a plane tangent to the WGS84 ellipse at our current lat/long
                v1, v2, v3 = calcVelocity(a,b,xlast,ylast,zlast,x,y,z, mpf(time.time() - velTimer))
                print("v1")
                print(v1)
                print("v2")
                print(v2)
                print("v3")
                print(v3)
                velTimer = time.time()
                velMagnitude = numpy.sqrt(v1**2 + v2**2 + v3**2)
                #print("velocity:")
                #print(velMagnitude)

                # no need for direction unit vector, because we only care about sign
                directionSign = initDirV1 * v1 + initDirV2 * v2 + initDirV3 * v3
                error = mpf(desiredValue) - actualValue
                #print("error:")
                #print(error)
                
                #phi = calcDirAndNormalProjection(a,b,xlast,ylast,zlast,x,y,z,pl1,pl2,pl3,lastError,error)
                phi = -1*math.atan2((points[0][0] - lastPoints[0][0]), (points[1][0] - lastPoints[1][0]))
                lastPoints = points
                #phiArr[phiIndex] = phi
                #phiIndex = phiIndex + 1
                #if(phiIndex >= numPhi):
                #    phiIndex = 0
                #phi = sum(phiArr)/len(phiArr)

                ######## Control Logic Start

                if(abs(phi) > 1.57):
                    phiNormalized = -1*sign(phi)*(abs(phi)-1.57)
                    error = -1*error
                else:
                    phiNormalized = phi

                #print("mode")
                #print(mode)
                #print("phiNorm")
                #print(phiNormalized)
                #print("deltaY")
                #print(lineStartY)
                #print("rowNum")
                #print(rowNum)
                #print("lastRowNum")
                #print(lastRowNum)
                controlType = 1
                if(mode == 1):
                    if (controlType == 0):
                        if(abs(error) > 2):
                            kPhiCorrected = 0
                        else:
                            kPhiCorrected = kPhi
                        steeringAngle = kX * error + kPhiCorrected * phiNormalized + kD * (error - lastError) + kDPhi * (phiNormalized - lastPhiNormalized)
                        if(steeringAngle > 0.3):
                            steeringAngle = 0.3
                        if(steeringAngle < -0.3):
                            steeringAngle = -0.3
                    elif(controlType == 1):
                        if(velMagnitude > 0):
                            kP = 0.938/(velMagnitude**2)
                            kD = 1.6/velMagnitude
                            print("velMagnitude")
                            print(velMagnitude)
                        else:
                            kP = 1
                            kD = 1
                            print("velocity is zero")
                        L = 1.72
                        Lh = 1.5 * velMagnitude
                        dE = error + Lh*numpy.sin(phiNormalized)
                        K = 0.3/L
                        if(phiNormalized != 0):
                            expont = -K*(numpy.sin(phiNormalized)*(kD*numpy.tan(phiNormalized) + kP*dE)/numpy.sin(phiNormalized) + Lh*((numpy.cos(phiNormalized))**4)*(kD*numpy.tan(phiNormalized) + kP*dE))
                        else:
                            expont = 0
                        print("expont")
                        print(expont)
                        #steeringAngle = -1*mp.atan(-K*L*(numpy.cos(phiNormalized))**3 * (1 - math.exp(expont))/(1 + math.exp(expont)))
                        steeringAngle = -1*mp.atan(-K*L*(numpy.cos(phiNormalized))**3 * -1 * mp.tanh(expont))
                        #top = -1 * L * numpy.sin(phiNormalized) * ((numpy.cos(phiNormalized))**3) * (kD * numpy.tan(phiNormalized) + kP*dE)
                        #bot = numpy.sin(phiNormalized) + Lh*((numpy.cos(phiNormalized))**4) * (kD * numpy.tan(phiNormalized) + kP*dE)
                        #steeringAngle = -1*mp.atan(top/bot)
                        #print("top")
                        #print(top)
                        #print("bottom")
                        #print(bot)
                        print("steering angle")
                        print(steeringAngle)
                elif(mode == 2):
                    if(travelDirection == 1 and (points[1][0] - lineStartY) > 0):
                        steeringAngle = 0.3
                    elif(travelDirection == 1 and (points[1][0] - lineStartY) < 0):
                        steeringAngle = -0.3
                    elif(travelDirection == 2 and (points[1][0] - lineStartY) > 0):
                        steeringAngle = -0.3
                    elif(travelDirection == 2 and (points[1][0] - lineStartY) < 0):
                        steeringAngle = 0.3

                    if(points[1][0] - lineStartY > 0 and abs(phi) > 1.7):
                        if(travelDirection == 1 and (rowNum - lastRowNum) == 1):
                            mode = 1
                            lineStartY = points[1][0]
                            lastRowNum = rowNum
                        elif(travelDirection == 2 and (rowNum - lastRowNum) == -1):
                            mode = 1
                            lineStartY = points[1][0]
                            lastRowNum = rowNum
                    elif(points[1][0] - lineStartY < 0 and abs(phi) < 1.44):
                        if(travelDirection == 1 and (rowNum - lastRowNum) == 1):
                            mode = 1
                            lineStartY = points[1][0]
                            lastRowNum = rowNum
                        elif(travelDirection == 2 and (rowNum - lastRowNum) == -1):
                            mode = 1
                            lineStartY = points[1][0]
                            lastRowNum = rowNum
        
                #~27.5 deg per steering wheel rev
                #~13.75 deg per motor rev
                #x deg at wheels x/27.5 * 360 deg at steering wheel
                #x deg at wheels x/27.5 * 360 * 2 deg at motor
                steeringAngleStr = (str(float(steeringAngle*180/3.1415))+'\n').encode('latin-1')
                print("steeringAngle")
                print(steeringAngleStr)
                print("steeringOffset")
                print(steeringOffset)
                print("phi")
                print(phi)

                arduino.write(steeringAngleStr)
				########## Control Logic End

                distFromPrev = math.sqrt((x-prevXPost)**2 + (y-prevYPost)**2 + (z-prevZPost)**2)

                steer = 0
                if (directionSign < 0 and error != 0):
                    # If we are going in the original direction, and error is positive
                    # then we are on the side of the line.  If the error is negative,
                    # then we are on the left side of the desired line.
                    if(error > 0):
                        steer = 1
                    else:
                        steer = -1

                lastError = error
                lastPhiNormalized = phiNormalized
                xlast = x
                ylast = y
                zlast = z
                
                #window.setTractorPos(float(error * -1))
                #window.setErrorString(float(error))
                self.errorOut.emit(float(error))
                self.steerOut.emit(steer)
                self.distOut.emit(float(distFromPrev))
                self.normVectOut.emit((mpf(pl1),mpf(pl2),mpf(pl3)))
                #App.processEvents()

            self.qualOut.emit(quality)
            while time.time() - mainTimer < (1/rate):
                pass

        print("At end of thread")
        keepRunning = False
        gps.close()
        arduino.close()
        self.finished.emit()

class CartesianPoint():
    def __init__(self, x, y, ID):
        self.x = x
        self.y = y
        self.parent = None
        self.upRight = None
        self.upLeft = None
        self.downRight = None
        self.downLeft = None
        self.ID = ID

    def getChild(self,type):
        if(type == 1):
            return self.upRight
        elif(type == 2):
            return self.upLeft
        elif(type == 3):
            return self.downLeft
        elif(type == 4):
            return self.downRight

    def addChild(self,type,child):
        if(type == 1):
            self.upRight = child
        elif(type == 2):
            self.upLeft = child
        elif(type == 3):
            self.downLeft = child
        elif(type == 4):
            self.downRight = child

class CartesianTree():
    def __init__(self):
        self.ID = 0
        self.root = CartesianPoint(0,0,self.ID)
        self.currPoints = []
        self.ID = 1

    def getQuadrant(self, parent, x, y):
        if(x > parent.x and y > parent.y):
            return 1
        elif(x < parent.x and y > parent.y):
            return 2
        elif(x < parent.x and y < parent.y):
            return 3
        elif(x > parent.x and y < parent.y):
            return 4
        elif(x > parent.x):
            if(not(parent.getChild(1) == None)):
                return 1
            else:
                return 4
        elif(x < parent.x):
            if(not(parent.getChild(2) == None)):
                return 2
            else:
                return 3
        elif(y > parent.y):
            if(not(parent.getChild(1) == None)):
                return 1
            else:
                return 2
        elif(y < parent.y):
            if(not(parent.getChild(3) == None)):
                return 3
            else:
                return 4
        else:
            return -1

    def addPoint(self, x, y):
        self.addPointRec(self.root, x, y)

    def addPointRec(self, parent, x, y):
        quad = self.getQuadrant(parent, x, y)
        if(parent.getChild(quad) == None):
            parent.addChild(quad,CartesianPoint(x,y,self.ID))
            self.ID = self.ID + 1
        else:
            self.addPointRec(parent.getChild(quad), x, y)

    def getPoints(self, xMax, xMin, yMax, yMin):
        self.currPoints = []
        self.getPointsRec(self.root, xMax, xMin, yMax, yMin)
        #This takes care of not returning an empty list.  This is OK because
        #the path drawing logic needs at least two IDs to draw anything
        if(len(self.currPoints) == 0):
            self.currPoints = [self.root]

        return self.currPoints

    def getPointsRec(self, parent, xMax, xMin, yMax, yMin):
        if(not(parent == None)):
            if(parent.x <= xMax and parent.x >= xMin and parent.y <= yMax and parent.y >= yMin):
                self.currPoints.append(parent)
                for i in range(1,5):
                    self.getPointsRec(parent.getChild(i), xMax, xMin, yMax, yMin)
            elif(parent.x >= xMax and parent.y >= yMax):
                self.getPointsRec(parent.getChild(3),xMax,xMin,yMax,yMin)
            elif(parent.x >= xMax and parent.y <= yMin):
                self.getPointsRec(parent.getChild(2),xMax,xMin,yMax,yMin)
            elif(parent.x <= xMin and parent.y >= yMax):
                self.getPointsRec(parent.getChild(4),xMax,xMin,yMax,yMin)
            elif(parent.x <= xMin and parent.y <= yMin):
                self.getPointsRec(parent.getChild(1),xMax,xMin,yMax,yMin)
            elif(parent.x >= xMax):
                self.getPointsRec(parent.getChild(2),xMax,xMin,yMax,yMin)
                self.getPointsRec(parent.getChild(3),xMax,xMin,yMax,yMin)
            elif(parent.x <= xMin):
                self.getPointsRec(parent.getChild(1),xMax,xMin,yMax,yMin)
                self.getPointsRec(parent.getChild(4),xMax,xMin,yMax,yMin)
            elif(parent.y >= yMax):
                self.getPointsRec(parent.getChild(3),xMax,xMin,yMax,yMin)
                self.getPointsRec(parent.getChild(4),xMax,xMin,yMax,yMin)
            elif(parent.y <= yMin):
                self.getPointsRec(parent.getChild(1),xMax,xMin,yMax,yMin)
                self.getPointsRec(parent.getChild(2),xMax,xMin,yMax,yMin)







keepRunning = True

class Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "PyQt5 QGraphicView"
        self.top = 200
        self.left = 500
        self.width = 750
        self.height = 400

        self.sceneModifier = 5
        self.sceneWidth = int(self.width/2)-self.sceneModifier
        self.sceneHeight = int(self.height)-self.sceneModifier


        self.lineWidth = 20
        self.thinLineWidth = 5
        self.veryThinLineWidth = 1
        self.trackLineWidth = 2
        self.oneFootPixels = int(self.width/8)
        self.error = 0.0
        self.pl1 = mpf(0)
        self.pl2 = mpf(0)
        self.pl3 = mpf(0)
        self.firstPoint = False
        self.secondPoint = False
        self.thirdPoint = False

        self.worldCoordsInit = False
        self.xMaxWorldDisplay = 0
        self.xMinWorldDisplay = 0
        self.yMaxWorldDisplay = 0
        self.yMinWorldDisplay = 0
        self.xMaxWorld = 0
        self.xMinWorld = 0
        self.yMaxWorld = 0
        self.yMinWorld = 0
        self.xMax = self.sceneWidth - 10
        self.xMin = 10
        self.yMax = self.sceneHeight - 10 + self.sceneModifier
        self.yMin = self.sceneModifier + 10
        #List to hold the track lines
        self.lines = []

        self.tree = CartesianTree()

        #Point to hold the points
        self.pointStart = WorldPoint(0,0,0,True)
        self.pointsIndex = self.pointStart
        self.redraw = False
        self.lastX = 0
        self.lastY = 0

        self.centerDrawMode = False

        self.implementWidth = 1.803 #implement width in meters

        self.scale = 20 #scale in meters
        self.scaleMax = 150
        self.scaleMin = 5

        self.InitWindow()

    def openFile(self):
        fileToOpen = QFileDialog.getOpenFileName(self, 'Open desired file',\
            'c:\\', "Text files (*.txt)")
        if(not(fileToOpen[0] == '')):
            f = open(fileToOpen[0],"r")
            pl1Temp = 0
            pl2Temp = 0
            pl3Temp = 0
            global pl1ToCalc
            global pl2ToCalc
            global pl3ToCalc
            global newplToCalc
            global quit
            global keepRunning
            line = f.readline()
            x = 0
            while((pl1Temp == 0 or pl2Temp == 0 or pl3Temp == 0) and not(line == '')):
                print(line)
                pos1 = line.find('pl1:')
                if(pos1>=0):
                    pl1Temp = mpf(line[pos1+4:].rstrip())
                    print('pl1')
                    print(pl1Temp)
                pos2 = line.find('pl2:')
                if(pos2>=0):
                    pl2Temp = mpf(line[pos2+4:].rstrip())
                    print('pl2')
                    print(pl2Temp)
                pos3 = line.find('pl3:')
                if(pos3>=0):
                    print('pl3')
                    pl3Temp = mpf(line[pos3+4:].rstrip())
                    print(pl3Temp)
                
                line = f.readline()
                x = x+1
                print(x)
            data_lock.acquire()
            pl1ToCalc = pl1Temp
            pl2ToCalc = pl2Temp
            pl3ToCalc = pl3Temp
            newplToCalc = True
            data_lock.release()
            f.close()

    def saveFile(self):
        fileToSave = QFileDialog.getSaveFileName(self, 'Save desired file',\
            '/', "Text files (*.txt)")
        print(fileToSave)
        if(not fileToSave == ('','')):
            f = open(fileToSave[0], "w")
            f.write('pl1:'+str(self.pl1)+'\n')
            f.write('pl2:'+str(self.pl2)+'\n')
            f.write('pl3:'+str(self.pl3)+'\n')

    def setNormVect(self, normVect):
        self.pl1 = normVect[0]
        self.pl2 = normVect[1]
        self.pl3 = normVect[2]

     
    def newPostButtonClicked(self):
        global newPost
        data_lock.acquire()
        newPost = True
        data_lock.release()

    def headlandTurnButtonClicked(self):
        global endRowTurnaround
        data_lock.acquire()
        endRowTurnaround = True
        data_lock.release()

    def pointButtonClicked(self):
        global pointButtonPushed
        if(not self.thirdPoint):
            point_lock.acquire()
            pointButtonPushed = True
            point_lock.release()
        if(not self.firstPoint):
            self.firstPoint = True
            self.pointButton.setText("Push At Point 2")
        elif(self.firstPoint and not self.secondPoint):
            self.pointButton.setText("Push At Point 3")
            self.secondPoint = True
        elif(self.secondPoint and not self.thirdPoint):
            self.pointButton.setText("Done")
            self.thirdPoint = True

    def InitWindow(self):
        self.setWindowIcon(QtGui.QIcon("icon.png"))
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width+20, self.height+20)

        self.errorReadout = QLabel(self)
        self.errorReadout.setText("0")
        font = self.errorReadout.font()
        font.setPointSize(40)
        self.errorReadout.setFont(font)
        self.errorReadout.setAlignment(Qt.AlignCenter)
        self.errorReadout.setGeometry(0,40,350, int(self.height/7))

        self.distReadout = QLabel(self)
        self.distReadout.setText("0")
        font = self.distReadout.font()
        font.setPointSize(40)
        self.distReadout.setFont(font)
        self.distReadout.setAlignment(Qt.AlignCenter)
        self.distReadout.setGeometry(0,int(self.height/7)+60,350, int(self.height/7))
        
        self.steerReadout = QLabel(self)
        self.steerReadout.setText("On Line")
        font = self.steerReadout.font()
        font.setPointSize(30)
        self.steerReadout.setFont(font)
        self.steerReadout.setAlignment(Qt.AlignCenter)
        self.steerReadout.setGeometry(0,int(self.height/7)*2 + 80,350, int(self.height/7))

        self.qualityReadout = QLabel(self)
        self.qualityReadout.setText("No Fix")
        font = self.qualityReadout.font()
        font.setPointSize(30)
        self.qualityReadout.setFont(font)
        self.qualityReadout.setAlignment(Qt.AlignCenter)
        self.qualityReadout.setGeometry(0,int(self.height/7)*4 + 130,350, int(self.height/7))

        self.nextPostButton = QPushButton(self)
        self.nextPostButton.setText("At Next Post Hole")
        self.nextPostButton.setGeometry(0,int(self.height/8)*3 + 100,350,int(self.height/8))
        self.nextPostButton.clicked.connect(self.newPostButtonClicked)

        self.pointButton = QPushButton(self)
        self.pointButton.setText("Push At Point 1")
        self.pointButton.setGeometry(0,int(self.height/8)*5 + 60,350, int(self.height/8))
        self.pointButton.clicked.connect(self.pointButtonClicked)

        self.zoomInButton = QPushButton(self)
        self.zoomInButton.setText("+")
        self.zoomInButton.setGeometry(int(self.width/2)-60,50,40,20)
        self.zoomInButton.clicked.connect(self.zoomIn)

        self.zoomOutButton = QPushButton(self)
        self.zoomOutButton.setText("-")
        self.zoomOutButton.setGeometry(int(self.width/2)-60,80,40,20)
        self.zoomOutButton.clicked.connect(self.zoomOut)

        self.viewModeButton = QPushButton(self)
        self.viewModeButton.setText("View Mode")
        self.viewModeButton.setGeometry(int(self.width/2) - 100, 110, 60, 30)
        self.viewModeButton.clicked.connect(self.toggleViewMode)

        self.headlandTurnButton = QPushButton(self)
        self.headlandTurnButton.setText("Headland Turn")
        self.headlandTurnButton.setGeometry(100, 110, 110, 30)
        self.headlandTurnButton.clicked.connect(self.headlandTurnButtonClicked)

        self.menuBar = QMenuBar(self)
        self.setMenuBar(self.menuBar)
        self.fileMenu = QMenu("&File",self)
        self.menuBar.addMenu(self.fileMenu)
        self.openAction = QAction("&Open",self)
        self.fileMenu.addAction(self.openAction)
        self.openAction.triggered.connect(self.openFile)

        self.saveAction = QAction("&Save",self)
        self.fileMenu.addAction(self.saveAction)
        self.saveAction.triggered.connect(self.saveFile)

        self.createGraphicView()
 
        self.show()
 
 
    def createGraphicView(self):
        self.scene  =QGraphicsScene()
        self.greenBrush = QBrush(Qt.green)
        self.grayBrush = QBrush(Qt.gray)
 
        self.greenPen = QPen(Qt.green)
        self.grayPen = QPen(Qt.gray)
        self.grayPen.setWidth(self.lineWidth)
        self.redPen = QPen(Qt.red)
        self.bluePen = QPen(Qt.blue)
        self.blueBrush = QBrush(Qt.blue)
        self.purplePen = QPen(Qt.magenta)
        self.purplePen.setWidth(self.trackLineWidth)
        self.purplePen.setCapStyle(Qt.FlatCap)

        graphicView = QGraphicsView(self.scene, self)
        #If any of these height modifiers are changed, we must also change the modifiers
        #in drawLine()
        graphicView.setGeometry(int(self.width/2),0,int(self.width/2),int(self.height)+5)
        graphicView.setSceneRect(0,0,self.sceneWidth,self.sceneHeight)
        self.initShapes()

    def initShapes(self):
        self.background = self.scene.addRect(0,0,int(self.width/2),int(self.height),self.greenPen,self.greenBrush) 
        self.trackLine = self.scene.addLine(0,0,0,int(self.height),self.grayPen)
        self.trackLine.setPos(int(self.width/4),0)
        self.grayPen.setWidth(self.thinLineWidth)
        self.xScaleLine = self.scene.addLine(0,0,int(self.width/2),0,self.grayPen)
        self.xScaleLine.setPos(0, int(self.height/2)-int(self.thinLineWidth/2))
        self.grayPen.setWidth(self.veryThinLineWidth)
        self.oneFootLineRight = self.scene.addLine(0,0,0,int(self.height/20))
        self.oneFootLineRight.setPos(3*self.oneFootPixels,\
            int(self.height/2)-int(self.height/40))
        self.oneFootLineLeft = self.scene.addLine(0,0,0,int(self.height/20))
        self.oneFootLineLeft.setPos(self.oneFootPixels,\
            int(self.height/2)-int(self.height/40))

        self.tractorWidth = int(self.width/40)
        self.tractorHeight = int(self.height/25)
        self.tractor = self.scene.addRect(0,0,self.tractorWidth,self.tractorHeight,self.bluePen,self.blueBrush)
        self.tractor.setPos(int(self.width/4)-int(self.tractorWidth/2),\
         int(self.height/2)-int(self.tractorHeight/2))

    def toggleViewMode(self):
        self.centerDrawMode = (not self.centerDrawMode)
        self.redraw = True

    def drawLine(self, x1, y1, x2, y2):
        #Note, this works in coordinates with (0,0) as the bottom left
        #corner, not the top left corner.
        #if(x1 < 0 or x2 < 0 or y1 < 0 or y2 < 0):
        #    print("All inputs must be positive")
        #else:
        lineLength = sqrt((x2-x1)**2 + (y2-y1)**2)
        self.lines.append(self.scene.addLine(0,0,int(x2-x1),int(y1-y2),self.purplePen))
        self.lines[len(self.lines)-1].setPos(x1, (self.height-y1))

    def setTractorPos(self, error):
        pixels = int(error * self.width * 0.125)
        self.tractor.setPos(int(self.width/4)-int(self.tractorWidth/2)+\
            pixels, int(self.height/2)-int(self.tractorHeight/2))

    def setErrorString(self, error):
        error = round(error, 3)
        self.errorReadout.setText(str(error))

    def setError(self, error):
        self.error = error

    def setSteer(self, steer):
        if(steer > 0):
            self.steerReadout.setText("Steer Right")
        elif(steer < 0):
            self.steerReadout.setText("Steer Left")
        else:
            self.steerReadout.setText("On Line")
        self.steerReadout.update()

    def setDist(self, dist):
        self.distReadout.setText(str(round(dist * 3.28,2)))

    def setQualityReadout(self, qual):
        if(qual == 1):
            self.qualityReadout.setText("GPS")
        elif(qual == 2):
            self.qualityReadout.setText("DGNSS")
        elif(qual == 4):
            self.qualityReadout.setText("RTK Fix")
        elif(qual == 5):
            self.qualityReadout.setText("RTK Float")
        else:
            self.qualityReadout.setText("No Fix")

    def zoomIn(self):
        if(self.centerDrawMode):
            self.scale = int(self.scale/1.5)

            if(self.scale > self.scaleMax):
                self.scale = self.scaleMax

    def zoomOut(self):
        if(self.centerDrawMode):
            self.scale = int(self.scale*1.5)

            if(self.scale < self.scaleMin):
                self.scale = self.scaleMin

    def setPoints(self, pointToDraw):
        #print("In set points")
        needToRedraw = False
        if(not self.worldCoordsInit):
            self.xMinWorld,self.yMinWorld,tmp = pointToDraw.get()
            #Default distance across = 20m
            self.xMaxWorld = self.xMinWorld + 20
            self.yMaxWorld = self.yMinWorld + 20

            self.tree.addPoint(self.xMinWorld,self.yMinWorld)

            self.xMaxWorldDisplay = self.xMaxWorld
            self.xMinWorldDisplay = self.xMinWorld
            self.yMaxWorldDisplay = self.yMaxWorld
            self.yMinWorldDisplay = self.yMinWorld

            self.worldCoordsInit = True
        else:
            xCoord, yCoord,tmp = pointToDraw.get()
            self.tree.addPoint(xCoord, yCoord)
            if(xCoord > self.xMaxWorld):
                self.xMaxWorld = xCoord
                needToRedraw = True
            elif(xCoord < self.xMinWorld):
                self.xMinWorld = xCoord
                needToRedraw = True

            if(yCoord > self.yMaxWorld):
                self.yMaxWorld = yCoord
                needToRedraw = True
            elif(yCoord < self.yMinWorld):
                self.yMinWorld = yCoord
                needToRedraw = True
            
            if(not self.centerDrawMode):

                if(needToRedraw):
                    self.redraw = True
                    needToRedraw = False     

                self.xMaxWorldDisplay = self.xMaxWorld
                self.xMinWorldDisplay = self.xMinWorld
                self.yMaxWorldDisplay = self.yMaxWorld
                self.yMinWorldDisplay = self.yMinWorld
                #print("Setting xMaxWorldDisplay to")
                #print(self.xMaxWorldDisplay)
            else:

                self.redraw = True
                               
                self.xMaxWorldDisplay = xCoord + self.scale
                self.xMinWorldDisplay = xCoord - self.scale
                self.yMaxWorldDisplay = yCoord + self.scale
                self.yMinWorldDisplay = yCoord - self.scale

        self.redraw = True


        self.drawPoints()

    #This function converts the pointToConvert values to the local coordinates
    #defined by the maximum and minimum world and scene values
    def convertToLocalCoords(self, pointToConvert):
        xWorld, yWorld, tmp = pointToConvert.get()
        '''
        print("xWorld")
        print(xWorld)
        print("yWOrld")
        print(yWorld)
        print("xmaxworld")
        print(self.xMaxWorld)
        print("xminworld")
        print(self.xMinWorld)
        print("xmaxworlddisplay")
        print(self.xMaxWorldDisplay)
        print("xminworlddisplay")
        print(self.xMinWorldDisplay)
        '''
        #Essentially, the equation gets the fraction of the world max/min, and
        #multiplies that by the length of the scene, plus the scene offset
        xScene = (xWorld - self.xMinWorldDisplay) * (self.xMax - self.xMin) / \
        (self.xMaxWorldDisplay - self.xMinWorldDisplay) + self.xMin

        yScene = (yWorld - self.yMinWorldDisplay) * (self.yMax - self.yMin) / \
        (self.yMaxWorldDisplay - self.yMinWorldDisplay) + self.yMin

        return xScene, yScene

    def convertToLocalCoords(self, xWorld, yWorld):
        #Essentially, the equation gets the fraction of the world max/min, and
        #multiplies that by the length of the scene, plus the scene offset
        xScene = (xWorld - self.xMinWorldDisplay) * (self.xMax - self.xMin) / \
        (self.xMaxWorldDisplay - self.xMinWorldDisplay) + self.xMin

        yScene = (yWorld - self.yMinWorldDisplay) * (self.yMax - self.yMin) / \
        (self.yMaxWorldDisplay - self.yMinWorldDisplay) + self.yMin

        return xScene, yScene

    #This function converts from scene coordinates to world coordinates
    #This is NOT USED
    def convertToWorldCoords(self, xScene, yScene):
        xWorld = (xScene - self.xMin)*(self.xMaxWorldDisplay - self.xMinWorldDisplay) / \
        (self.xMax - self.xMin) + self.xMinWorldDisplay

        yWorld = (yScene - self.yMin) * (self.yMaxWorldDisplay - self.yMinWorldDisplay) / \
        (self.yMax - self.yMin) + self.yMinWorldDisplay

        return xWorld, yWorld

    def drawPoints(self):
        if(self.redraw):

            if(self.xMaxWorldDisplay == self.xMinWorldDisplay):
                self.trackLineWidth = 1
            else:
                self.trackLineWidth = int(abs(self.implementWidth * (self.xMax - self.xMin) / (self.xMaxWorldDisplay - self.xMinWorldDisplay)))
            if(self.trackLineWidth == 0):
                self.trackLineWidth = 1
            if(self.trackLineWidth == 1):
                self.purplePen.setColor(Qt.red)
            else:
                self.purplePen.setColor(Qt.magenta)
            #print("Width: ")
            #print(self.trackLineWidth)

            self.purplePen.setWidth(self.trackLineWidth)

            for currentLine in self.lines:
                self.scene.removeItem(currentLine)

            self.lines = []


            #print("xMax")
            #print(self.xMaxWorldDisplay)
            #print("xMin")
            #print(self.xMinWorldDisplay)
            #print("yMax")
            #print(self.yMaxWorldDisplay)
            #print("yMin")
            #print(self.yMinWorldDisplay)
            linesInScene = self.tree.getPoints(self.xMaxWorldDisplay,self.xMinWorldDisplay, \
                self.yMaxWorldDisplay, self.yMinWorldDisplay)

            linesInScene.sort(key=lambda x: x.ID)

            lastID = -2
            lastLine = None
            if(len(linesInScene) > 150):
                linesInScene = linesInScene[len(linesInScene)-150:]
            for line in linesInScene:
                #print("ID")
                #print(line.ID)
                #print("Last ID")
                #print(lastID)
                if(line.ID-lastID == 1):
                    lastXLocal, lastYLocal = self.convertToLocalCoords(lastLine.x, lastLine.y)
                    xLocal, yLocal = self.convertToLocalCoords(line.x, line.y)
                    #print("xLocal")
                    #print(xLocal)
                    #print("yLocal")
                    #print(yLocal)
                    self.drawLine(lastXLocal, lastYLocal, xLocal, yLocal)
                lastID = line.ID
                lastLine = line


            '''
            self.pointsIndex, hasNext = self.pointStart.iterateSafe()
            firstRedrawPoint = True
            while(hasNext):
                if(firstRedrawPoint):
                    self.lastX, self.lastY = self.convertToLocalCoords(self.pointsIndex)
                    firstRedrawPoint = False
                else:
                    currX, currY = self.convertToLocalCoords(self.pointsIndex)
                    if(self.centerDrawMode):
            '''
            '''
                        if((currX >= self.xMin and currX <= self.xMax and \
                            currY >= self.yMin and currY <= self.yMax) or \
                            (self.lastX >= self.xMin and self.lastX <= self.xMax and \
                            self.lastY >= self.yMin and self.lastY <= self.yMax)):

                            self.drawLine(self.lastX, self.lastY, currX, currY)
                        elif(doesIntersectRect(self.xMin,self.yMin,self,xMin,self.yMax,\
                            self.xMax,self.yMax,self.xMax,self.yMin,self.lastX,self.lastY,\
                            currX,currY)):
            '''
            '''
                        if(not(currX == self.lastX or currY == self.lastY)):
                            self.drawLine(self.lastX, self.lastY, currX, currY)
                    else:
                        if(not(currX == self.lastX or currY == self.lastY)):
                            self.drawLine(self.lastX, self.lastY, currX, currY)
                    self.lastX = currX
                    self.lastY = currY
                self.pointsIndex, hasNext = self.pointsIndex.iterateSafe()
            '''
            self.redraw = False
        else:
            currX, currY = self.convertToLocalCoords(self.pointsIndex)
            if(self.lastX == None or self.lastY == None):
                self.lastX = currX
                self.lastY = currY
            else:
                '''
                print("Lastx")
                print(self.lastX)
                print("Lasty")
                print(self.lastY)
                print("CurrX")
                print(currX)
                print("CurrY")
                print(currY)
                '''
                if(not(currX == self.lastX or currY == self.lastY)):
                    self.drawLine(self.lastX, self.lastY, currX, currY)

                self.lastX = currX
                self.lastY = currY

    def runMainTask(self):
        self.thread = QThread()
        self.worker = Worker()
        self.worker.moveToThread(self.thread)

        self.thread.started.connect(self.worker.runMainLoop)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.worker.errorOut.connect(self.setError)
        self.worker.steerOut.connect(self.setSteer)
        self.worker.distOut.connect(self.setDist)
        self.worker.normVectOut.connect(self.setNormVect)
        self.worker.qualOut.connect(self.setQualityReadout)
        self.worker.localPointOut.connect(self.setPoints)

        self.thread.start()

    def runSaveTask(self):
        self.saveThread = QThread()
        self.saveWorker = fileWorker()
        self.saveWorker.moveToThread(self.saveThread)

        self.saveThread.started.connect(self.saveWorker.runFileWorker)
        self.saveWorker.finished.connect(self.saveThread.quit)
        self.saveWorker.finished.connect(self.saveWorker.deleteLater)
        self.saveThread.finished.connect(self.saveThread.deleteLater)

        self.saveThread.start()

    def closeEvent(self, event):
        self.thread.requestInterruption()
        self.saveThread.requestInterruption()


App = QApplication(sys.argv)
window = Window()
print(window.centerDrawMode)
window.runMainTask()
#window.runSaveTask()

#print(doesIntersect(0,0,1,5,-1,1,3,2))
#print(doesIntersect(0,0,0,5,0,0,-1,1))


while(keepRunning):
    window.setErrorString(float(window.error * 3.28))
    window.setTractorPos(float(window.error * -1))
    App.processEvents()

App.quit()

#sys.exit(App.exec_())