#!/usr/bin/env python3 -i

from PyQt5 import QtGui
from PyQt5.QtCore import QObject, pyqtSignal, QThread

from PyQt5.QtWidgets import QApplication, QMainWindow, QGraphicsScene,\
 QGraphicsView, QGraphicsItem, QGraphicsRectItem, QLabel, QVBoxLayout, QWidget,\
 QPushButton, QMenuBar, QMenu, QAction, QFileDialog
from PyQt5.QtGui import QPen, QBrush
from PyQt5.Qt import Qt
 
from mpmath import *
mp.dps = 20
print(mp)

import sys

import time
import serial
import numpy
import math

from threading import Lock
from queue import Queue

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
fileQueue = Queue()

#Add x,y,z,lat,longi to the fileQueue to put
#the values into the file
def addToQueue(x, y, z, lat, longi):
    global fileQueue
    print("Adding to queue")
    fileQueue.put((x,y,z,lat,longi))

class fileWorker(QObject):
    global fileQueue
    finished = pyqtSignal()

    def runFileWorker(self):
        run = True
        self.openFileWrite('saveFile.txt')
        print("Opened file")
        while(run):
            if(QThread.currentThread().isInterruptionRequested()):
                run = False
                print("Interrupted")

            if(not fileQueue.empty()):
                data = fileQueue.get(False)
                print("Committing1!")

                if(not data == None):
                    print("Committing2!")
                    self.commitToFile(data)
            
        self.closeFile()

    def openFileRead(self, name):
        if(not name[-4:] == '.txt'):
            raise ReferenceError("Filename must end in \'.txt\'")
        self.saveFile = open(name, "r")
        print("Opened file read")

    def openFileWrite(self, name):
        if(not name[-4:] == '.txt'):
            raise ReferenceError("Filename must end in \'.txt\'")
        self.saveFile = open(name, "w")
        print("Opened file write")

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

    def closeFile(self):
        self.saveFile.close()
        print("Closed file")


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

def calcDist(x,y,z,v1,v2,v3):
    if(not(v1 == 0 and v2 == 0 and v3 == 0)):
        dist = mpf(-(x * v1 + y * v2 + z * v3)/mpf(mp.sqrt(mp.power(v1,2) +\
         mp.power(v2,2) + mp.power(v3,2))))
    else:
        dist = 0
    return dist

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
                    latOut = mpf(latDeg + (latMinutes / 60))
                else:
                    latDMS = mpf(latStr)
                    latDeg = mp.floor(latDMS) * 0.01
                    latMinutes = latDMS - (latDeg * 100)
                    latOut = mpf(latDeg + (latMinutes / 60))

                if(longDir == 'W'):
                    longDMS = mpf(longStr)
                    longDeg = mp.floor(longDMS) * -0.01
                    longMinutes = longDMS - (longDeg * 100)
                    longOut = mpf(longDeg + (longMinutes / 60))
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

def updatePosition(gps):
    inputNMEA = getNMEA(gps)
    lat,longi,qual,alt,geoidSep = parseNMEA(inputNMEA)
    if (qual >= 0):
        height = alt + geoidSep
        x,y,z = getECEF(lat, longi, height)
    return mpf(x),mpf(y),mpf(z),qual

def calcPlaneThroughOrigin(x1,y1,z1,x2,y2,z2):
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

class Worker(QObject):
    finished = pyqtSignal()
    errorOut = pyqtSignal(float)
    steerOut = pyqtSignal(int)
    distOut = pyqtSignal(float)
    normVectOut = pyqtSignal(tuple)
    qualOut = pyqtSignal(int)
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
        prevXPost = mpf(0)
        prevYPost = mpf(0)
        prevZPost = mpf(0)
        v1 = 0
        v2 = 0
        v3 = 0
        pl1 = mpf(0.0)
        pl2 = mpf(0.0)
        pl3 = mpf(0.0)
        offset = 1.52 # 5ft in meters
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
        global newPost
        global pl1ToCalc
        global pl2ToCalc
        global pl3ToCalc
        global newplToCalc
        global keepRunning
        global pointButtonPushed

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
                if(serialNumber > 5):
                    raise serial.SerialException('Unable to connect to GPS.  Check the serial connection.')
                print('trying /dev/ttyACM'+str(serialNumber))

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
                        pl1, pl2, pl3 = calcPlaneThroughOrigin(x1,y1,z1,x2,y2,z2)
                        addToQueue(x2,y2,z2,0,0)
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
                    print("Updating")
                    x, y, z, quality = updatePosition(gps)

                point_lock.release()

            if(thirdPos):
                row = 0
                error = 0
                effort = 0 #- for turning left + for turning right
                kp = 0 #unused for now
                x, y, z, quality = updatePosition(gps)
                addToQueue(x,y,z,0,0)
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
                data_lock.release()

                if(first == 1):
                    xlast = x
                    ylast = y
                    zlast = z
                    first = 0
                    prevXPost = x
                    prevYPost = y
                    prevZPost = z

                # this is a simple calculation for our desired distance
                # from the original plane
                desiredValue = sideOfLine * offset * rowNum

                # this is our actual distance from the original plane along
                # a normal to the plane
                actualValue = calcDist(x,y,z,pl1,pl2,pl3)

                # this is our 3D earth-referenced velocity vector projected onto
                # a plane tangent to the WGS84 ellipse at our current lat/long
                v1, v2, v3 = calcVelocity(a,b,xlast,ylast,zlast,x,y,z, mpf(time.time() - velTimer))
                velTimer = time.time()
                xlast = x
                ylast = y
                zlast = z

                # no need for direction unit vector, because we only care about sign
                directionSign = initDirV1 * v1 + initDirV2 * v2 + initDirV3 * v3
                error = mpf(desiredValue) - actualValue
                distFromPrev = math.sqrt((x-prevXPost)**2 + (y-prevYPost)**2 + (z-prevZPost)**2)
                #print("Error: " + str(actualValue))
                #print("x: " + str(x) + "y: " + str(y) + "z: " + str(z))
                print("pl1: " + str(pl1) + "pl2: " + str(pl2) + "pl3: " + str(pl3))
                steer = 0
                if (directionSign < 0 and error != 0):
                    # If we are going in the original direction, and error is positive
                    # then we are on the side of the line.  If the error is negative,
                    # then we are on the left side of the desired line.
                    if(error > 0):
                        steer = 1
                    else:
                        steer = -1
                
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
        self.finished.emit()


keepRunning = True

class Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "PyQt5 QGraphicView"
        self.top = 200
        self.left = 500
        self.width = 750
        self.height = 400
        self.lineWidth = 20
        self.thinLineWidth = 5
        self.veryThinLineWidth = 1
        self.oneFootPixels = int(self.width/8)
        self.error = 0.0
        self.pl1 = mpf(0)
        self.pl2 = mpf(0)
        self.pl3 = mpf(0)
        self.firstPoint = False
        self.secondPoint = False
        self.thirdPoint = False
 
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

        graphicView = QGraphicsView(self.scene, self)
        graphicView.setGeometry(int(self.width/2),0,int(self.width/2),int(self.height)+5)
        graphicView.setSceneRect(0,0,int(self.width/2)-5,int(self.height)-5)
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
window.runMainTask()
window.runSaveTask()

while(keepRunning):
    window.setErrorString(float(window.error * 3.28))
    window.setTractorPos(float(window.error * -1))
    App.processEvents()

App.quit()

#sys.exit(App.exec_())


