from PyQt5 import QtGui
from PyQt5.QtCore import QObject, pyqtSignal, QThread

from PyQt5.QtWidgets import QApplication, QMainWindow, QGraphicsScene,\
 QGraphicsView, QGraphicsItem, QGraphicsRectItem, QLabel, QVBoxLayout, QWidget,\
 QPushButton, QMenuBar, QMenu, QAction, QFileDialog
import PyQt5.QtWidgets as QtWidgets
import PyQt5.QtCore as QtCore
from PyQt5.QtGui import QPen, QBrush
from PyQt5.Qt import Qt
import serial
import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import socket
import sys
from threading import Lock
from mpmath import *
mp.dps = 30
print(mp)

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
xs = []
ts = []

data_lock = Lock()
network_lock = Lock()


#desiredWheelAngle = 20000
global wheelAngle
global runKinematics
global setXVal
global connect
global connectAddr
global resetFlag

runKinematics = True
setXVal = 0
connect = False
connectAddr = "192.168.86.33"
resetFlag = False

lastTimeMain = time.time()

wheelAngle = np.int16(0)


#stepperDriver = serial.Serial(port = ('COM3'),baudrate=115200,\
#bytesize=8, timeout=10, stopbits=serial.STOPBITS_ONE, write_timeout=0)


def importWheelAngle(prevAngle):
	stepperDriver.flushInput()
	data = stepperDriver.readline()
	#print(data)
	try:
		return int(data)
	except:
		return prevAngle
		print("prev ANGLE")


def exportDesiredWheelAngle(desiredWheelAngle):
	stepperDriver.write(bytes(str(desiredWheelAngle)+'\n','utf-8'))

class controller:
	def __init__(self,kX,kPhi,kD,kDPhi):
		self.kX = kX
		self.kPhi = kPhi
		self.kD = kD
		self.kDPhi = kDPhi
		self.error = 0
		self.setPoint = 0
		self.lastError = 0
		self.phi = 0
		self.lastPhi = 0
		self.currTime = 0
		self.lastTime = 0
		self.steeringAngle = 0

	def runController(self, phi, x):
		self.phi = phi
		self.error = x - self.setPoint
		self.currTime = time.time()

		self.steeringAngle = self.error * self.kX + self.phi * self.kPhi + \
		(self.error - self.lastError) * self.kD / (self.currTime - self.lastTime) + \
		(self.phi - self.lastPhi) * self.kDPhi / (self.currTime - self.lastTime)
		self.lastTime = self.currTime

		if(self.steeringAngle > 0.3):
			self.steeringAngle = 0.3
		elif(self.steeringAngle < -0.3):
			self.steeringAngle = -0.3



class position:
	def __init__(self, l, d, vel):
		self.l = l
		self.d = d
		self.vel = vel
		self.xVal = 0
		self.yVal = 0
		self.phiVal = 0
		self.lastTime = time.time()
		self.currTime = time.time()
		self.a = 6378137
		self.b = mpf(6356752.3142)
		self.initLat = 42.03
		self.initLong = -93.6
		self.initAlt = 250
		self.esq = 1 - ((self.b**2)/(self.a**2))
		self.ePrime2 = ((self.a**2)/(self.b**2)) - 1
		self.initX,self.initY,self.initZ = self.getECEF(self.initLat,self.initLong,self.initAlt)
		self.X = self.initX
		self.Y = self.initY
		self.Z = self.initZ
		self.normBasis = self.getEllipseNormal(self.initX,self.initY,self.initZ)
		self.planeNorm = self.getPlaneNormal(self.initX,self.initY,self.initZ)
		self.transform = [[0,0,0],[0,0,0],[0,0,0]]
		self.getBasis()

	def reset(self):
		self.xVal = 0
		self.phiVal = 0

	def setSpeed(self, speed):
		self.vel = speed

	def setX(self, x):
		self.xVal = x

	def integrate(self, dxDt, dyDt, dPhiDt):
		self.currTime = time.time()
		self.xVal = self.xVal + (self.currTime - self.lastTime) * dxDt
		self.phiVal = self.phiVal + (self.currTime - self.lastTime) * dPhiDt	
		self.yVal = self.yVal + (self.currTime - self.lastTime) * dyDt
		self.lastTime = self.currTime

	def runKinematics(self, thetaInput, shouldRun):
		if(shouldRun):
			dxDt = self.vel * self.d * np.tan(thetaInput) * np.cos(self.phiVal) / self.l - \
			np.sin(self.phiVal)*self.vel	
			dyDt = self.vel * self.d * np.tan(thetaInput) * np.sin(self.phiVal) / self.l + \
			np.cos(self.phiVal)*self.vel	
			dPhiDt = -self.vel * np.tan(thetaInput)/self.l
		else:
			dxDt = 0
			dyDt = 0
			dPhiDt = 0
		self.integrate(dxDt, dyDt, dPhiDt)

	def getECEF(self,lat,longi,height):
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

	def getEllipseNormal(self,x,y,z):
		t = mpf(math.sqrt((self.a**2 * self.b**2)/(x**2 * self.b**2 + y**2 * self.b**2 + z**2 * self.a*2)))
		xEllipse = x * t
		yEllipse = y * t
		zEllipse = z * t
		n1 = mpf(2 * xEllipse / (self.a**2))
		n2 = mpf(2 * yEllipse / (self.a**2))
		n3 = mpf(2 * zEllipse / (self.b**2))
		return np.divide([n1,n2,n3],mpf(sqrt(n1**2 + n2**2 + n3**2)))

	def getPlaneNormal(self,x,y,z):
		v1 = [0,0,100000]
		v2 = [x,y,z]
		norm = np.cross(v1, v2)
		return norm

	def getBasis(self):
		zBasis = np.divide(self.normBasis,mpf(sqrt(self.normBasis[0]**2+self.normBasis[1]**2+self.normBasis[2]**2)))
		xBasis = np.divide(self.planeNorm,mpf(sqrt(self.planeNorm[0]**2+self.planeNorm[1]**2+self.planeNorm[2]**2)))
		yBasis = np.cross(zBasis,xBasis)
		self.transform = [[xBasis[0],yBasis[0],zBasis[0]],[xBasis[1],yBasis[1],\
		zBasis[1]],[xBasis[2],yBasis[2],zBasis[2]]]

	def updateXYZ(self):
		worldVector = np.matmul(self.transform,[[self.xVal],[self.yVal],[0]]) + [[self.initX],[self.initY],[self.initZ]]
		self.X = worldVector[0][0]
		self.Y = worldVector[1][0]
		self.Z = worldVector[2][0]

	def convertToLatLongAlt(self):
		self.updateXYZ()
		p = sqrt(self.X**2 + self.Y**2)
		F = 54*(self.b**2)*(self.Z**2)
		G = p**2 + (1-self.esq)*self.Z**2 - self.esq*(self.a**2-self.b**2)
		c = self.esq**2 * F * p**2 / (G**3)
		s = (1 + c + sqrt(c**2 + 2*c))**(1/3)
		k = s + 1 + 1/s
		P = F/(3 * k**2 * G**2)
		Q = sqrt(1 + 2 * self.esq**2 * P)
		r0 = ((-1*P*self.esq*p)/(1+Q)) + sqrt(0.5 * self.a**2 * (1 + 1/Q) - ((P * (1-self.esq) * (self.Z**2))/(Q*(1+Q))) - 0.5*P*p**2)
		U = sqrt((p - self.esq*r0)**2 + self.Z**2)
		V = sqrt((p - self.esq*r0)**2 + (1-self.esq)*(self.Z**2))
		z0 = self.b**2 * self.Z / (self.a * V)
		h = U * (1 - (self.b**2 / (self.a * V)))
		lat = mp.atan((self.Z + self.ePrime2*z0)/p)
		longi = mp.atan2(self.Y, self.X)
		return mp.degrees(lat),mp.degrees(longi),h

tractorPos = position(1.72, 0.4, 10)
#tractorPos = position(1.72, 0.4, 0.5)


tractorController = controller(-0.12, 0.6, -0.05, 0.1)


class Ui_MainWindow(QMainWindow):
	def __init__(self):
		super().__init__()
		self.setupUi(self)

	def setupUi(self, MainWindow):
		MainWindow.setObjectName("MainWindow")
		MainWindow.resize(800, 600)
		self.centralwidget = QWidget(MainWindow)
		self.centralwidget.setObjectName("centralwidget")
		self.ipAddress = QtWidgets.QPlainTextEdit(self.centralwidget)
		self.ipAddress.setGeometry(QtCore.QRect(10, 10, 131, 31))
		self.ipAddress.setObjectName("ipAddress")
		self.connectButton = QtWidgets.QPushButton(self.centralwidget)
		self.connectButton.setGeometry(QtCore.QRect(30, 50, 93, 28))
		self.connectButton.setObjectName("connectButton")
		self.runButton = QtWidgets.QPushButton(self.centralwidget)
		self.runButton.setGeometry(QtCore.QRect(10, 490, 171, 41))
		self.runButton.setObjectName("runButton")
		self.resetButton = QtWidgets.QPushButton(self.centralwidget)
		self.resetButton.setGeometry(QtCore.QRect(10, 437, 171, 41))
		self.resetButton.setObjectName("resetButton")
		self.xValue = QtWidgets.QPlainTextEdit(self.centralwidget)
		self.xValue.setGeometry(QtCore.QRect(10, 110, 131, 31))
		self.xValue.setObjectName("xValue")
		self.setXButton = QtWidgets.QPushButton(self.centralwidget)
		self.setXButton.setGeometry(QtCore.QRect(30, 150, 93, 28))
		self.setXButton.setObjectName("setXButton")
		MainWindow.setCentralWidget(self.centralwidget)
		self.menubar = QtWidgets.QMenuBar(MainWindow)
		self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 26))
		self.menubar.setObjectName("menubar")
		MainWindow.setMenuBar(self.menubar)
		self.statusbar = QtWidgets.QStatusBar(MainWindow)
		self.statusbar.setObjectName("statusbar")
		MainWindow.setStatusBar(self.statusbar)

		self.retranslateUi(MainWindow)
		QtCore.QMetaObject.connectSlotsByName(MainWindow)

		self.setCallbacks()

		self.show()

	def setCallbacks(self):
		self.connectButton.clicked.connect(self.connectButtonClicked)
		self.setXButton.clicked.connect(self.setXButtonClicked)
		self.runButton.clicked.connect(self.runButtonClicked)
		self.resetButton.clicked.connect(self.resetButtonClicked)

	def connectButtonClicked(self):
		global connect
		global connectAddr
		data_lock.acquire()
		connect = 1
		connectAddr = self.ipAddress.toPlainText()
		data_lock.release()

	def setXButtonClicked(self):
		global setXVal
		data_lock.acquire()
		setXVal = float(self.xValue.toPlainText())
		data_lock.release()

	def runButtonClicked(self):
		global runKinematics
		data_lock.acquire()
		runKinematics = not(runKinematics)
		data_lock.release()

	def resetButtonClicked(self):
		global resetFlag
		data_lock.acquire()
		resetFlag = True
		data_lock.release()


	def retranslateUi(self, MainWindow):
		_translate = QtCore.QCoreApplication.translate
		MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
		self.ipAddress.setPlainText(_translate("MainWindow", "192.168.86.39"))
		self.connectButton.setText(_translate("MainWindow", "Connect"))
		self.runButton.setText(_translate("MainWindow", "Run"))
		self.xValue.setPlainText(_translate("MainWindow", "0"))
		self.setXButton.setText(_translate("MainWindow", "Set X"))
		self.resetButton.setText(_translate("MainWindow", "Reset"))

	def runMainTask(self):
		self.thread = QThread()
		self.worker = Worker()
		self.worker.moveToThread(self.thread)

		self.thread.started.connect(self.worker.runMainLoop)
		#self.worker.finished.connect(self.thread.quit)
		#self.worker.finished.connect(self.worker.deleteLater)
		#self.thread.finished.connect(self.thread.deleteLater)

		self.thread.start()

	def runNetworkingTask(self):
		self.networkThread = QThread()
		self.networkWorker = networkWorker()
		self.networkWorker.moveToThread(self.networkThread)

		self.networkThread.started.connect(self.networkWorker.runMainLoop)

		self.networkThread.start()

class networkWorker(QObject):
	def runMainLoop(self):
		global wheelAngle
		wheelAngleSock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
		myAddr = ('192.168.86.37',5005)
		wheelAngleSock.bind(myAddr)

		while(True):
			wheelAngleBytes = wheelAngleSock.recvfrom(15)
			wheelAngleStr = wheelAngleBytes[0].decode('latin-1').strip()
			if(not(wheelAngleStr == '')):
				last = wheelAngle
				network_lock.acquire()
				try:
					wheelAngle = np.radians(float(wheelAngleStr))
				except:
					wheelAngle = last
				print(wheelAngle)
				network_lock.release()


class Worker(QObject):

	def runMainLoop(self):
		self.lastTimeMain = time.time()
		global runKinematics
		global setXVal
		global connect
		global connectAddr
		global wheelAngle
		global resetFlag
		data_lock.acquire()
		lastSetX = setXVal
		localConnectAddr = connectAddr
		data_lock.release()
		outSock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
		serverAddress = (localConnectAddr,9999)

		while(True):
			data_lock.acquire()
			localRunKin = runKinematics
			localSetXVal = setXVal
			localConnect = connect
			localConnectAddr = connectAddr
			localResetFlag = resetFlag
			resetFlag = False
			data_lock.release()

			network_lock.acquire()
			if(localResetFlag):
				wheelAngle = 0
				localWheelAngle = 0
			else:
				localWheelAngle = wheelAngle
			network_lock.release()
			#print('wheelAngle')
			#print(localWheelAngle)

			if(not(localSetXVal == lastSetX)):
				tractorPos.setX(localSetXVal)

			if(localConnect):
				localConnect = not(localConnect)
				serverAddress = (localConnectAddr,9999)

			if(localResetFlag):
				tractorPos.reset()
				localResetFlag = False
			tractorPos.runKinematics(wheelAngle,localRunKin)
			curLat,curLongi,curAlt = tractorPos.convertToLatLongAlt()
			self.currTimeMain = time.time()
			if self.currTimeMain-self.lastTimeMain > 0.1:

				#exportDesiredWheelAngle(tractorController.steeringAngle*18000/3.1415)
				outStr = "$GNGGA,-," + str(curLat) + ",N," + str(abs(curLongi)) + ",W,4," + "0,0," + str(curAlt) + ",0,0,0,0,$"
				numBytes = outSock.sendto(bytes(outStr,'utf-8'), serverAddress)
				#print("sent " + str(numBytes) + " bytes.")
				#print(outStr)
				#print('x')
				#print(tractorPos.X)
				#print('y')
				#print(tractorPos.Y)
				#print('z')
				#print(tractorPos.Z)
				self.lastTimeMain = self.currTimeMain
			lastSetX = localSetXVal



App = QApplication(sys.argv)
window = Ui_MainWindow()
window.setupUi(window)
window.runMainTask()
window.runNetworkingTask()

keepRunning = 1
while(keepRunning):
	App.processEvents()

App.quit()

