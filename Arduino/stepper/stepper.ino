#include <math.h>

#define dirPin 12
#define stepPin 13
#define pulsePerRev 200
#define accel 1.5 // rev per second^2 0.5
#define accelPulse (accel*pulsePerRev)
#define clockSpeed 16000000
#define prescalar 1024
#define periodScalar 10000
#define startupRPS 0.2
#define secPerCount 0.000064
#define cruiseSpeed 3 // 3rps
#define degPerPulse 1.8

float targetPeriod;

float finalPeriod;

float velocity = 0;
float lastVelocity = 0;
float rpsModify = 0;

float targetCyclesPerPulse;
unsigned int targetCount;
unsigned int targetCountMultiplier;

unsigned long finalCyclesPerPulse;
unsigned int finalCount;
unsigned int finalCountMultiplier;

unsigned int deltaPeriod;


char currentMultiplier = 0;

int startRPM = 0;

bool inInterrupt = 0;
bool motorOn = 0;
bool accelDirection = 1; // 1 = positive(CCW) 0 = negative(CW)
bool velocityDirection = 1; // 1 = positive(CCW) 0 = negative(CW)

float lastRPS = 0;

float currentPos = 0;
int targetPos = 0;
int lastTargetPos = 0;

int currentPulsePerSec = 0;

int leadPulse = 0;

int posState = -1;

bool stepPinState = 0;
bool dirPinState = 0;

const byte numInputChars = 32;
char receivedChars[numInputChars];
char endMarker = '\n';
byte inputPos = 0;

float currentSteeringAngle = 0;
float desiredSteeringAngle = 0;

unsigned long lastTime = 0;
int setSpeedPeriodMillis = 50;

void setup() {
  // put your setup code here, to run once:
  pinMode(dirPin, OUTPUT);
  pinMode(stepPin, OUTPUT);

  digitalWrite(dirPin, dirPinState);
  
  
  // set TCCR1A to 00000000
  TCCR1A = 0;
  
  // set OC1A to desired number (targetCount)
  OCR1A = 65535;
  
  // set TCNT1 to 0
  TCNT1 = 0;
  
  // set TCCR1B to 00000101
  TCCR1B = 5; 

  // set TIMSK1 to 00000010
  TIMSK1 = 2;

  Serial.begin(9600);

  sei();
}

int sign(float input) {
  if(input >= 0) {
    return 1;
  } else {
    return -1;
  }
}

void loop() {

  //Verify the current steering angle and the desired steering angle refer
  //to the same thing - either motor or axle angle.

  currentSteeringAngle = (currentSteeringAngle + (230-analogRead(0))*27.5/79)/2;
  
  if(Serial.available() > 0) {
    //targetPos = Serial.parseFloat(SKIP_ALL);
    char rc = Serial.read();

    if (rc != endMarker) {
      receivedChars[inputPos] = rc;
      inputPos++;
      if(inputPos >= numInputChars) {
        inputPos = numInputChars - 1;
      }
    }
    else {
      receivedChars[inputPos] = '\0';
      inputPos = 0;
      ///targetPos = atoi(receivedChars);
      desiredSteeringAngle = atof(receivedChars);
    }
  }

  //if(millis() - lastTime > setSpeedPeriodMillis){
    //setSpeedRPS(desiredSteeringAngle);
    float error = (desiredSteeringAngle-currentSteeringAngle);
    /*Serial.println(error);
    Serial.print("Motor On: ");
    Serial.println(motorOn);
    Serial.print("Final cycles per pulse: ");
    Serial.println(finalCyclesPerPulse);
    Serial.print("Accel Direction: ");
    Serial.println(accelDirection);
    Serial.print("Target cycles per pulse: ");
    Serial.println(targetCyclesPerPulse);*/
    if(abs(error) > 0.5) {
      setSpeedRPS(error/10);
    } else {
      setSpeedRPS(0);
    }
    /*
    Serial.print("OCR1A: ");
    Serial.println(OCR1A);*/
    //Serial.print("setSpeedRPS");
    //Serial.println((desiredSteeringAngle-currentSteeringAngle)/3);
    //lastTime = millis();
  //}

  //Serial.println((desiredSteeringAngle-currentSteeringAngle)/3);
  manageDirection();
  //managePos();
}

void managePos() {
  if(lastTargetPos != targetPos) {
    posState = 0;
  }
  
  switch (posState) {
    case 0:
      if(currentPos < targetPos){
        setSpeedRPS(cruiseSpeed);
      } else if(currentPos > targetPos) {
        setSpeedRPS(-1*cruiseSpeed);
        //Serial.println("SET SPEED");
        //Serial.print("velocity: ");
        //Serial.println(velocity);
      }
      posState++;
    case 1:
      if((currentPos < targetPos) && dirPinState == 0) {
        currentPulsePerSec = 1/(OCR1A * 2 * secPerCount);
        leadPulse = pow(currentPulsePerSec, 2) / (2 * accelPulse);
        //Serial.print("leadPulse: ");
        //Serial.println(leadPulse);
        if(currentPos > targetPos - (degPerPulse * leadPulse + 36)) {
          setSpeedRPS(0.3);
          posState++;
        }
      } else if((currentPos > targetPos) && dirPinState == 1) {
        currentPulsePerSec = 1/(OCR1A * 2 * secPerCount);
        leadPulse = pow(currentPulsePerSec, 2) / (2 * accelPulse);
        //Serial.print("leadPulse: ");
        //Serial.println(leadPulse);
        if(currentPos < targetPos + (degPerPulse * leadPulse + 36)) {
          setSpeedRPS(-0.3);
          posState++;
        }
      }
  }
  lastTargetPos = targetPos;
}

void manageDirection() {
  if((velocity < 0 && !digitalRead(dirPin)) || (velocity > 0 && digitalRead(dirPin))) {
    if(targetCyclesPerPulse > 194) {
      //Serial.println("SWITCHING DIRECTION!!");
      dirPinState = !dirPinState;
      digitalWrite(dirPin, dirPinState);
      lastVelocity = -1 * velocity;
      setSpeedRPS(velocity);
    }
  }
}

void setSpeedRPS(float rps) {
  velocity = rps;
  if(!(abs(lastVelocity - velocity) < 0.09)){
    //digitalRead(dirPin) == 1 -> CCW rotation
    if(motorOn && ((rps < -0.01 && !digitalRead(dirPin)) || (rps > 0.01 && digitalRead(dirPin)))) {
      rpsModify = 0.1;
    } else {
      rpsModify = abs(rps);
    }
    accelDirection = !(abs(rpsModify) > abs(lastRPS));
   
    if(rpsModify <= 0.01) {
      OCR1A = 300;
      motorOn = 0;
      lastRPS = 0;
    } else {
      float currRPS = 1/(secPerCount*targetCyclesPerPulse*2*pulsePerRev);
      if(currRPS < 0.1 || !motorOn) {
        TCNT1 = 0;
        dirPinState = (velocity < 0);
        digitalWrite(dirPin, dirPinState);
        finalPeriod = periodScalar/(pulsePerRev * rpsModify);
        finalCyclesPerPulse = (finalPeriod/2) * clockSpeed / prescalar / periodScalar;
        finalCount = finalCyclesPerPulse % 65535;
        finalCountMultiplier = (int)(finalCyclesPerPulse/65535);
        targetPeriod = periodScalar/(pulsePerRev * startupRPS); 
        targetCyclesPerPulse = (targetPeriod/2) * clockSpeed / prescalar / periodScalar;
        targetCount = (int)targetCyclesPerPulse % 65535;
        targetCountMultiplier = (int)(targetCyclesPerPulse/65535);
        OCR1A = targetCount;  
      } else {
        finalPeriod = periodScalar/(pulsePerRev * rpsModify);
        finalCyclesPerPulse = (finalPeriod/2) * clockSpeed / prescalar / periodScalar;
        finalCount = finalCyclesPerPulse % 65535;
        finalCountMultiplier = (int)(finalCyclesPerPulse/65535);  
      }
      //setPeriod(targetPeriod);
      motorOn = 1;
      lastRPS = rpsModify;
    }
    lastVelocity = velocity;
  }
}

void setPeriod(int period) {  
  targetCyclesPerPulse = (targetPeriod/2) * clockSpeed / prescalar / periodScalar;
  targetCount = (int)targetCyclesPerPulse % 65535;
  targetCountMultiplier = (int)(targetCyclesPerPulse/65535);
  OCR1A = targetCount;

  finalCyclesPerPulse = (finalPeriod/2) * clockSpeed / prescalar / periodScalar;
  finalCount = finalCyclesPerPulse % 65535;
  finalCountMultiplier = (int)(finalCyclesPerPulse/65535);

  //Serial.println("Set Timer!");
  //Serial.println(OCR1A);
}

ISR(TIMER1_COMPA_vect) {
  if(currentMultiplier == targetCountMultiplier) {
    currentMultiplier = 0;
    
    // set TCNT1 to 0
    TCNT1 = 0;

    if(motorOn) {
      /*if((currentPos >= targetPos) && velocity > 0) {
        velocity = 0;
        motorOn = 0;
        lastRPS = 0;
        OCR1A = 10;
      } else if((currentPos <= targetPos) && velocity < 0) {
        velocity = 0;
        motorOn = 0;
        lastRPS = 0;
        OCR1A = 10;
      } else {*/

        stepPinState = (PORTB & ~(B11000000)) >> 5; //digitalRead(stepPin);
        if(!stepPinState) {
          if(dirPinState == 1) {
            currentPos -= 1.8; //CCW
          } else {
            currentPos += 1.8; //CW
          }
        }
        //digitalWrite(stepPin, !stepPinState);
        PORTB ^= B00100000;

        // Accelerate
        if(!accelDirection) {
            if(targetCyclesPerPulse > finalCyclesPerPulse && !stepPinState) {
              
              //targetCyclesPerPulse = (1/ (1/ (secPerCount*2 * targetCyclesPerPulse) + (targetCyclesPerPulse * secPerCount*2 * accelPulse)))/(secPerCount*2);
              //targetCyclesPerPulse = (1/ (1/(secPerCount*targetCyclesPerPulse) + (accelPulse * targetCyclesPerPulse * 2 * secPerCount))) / (secPerCount) * 2;
              //targetCyclesPerPulse = ((1/(targetCyclesPerPulse * secPerCount)) + accelPulse * targetCyclesPerPulse * secPerCount * 2) * 2 / secPerCount;
              //targetCyclesPerPulse = 1/((1/(targetCyclesPerPulse * secPerCount) + (accelPulse * targetCyclesPerPulse * 2)) * secPerCount);
              targetCyclesPerPulse = 0.5/((1/(targetCyclesPerPulse*secPerCount*2) + (accelPulse*secPerCount*targetCyclesPerPulse*2))*secPerCount);
              
              if(targetCyclesPerPulse < finalCyclesPerPulse) {
                targetCyclesPerPulse = finalCyclesPerPulse;
              }
              targetCount = (int)targetCyclesPerPulse % 65535;
              targetCountMultiplier = (int)(targetCyclesPerPulse/65535);
              OCR1A = targetCount;
            }
        } 
        // Decelerate
        else if(accelDirection) {
            if(targetCyclesPerPulse < finalCyclesPerPulse && !stepPinState) {
              
              //targetCyclesPerPulse = (1/ (1/ (secPerCount*2 * targetCyclesPerPulse) - (targetCyclesPerPulse * secPerCount*2 * accelPulse)))/(secPerCount*2);
              //targetCyclesPerPulse = (1/ (1/(secPerCount*targetCyclesPerPulse) - (accelPulse * targetCyclesPerPulse * 2 * secPerCount))) / (secPerCount);
              targetCyclesPerPulse = 0.5/((1/(targetCyclesPerPulse*secPerCount*2) - (accelPulse*secPerCount*targetCyclesPerPulse*2))*secPerCount);
              
              
              if(targetCyclesPerPulse > finalCyclesPerPulse) {
                targetCyclesPerPulse = finalCyclesPerPulse;
              }
              targetCount = (int)targetCyclesPerPulse % 65535;
              targetCountMultiplier = (int)(targetCyclesPerPulse/65535);
              OCR1A = targetCount;
          }
        }
      //}      
    }

    
  } else {
    currentMultiplier++;
  }
}
