#include <ICM_20948.h>
#include <Wire.h>

#define AD0_VAL 1

ICM_20948_I2C ICM;

#define dps500ScaleFactor 0.0152671755

float offset;

float calibrate() {
  Wire.flush();
  float calibrationOffset = 0;
  int calibrationTrials = 120;
  for(int i = 0; i < calibrationTrials; i++) {
    while(!ICM.dataReady()){}
    ICM.getAGMT();
    float xRate = ICM.agmt.gyr.axes.x / 131;
    Serial.print(xRate);
    if(i > 19){
      calibrationOffset += xRate;
    }
    
  }
  Serial.print("Calibration Offset: ");
  Serial.println(calibrationOffset);
  return calibrationOffset/(calibrationTrials - 20);
}

void setup() {
  // put your setup code here, to run once:
  Serial.begin(9600);
  while(!Serial)
  {};

  Wire.begin();
  Wire.setClock(400000);
  ICM.begin(Wire, AD0_VAL);

  Serial.print("Initalization of the sensor returned: ");
  Serial.println(ICM.statusString());

  bool initialized = false;
  
  while(ICM.status != ICM_20948_Stat_Ok)
  {
    Serial.println("Trying again...");
    delay(500);
  }
  
  initialized = true;

  ICM_20948_fss_t fss;

  offset = calibrate();
  /*fss.a = gpm2;
  fss.g = dps500;
  
  ICM.setFullScale(ICM_20948_Internal_Gyr, fss);
  */

}

/*typedef struct sensorData{
  
}*/
unsigned long StartTime;
unsigned long EndTime;
bool first = true;
float rotation = 0;
float xRate = 0;
long xAccel = 0;
long yAccel = 0;
long zAccel = 0;
float xRateFiltered = 0;

const int filterBufferSize = 6;
float filterBuffer[filterBufferSize];
int bufferIndex = 0;
float meanSum = 0;

// This function returns the oldest index and wraps
// around at the end of the buffer
int getOldestIndex(int index, int bufferSize) {
  if(index+1 >= bufferSize) {
    return 0;
  } else {
    return index + 1;
  }
}

void loop() {
  if(ICM.dataReady()){
    ICM.getAGMT();
    long accelDivisorSize = 1000000000000;
    xRate = (ICM.agmt.gyr.axes.x) / 131;
    xAccel = int(ICM.agmt.acc.axes.x);
    long yAccelPre = ICM.agmt.acc.axes.y;
    yAccel = long(yAccelPre * 16384 * 981);
    zAccel = int(ICM.agmt.acc.axes.z * 16384 * 981);
    Serial.print(xAccel);
    Serial.print(" ");
    Serial.print(yAccelPre);
    Serial.print(" ");
    Serial.print(float(yAccel/accelDivisorSize));
    Serial.print(" ");
    Serial.print(zAccel);
    Serial.print(" ");
    float totalAccel = sqrt(pow(xAccel, 2) + pow(yAccel, 2) + pow(zAccel, 2));

    // On initialization, just use the first value
    // for everything in the buffer
    if(first) {
      Serial.print("Offset: ");
      Serial.println(offset);
      Wire.flush();
      for(int i = 0; i < filterBufferSize; i++){
        filterBuffer[i] = xRate;
      }
      meanSum = xRate * filterBufferSize;
    }

    int oldestVal = filterBuffer[bufferIndex];
    filterBuffer[bufferIndex] = xRate;
 
    meanSum = meanSum - oldestVal;
    meanSum = meanSum + xRate;
    bufferIndex++;
    
    // Wrap the index around at max array size
    if(bufferIndex >= filterBufferSize) {
      bufferIndex = 0;
    }

    xRateFiltered = (meanSum / filterBufferSize) - offset;

    if(!first){
      EndTime = millis();
      if(abs(xRateFiltered) > 1){
        rotation += xRateFiltered * (EndTime - StartTime) * 0.001;
      }
      StartTime = EndTime;
    }
    first = false;

    Serial.print(millis());
    Serial.print(",");
    Serial.println(totalAccel);
  }
}
