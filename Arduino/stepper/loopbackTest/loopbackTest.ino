void setup() {
  // put your setup code here, to run once:
  Serial.begin(115200);
}


const byte numInputChars = 32;
char receivedChars[numInputChars];
char endMarker = '\n';
byte inputPos = 0;
int targetPos = 0;

void loop() {
  // put your main code here, to run repeatedly:
  /*if(Serial.available() > 0) {
    char x = Serial.read();
    delay(100);
    Serial.println(x);
  }*/
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
      targetPos = atoi(receivedChars);
      Serial.println(targetPos);
    }  

  }
}
