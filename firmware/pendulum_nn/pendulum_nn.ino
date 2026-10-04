#include "AccelStepper.h"
#include "nn_weights.h"

// Pin definitions
#define DIR_PIN 25
#define PUL_PIN 26
#define ENC_X_A_PIN 27
#define ENC_X_B_PIN 14
#define ENC_TH_A_PIN 16
#define ENC_TH_B_PIN 17
#define SWITCH_PIN 0
#define POTMET_PIN 34

// Constants
const double PULLEY_DIAMETER_MM = 37.8;
const double PULLEY_LENGTH_MM = (PULLEY_DIAMETER_MM / 2) * PI * 2;
const double NORMAL_STEPS_PER_REV = 200;
const double MICROSTEPPING_RES = NORMAL_STEPS_PER_REV * 16;
const double MM_PER_STEP = PULLEY_LENGTH_MM / MICROSTEPPING_RES;
const double MM_PER_ENC = PULLEY_LENGTH_MM / 600;
const double RAD_PER_ENC = 0.00306796157;
const double DT_SEC = 0.02;

// Encoder counts
volatile int countEncTheta = 0;
volatile int countEncPos = 0;

// Stepper motor
AccelStepper stepper(AccelStepper::DRIVER, PUL_PIN, DIR_PIN);

// State variables
int printFlag = 1;
int state = 0;
unsigned int currentTimeMicros = 0;
bool switchIsPressed = false;

// Potentiometer variables
int rawPotentiometerValue = 0;
double potentiometerValue = 0;

// Wait counter
unsigned int waitCounter = 0;

// Low-pass filter class
class LowPassFilter {
public:
    LowPassFilter(double alpha, double beta)
        : alpha(alpha), beta(beta), yk_0(0), yk_1(0), uk_1(0) {}

    double applyFilter(double uk_0) {
        yk_0 = alpha * uk_1 + beta * yk_1;
        yk_1 = yk_0;
        uk_1 = uk_0;
        return yk_0;
    }
private:
    double alpha, beta;
    double yk_0, yk_1, uk_1;
};

LowPassFilter angularVelocityFilter(0.8647, 0.1353);
LowPassFilter cartVelocityFilter(0.8647, 0.1353);

// Stare citita O SINGURA DATA per ciclu
double x1_m    = 0;   // pozitie [m]
double x2_ms   = 0;   // viteza carucior [m/s]
double x3_rad  = 0;   // unghi [rad]
double x4_rads = 0;   // viteza unghiulara [rad/s]

double command_ms2 = 0;   // comanda NN [m/s^2]

double prevAngle = 0;
double prevPos   = 0;

// ================================================================
// Forward pass NN
// ================================================================
static double _h0[NN_M];
static double _h1[NN_M];

double nn_evaluate(double x1, double x2, double x3, double x4) {
    double input[NN_NX] = {x1, x2, x3, x4};

    double *h_in  = _h0;
    double *h_out = _h1;

    // Normalizare intrare
    for (int j = 0; j < NN_NX; j++) {
        h_in[j] = (input[j] - NN_X_MEAN[j]) / NN_X_STD[j];
    }

    // Straturi ascunse cu ReLU
    for (int li = 0; li < NN_N_FC - 1; li++) {
        int rows = NN_ROWS[li];
        int cols = NN_COLS[li];
        const double *W = NN_W[li];
        const double *b = NN_B[li];

        for (int i = 0; i < rows; i++) {
            double sum = b[i];
            for (int j = 0; j < cols; j++) {
                sum += W[i * cols + j] * h_in[j];
            }
            h_out[i] = sum > 0.0 ? sum : 0.0;
        }
        double *tmp = h_in; h_in = h_out; h_out = tmp;
    }

    // Strat de iesire (fara ReLU)
    {
        int li   = NN_N_FC - 1;
        int cols = NN_COLS[li];
        const double *W = NN_W[li];
        const double *b = NN_B[li];
        double sum = b[0];
        for (int j = 0; j < cols; j++) {
            sum += W[j] * h_in[j];
        }
        // Denormalizare iesire
        return NN_U_MEAN + NN_U_STD * sum;
    }
}

// ================================================================
// Function prototypes
// ================================================================
int    cmToSteps(double l_cm);
double getPosInCm();
double getAngleInRad_raw();
void   readAllStates();
void   printAllStates();
void   waitAndGoToNextStep(int nextState);
float  mapFloat(float value, float fromLow, float fromHigh,
                float toLow, float toHigh);
double readPotentiometer();
void   updateStepperSpeed(void *parameter);

// ================================================================
// ISR
// ================================================================
void IRAM_ATTR isrEncPos() {
    if (digitalRead(ENC_X_B_PIN)) countEncPos++;
    else                           countEncPos--;
}

void IRAM_ATTR isrEncTheta() {
    if (digitalRead(ENC_TH_B_PIN)) countEncTheta++;
    else                            countEncTheta--;
}

// ================================================================
// Setup
// ================================================================
void setup() {
    Serial.begin(115200);

    pinMode(DIR_PIN, OUTPUT);
    pinMode(PUL_PIN, OUTPUT);
    pinMode(POTMET_PIN, INPUT);
    pinMode(SWITCH_PIN, INPUT_PULLUP);

    pinMode(ENC_TH_A_PIN, INPUT);
    pinMode(ENC_TH_B_PIN, INPUT);
    attachInterrupt(digitalPinToInterrupt(ENC_TH_A_PIN), isrEncTheta, RISING);

    pinMode(ENC_X_A_PIN, INPUT);
    pinMode(ENC_X_B_PIN, INPUT);
    attachInterrupt(digitalPinToInterrupt(ENC_X_A_PIN), isrEncPos, RISING);

    stepper.setMaxSpeed(cmToSteps(70));
    stepper.setAcceleration(cmToSteps(20));

    countEncTheta = 0;
    countEncPos   = 0;

    delay(1000);

    xTaskCreatePinnedToCore(updateStepperSpeed, "updateStepperSpeed",
                            5000, NULL, 1, NULL, 0);
}

// ================================================================
// Loop
// ================================================================
void loop() {
    if (micros() - currentTimeMicros > 20000) {
        currentTimeMicros = micros();

        switchIsPressed = !digitalRead(SWITCH_PIN);
        readPotentiometer();

        // Citire stare O SINGURA DATA per ciclu
        readAllStates();

        switch (state) {

            case 0:
                stepper.setSpeed(cmToSteps(-10));
                if (switchIsPressed) {
                    stepper.setSpeed(0);
                    state = 1;
                    countEncPos = 0;
                }
                break;

            case 1:
                if (getPosInCm() < 20) {
                    stepper.setSpeed(cmToSteps(30));
                } else {
                    stepper.setSpeed(0);
                    waitAndGoToNextStep(2);
                }
                break;

            case 2:
                countEncPos = 0;
                if (abs(x3_rad) < 0.1) {
                    state = 3;
                }
                break;

            case 3: {
                double acc_ms2  = nn_evaluate(x1_m, x2_ms, x3_rad, x4_rads);
                command_ms2     = acc_ms2;
                double acc_cms2 = acc_ms2 * 100.0;
                int newSpeed    = stepper.speed() + cmToSteps(acc_cms2 * DT_SEC);
                stepper.setSpeed(newSpeed);

                // Siguranta
                if (abs(x3_rad) > 0.5 ||
                    x1_m < -0.17 || x1_m > 0.17) {
                    stepper.setSpeed(0);
                    state = 0;
                }
                break;
            }
        }

        if (printFlag) {
            printAllStates();
        }
    }
}

// ================================================================
// Citire stare completa — apelata O SINGURA DATA per ciclu
// ================================================================
void readAllStates() {
    double posInCm = MM_PER_ENC * countEncPos / 10.0;
    double angle   = countEncTheta * RAD_PER_ENC + potentiometerValue - 3.14;

    double rawCartVel = (posInCm - prevPos)   / DT_SEC;
    double rawAngVel  = (angle   - prevAngle) / DT_SEC;

    prevPos   = posInCm;
    prevAngle = angle;

    x1_m    = posInCm * 0.01;
    x2_ms   = cartVelocityFilter.applyFilter(rawCartVel) * 0.01;
    x3_rad  = angle;
    x4_rads = angularVelocityFilter.applyFilter(rawAngVel);
}

// ================================================================
// Functii auxiliare
// ================================================================
double getPosInCm() {
    return MM_PER_ENC * countEncPos / 10.0;
}

double readPotentiometer() {
    rawPotentiometerValue = analogRead(POTMET_PIN);
    potentiometerValue = mapFloat(rawPotentiometerValue, 0, 4096, -0.1, 0.1);
    return potentiometerValue;
}

float mapFloat(float value, float fromLow, float fromHigh,
               float toLow, float toHigh) {
    return (value - fromLow) * (toHigh - toLow) /
           (fromHigh - fromLow) + toLow;
}

int cmToSteps(double l_cm) {
    return -static_cast<int>((l_cm * 10) / MM_PER_STEP);
}

void waitAndGoToNextStep(int nextState) {
    waitCounter++;
    if (waitCounter == 25) {
        state = nextState;
        waitCounter = 0;
    }
}

void printAllStates() {
    Serial.print(x1_m);        Serial.print("\t");
    Serial.print(x2_ms);       Serial.print("\t");
    Serial.print(x3_rad);      Serial.print("\t");
    Serial.print(x4_rads);     Serial.print("\t");
    Serial.println(command_ms2);
}

void updateStepperSpeed(void *parameter) {
    disableCore0WDT();
    while (true) {
        stepper.runSpeed();
    }
}