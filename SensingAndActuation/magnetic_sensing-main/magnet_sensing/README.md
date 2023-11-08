# Matlab code for magnet sensing

## File Names

#### UnbiasMagnetSensing.m
Input: Filename without magnet, filename with magnet
Output: Excel sheet of Earth's Magnetic Field and Unbias magnet sensing per sensor
Will also output plots

#### SensorDataRead.m
Processes data from lab view, cleans data to use in UnibasMAgnetSensing.m
input nmust be a table with column names: Time, X1, X2, Y1,Y2,Z1, Z2, Sensor, in any order

#### average_XYZ.m
function used in SensorDataRead.m to calculate average values per column

