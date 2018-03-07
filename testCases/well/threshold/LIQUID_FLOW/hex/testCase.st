#SOURCE_TERM
 $PCS_TYPE
  LIQUID_FLOW
 $PRIMARY_VARIABLE
  PRESSURE1
 $GEO_TYPE
  POLYLINE PLY_WARMWELL
 $THRESHOLD_GEOMETRY
  POINT POINT13
 $DIS_TYPE
  CONSTANT_NEUMANN 1e-4 ; 6.25E-03 ;[m^3/s] 19/20 von 30 m^2/h (da 20 Layer und 2 Grenzpunkte, die zusammen einen Punkt ausmachen)
 $THRESHOLD
  1 HEAT_TRANSPORT 350 1 2 2 ; type [1: lower, 2: upper], process name, threshold value, scheme [0: explicit, 1: implicit], delta [only if scheme is implicit], verbosity [0, 1, 2]
 $FCT_TYPE
  IN_OUT
 
#STOP
