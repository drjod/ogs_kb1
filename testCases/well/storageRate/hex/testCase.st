
#SOURCE_TERM
 $PCS_TYPE
  LIQUID_FLOW
 $PRIMARY_VARIABLE
  PRESSURE1
 $GEO_TYPE
  POINT POINT24
 $STORAGE_RATE_GEOMETRY
  POINT POINT0
  POINT POINT2
; $THRESHOLD_GEOMETRY
;  POINT POINT11
  $DIS_TYPE
   CONSTANT 1e7
;  $THRESHOLD
;  type [1: lower, 2: upper], process name, threshold value, scheme [0: explicit, 1: implicit], delta [only if scheme is implicit], verbosity [0, 1, 2]
;   2 HEAT_TRANSPORT 310. 1 2. 2
  $STORAGE_RATE
; process name, maximum absolute flow rate, verbosity [0, 1, 2]
   HEAT_TRANSPORT  1.e-1 2
  $TIME_CONTROLLED_ACTIVE
   FUNCTION WARMWELL

#SOURCE_TERM
 $PCS_TYPE
  LIQUID_FLOW
 $PRIMARY_VARIABLE
  PRESSURE1
 $GEO_TYPE
  POINT POINT25
 $STORAGE_RATE_GEOMETRY
  POINT POINT0
  POINT POINT2
; $THRESHOLD_GEOMETRY
;  POINT POINT11
 $DIS_TYPE
  CONSTANT -1e7
; $THRESHOLD
; type [1: lower, 2: upper], process name, threshold value, scheme [0: explicit, 1: implicit], delta [only if scheme is implicit], verbosity [0, 1, 2]
;  2 HEAT_TRANSPORT 310. 1 2. 2
 $STORAGE_RATE
; process name, maximum absolute flow rate, verbosity [0, 1, 2]
  HEAT_TRANSPORT  1.e-1 2
 $TIME_CONTROLLED_ACTIVE
  FUNCTION WARMWELL



#SOURCE_TERM
 $PCS_TYPE
  LIQUID_FLOW
 $PRIMARY_VARIABLE
  PRESSURE1
 $GEO_TYPE
  POINT POINT24
 $STORAGE_RATE_GEOMETRY
  POINT POINT0
  POINT POINT2
; $THRESHOLD_GEOMETRY
;  POINT POINT11
  $DIS_TYPE
   CONSTANT -1e7
;  $THRESHOLD
;  type [1: lower, 2: upper], process name, threshold value, scheme [0: explicit, 1: implicit], delta [only if scheme is implicit], verbosity [0, 1, 2]
;   2 HEAT_TRANSPORT 330. 1 2. 2
  $STORAGE_RATE
; process name, maximum absolute flow rate, verbosity [0, 1, 2]
   HEAT_TRANSPORT  1.e-1 2
  $TIME_CONTROLLED_ACTIVE
   FUNCTION COLDWELL

#SOURCE_TERM
 $PCS_TYPE
  LIQUID_FLOW
 $PRIMARY_VARIABLE
  PRESSURE1
 $GEO_TYPE
  POINT POINT25
 $STORAGE_RATE_GEOMETRY
  POINT POINT0
  POINT POINT2
; $THRESHOLD_GEOMETRY
;  POINT POINT11
 $DIS_TYPE
  CONSTANT 1e7
; $THRESHOLD
; type [1: lower, 2: upper], process name, threshold value, scheme [0: explicit, 1: implicit], delta [only if scheme is implicit], verbosity [0, 1, 2]
;  2 HEAT_TRANSPORT 330. 1 2. 2
 $STORAGE_RATE
; process name, maximum absolute flow rate, verbosity [0, 1, 2]
  HEAT_TRANSPORT  1.e-1 2
 $TIME_CONTROLLED_ACTIVE
  FUNCTION COLDWELL

#STOP
