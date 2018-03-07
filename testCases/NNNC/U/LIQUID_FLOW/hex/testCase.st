#SOURCE_TERM
 $PCS_TYPE
  LIQUID_FLOW
 $PRIMARY_VARIABLE
  PRESSURE1
 $GEO_TYPE
  SURFACE SURF_TOP_LEFT
 $DIS_TYPE
  CONSTANT_NEUMANN -.1
 $FCT_TYPE
  Time_ST_curve


#SOURCE_TERM
 $PCS_TYPE
  LIQUID_FLOW
 $PRIMARY_VARIABLE 
  PRESSURE1
 $GEO_TYPE
  SURFACE SURF_BOTTOM_LEFT
 $DIS_TYPE
  CONSTANT_NEUMANN      1.0
 $CONNECTED_GEOMETRY
  SURFACE SURF_BOTTOM_RIGHT
 $CONNECT_PARAMETERS 
  1e3 1  ; exchange coefficient   verbosity level = 0, 1
 $CONNECT_MODE
  0 ; mode  ; 0: symmetric, 1: non-symmetric (downwind fixed), 2 variable (dependent on velocity in reference element)
  ;  2  ref_element_number   n_ref_x, n_ref_y, n_ref_z minimum_velocity_abs 


#SOURCE_TERM
 $PCS_TYPE
  HEAT_TRANSPORT
 $PRIMARY_VARIABLE 
  TEMPERATURE1
 $GEO_TYPE
  SURFACE SURF_BOTTOM_LEFT
 $DIS_TYPE
  CONSTANT_NEUMANN      1.0
 $CONNECTED_GEOMETRY
  SURFACE SURF_BOTTOM_RIGHT
 $CONNECT_PARAMETERS 
  1e6 1  ; exchange coefficient   verbosity level = 0, 1
 $CONNECT_MODE
  2 226 0 0 1 1.e-10
  ; mode  ; 0: symmetric, 1: non-symmetric (downwind fixed), 2 variable (dependent on velocity in reference element)
  ;  2  ref_element_number   n_ref_x, n_ref_y, n_ref_z minimum_velocity_abs 

#STOP

$EVERYONE_WITH_EVERYONE
#STOP








