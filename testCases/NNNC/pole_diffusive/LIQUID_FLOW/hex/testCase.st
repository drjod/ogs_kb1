GeoSys-ST: Source Terms ------------------------------------------------

#SOURCE_TERM
 $PCS_TYPE
  MASS_TRANSPORT
 $PRIMARY_VARIABLE 
  Tracer
 $GEO_TYPE
  SURFACE SURF_3   ;;POLYLINE PLY_2
 $DIS_TYPE
  CONSTANT_NEUMANN      2.0  ; // 2000  
 $FCT_TYPE
  Time_ST_curve

#SOURCE_TERM
 $PCS_TYPE
  HEAT_TRANSPORT
 $PRIMARY_VARIABLE 
  TEMPERATURE1
 $GEO_TYPE
  SURFACE SURF_3   
 $DIS_TYPE
  CONSTANT_NEUMANN      2000.0  ; // 2000  
 $FCT_TYPE
  Time_ST_curve

#SOURCE_TERM
 $PCS_TYPE
  MASS_TRANSPORT
 $PRIMARY_VARIABLE 
  Tracer
 $GEO_TYPE
  SURFACE SURF_1  
 $DIS_TYPE
  CONSTANT_NEUMANN      1.0
; $COPY_VALUE
; POINT P_OUT 0
 $CONNECTED_GEOMETRY
  SURFACE SURF_2  
 $CONNECT_PARAMETERS 
  1000 1  ; exchange coefficient   verbosity level = 0, 1
 $CONNECT_MODE
  0 ; mode  ; 0: symmetric, 1: non-symmetric (downwind fixed), 2 non-symmetric (downwind)
  ;  2  ref_element_number   n_ref_x, n_ref_y, n_ref_z minimum_velocity_abs 
;$EVERYONE_WITH_EVERYONE

#SOURCE_TERM
 $PCS_TYPE
  HEAT_TRANSPORT
 $PRIMARY_VARIABLE 
  TEMPERATURE1
 $GEO_TYPE
  SURFACE SURF_1  
 $DIS_TYPE
  CONSTANT_NEUMANN      1.0
; $COPY_VALUE
; POINT P_OUT 0
 $CONNECTED_GEOMETRY
  SURFACE SURF_2  
 $CONNECT_PARAMETERS 
  1000 1  ; exchange coefficient   verbosity level = 0, 1
 $CONNECT_MODE
  0 ; mode  ; 0: symmetric, 1: non-symmetric (downwind fixed), 2 non-symmetric (downwind)
  ;  2  ref_element_number   n_ref_x, n_ref_y, n_ref_z minimum_velocity_abs 
;$EVERYONE_WITH_EVERYONE
#STOP









#SOURCE_TERM
 $PCS_TYPE
  MASS_TRANSPORT
 $PRIMARY_VARIABLE 
  Tracer
 $GEO_TYPE
  SURFACE SURF_1   ;;POLYLINE PLY_2
 $DIS_TYPE
  CONSTANT_NEUMANN      1.0
; $COPY_VALUE
; POINT P_OUT 0
 $CONNECTED_GEOMETRY
  SURFACE SURF_2  1000 1 0 ; alpha verbosity geomodel


#STOP


#SOURCE_TERM
 $PCS_TYPE
  LIQUID_FLOW
 $PRIMARY_VARIABLE 
  PRESSURE1
 $GEO_TYPE
  SURFACE SURF_2   ;;POLYLINE PLY_2
 $DIS_TYPE
  CONSTANT_NEUMANN      -5e-4  ; // 5e-4
