; $ID:	ATLANTIS_MAIN_V1.PRO,	2021-09-16-17,	USER-KJWH	$
  PRO ATLANTIS_MAIN_V1, SUBAREA_EXTRACT=SUBAREA_EXTRACT, PLOT_DATA=PLOT_DATA, HIRATA_RATIO=HIRATA_RATIO

;+
; NAME:
;   ATLANTIS_MAIN_V1
;
; PURPOSE:
;   The MAIN processing program for creating data inputs for the Atlantis modeling project (2020)
;
; CATEGORY:
;   MAIN
;
; CALLING SEQUENCE:
;   ATLANTIS_MAIN
;
; REQUIRED INPUTS:
;   None
;
; OPTIONAL INPUTS:
;   SUBAREA_EXTRACT....... The SWITCH information to initiate the SUBAREA_EXTRACT block of code
;   PLOT_DATA............. The SWITCH information to initiate the PLOT_DATA block of code
;   HIRATA_RATIO.......... The SWITCH information to initiate the HIRATA_RATIO block of code
;
; KEYWORD PARAMETERS:
;  None
;
; OUTPUTS:
;   This procedure runs various blocks of code for the ATLANTIS project
;
; OPTIONAL OUTPUTS:
;   None
;
; COMMON BLOCKS: 
;   None
;
; SIDE EFFECTS:  
;   None
;
; RESTRICTIONS:  
;   None
;
; EXAMPLE:
; 
;
; NOTES:
;   
;   
; COPYRIGHT: 
; Copyright (C) 2020, Department of Commerce, National Oceanic and Atmospheric Administration, National Marine Fisheries Service,
;   Northeast Fisheries Science Center, Narragansett Laboratory.
;   This software may be used, copied, or redistributed as long as it is not sold and this copyright notice is reproduced on each copy made.
;   This routine is provided AS IS without any express or implied warranties whatsoever.
;
; AUTHOR:
;   This program was written on July 16, 2020 by Kimberly J. W. Hyde, Northeast Fisheries Science Center | NOAA Fisheries | U.S. Department of Commerce, 28 Tarzwell Dr, Narragansett, RI 02882
;    
; MODIFICATION HISTORY:
;   Jul 16, 2020 - KJWH: Initial code written
;   Aug 10, 2020 - KJWH: Added program blocks to run different chunks of code
;                        Added SUBAREA_EXTRACT block to extract the data for the Atlantis subregions
;                        In SUBAREA_EXTRACT, now looping on year to speed up the data extractions
;   Sep 16, 2021 - KJWH: Added NETCDFS block to create D8 and DOY netcdf files
;   Feb 01, 2023 - KJWH: Renamed to ATLANTIS_MAIN_V1 and created a new ATLANTIS_MAIN to work with the new and updated stacked files etc.
;-
; ****************************************************************************************************
  ROUTINE_NAME = 'ATLANTIS_MAIN'
  COMPILE_OPT IDL2
  LUN = []
  SL = PATH_SEP()
  
  PRINT, 'This program is for the first iteration of Atlantis data.  New versions are created in ATLANTIS_MAIN'
  STOP
  
  VERSION = 'VER_2'
  
  DIR_MAIN = !S.ATLANTIS
  DIR_DATA = DIR_MAIN + 'DATA' + SL + VERSION + SL & DIR_TEST, DIR_DATA
  DIR_COMP = DIR_MAIN + 'COMPARE_PLOTS' + SL + VERSION + SL & DIR_TEST, DIR_COMP
  SUBAREA = 'ATLANTIS_NEUS'
  
  OCCCI_VERSION = '4.2'
  
  IF NONE(NETCDFS)         THEN NETCDFS         = ''
  IF NONE(SUBAREA_EXTRACT) THEN SUBAREA_EXTRACT = ''
  IF NONE(PLOT_DATA)       THEN PLOT_DATA       = ''
  IF NONE(HIRATA_RATIO)    THEN HIRATA_RATIO    = ''
  
  
; **************************************************************************************************
  IF KEYWORD_SET(NETCDFS) THEN BEGIN
; **************************************************************************************************
    SNAME = 'NETCDFS'
    SWITCHES,NETCDFS,OVERWRITE=OVERWRITE,VERBOSE=VERBOSE,INIT=INIT,R_FILES=R_FILES,R_DATASETS=R_DATASETS,R_MAPS=R_MAPS,DMAPS=D_MAPS,DPRODS=D_PRODS,DATERANGE=DATERANGE,DATASETS=DATASETS
    IF NONE(DATASETS) THEN DATASET = 'OCCCI' ELSE DATASET=DATASETS
    IF DATERANGE EQ [] THEN DATERANGE = GET_DATERANGE(SENSOR_DATES(DATASET))
    
    PLUN, LUN, 'Starting ' + SNAME + '...', 1
    PERIODS = ['D8','DOY']
    MAP_OUT = 'NESGRID4'
    IMG_MAP = 'NES'
        
    BFILES = []
    FOR PER=0, N_ELEMENTS(PERIODS)-1 DO BEGIN
      PERIOD = PERIODS[PER]
      EFILES = []
      
      CASE PERIOD OF
        'DOY': PRODS = ['CHLOR_A-CCI','MICRO-HIRATA','DIATOM-HIRATA','DINOFLAGELLATE-HIRATA']
        'D8':  PRODS = ['NANOPICO-BREWINSST_NES','CHLOR_A-CCI','MICRO-BREWINSST_NES','NANO-BREWINSST_NES','PICO-BREWINSST_NES']
      ENDCASE
    
      FOR PR=0, N_ELEMENTS(PRODS)-1 DO BEGIN
        APROD = PRODS[PR]
        
        FILES = GET_FILES('OCCCI',PERIOD=PERIOD,PRODS=APROD,DATERANGE=DR,COUNT=CT,VERSION=OCCCI_VERSION)
        IF CT EQ 0 THEN CONTINUE
        FP = PARSE_IT(FILES[0],/ALL)
        
        DIR_CDF = REPLACE(FP.DIR,[FP.MAP,'STATS'],[MAP_OUT,'NETCDF'])  & DIR_TEST, DIR_CDF
    ;    WRITE_NETCDF, FILES, DIR_OUT=DIR_CDF, MAP_OUT=MAP_OUT,/REVERSE_FILES

        DIR_PNG = REPLACE(FP.DIR,[FP.MAP,'STATS'],[IMG_MAP,'PNG']) & DIR_TEST, DIR_PNG
        PRODS_2PNG, FILES, DIR_OUT=DIR_PNG, /ADD_CB, /ADD_NAME, BUFFER=1, MAPP=IMG_MAP

      ENDFOR ; PRODS
    ENDFOR ; PERIODS  
  ENDIF ; NETCDFS  

  
; **************************************************************************************************    
  IF KEYWORD_SET(SUBAREA_EXTRACT) THEN BEGIN
; **************************************************************************************************
    SNAME = 'SUBAREA_EXTRACT'
    SWITCHES,SUBAREA_EXTRACT,OVERWRITE=OVERWRITE,VERBOSE=VERBOSE,INIT=INIT,R_FILES=R_FILES,R_DATASETS=R_DATASETS,R_MAPS=R_MAPS,DMAPS=D_MAPS,DPRODS=D_PRODS,DATERANGE=DATERANGE,DATASETS=DATASETS
    IF NONE(DATASETS) THEN DATASET = 'OCCCI' ELSE DATASET=DATASETS
    IF DATERANGE EQ [] THEN DATERANGE = GET_DATERANGE(SENSOR_DATES(DATASET))
    
    PLUN, LUN, 'Starting ' + SNAME + '...', 1
    PERIODS = ['D8','DOY']
    IF NONE(D_PRODS) THEN PRODS = ['ZEU','PICO-BREWINSST_NES','CHLOR_A-OCI','ZEU','PPD-VGPM2','MICRO-BREWINSST_NES','NANO-BREWINSST_NES'] ELSE PRODS = D_PRODS
    IF KEYWORD_SET(R_PRODS) THEN PRODS = REVERSE(PRODS)
        
    BFILES = []
    FOR PER=0, N_ELEMENTS(PERIODS)-1 DO BEGIN
      PERIOD = PERIODS[PER]
      EFILES = []
      DATFILE = DIR_DATA + PERIOD + '-' + DATASET + '-' + SUBAREA + '-' + VERSION + '.SAV'
      IF PERIOD NE 'DOY' THEN BEGIN
        YEARS = YEAR_RANGE(DATE_2YEAR(DATERANGE[0]),DATE_2YEAR(DATERANGE[1])) 
        IF KEYWORD_SET(R_FILES) THEN YEARS = REVERSE(YEARS)
      ENDIF ELSE YEARS=''
      FOR PR=0, N_ELEMENTS(PRODS)-1 DO BEGIN
        APROD = PRODS[PR]
        SAV = DIR_DATA + PERIOD + '-' + DATASET + '-' + SUBAREA + '-' + APROD + '.SAV' 
        FOR Y=0, N_ELEMENTS(YEARS)-1 DO BEGIN
          IF YEARS[Y] EQ '' THEN DR = DATERANGE ELSE DR = YEARS[Y]
          FILES = GET_FILES('OCCCI',PERIOD=PERIOD,PRODS=APROD,DATERANGE=DR,COUNT=CT)
          IF CT EQ 0 THEN CONTINUE
          SUBAREAS_EXTRACT, FILES, SHP_NAME=SUBAREA, INIT=INIT, VERBOSE=VERBOSE, DIR_OUT=DIR_DATA, STRUCT=STR, SAVEFILE=SAV, OUTPUT_STATS=OSTATS, BAD_FILES=BAD_FILES, /SKIP_BAD
          IF ANY(BAD_FILES) THEN BFILES = [BFILES,BAD_FILES]  
        ENDFOR
        EFILES = [EFILES,SAV]          
        IF PR EQ 0 THEN STRUCT = STR ELSE STRUCT = STRUCT_CONCAT(STRUCT,STR)
      ENDFOR ; DSETS
    
      IF FILE_MAKE(EFILES,DATFILE,OVERWRITE=OVERWRITE) EQ 1 THEN BEGIN
        SAVE, STRUCT, FILENAME=DATFILE ; ===> SAVE THE MERGED DATAFILE
        SAVE_2CSV, DATFILE
      ENDIF
  
    ENDFOR
    IF N_ELEMENTS(BFILES) GT 0 THEN BEGIN
      LI, BFILES
      STOP
      FILE_DELETE, BFILES
    ENDIF
   
  ENDIF ; SUBAREA_EXTRACT


; **************************************************************************************************
  IF KEYWORD_SET(PLOT_DATA) THEN BEGIN
; **************************************************************************************************
    SNAME = 'PLOT_DATA'
    SWITCHES,PLOT_DATA,OVERWRITE=OVERWRITE,VERBOSE=VERBOSE,INIT=INIT,R_FILES=R_FILES,R_DATASETS=R_DATASETS,R_MAPS=R_MAPS,DMAPS=D_MAPS,DPRODS=D_PRODS,DATERANGE=DATERANGE,DATASETS=DATASETS
    IF NONE(DATASETS) THEN DATASET = 'OCCCI' ELSE DATASET=DATASETS
    IF DATERANGE EQ [] THEN DATERANGE = GET_DATERANGE(SENSOR_DATES(DATASET))
    BUFFER=1
    PLTPROD='CHLOR_A' & YRANGE=[0,3]
    DATERANGE = ['19970904','20191231']

    COMPARE_SAT_PRODS, ['MICRO-BREWINSST_NES'],SENSORS='OCCCI',PERIODS='M',PROD_VERSION=['','PVER_1'],SHPFILES=SUBAREA,DIR_OUT=DIR_COMP,DATERANGE=DATERANGE
    COMPARE_SAT_PRODS, ['NANOPICO-BREWINSST_NES'],SENSORS='OCCCI',PERIODS='M',PROD_VERSION=['','PVER_1'],SHPFILES=SUBAREA,DIR_OUT=DIR_COMP,DATERANGE=DATERANGE
    COMPARE_SAT_PRODS, ['MICRO-BREWINSST_NES'],SENSORS='OCCCI',PERIODS='D8',PROD_VERSION=['','PVER_1'],SHPFILES=SUBAREA,DIR_OUT=DIR_COMP,DATERANGE=DATERANGE
    COMPARE_SAT_PRODS, ['MICRO-BREWINSST_NES'],SENSORS='OCCCI',PERIODS='DOY',PROD_VERSION=['','PVER_1'],SUBAREAS=SUBAREA,DIR_OUT=DIR_COMP,DATERANGE=DATERANGE

 stop   
    COMPARE_SAT_PRODS, ['MICRO-BREWINSST_NES','MICRO-HIRATA'],SENSORS='OCCCI',PERIODS='DOY',SHPFILES=SUBAREA,DIR_OUT=DIR_COMP,DATERANGE=DATERANGE
    COMPARE_SAT_PRODS, ['MICRO-BREWINSST_NES','MICRO-HIRATA_NES','MICRO-HIRATA'],SENSORS='OCCCI',PERIODS='DOY',SHPFILES=SUBAREA,DIR_OUT=DIR_COMP+SL,DATERANGE=DATERANGE
    COMPARE_SAT_PRODS, COMBO=['OCCCI,MICRO-HIRATA;OCCCI,DIATOM-HIRATA'], PERIODS='DOY', SHPFILES=SUBAREA, DIR_OUT=DIR_COMP, BUFFER=BUFFER,PLTPROD=PLTPROD,YRANGE=YRANGE,DATERANGE=DATERANGE
    COMPARE_SAT_PRODS, COMBO=['OCCCI,MICRO-HIRATA;OCCCI,DIATOM-HIRATA;OCCCI,DINOFLAGELLATE-HIRATA'], PERIODS='DOY', SHPFILES=SUBAREA, DIR_OUT=DIR_COMP,PLTPROD=PLTPROD,BUFFER=BUFFER,YRANGE=YRANGE,DATERANGE=DATERANGE
    COMPARE_SAT_PRODS, COMBO=['OCCCI,MICRO-BREWINSST_NES;OCCCI,MICRO-HIRATA;OCCCI,DIATOM-HIRATA'], PERIODS='DOY', SHPFILES=SUBAREA, DIR_OUT=DIR_COMP,PLTPROD=PLTPROD,YRANGE=YRANGE,DATERANGE=DATERANGE

stop
  ENDIF ; PLOT_DATA
  
; **************************************************************************************************
  IF KEYWORD_SET(HIRATA_RATIO) THEN BEGIN
; **************************************************************************************************
    SNAME = 'HIRATA_RATIO'
    SWITCHES,HIRATA_RATIO,OVERWRITE=OVERWRITE,VERBOSE=VERBOSE,INIT=INIT,R_FILES=R_FILES,R_DATASETS=R_DATASETS,R_MAPS=R_MAPS,DMAPS=D_MAPS,DPRODS=D_PRODS,DATERANGE=DATERANGE,DATASETS=DATASETS
    IF NONE(DATASETS) THEN DATASET = 'OCCCI' ELSE DATASET=DATASETS
    IF DATERANGE EQ [] THEN DATERANGE = GET_DATERANGE(SENSOR_DATES(DATASET))
    
    ENDIF ; HIRATA_RATIO  

END ; End of ATLANTIS_MAIN
