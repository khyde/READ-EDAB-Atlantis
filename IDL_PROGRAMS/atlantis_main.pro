; $ID:	ATLANTIS_MAIN.PRO,	2023-02-01-13,	USER-KJWH	$
  PRO ATLANTIS_MAIN

;+
; NAME:
;   ATLANTIS_MAIN
;
; PURPOSE:
;   The MAIN processing program for creating data inputs for the Atlantis modeling project (2020)
;
; PROJECT:
;   ATLANTIS
;
; CALLING SEQUENCE:
;   ATLANTIS_MAIN,$Parameter1$, $Parameter2$, $Keyword=Keyword$, ....
;
; REQUIRED INPUTS:
;   Parm1.......... Describe the positional input parameters here. 
;
; OPTIONAL INPUTS:
;   Parm2.......... Describe optional inputs here. If none, delete this section.
;
; KEYWORD PARAMETERS:
;   KEY1........... Document keyword parameters like this. Note that the keyword is shown in ALL CAPS!
;
; OUTPUTS:
;   OUTPUT.......... Decribe the output of this program or function
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
;   $Citations or any other useful notes$
;   
; COPYRIGHT: 
; Copyright (C) 2023, Department of Commerce, National Oceanic and Atmospheric Administration, National Marine Fisheries Service,
;   Northeast Fisheries Science Center, Narragansett Laboratory.
;   This software may be used, copied, or redistributed as long as it is not sold and this copyright notice is reproduced on each copy made.
;   This routine is provided AS IS without any express or implied warranties whatsoever.
;
; AUTHOR:
;   This program was written on February 01, 2023 by Kimberly J. W. Hyde, Northeast Fisheries Science Center | NOAA Fisheries | U.S. Department of Commerce, 28 Tarzwell Dr, Narragansett, RI 02882
;    
; MODIFICATION HISTORY:
;   Feb 01, 2023 - KJWH: Initial code written
;-
; ****************************************************************************************************
  ROUTINE_NAME = 'ATLANTIS_MAIN'
  COMPILE_OPT IDL2
  SL = PATH_SEP()
  
  IF ~N_ELEMENTS(VERSION) THEN VERSION = 'V3'
  
  DIR_OUT = !S.ATLANTIS + VERSION + SL 
  DIR_DATA = DIR_OUT + 'DATA' + SL
  DIR_COMP = DIR_OUT + 'COMPARE_PLOTS' + SL
  DIR_TEST, [DIR_OUT,DIR_DATA,DIR_COMP]

  CASE VERSION OF
    'V3': BEGIN
      DATASETS = 'OCCCI'  
      PERIODS = ['D8','DOY','M']
      PRODS = ['CHLOR_A-CCI','PPD-VGPM2','PSC-TURNER','PSC_FDIATOM-HIRATA','ZEU-VGPM2']
      DO_STATS = 0
      SUBAREA_EXTRACT = 0
      COMPARE_DATA = 1
      COMPARE_PRODS = []
      COMPARE_PERIODS = []
      DATERANGE = ['19971001','20220630']
      SHPFILE = 'ATLANTIS_NEUS'
    END  
    
  ENDCASE
  
  FOR D=0, N_ELEMENTS(DATASETS)-1 DO BEGIN
    DSET = DATASETS[D]
    
    ; ===> First make sure that the stats are current
    IF KEYWORD_SET(DO_STATS) THEN STACKED_STATS_WRAPPER, DSET, PRODS=PRODS, PERIODS=PERIODS, OUTSTATS=OUTSTATS, DATERANGE=DATERANGE, OVERWRITE=OVERWRITE

    IF KEYWORD_SET(SUBAREA_EXTRACT) THEN BEGIN
      FOR R=0, N_ELEMENTS(PRODS)-1 DO BEGIN
        APROD = PRODS[R]
        CASE APROD OF
          'PSC-TURNER': EPRODS='PSC_'+['MICRO','NANO','PICO','NANOPICO','FMICRO','FNANO','FPICO','FNANOPICO'] + '-TURNER'        
          ELSE: EPRODS = APROD
        ENDCASE
        
        FOR E=0, N_ELEMENTS(EPRODS)-1 DO BEGIN
          EPROD = EPRODS[E]        
          EFILES = []
          FOR I=0, N_ELEMENTS(PERIODS)-1 DO BEGIN
            APER = PERIODS[I]
            SAV = DIR_DATA + APER + '-' + DSET + '-' + SHPFILE + '-' + EPROD + '.SAV'
  
            FILES = GET_FILES(DSET,PRODS=EPROD,FILE_TYPE='STACKED_STATS',PERIOD=APER,COUNT=CT)
            IF CT EQ 0 THEN CONTINUE
            SUBAREAS_EXTRACT, FILES, SHP_NAME=SHPFILE, SUBAREAS=NAMES, VERBOSE=VERBOSE, DIR_OUT=DIR_EXTRACT, STRUCT=STR, SAVEFILE=SAV, EXTRACT_STAT=ESTAT, /ADD_DIR
            
            EFILES = [EFILES,SAV]
            IF I EQ 0 THEN STRUCT = STR ELSE STRUCT = STRUCT_CONCAT(STRUCT,STR)
          ENDFOR ; PERIODS
        ENDFOR ; EPRODS
      ENDFOR ; PRODS    

      IF FILE_MAKE(EFILES,DATFILE,OVERWRITE=OVERWRITE) EQ 1 THEN BEGIN
        SAVE, STRUCT, FILENAME=DATFILE ; ===> SAVE THE MERGED DATAFILE
        SAVE_2CSV, DATFILE
      ENDIF
        
    ENDIF ; KEYWORD_SET(SUBAREA_EXTRACT
      
    IF KEYWORD_SET(COMPARE_DATA) THEN BEGIN
      BUFFER=0
      PLTPROD='CHLOR_A' & YRANGE=[0,3]
      DATERANGE = ['19970904','20221231']
      PERIODS = ['M','DOY']
      SHAPES = ['NES_EPU_NOESTUARIES',SHPFILE]
      PRODS = ['CHLOR_A-CCI','PSC_'+['MICRO','NANO','PICO']+'-TURNER','PSC_FDIATOM-HIRATA']
      VCOMBOS = LIST(['4.2','5.0'],['4.2','6.0'],['5.0','6.0'],['4.2','5.0','6.0'])
      
      FOR N=0, N_ELEMENTS(PRODS)-1 DO BEGIN
        PLOTPROD = VALIDS('PRODS',PRODS[N])
        CASE PRODS[N] OF
          'CHLOR_A-CCI': BPROD='CHLOR_A-CCI'
          'PSC_MICRO-TURNER': BPROD='MICRO-BREWINSST_NES'
          'PSC_NANO-TURNER': BPROD='NANO-BREWINSST_NES'
          'PSC_PICO-TURNER': BPROD='PICO-BREWINSST_NES'
          'PSC_FDIATOM-HIRATA': BPROD='DIATOM_PERCENTAGE-HIRATA'
        ENDCASE
        FOR C=0, N_ELEMENTS(VCOMBOS)-1 DO BEGIN
          VCMB = VCOMBOS[C]
          CMBS = []
          FOR V=0, N_ELEMENTS(VCMB)-1 DO BEGIN
            CASE VCMB[V] OF
              '4.2': CMBS = [CMBS,'OCCCI,V4.2,STATS,'+BPROD+', ']
              '5.0': CMBS = [CMBS,'OCCCI,V5.0,STATS,'+BPROD+', ']
              '6.0': CMBS = [CMBS,'OCCCI,V6.0,STACKED_STATS,'+PRODS[N]+', ']
            ENDCASE
            CMB = STRJOIN(CMBS,';')
          ENDFOR

          FOR R=0, N_ELEMENTS(PERIODS)-1 DO BEGIN
            FOR S=0, N_ELEMENTS(SHAPES)-1 DO BEGIN
              COMPARE_SAT_PRODS, COMBO=CMB,PLTPROD=PLOTPROD,PERIODS=PERIODS[R],SHPFILES=SHAPES[S],DIR_OUT=DIR_COMP,DATERANGE=DATERANGE, BUFFER=BUFFER     
            ENDFOR ; SHAPES
          ENDFOR ; PERIODS
        ENDFOR  ; VCOMBOS
      ENDFOR ; PRODS  
        
    ENDIF ; COMPARE_DATA
  ENDFOR ; DATASETS


END ; ***************** End of ATLANTIS_MAIN *****************
