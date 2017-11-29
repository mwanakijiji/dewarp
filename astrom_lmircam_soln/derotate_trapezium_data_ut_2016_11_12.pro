PRO DEROTATE_TRAPEZIUM_DATA_UT_2016_11_12
; this makes simple derotations of Trapezium data that was taken for LMIRCam dewarping

; created 13 Dec 2016 by E.S.

stem = ' ~/../../media/unasemaje/Seagate Expansion Drive/lbti_data_reduction/lmircam_astrometry/ut_2016_11_12/asterism/processedData/'

; COMMENT OUT THE BELOW LOOP TO WRITE OUT INDIVIDUAL DEROTATED IMAGES
FOR frameNum = 1892,2251 DO BEGIN
  PRINT, 'Derotating frame'
  PRINT, frameNum  
  PRINT, '-----------'
  array = READFITS(stem+'step02_dewarped/lm_161112_'+STRING(frameNum,format="(i05)")+'.fits', header)
  pa = GET_PARALLACTIC_ANGLE(stem+'step02_dewarped/lm_161112_'+STRING(frameNum,format="(i05)")+'.fits')  
  derotatedArray = ROT(array, -pa, /INTERP)  
  WRITEFITS, stem+'step03_derotate/lm_161112_'+STRING(frameNum,format="(i05)")+'.fits', derotatedArray, header  
ENDFOR


FOR ditherPos = 0,17 DO BEGIN
  ditherCube = MAKE_ARRAY(2048,2048,20,/FLOAT,VALUE=0.)
  FOR sliceNum = 0,19 DO BEGIN
    frameNum = 1892+ditherPos*20+sliceNum
    PRINT, 'Frame number'
    PRINT, frameNum
    PRINT, 'Dither position'
    PRINT, ditherPos
    PRINT, '-----------'
    ditherCube[*,*,sliceNum] = READFITS(stem+'step03_derotate/lm_161112_'+STRING(frameNum,format="(i05)")+'.fits')
  ENDFOR  
  ditherMedian = MEDIAN(ditherCube, DIMENSION=3)
  WRITEFITS, stem+'step04_ditherMedians/median_dither_'+STRING(ditherPos, FORMAT="(i02)")+'.fits', ditherMedian, header
ENDFOR

END