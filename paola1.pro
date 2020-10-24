
;;;;;;;;;;;;;;;
; MAIN FUNCTION
;;;;;;;;;;;;;;;
FUNCTION PAOLA1,MODE,DIM,TSC,W0,L0,ZA,HEIGHT,DISTCN2,P09,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20,P21,P22,$
                DISPERSION=DISPERSION,$
                FWHM_ANALYSIS=FWHM_ANALYSIS,$
                INFO=INFO,PSF=PSF,OTF=OTF,SF=SF,PSD=PSD,ONLY_PSD=ONLY_PSD,$
                ANTI_ALIAS=ANTI_ALIAS,SCINTILLATION=SCINTILLATION,$
                FITSCODE=FITSCODE,SINGLE=SINGLE,LOGCODE=LOGCODE,X_PSF=X_PSF,Y_PSF=Y_PSF,PRECISION=PRECISION,$
                OPTIMIZE_WFS_INT=OPTIMIZE_WFS_INT,OPTIMIZE_LOOP_GAIN=OPTIMIZE_LOOP_GAIN,OPTIMIZE_ALL=OPTIMIZE_ALL,$
                POST_TIPTILT=POST_TIPTILT,TILT_ANGLE_STEP=TILT_ANGLE_STEP,FRFFT=FRFFT

  COMMON OPT_BLOCK,block
  ;CHECKING GLAO MODE
  ;------------------
  if  (strtrim(MODE,2)) eq 'glao' then begin
    if (size(HEIGHT))(0) ne 1 then message,'LAYERS ALTITUDES (HEIGHT) MUST BE A VECTOR IN GLAO MODE'
    if (size(DISTCN2))(0) ne 1 then message,'LAYERS CN2DH DISTRIBUTION (DISTCN2) MUST BE A VECTOR IN GLAO MODE'
    if n_params() eq 13 then begin
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS',P09,'structure',{dim:[1,3],tags:['dm_height','dmtf','actpitch']},'no','free','free'
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS.DM_HEIGHT',P09.DM_HEIGHT,'real',0,'no','free','free'
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS.ACTPITCH',P09.ACTPITCH,'real',0,'yes','++','free'
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS.DMTF',P09.DMTF,'real',2,'yes','free','free'
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P10,'structure',{dim:[1,1],tags:['wfs_pitch']},'no','free','free'
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.WFS_PITCH',P10.WFS_PITCH,'real',0,'yes','++','free'
      VALID_INPUT,'PAOLA.PRO','SO_ANG',P11,'real',0,'no','0+','free'
      VALID_INPUT,'PAOLA.PRO','SO_ORI',P12,'real',0,'no','free','free'
      VALID_INPUT,'PAOLA.PRO','GLAO_WFS',P13,'structure',{dim:[1,2],tags:['type','ang']},'no','free','free'
      VALID_INPUT,'PAOLA.PRO','GLAO_WFS.TYPE',P13.TYPE,'string',0,'no',['full','edge','FULL','EDGE'],'free'
      VALID_INPUT,'PAOLA.PRO','GLAO_WFS.ANG',P13.ANG,'real',0,'no','0+','free'
    endif else begin
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS',P10,'structure',{dim:[1,3],tags:['dm_height','dmtf','actpitch']},'no','free','free'
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS.DM_HEIGHT',P10.DM_HEIGHT,'real',0,'no','free','free'
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS.DMTF',P10.DMTF,'real',2,'yes','free','free'
      VALID_INPUT,'PAOLA.PRO','DM_PARAMS.ACTPITCH',P10.ACTPITCH,'real',0,'yes','++','free'
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P11,'structure',[1,-1],'no','free','free'
      if n_params() eq 20 and n_tags(P11) ne 1 then message,'WFS_PARAMS STRUCTURE VARIABLE INPUT MUST HAVE 1 COMPONENTS'
      if n_params() eq 22 and n_tags(P11) ne 6 and n_tags(P11) ne 9 then message,'WFS_PARAMS STRUCTURE VARIABLE INPUT MUST HAVE 6 OR 9 COMPONENTS'
      if n_tags(P11) eq 1 then VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P11,'structure',{dim:[1,1],tags:['wfs_pitch']},'no','free','free'
      if n_tags(P11) eq 6 then VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P11,'structure',{dim:[1,6],tags:['mirvec','nblenses','extrafilter','wfs_pitch','wfs_ron','algorithm']},'no','free','free'
      if n_tags(P11) eq 9 then VALID_INPUT,'PAOLA.PRO','WFS_PARAMS',P11,'structure',{dim:[1,9],$
        tags:['mirvec','nblenses','extrafilter','wfs_pitch','wfs_pxfov','wfs_pxsize','wfs_jitter','wfs_ron','algorithm']},'no','free','free'
      VALID_INPUT,'PAOLA.PRO','WFS_PARAMS.WFS_PITCH',P11.WFS_PITCH,'real',0,'yes','++','free'
      VALID_INPUT,'PAOLA.PRO','SO_ANG',P12,'real',0,'no','0+','free'
      VALID_INPUT,'PAOLA.PRO','SO_ORI',P13,'real',0,'no','free','free'
      VALID_INPUT,'PAOLA.PRO','WFS_INT',P14,'real',0,'yes','0+','free'
      VALID_INPUT,'PAOLA.PRO','LAG',P15,'real',0,'no','0+','free'
      if P14 eq 0 and P15 ne 0 then message,'GS TIME LAG MUST BE 0 IF WFS INTEGRATION TIME IS 0'
      VALID_INPUT,'PAOLA.PRO','LOOP_MODE',P16,'string',0,'no',['open','closed','OPEN','CLOSED'],'free'
      VALID_INPUT,'PAOLA.PRO','LOOP_GAIN',P17,'real',0,'no',[0,4.9348],0
      VALID_INPUT,'PAOLA.PRO','GLAO_WFS',P18,'structure',{dim:[1,2],tags:['type','ang']},'no','free','free'
      VALID_INPUT,'PAOLA.PRO','GLAO_WFS.TYPE',P18.TYPE,'string',0,'no',['star','STAR'],'free'
      VALID_INPUT,'PAOLA.PRO','GLAO_WFS.ANG',P18.ANG,'real',[2,2,-1],'no','free','free'
      VALID_INPUT,'PAOLA.PRO','GS_WEIGHT',P19,'real',1,'yes','++','free'
      if n_elements(P19) gt 1 then if n_elements(P19) ne n_elements(P18.ANG(0,*)) then message,'GUIDE STARS NUMBERS IN GS_WEIGHT AND GLAO_WFS.ANG ARE NOT THE SAME'
      if n_params() eq 20 then begin
        VALID_INPUT,'PAOLA.PRO','WFS_NEA',P20,'real',1,'no','0+','free'
        if n_elements(P20) ne n_elements(P18.ANG(0,*)) then $
          message,'NUMBER OF GUIDE STARS IN GLAO_WFS.ANG AND IN THE WFS_NEA ARRAY ARE NOT THE SAME - CHECK THE NUMBER OF ELEMENTS IN WFS_NEA'
      endif

    endelse
  endif
  ;GLAO MODE
  ;---------
  if strlowcase(strtrim(MODE,2)) eq 'glao' then begin
;    if keyword_set(OPTIMIZE_WFS_INT) then begin ; WFS integration time optimization
;      if keyword_set(INFO) then print,'...GLAO WFS INTEGRATION TIME OPTIMIZATION...'
;      if keyword_set(INFO) then print,''
;      if keyword_set(INFO) then print,'    WFS INT [MS]     WFE SERVO+NOISE [NM] '
;      nstar=n_elements(reform(GLAO_WFS.ANG(0,*)))
;      if n_elements(GS_WEIGHT) eq 1 then gsw=dblarr(nstar)+1.d/nstar
;      if n_elements(GS_WEIGHT) ne 1 then gsw=GS_WEIGHT
;      block={n_psd:DIM.N_LF,fpcoo:fpcoolf,dfplam:DIM.DFP_LF,dext:TSC.DEXTMAX,lambda:DIM.LAMBDA,height:newheight,$
;        wind:newwind,r0500:r0500,r0500_i:r0500_i,l0:double(L0),za:ZA,mirvec:MIRVEC,nblenses:NBLENSES,$
;        extrafilter:EXTRAFILTER,wfs_pitch:WFS_PITCH,wfs_ron:WFS_RON,$
;        dmtf:FTMD,dmh:DM_HEIGHT,ang:GS_ANG,ori:GS_ORI,wfs_int:gsint,lag:LAG,$
;        algor:ALGORITHM,ngs_mag:NGS_MAG,filter:FILTER,ngs_tem:NGS_TEM,gs_weight:gsw,glao_wfs:GLAO_WFS,$
;        prec:PRECISION,prt:keyword_set(INFO),alias:1-keyword_set(ANTI_ALIAS)}
;      if strlowcase(ALGORITHM) eq 'cg' then block=create_struct(block,'wfs_fov',WFS_PXFOV,'wfs_pxs',WFS_PXSIZE,'wfs_jit',WFS_JITTER)
;      ;      if gsint lt 0 then minF_parabolic,1e-4,timescal,1e4,gsint,fmin,FUNC_NAME='OPTIMIZE_GLAO',TOLERANCE=1e-5
;      ;      if gsint gt 0 then minF_parabolic,1e-4,gsint,1e4,gsint,fmin,FUNC_NAME='OPTIMIZE_GLAO',TOLERANCE=1e-5
;      if keyword_set(INFO) then print,'======================================='
;      if keyword_set(INFO) then printf,-1,format='("OPTIMIZED WFS_INT     [MS]",d13.6)',gsint
;      ;      if keyword_set(INFO) then printf,-1,format='("OPTIMIZED VARIANCE [RAD^2]",d13.6)',fmin(0)
;      if keyword_set(INFO) then print,'================================================'
;      if keyword_set(INFO) then print,''
;      if DIM.DFP_LF gt 0.1/((gsint+2*LAG)*1e-3*meanwind) then begin
;        print,'################################################'
;        print,'WARNING WARNING WARNING WARNING WARNING WARNING'
;        print,'################################################'
;        print,''
;        print,'WFS integration time is now too large for the'
;        print,'current servo-lag power spectrum: numerical'
;        print,'undersampling will occur in the calculation'
;        print,'of the servo-lag PSD, and the servo-lag'
;        print,'variance will be strongly underestimated.'
;        print,'To avoid this:'
;        print,'Re-run the function PIXMATSIZE.PRO with the'
;        print,'optimized WFS integration time value.'
;        print,''
;        print,'COMPUTATION IS CONTINUING ANYWAY'
;        print,''
;        print,'################################################'
;        print,''
;        warning='yes'
;      endif
;    endif
    ;==============================================================================
    ;=========================END OF ARGUMENTS CHECK===============================
    ;==============================================================================

    if not keyword_set(PRECISION) then PRECISION='double'

    ;SETTING VARIABLES ACCORDING TO MODE
    ;-----------------------------------
    if strlowcase(strtrim(MODE,2)) eq 'seli' and n_params() eq 9 then WIND=double(P09)
    if strlowcase(strtrim(MODE,2)) eq 'ngs' then begin
      WIND=double(P09)
      DM_HEIGHT=double(P10.DM_HEIGHT)
      ACTPITCH=double(P10.ACTPITCH)
      DMTF=double(P10.DMTF)
      WFS_PITCH=double(P11.WFS_PITCH)
      GS_ANG=double(P12)
      GS_ORI=double(P13)
      WFS_INT=double(P14)
      LAG=double(P15)
      LOOP_MODE=P16
      LOOP_GAIN=double(P17)
      if n_params() eq 20 then begin
        MIRVEC=P11.MIRVEC
        NBLENSES=P11.NBLENSES
        EXTRAFILTER=double(P11.EXTRAFILTER)
        WFS_RON=double(P11.WFS_RON)
        ALGORITHM=P11.ALGORITHM
        if strlowcase(ALGORITHM) eq 'cg' then begin
          WFS_JITTER=double(P11.WFS_JITTER)
          WFS_PXFOV=long(P11.WFS_PXFOV)
          WFS_PXSIZE=double(P11.WFS_PXSIZE)
        endif
        NGS_MAG=double(P18)
        FILTER=strlowcase(strtrim(P19,2))
        NGS_TEM=double(P20)
      endif else WFS_NEA=double(P18)
    endif
    if strlowcase(strtrim(MODE,2)) eq 'glao' then begin
      if n_params() eq 13 then begin
        DM_HEIGHT=double(P09.DM_HEIGHT)
        ACTPITCH=double(P09.ACTPITCH)
        DMTF=double(P09.DMTF)
        WFS_PITCH=double(P10.WFS_PITCH)
        GS_ANG=double(P11)
        GS_ORI=double(P12)
        GLAO_WFS=P13
      endif
      if n_params() eq 20 or n_params() eq 22 then begin
        WIND=double(P09)
        DM_HEIGHT=double(P10.DM_HEIGHT)
        ACTPITCH=double(P10.ACTPITCH)
        DMTF=double(P10.DMTF)
        WFS_PITCH=double(P11.WFS_PITCH)
        GS_ANG=double(P12)
        GS_ORI=double(P13)
        WFS_INT=double(P14)
        LAG=double(P15)
        LOOP_MODE=P16
        if strlowcase(LOOP_MODE) eq 'closed' then message,'CLOSED LOOP NOT IMPLEMENTED IN GLAO MODE, YET. SORRY. SET LOOP MODE TO ''OPEN'' FOR NOW.'
        LOOP_GAIN=double(P17)
        GLAO_WFS=P18
        GS_WEIGHT=double(P19)
      endif
      if n_params() eq 20 then WFS_NEA=double(P20)
      if n_params() eq 22 then begin
        MIRVEC=P11.MIRVEC
        NBLENSES=P11.NBLENSES
        EXTRAFILTER=double(P11.EXTRAFILTER)
        WFS_RON=double(P11.WFS_RON)
        ALGORITHM=P11.ALGORITHM
        if strlowcase(ALGORITHM) eq 'cg' then begin
          WFS_JITTER=double(P11.WFS_JITTER)
          WFS_PXFOV=long(P11.WFS_PXFOV)
          WFS_PXSIZE=double(P11.WFS_PXSIZE)
        endif
        NGS_MAG=double(P20)
        FILTER=strlowcase(strtrim(P21,2))
        NGS_TEM=double(P22)
      endif
    endif

    ;HANDLING OUTER SCALE
    ;--------------------
    if L0 gt 10000 and keyword_set(INFO) then print,'***********WARNING: OUTER SCALE LARGER THAN 10 KM. CONSIDERED INFINITE.'
    if L0 gt 10000 then L0=-1

    ;SEPARE WFS_INT AND LOOP_GAIN INPUT FROM OUTSIDE WORLD,
    ;BECAUSE OPTIMIZATION CHANGES INPUT VALUES
    ;---------------------------------------------------
    if size(WFS_INT,/type) ne 0 then gsint=double(WFS_INT)
    if size(LOOP_GAIN,/type) ne 0 then loopgain=double(LOOP_GAIN)
    if strlowcase(strtrim(MODE,2)) eq 'glao' then if strlowcase(GLAO_WFS.TYPE) ne 'star' then gsint=0

    ;FRIED PARAMETER VERSUS WAVELENGTH
    ;---------------------------------
    rad2asec=3600.d*180.d/!dpi
    asec2rad=1.d/rad2asec
    W0rad=W0*asec2rad
    r0500=0.98d*0.5d-6/W0rad*cos(double(ZA)/180*!dpi)^(3.d/5) ; Fried's r0 @ 500 nm @ z=ZA
    r0LAM=r0500*(DIM.LAMBDA/0.5)^(1.2d)

    ;SEEING ANGLE ATTENUATION DUE TO OUTER SCALE
    ;-------------------------------------------
    if L0 eq -1 then aos=1
    if L0 ne -1 then aos=(ATTOS_FWHM(double(TSC.DEXTMAX)/r0LAM,double(TSC.DEXTMAX)/L0))(0)

    ;REMOVING LAYERS FOR WHICH CN2 = 0 AND ADAPTING LAYERS HEIGHT TO ZENITH ANGLE
    ;----------------------------------------------------------------------------
    if n_params() gt 6 then begin
      tmp=DISTCN2/total(DISTCN2) ; normalization of cn2 profile distribution
      w=where(tmp gt 0)
      if w(0) ne -1 then begin
        newdistcn2=tmp(w)
        newheight=HEIGHT(w)/cos(double(ZA)/180*!dpi)
        if not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then begin
          newwindx=(reform(WIND(*,0)))(w)
          newwindy=(reform(WIND(*,1)))(w)
          newwind=[[newwindx],[newwindy]]
        endif
      endif else begin
        newdistcn2=tmp
        newheight=HEIGHT/cos(double(ZA)/180*!dpi)
        if not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then newwind=WIND
      endelse
    endif

    ;FRIED PARAMETER AND Cn2*dh FOR EACH LAYER
    ;-----------------------------------------
    if n_params() gt 6 then begin
      cn2dh_i=(500d-9/2/!dpi)^2/0.423d*r0500^(-5.d/3)*newdistcn2
      r0500_i=r0500*newdistcn2^(-3.d/5)
    endif

    ;LAYERS MEAN ALTITUDE AND ISOPLANATIC ANGLE
    ;------------------------------------------
    if n_params() gt 6 then begin
      meanalti=total(newdistcn2*abs(newheight)^(5.d/3))^(3.d/5)
      if meanalti ne 0 then anisoang=0.314d/meanalti*r0500*(DIM.LAMBDA/0.5d)^(1.2d)*rad2asec
      if meanalti eq 0 then anisoang=0.d
    endif

    ;LAYERS MEAN VELOCITY AND TURBULENT PHASE LIFE TIME
    ;--------------------------------------------------
    if n_params() gt 6 then if not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then begin
      if (size(newwind))(0) ne 0 then begin
        meanwind=total(newdistcn2*sqrt(newwind(*,0)^2+newwind(*,1)^2)^(5.d/3))^(3.d/5)
        if meanwind ne 0 then timescal=0.314d/meanwind*r0500*(DIM.LAMBDA/0.5d)^(1.2d)*1d3
        if meanwind eq 0 then timescal=0
      endif else begin
        meanwind=newwind
        if meanwind ne 0 then timescal=0.314d/meanwind*r0500*(DIM.LAMBDA/0.5d)^(1.2d)*1d3
        if meanwind eq 0 then timescal=0
      endelse
    endif

    ;ACTUATOR AND WFS PITCH
    ;----------------------
    if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
      ACTPITCHnom=r0LAM ; NOMINAL PITCH = r0 @ <lambda>
      if ACTPITCH eq -1 then ACTPITCH=ACTPITCHnom ; THE DEFAULT VALUE (-1) IS SET TO THE NOMINAL VALUE ACTPITCHnom
      dm_nactnom=TSC.DEXTMAX/ACTPITCHnom
      dm_nact=double(TSC.DEXTMAX)/ACTPITCH
      if WFS_PITCH eq -1 then WFS_PITCH=ACTPITCH
    endif

    ;WFS OPTICAL BANDWIDTH
    ;---------------------
    if size(WFS_BDWTH,/type) ne 0 then if WFS_BDWTH eq -1 then WFS_BDWTH=[0.332,1.000]

    ;PUPIL PLANE COORDINATE RADIUS [m]
    ;---------------------------------
    if strlowcase(PRECISION) eq 'single' then xpcoohf=(COOGRID(DIM.N_OTF,DIM.N_OTF,SCALE=DIM.N_OTF/2*DIM.DXP,/FT,/RADIUS,/SINGLE)).r
    if strlowcase(PRECISION) eq 'double' then xpcoohf=(COOGRID(DIM.N_OTF,DIM.N_OTF,SCALE=DIM.N_OTF/2*DIM.DXP,/FT,/RADIUS)).r

    ;PUPIL PLANE SPATIAL FREQUENCY [1/m]
    ;-----------------------------------
    if strlowcase(PRECISION) eq 'single' and strlowcase(strtrim(MODE,2)) ne 'seli' then fpcoolf=(COOGRID(DIM.N_LF,DIM.N_LF,SCALE=(DIM.N_LF-1)/2*DIM.DFP_LF,/COO_X,/SINGLE)).x
    if strlowcase(PRECISION) eq 'double' and strlowcase(strtrim(MODE,2)) ne 'seli' then fpcoolf=(COOGRID(DIM.N_LF,DIM.N_LF,SCALE=(DIM.N_LF-1)/2*DIM.DFP_LF,/COO_X)).x
    if strlowcase(PRECISION) eq 'single' then fpcoohf=(COOGRID(DIM.N_OTF,DIM.N_OTF,SCALE=DIM.N_OTF/2*DIM.DFP,/FT,/COO_X,/SINGLE)).x
    if strlowcase(PRECISION) eq 'double' then fpcoohf=(COOGRID(DIM.N_OTF,DIM.N_OTF,SCALE=DIM.N_OTF/2*DIM.DFP,/FT,/COO_X)).x

    ;DM SPATIAL TRANSFER FUNCTION
    ;----------------------------
    if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
      if strlowcase(PRECISION) eq 'single' then begin ; the 1+1e-3 is there to avoid a boundary issue
        if (size(DMTF))(0) eq 0 then FTMD=float(abs(fpcoolf) le (1+1e-3)*0.5d/ACTPITCH and abs(rotate(fpcoolf,1)) le (1+1e-3)*0.5d/ACTPITCH)
        if (size(DMTF))(0) eq 2 then FTMD=float(DMTF)
      endif else begin
        if (size(DMTF))(0) eq 0 then FTMD=double(abs(fpcoolf) le (1+1e-3)*0.5d/ACTPITCH and abs(rotate(fpcoolf,1)) le (1+1e-3)*0.5d/ACTPITCH)
        if (size(DMTF))(0) eq 2 then FTMD=DMTF
      endelse
    endif

    ;VARIANCE OF PHASE COMPONENTS
    ;----------------------------
    var=dblarr(4)

    ;DISPERSION DEFAULT VALUE
    ;------------------------
    if not keyword_set(DISPERSION) then DISPERSION=1

    ;==========================================================================
    ;======================END OF SETTINGS SECTION=============================
    ;==========================================================================

    ; Fitting error PSD
    tmp=PSD_HFERR_NGS(DIM.N_OTF,fpcoohf,r0LAM,double(L0),0,WFS_PITCH,PRECISION) ; fitting error PSD
    Whf=tmp
    var(0)=0.2313d*(WFS_PITCH/r0LAM)^(5.d/3)
    if strlowcase(PRECISION) eq 'single' then Wlf=fltarr(DIM.N_LF,DIM.N_LF)
    if strlowcase(PRECISION) eq 'double' then Wlf=dblarr(DIM.N_LF,DIM.N_LF)
    ; WFS aliasing PSD
    if not keyword_set(ANTI_ALIAS) and strlowcase(GLAO_WFS.TYPE) eq 'star' then begin
      if size(gsint,/type) eq 0 then begin
        tmp=DISPERSION(0)^2*PSD_ALIAS_NGS_GLAO(DIM.N_LF,fpcoolf,DIM.LAMBDA,newheight,[[0],[0]],r0500,double(L0),0,WFS_PITCH,FTMD,DM_HEIGHT,0,GS_WEIGHT,GLAO_WFS,PRECISION)
      endif else tmp=DISPERSION(0)^2*PSD_ALIAS_NGS_GLAO(DIM.N_LF,fpcoolf,DIM.LAMBDA,newheight,newwind,r0500_i,double(L0),0,WFS_PITCH,FTMD,DM_HEIGHT,gsint,GS_WEIGHT,GLAO_WFS,PRECISION)
      if keyword_set(PSD) or keyword_set(ONLY_PSD) then Walias=tmp
      Wlf=Wlf+tmp
      var(1)=total(tmp)*DIM.DFP_LF^2
    endif
    ; aniso-servo PSD
    if strlowcase(GLAO_WFS.TYPE) eq 'star' then begin
      if n_params() eq 13 then tmp=PSD_ANISO_SERVO_GLAO_STAR(DIM.N_LF,fpcoolf,DIM.LAMBDA,newheight,dblarr(n_elements(newheight),2),r0500_i,double(L0),WFS_PITCH,FTMD,$
        DM_HEIGHT,GS_ANG,GS_ORI,0,0,GS_WEIGHT,GLAO_WFS,PRECISION)
      if n_params() ne 13 then tmp=PSD_ANISO_SERVO_GLAO_STAR(DIM.N_LF,fpcoolf,DIM.LAMBDA,newheight,newwind,r0500_i,double(L0),WFS_PITCH,FTMD,$
        DM_HEIGHT,GS_ANG,GS_ORI,gsint,LAG,GS_WEIGHT,GLAO_WFS,PRECISION)
    endif
    if keyword_set(PSD) or keyword_set(ONLY_PSD) then Wanisoservo=tmp
    Wlf=Wlf+tmp
    var(2)=total(tmp)*DIM.DFP_LF^2
    ; WFS noise PSD
    if strlowcase(GLAO_WFS.TYPE) eq 'star' then begin
      tmp=dblarr(DIM.N_LF,DIM.N_LF)
      if size(NGS_MAG,/type) ne 0 then nstar=n_elements(NGS_MAG)
      if size(WFS_NEA,/type) ne 0 then nstar=n_elements(WFS_NEA)
      if n_elements(GS_WEIGHT) eq 1 then gsw=dblarr(nstar)+1.d/nstar
      if n_elements(GS_WEIGHT) ne 1 then gsw=GS_WEIGHT
      ; WFS noise PSD from the GS magnitudes and spectrums and RON
      if size(NGS_MAG,/type) ne 0 then begin
        tmp=NBPHOTONS(WFS_PITCH^2,MIRVEC,NBLENSES,ZA,NGS_MAG,FILTER,NGS_TEM,gsint,EXTRAFILTER)
        wfsnph=tmp.nph
        wfslam=tmp.lam
        wfsbdw=tmp.bdw
        wfstau=tmp.aot
        if strlowcase(ALGORITHM) eq '4q' then nea=(SH_NEA(wfsnph,r0500,wfslam,WFS_PITCH,WFS_RON,ALGORITHM)).nea
        if strlowcase(ALGORITHM) eq 'cg' then nea=(SH_NEA(wfsnph,r0500,wfslam,WFS_PITCH,WFS_RON,ALGORITHM,WFS_JITTER,WFS_PXFOV,WFS_PXSIZE)).nea
        
        tmp=PSD_NOISE_NGS_GLAO(DIM.N_LF,fpcoolf,WFS_PITCH,FTMD,gsw,nea,DIM.LAMBDA,PRECISION)
        if keyword_set(PSD) or keyword_set(ONLY_PSD) then Wnoise=tmp
        Wlf=Wlf+tmp
        var(3)=total(tmp)*DIM.DFP_LF^2
      endif
      ; WFS noise PSD from GS Noise Equivalent Angle
      if size(WFS_NEA,/type) ne 0 then begin
        tmp=PSD_NOISE_NGS_GLAO(DIM.N_LF,fpcoolf,WFS_PITCH,FTMD,gsw,WFS_NEA,DIM.LAMBDA,PRECISION)
        if keyword_set(PSD) or keyword_set(ONLY_PSD) then Wnoise=tmp
        Wlf=Wlf+tmp
        var(3)=total(tmp)*DIM.DFP_LF^2
      endif
    endif
  endif
print,'######################################'
print,var


  ;ADAPTIVE OPTICS STRUCTURE FUNCTION
  ;----------------------------------
  ;it is computed from the Fourier transform of the phase PSD
  if strlowcase(strtrim(MODE,2)) ne 'seli' then if not keyword_set(ONLY_PSD) then begin
    if strlowcase(PRECISION) eq 'single' then begin
      ;low frequency structure function
      varLF=float(total(Wlf)*DIM.DFP_LF^2)
      if varLF gt 0 then begin
        if not keyword_set(FRFFT) then begin
          dimPSD=DIM.N_LF_PADDED>DIM.N_OTF ; we take the largest dimension
          tmp=fltarr(dimPSD,dimPSD)
          tmp(dimPSD/2-(DIM.N_LF-1)/2:dimPSD/2+(DIM.N_LF-1)/2,dimPSD/2-(DIM.N_LF-1)/2:dimPSD/2+(DIM.N_LF-1)/2)=float(Wlf) ; we insert the LF PSD into the matrix tmp
          tmp=float(MATHFT(tmp,ic=dimPSD/2,jc=dimPSD/2))
        endif else begin
          dimPSD=2*DIM.N_LF>DIM.N_OTF ; we take the largest dimension
          tmp=fltarr(dimPSD,dimPSD)
          tmp(dimPSD/2-(DIM.N_LF-1)/2:dimPSD/2+(DIM.N_LF-1)/2,dimPSD/2-(DIM.N_LF-1)/2:dimPSD/2+(DIM.N_LF-1)/2)=float(Wlf) ; we insert the LF PSD into the matrix tmp
          tmp=float(FRACTIONAL_FFT(tmp,dimPSD,DIM.N_LF_PADDED,dimPSD/2,dimPSD/2,dimPSD/2,dimPSD/2))
        endelse
        tmp=tmp(dimPSD/2-DIM.N_OTF/2:dimPSD/2+DIM.N_OTF/2-1,dimPSD/2-DIM.N_OTF/2:dimPSD/2+DIM.N_OTF/2-1)/tmp(dimPSD/2,dimPSD/2)*varLF ; we keep the DIM.N_OTF part and renormalize with variance
        dphiLF=2*(varLF-tmp)
        dphiLF=dphiLF-min(dphiLF) ; in principle min(SF)=0, but there are always round-off errors
      endif else dphiLF=0
      ;high frequency structure function
      tmp=float(MATHFT(Whf,dx=DIM.DXP,ic=DIM.N_OTF/2,jc=DIM.N_OTF/2,/inverse))
      tmp=tmp/tmp(DIM.N_OTF/2,DIM.N_OTF/2)*var(0)
      dphiHF=2*(float(var(0))-tmp)
    endif else begin
      ;low frequency structure function
      varLF=total(Wlf)*DIM.DFP_LF^2
      if varLF gt 0 then begin
        if not keyword_set(FRFFT) then begin
          dimPSD=DIM.N_LF_PADDED>DIM.N_OTF ; we take the largest dimension
          tmp=dblarr(dimPSD,dimPSD)
          tmp(dimPSD/2-(DIM.N_LF-1)/2:dimPSD/2+(DIM.N_LF-1)/2,dimPSD/2-(DIM.N_LF-1)/2:dimPSD/2+(DIM.N_LF-1)/2)=Wlf ; we insert the LF PSD into the matrix tmp
          tmp=double(MATHFT(tmp,ic=dimPSD/2,jc=dimPSD/2))
        endif else begin
          dimPSD=2*DIM.N_LF>DIM.N_OTF ; we take the largest dimension
          tmp=dblarr(dimPSD,dimPSD)
          tmp(dimPSD/2-(DIM.N_LF-1)/2:dimPSD/2+(DIM.N_LF-1)/2,dimPSD/2-(DIM.N_LF-1)/2:dimPSD/2+(DIM.N_LF-1)/2)=Wlf ; we insert the LF PSD into the matrix tmp
          tmp=double(FRACTIONAL_FFT(tmp,dimPSD,DIM.N_LF_PADDED,dimPSD/2,dimPSD/2,dimPSD/2,dimPSD/2))
        endelse
        tmp=tmp(dimPSD/2-DIM.N_OTF/2:dimPSD/2+DIM.N_OTF/2-1,dimPSD/2-DIM.N_OTF/2:dimPSD/2+DIM.N_OTF/2-1)/tmp(dimPSD/2,dimPSD/2)*varLF ; we keep the DIM.N_OTF part and renormalize with variance
        dphiLF=2*(varLF-tmp)
        dphiLF=dphiLF-min(dphiLF)
      endif else dphiLF=0
      ;high frequency structure function
      tmp=double(MATHFT(Whf,dx=DIM.DXP,ic=DIM.N_OTF/2,jc=DIM.N_OTF/2,/inverse))
      tmp=tmp/tmp(DIM.N_OTF/2,DIM.N_OTF/2)*var(0)
      dphiHF=2*(var(0)-tmp)
    endelse
    dphi=dphiLF+dphiHF
    if keyword_set(SCINTILLATION) then dphi=dphi+dphiamp
  endif

  ;INDEPENDANT POST TIP-TILT CORRECTION [OPTION]
  ;---------------------------------------------
  if keyword_set(POST_TIPTILT) then begin

    ;first task is to compute the G-tilt covariance Zernike matrix
    ;note that we compute both tilt components because
    ;the covariance in x and y are not necessarily the same
    ;for instance <a2*a8> can be different from <a3*a7>
    covmat=dblarr(8,8)
    ;we start with the low spatial frequency part
    if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
      Q02=ZERFT( 2,fpcoolf*TSC.DEXTMAX*0.5,transpose(fpcoolf)*TSC.DEXTMAX*0.5) ; and Q03=transpose(Q02)
      Q08=ZERFT( 8,fpcoolf*TSC.DEXTMAX*0.5,transpose(fpcoolf)*TSC.DEXTMAX*0.5) ;     Q07=transpose(Q08)
      Q16=ZERFT(16,fpcoolf*TSC.DEXTMAX*0.5,transpose(fpcoolf)*TSC.DEXTMAX*0.5) ;     Q17=transpose(Q16)
      Q30=ZERFT(30,fpcoolf*TSC.DEXTMAX*0.5,transpose(fpcoolf)*TSC.DEXTMAX*0.5) ;     Q29=transpose(Q30)
      ;covariance matrix elements
      covmat(0,0)=total(Q02*Wlf*conj(Q02))*DIM.DFP_LF^2
      covmat(1,1)=total(transpose(Q02)*Wlf*conj(transpose(Q02)))*DIM.DFP_LF^2
      covmat(0,3)=total(Q02*Wlf*conj(Q08))*DIM.DFP_LF^2
      covmat(3,3)=total(Q08*Wlf*conj(Q08))*DIM.DFP_LF^2
      covmat(1,2)=total(transpose(Q02)*Wlf*conj(transpose(Q08)))*DIM.DFP_LF^2
      covmat(2,2)=total(transpose(Q08)*Wlf*conj(transpose(Q08)))*DIM.DFP_LF^2
      covmat(0,4)=total(Q02*Wlf*conj(Q16))*DIM.DFP_LF^2
      covmat(3,4)=total(Q08*Wlf*conj(Q16))*DIM.DFP_LF^2
      covmat(4,4)=total(Q16*Wlf*conj(Q16))*DIM.DFP_LF^2
      covmat(1,5)=total(transpose(Q02)*Wlf*conj(transpose(Q16)))*DIM.DFP_LF^2
      covmat(2,5)=total(transpose(Q08)*Wlf*conj(transpose(Q16)))*DIM.DFP_LF^2
      covmat(5,5)=total(transpose(Q16)*Wlf*conj(transpose(Q16)))*DIM.DFP_LF^2
      covmat(0,7)=total(Q02*Wlf*conj(Q30))*DIM.DFP_LF^2
      covmat(3,7)=total(Q08*Wlf*conj(Q30))*DIM.DFP_LF^2
      covmat(4,7)=total(Q16*Wlf*conj(Q30))*DIM.DFP_LF^2
      covmat(7,7)=total(Q30*Wlf*conj(Q30))*DIM.DFP_LF^2
      covmat(1,6)=total(transpose(Q02)*Wlf*conj(transpose(Q30)))*DIM.DFP_LF^2
      covmat(2,6)=total(transpose(Q08)*Wlf*conj(transpose(Q30)))*DIM.DFP_LF^2
      covmat(5,6)=total(transpose(Q16)*Wlf*conj(transpose(Q30)))*DIM.DFP_LF^2
      covmat(6,6)=total(transpose(Q30)*Wlf*conj(transpose(Q30)))*DIM.DFP_LF^2
    endif
    ;and now the high spatial frequency part
    Q02=ZERFT( 2,fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5) ; and Q03=transpose(Q02)
    Q08=ZERFT( 8,fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5) ;     Q07=transpose(Q08)
    Q16=ZERFT(16,fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5) ;     Q17=transpose(Q16)
    Q30=ZERFT(30,fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5) ;     Q29=transpose(Q30)
    ;covariance matrix elements HF added to LF ones
    if strlowcase(strtrim(MODE,2)) eq 'seli' then Whf=Watm
    covmat(0,0)=covmat(0,0)+total(Q02*Whf*conj(Q02))*DIM.DFP^2
    covmat(1,1)=covmat(1,1)+total(transpose(Q02)*Whf*conj(transpose(Q02)))*DIM.DFP^2
    covmat(0,3)=covmat(0,3)+total(Q02*Whf*conj(Q08))*DIM.DFP^2
    covmat(3,3)=covmat(3,3)+total(Q08*Whf*conj(Q08))*DIM.DFP^2
    covmat(1,2)=covmat(1,2)+total(transpose(Q02)*Whf*conj(transpose(Q08)))*DIM.DFP^2
    covmat(2,2)=covmat(2,2)+total(transpose(Q08)*Whf*conj(transpose(Q08)))*DIM.DFP^2
    covmat(0,4)=covmat(0,4)+total(Q02*Whf*conj(Q16))*DIM.DFP^2
    covmat(3,4)=covmat(3,4)+total(Q08*Whf*conj(Q16))*DIM.DFP^2
    covmat(4,4)=covmat(4,4)+total(Q16*Whf*conj(Q16))*DIM.DFP^2
    covmat(1,5)=covmat(1,5)+total(transpose(Q02)*Whf*conj(transpose(Q16)))*DIM.DFP^2
    covmat(2,5)=covmat(2,5)+total(transpose(Q08)*Whf*conj(transpose(Q16)))*DIM.DFP^2
    covmat(5,5)=covmat(5,5)+total(transpose(Q16)*Whf*conj(transpose(Q16)))*DIM.DFP^2
    covmat(0,7)=covmat(0,7)+total(Q02*Whf*conj(Q30))*DIM.DFP^2
    covmat(3,7)=covmat(3,7)+total(Q08*Whf*conj(Q30))*DIM.DFP^2
    covmat(4,7)=covmat(4,7)+total(Q16*Whf*conj(Q30))*DIM.DFP^2
    covmat(7,7)=covmat(7,7)+total(Q30*Whf*conj(Q30))*DIM.DFP^2
    covmat(1,6)=covmat(1,6)+total(transpose(Q02)*Whf*conj(transpose(Q30)))*DIM.DFP^2
    covmat(2,6)=covmat(2,6)+total(transpose(Q08)*Whf*conj(transpose(Q30)))*DIM.DFP^2
    covmat(5,6)=covmat(5,6)+total(transpose(Q16)*Whf*conj(transpose(Q30)))*DIM.DFP^2
    covmat(6,6)=covmat(6,6)+total(transpose(Q30)*Whf*conj(transpose(Q30)))*DIM.DFP^2
    ;building the other triangle of the covariance matrix
    tmp=covmat+transpose(covmat) ; the covariance matrix is symetric / diagonal
    for i=0,7 do tmp(i,i)=0.5*tmp(i,i)
    covmat=tmp

    ;here we build the Fourier transform of the Zernike products that
    ;appears in the G-tilt structure function formula
    a0202=ZPZQTOZJ(02,02) & Q0202=0
    for i=0,n_elements(a0202.j)-1 do Q0202=Q0202+a0202.aj(i)*ZERFT(a0202.j(i),fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5)
    a0208=ZPZQTOZJ(02,08) & Q0208=0
    for i=0,n_elements(a0208.j)-1 do Q0208=Q0208+a0208.aj(i)*ZERFT(a0208.j(i),fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5)
    a0216=ZPZQTOZJ(02,16) & Q0216=0
    for i=0,n_elements(a0216.j)-1 do Q0216=Q0216+a0216.aj(i)*ZERFT(a0216.j(i),fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5)
    a0230=ZPZQTOZJ(02,30) & Q0230=0
    for i=0,n_elements(a0230.j)-1 do Q0230=Q0230+a0230.aj(i)*ZERFT(a0230.j(i),fpcoohf*TSC.DEXTMAX*0.5,transpose(fpcoohf)*TSC.DEXTMAX*0.5)

    ;the alpha coefficients - see theory
    if keyword_set(TILT_ANGLE_STEP) then begin ; this is with quantification of the tilt correction, by TILT_ANGLE_STEP steps
      covxGxG=covmat(0,0)+2*sqrt(2)*covmat(0,3)+2*covmat(3,3)
      covyGyG=covmat(1,1)+2*sqrt(2)*covmat(1,2)+2*covmat(2,2)
      facX=QUANTCOVARIANCEGTILT(TILT_ANGLE_STEP,TSC.DEXTMAX,DIM.LAMBDA,TSC.DEXTMAX/r0LAM)
      facY=QUANTCOVARIANCEGTILT(TILT_ANGLE_STEP,TSC.DEXTMAX,DIM.LAMBDA,TSC.DEXTMAX/r0LAM)
      alpha0202=2*facX.f1(0)*(covmat(0,0)+sqrt(2)*covmat(0,3))-facX.f2(0)*covxGxG ; alpha[2,2]
      alpha0303=2*facY.f1(0)*(covmat(1,1)+sqrt(2)*covmat(1,2))-facY.f2(0)*covyGyG ; alpha[3,3]
      alpha0208=facX.f3(0)*(covmat(0,3)+sqrt(2)*covmat(3,3)) ; alpha[2,8]
      alpha0307=facY.f3(0)*(covmat(1,2)+sqrt(2)*covmat(2,2)) ; alpha[3,7]
    endif else begin
      alpha0202=covmat(0,0)-2*covmat(3,3)-3*covmat(4,4)-4*covmat(7,7)-2*sqrt(6)*covmat(3,4)-2*sqrt(8)*covmat(3,7)-2*sqrt(12)*covmat(4,7) ; alpha[2,2]
      alpha0208=covmat(0,3)+sqrt(2)*covmat(3,3)+sqrt(3)*covmat(3,4)+2*covmat(3,7) ; alpha[2,8]
      alpha0303=covmat(1,1)-2*covmat(2,2)-3*covmat(5,5)-4*covmat(6,6)-2*sqrt(6)*covmat(2,5)-2*sqrt(8)*covmat(2,6)-2*sqrt(12)*covmat(5,6) ; alpha[3,3]
      alpha0307=covmat(1,2)+sqrt(2)*covmat(2,2)+sqrt(3)*covmat(2,5)+2*covmat(2,6) ; alpha[3,7]
    endelse
    alpha0216=covmat(0,4)+sqrt(2)*covmat(3,4)+sqrt(3)*covmat(4,4)+2*covmat(4,7) ; alpha[2,16]
    alpha0230=covmat(0,7)+sqrt(2)*covmat(3,7)+sqrt(3)*covmat(4,7)+2*covmat(7,7) ; alpha[2,30]
    alpha0317=covmat(1,5)+sqrt(2)*covmat(2,5)+sqrt(3)*covmat(5,5)+2*covmat(5,6) ; alpha[3,17]
    alpha0329=covmat(1,6)+sqrt(2)*covmat(2,6)+sqrt(3)*covmat(5,6)+2*covmat(6,6) ; alpha[3,29]

    ;sum of the Fourier transform of the products Uij*Ap
    ;where we multiply the Qs by !dpi*(DIAMETER*0.5)^2
    ;to go back in the regular space (not-normalized)
    pupilFT=DISCFT(fpcoohf,transpose(fpcoohf),TSC.DEXTMAX,0) ; pupil Fourier transform
    tmp=alpha0202*2*          double(pupilFT*conj(Q0202)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q02)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2) +$
      alpha0303*2*transpose(double(pupilFT*conj(Q0202)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q02)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2))+$
      alpha0208*2*          double(pupilFT*conj(Q0208)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q08)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2) +$
      alpha0307*2*transpose(double(pupilFT*conj(Q0208)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q08)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2))+$
      alpha0216*2*          double(pupilFT*conj(Q0216)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q16)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2) +$
      alpha0317*2*transpose(double(pupilFT*conj(Q0216)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q16)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2))+$
      alpha0230*2*          double(pupilFT*conj(Q0230)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q30)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2) +$
      alpha0329*2*transpose(double(pupilFT*conj(Q0230)*!dpi*(TSC.DEXTMAX*0.5)^2-Q02*conj(Q30)*(!dpi*(TSC.DEXTMAX*0.5)^2)^2))

    ;inverse FT to get the average structure function
    tmp=double(MATHFT(tmp-mean(tmp),dx=DIM.DXP,ic=DIM.N_OTF/2,jc=DIM.N_OTF/2,/INVERSE))
    sfGtilt=tmp-tmp(DIM.N_OTF/2,DIM.N_OTF/2)

    ;now, we divide with the pupil auto-correlation to follow the
    ;definition of the pupil averaged structure function
    pupilAC=TSCOTF(xpcoohf/TSC.DEXTMAX,0)*!dpi*(TSC.DEXTMAX*0.5)^2 ; pupil autocorrelation (OTF times pupil surface, and we do not need to take into account the actual pupil shape here)
    mask=xpcoohf lt TSC.DEXTMAX-DIM.DXP
    sfGtilt(where(mask eq 1))=sfGtilt(where(mask eq 1))/pupilAC(where(mask eq 1))
    sfGtilt(where(mask eq 0))=0

    ;and now we remove the G-tilt SF from the overall SF
    ;to model the G-tilt correction
    dphi=(dphi-sfGtilt)*mask
    tmp=exp(-0.5*dphi)
    tmp(where(1-mask))=max(tmp)
    coomin,abs(tmp),pos,/silent
    kk=5.d/(6*xpcoohf(pos.i,pos.j)^6)
    filterOTF=exp(-xpcoohf^8*kk)

  endif

  ;FINAL MONOCHROMATIC OTF & PSF (ATM+TSC+AO)
  ;------------------------------------------
  if not keyword_set(ONLY_PSD) then begin
    tmp=dphi
    if strlowcase(PRECISION) eq 'single' then begin
      tmp=alog(complex(TSC.OTF))/alog(2.0)-0.5*tmp/alog(2.0) ; LOG2(OTF_TSC*OTF_AO)
      w=where(float(tmp) lt -1022)
      if w(0) ne -1 then tmp(w)=-1022 ; removing values of OTF that are too close to zero (to avoid underflow message)
      tototf=2.0^tmp ; final OTF
    endif else begin
      tmp=alog(TSC.OTF)/alog(2.0)-0.5*tmp/alog(2.0) ; LOG2(OTF_TSC*OTF_AO)
      w=where(double(tmp) lt -1022)
      if w(0) ne -1 then tmp(w)=-1022 ; removing values of OTF that are too close to zero (to avoid underflow message)
      tototf=2.d^tmp ; final OTF
      
      
    endelse
    if size(filterOTF,/type) ne 0 then tototf=tototf*filterOTF
    strehl=total(abs(tototf))/total(abs(TSC.OTF))*TSC.STREHL ; overall Strehl ratio
    cutfreq=sqrt(total(abs(tototf))*DIM.DFF^2/!dpi)*2*DIM.LAMBDA*1d-6/TSC.DEXTMAX<1 ; relative angular cutoff frequency
    if strlowcase(strtrim(MODE,2)) ne 'seli' then co2halo=DIM.LAMBDA*1d-6/2/max([ACTPITCH,WFS_PITCH])*rad2asec
    if keyword_set(PSF) or keyword_set(X_PSF) or keyword_set(Y_PSF) or keyword_set(EE_ANALYSIS_DISC) or $
      keyword_set(EE_ANALYSIS_SQUA) or keyword_set(EE_ANALYSIS_SLIT) then begin
      if strlowcase(PRECISION) eq 'single' then tmp=complexarr(DIM.N_PSF_USR,DIM.N_PSF_USR)
      if strlowcase(PRECISION) eq 'double' then tmp=dcomplexarr(DIM.N_PSF_USR,DIM.N_PSF_USR)
      tmp(DIM.N_PSF_USR/2-DIM.N_OTF/2:DIM.N_PSF_USR/2+DIM.N_OTF/2-1,DIM.N_PSF_USR/2-DIM.N_OTF/2:DIM.N_PSF_USR/2+DIM.N_OTF/2-1)=tototf
      if keyword_set(X_PSF) then begin
        if strlowcase(PRECISION) eq 'single' then begin
          totpsfx=float(strehl*NORMAL(abs(MATHFT(reform(total(tmp,2))*DIM.DFF,DX=DIM.DXF_USR,IC=DIM.N_PSF_USR/2,/INVERSE,/SINGLE))))
        endif else begin
          totpsfx=strehl*NORMAL(abs(MATHFT(reform(total(tmp,2))*DIM.DFF,DX=DIM.DXF_USR,IC=DIM.N_PSF_USR/2,/INVERSE)))
        endelse
        totpsfx=congrid(totpsfx,DIM.N_PSF_USR,cubic=-0.4,/center)
      endif
      if keyword_set(Y_PSF) then begin
        if strlowcase(PRECISION) eq 'single' then begin
          totpsfy=float(strehl*NORMAL(abs(MATHFT(reform(total(tmp,1))*DIM.DFF,DX=DIM.DXF_USR,IC=DIM.N_PSF_USR/2,/INVERSE,/SINGLE))))
        endif else begin
          totpsfy=strehl*NORMAL(abs(MATHFT(reform(total(tmp,1))*DIM.DFF,DX=DIM.DXF_USR,IC=DIM.N_PSF_USR/2,/INVERSE)))
        endelse
        totpsfy=congrid(totpsfy,DIM.N_PSF_USR,cubic=-0.4,/center)
      endif
      if not keyword_set(X_PSF) and not keyword_set(Y_PSF) then begin
        if strlowcase(PRECISION) eq 'single' then begin
          totpsf=float(strehl*NORMAL(abs(MATHFT(tmp,DX=DIM.DXF_USR,IC=DIM.N_PSF_USR/2,JC=DIM.N_PSF_USR/2,/INVERSE,/SINGLE))))
        endif else begin
          totpsf=strehl*NORMAL(abs(MATHFT(tmp,DX=DIM.DXF_USR,IC=DIM.N_PSF_USR/2,JC=DIM.N_PSF_USR/2,/INVERSE)))
        endelse
      endif
    endif
  endif
  

  ;==========================================================================
  ;===============END OF PSD, OTF & PSF CALCULATION SECTION==================
  ;==========================================================================

  ;===================================
  ;SAVING/PRINT RESULTS OF CALCULATION
  ;===================================
  ;
  ;PRINT THE RESULTS ON THE SCREEN OR INTO A LOG FILE [OPTIONS]
  ;------------------------------------------------------------
  if keyword_set(INFO) or keyword_set(LOGCODE) then begin
    if keyword_set(INFO) then unit=-1
    if keyword_set(LOGCODE) then openw,unit,'paola_'+LOGCODE+'.log',/get_lun
    aff:

    printf,unit,'TELESCOPE & FOCAL PLANE PARAMETERS'
    printf,unit,'----------------------------------'
    printf,unit,format='("            TELESCOPE STREHL    [1]",d13.6)',TSC.STREHL
    printf,unit,format='("             MIRROR DIAMETER    [m]",d13.6)',TSC.DEXTMAX
    printf,unit,format='("             MIRROR  SURFACE  [m^2]",d13.6)',TSC.SURF
    printf,unit,format='("       ANGULAR FIELD OF VIEW [asec]",d13.6)',DIM.FOV_PSF
    printf,unit,format='("             PSF MATRIX SIZE   [px]",i13)',DIM.N_PSF_USR
    printf,unit,format='("             OTF MATRIX SIZE   [px]",i13)',DIM.N_OTF
    printf,unit,format='("           PUPIL MATRIX SIZE   [px]",i13)',DIM.N_OTF
    printf,unit,format='("          IMAGING WAVELENGTH [mu-m]",d13.6)',DIM.LAMBDA
    printf,unit,format='("         ANGULAR PIXEL SCALE  [mas]",d13.6)',DIM.DXF_USR*1d3
    printf,unit,format='("        PSF THEORETICAL FWHM  [mas]",d13.6)',DIM.LAMBDA*1e-6/TSC.DEXTMAX*rad2asec*1e3
    printf,unit,''

    printf,unit,'ATMOSPHERIC PARAMETERS'
    printf,unit,'----------------------'
    printf,unit,format='("   SEEING AT ZENITH @ 500 NM [asec]",d13.6)',W0
    if L0 eq -1 then printf,unit,         '                 OUTER SCALE            INFINITE'
    if L0 ne -1 then printf,unit,format='("                 OUTER SCALE L0 [m]",d13.6)',double(L0)
    printf,unit,format='("                ZENITH ANGLE  [deg]",d13.6)',ZA
    printf,unit,format='("APPROX. SEEING(LAMBDA,L0,ZA) [asec]",d13.6)',0.98*DIM.LAMBDA*1e-6/r0LAM*rad2asec*aos
    printf,unit,format='("            r0 @ 500 nm @ ZA    [m]",d13.6)',r0500
    printf,unit,format='("            r0 @ LAMBDA @ ZA    [m]",d13.6)',r0LAM
    if n_params() ge 9 then printf,unit,format='("           <LAYERS ALTITUDE>    [m]",d13.6)',meanalti
    if n_params() ge 9 then printf,unit,format='("  ISOPLANATIC ANGLE @ LAMBDA [asec]",d13.6)',anisoang
    if n_params() ge 9 and not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then begin
      printf,unit,format='("           <LAYERS VELOCITY>  [m/s]",d13.6)',meanwind
      printf,unit,format='("   PHASE TIME SCALE @ LAMBDA   [ms]",d13.6)',timescal
    endif
    if keyword_set(SCINTILLATION) then printf,unit,format='("            SCINTILLATION INDEX [-]",e13.6)',sci_index
    printf,unit,''

    printf,unit,'SCIENCE INSTRUMENT MODE'
    printf,unit,'-----------------------'
    if strtrim(TSC.INST,2) eq 'IMAGER' then  printf,unit,'                                          IMAGER'
    if strtrim(TSC.INST,2) eq 'SPECTRO' then printf,unit,'                                    SPECTROSCOPY'
    printf,unit,''

    printf,unit,'AO CORRECTION MODE'
    printf,unit,'------------------'
    if strlowcase(strtrim(MODE,2)) eq 'seli' then printf,unit,'                           SEEING LIMITED, NO AO'
    if strlowcase(strtrim(MODE,2)) eq 'ngs' then printf,unit,'                       SINGLE NATURAL GUIDE STAR'
    if strlowcase(strtrim(MODE,2)) eq 'glao' then begin
      printf,unit,'       GROUND LAYER ADAPTIVE OPTICS'
      if strlowcase(GLAO_WFS.TYPE) eq 'full' then printf,unit,'             FULL WFS FOV AVERAGING'
      if strlowcase(GLAO_WFS.TYPE) eq 'edge' then printf,unit,'          EDGE OF WFS FOV AVERAGING'
      if strlowcase(GLAO_WFS.TYPE) eq 'star' then printf,unit,' GUIDE STAR CONSTELLATION AVERAGING'
    endif
    if keyword_set(POST_TIPTILT) and strlowcase(strtrim(MODE,2)) ne 'seli' then printf,unit,'                  WITH POST AO G-TILT CORRECTION'
    if keyword_set(POST_TIPTILT) and strlowcase(strtrim(MODE,2)) eq 'seli' then printf,unit,'                          WITH G-TILT CORRECTION'
    if keyword_set(TILT_ANGLE_STEP) then printf,unit,format='("   G-TILT STEP CORRECTION ANGLE [-]",e13.6)',TILT_ANGLE_STEP
    printf,unit,''

    if strlowcase(strtrim(MODE,2)) eq 'ngs' then begin
      printf,unit,'NATURAL GUIDE STAR PARAMETERS'
      printf,unit,'-----------------------------'
      if n_params() eq 18 then if WFS_NEA eq 0 then printf,unit,'                  BRIGHT STAR MODE, NO WFS NOISE'
      if n_params() eq 18 then if WFS_NEA gt 0 then printf,unit,'      NGS MAGNITUDE NOT GIVEN, SEE WFS NEA INPUT'
      if n_params() eq 20 then printf,unit,format='("                      MAGNITUDE [1]",d13.6)',NGS_MAG
      if n_params() eq 20 then printf,unit,         '                        FILTER     '+strupcase(FILTER)
      if n_params() eq 20 then printf,unit,format='("           BLK-BODY TEMPERATURE [K]",d13.6)',NGS_TEM
      printf,unit,''
    endif

    if strlowcase(strtrim(MODE,2)) eq 'glao' then if strlowcase(GLAO_WFS.TYPE) eq 'star' and n_params() eq 22 then begin
      printf,unit,'NATURAL GUIDE STARS PARAMETERS'
      printf,unit,'------------------------------'
      printf,unit,format='("                 NGS MAGNITUDES [1]",'+strtrim(string(n_elements(NGS_MAG)),2)+'d13.6)',NGS_MAG
      printf,unit,format='("FILTER ASSOCIATED TO MAGNITUDES    ",'+strtrim(string(n_elements(FILTER)),2)+'a13)',strupcase(strtrim(FILTER,2))
      printf,unit,format='("          BLK-BODY TEMPERATURES [K]",'+strtrim(string(n_elements(NGS_TEM)),2)+'d13.6)',NGS_TEM
      if max(wfsnph) le 99999.0 then printf,unit,format='("   NUMBER PHOTONS/FRAME/LENSLET [1]",'+strtrim(string(n_elements(wfsnph)),2)+'d13.6)',wfsnph
      if max(wfsnph) gt 99999.0 then printf,unit,format='("   NUMBER PHOTONS/FRAME/LENSLET [1]",'+strtrim(string(n_elements(wfsnph)),2)+'e13.6)',wfsnph
      printf,unit,''
    endif

    if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
      printf,unit,'SCIENCE TARGET POSITION / OPTICAL AXIS'
      printf,unit,'--------------------------------------'
      printf,unit,format='("              OFF-AXIS ANGLE [asec]",d13.6)',GS_ANG
      printf,unit,format='("              OFF-AXIS ANGLE [amin]",d13.6)',GS_ANG/60.d
      if GS_ANG gt 0 then printf,unit,format='(" ORIENTATION/X-AXIS (->WEST)  [deg]",d13.6)',GS_ORI
      printf,unit,''
    endif

    if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
      printf,unit,'LOOP PARAMETERS'
      printf,unit,'---------------'
      if not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then begin
        printf,unit,format='("              INTEGRATION TIME [ms]",d13.6)',gsint
        printf,unit,format='("       CONTROL SYSTEM TIME LAG [ms]",d13.6)',LAG
        printf,unit,         '                          LOOP MODE     '+strupcase(LOOP_MODE)
        if strlowcase(strtrim(LOOP_MODE,2)) eq 'closed' then begin
          printf,unit,format='("                     LOOP GAIN  [1]",d13.6)',loopgain
          printf,unit,         '               IS THE LOOP STABLE ?     '+strupcase(stability)
        endif
      endif
      printf,unit,''
    endif

    if strlowcase(strtrim(MODE,2)) ne 'seli' then begin
      printf,unit,'WAVEFRONT SENSOR PARAMETERS'
      printf,unit,'---------------------------'
      printf,unit,'              WAVEFRONT SENSOR TYPE           SH'
      if keyword_set(ANTI_ALIAS) then printf,unit,       '   WFS ANTI-ALIASING SPATIAL FILTER           ON'
      if not keyword_set(ANTI_ALIAS) then printf,unit,   '   WFS ANTI-ALIASING SPATIAL FILTER          OFF'
      printf,unit,format='("   LENSLET PITCH IN M1 PLANE    [m]",d13.6)',WFS_PITCH
      printf,unit,format='(" # OF LENSLETS / M1 DIAMETER    [1]",d13.6)',TSC.DEXTMAX/WFS_PITCH
      printf,unit,format='("  # OF LENSLETS / M1 SURFACE    [1]",d13.6)',TSC.SURF/WFS_PITCH^2
      if keyword_set(DISPERSION) then printf,unit,format='(" WFS TO DM DISPERSION FACTOR    [1]",d13.6)',DISPERSION(0)
      if not(strlowcase(strtrim(MODE,2)) eq 'glao' and n_params() eq 13) then begin
        if size(WFS_RON,/type) ne 0 then printf,unit, format='("          DETECTOR READNOISE [e/px]",d13.6)',WFS_RON
        if size(wfslam,/type) ne 0 then printf,unit,  format='("          CENTRAL WAVELENGTH   [mu]",'+strtrim(string(n_elements(wfslam)),2)+'d13.6)',wfslam
        if size(wfstau,/type) ne 0 then printf,unit,  format='("           OPTICAL BANDWIDTH   [mu]",'+strtrim(string(n_elements(wfsbdw)),2)+'d13.6)',wfsbdw
        if size(wfstau,/type) ne 0 then printf,unit,  format='("AVERAGE OPTICAL TRANSMISSION    [1]",'+strtrim(string(n_elements(wfstau)),2)+'d13.6)',wfstau
        if size(wfsnph,/type) ne 0 then if max(wfsnph) le 99999.0 then printf,unit,format='("NUMBER PHOTONS/FRAME/LENSLET    [1]",'+strtrim(string(n_elements(wfsnph)),2)+'d13.6)',wfsnph
        if size(wfsnph,/type) ne 0 then if max(wfsnph) gt 99999.0 then printf,unit,format='("NUMBER PHOTONS/FRAME/LENSLET    [1]",'+strtrim(string(n_elements(wfsnph)),2)+'e13.6)',wfsnph
        if size(nea,/type) ne 0 then printf,unit,format='("      NOISE EQUIVALENT ANGLE [asec]",'+strtrim(string(n_elements(nea)),2)+'d13.6)',nea
        if size(WFS_NEA,/type) ne 0 then printf,unit,format='("      NOISE EQUIVALENT ANGLE [asec]",d13.6)',WFS_NEA
      endif
      if strlowcase(strtrim(MODE,2)) eq 'glao' then begin
        if strlowcase(GLAO_WFS.TYPE) ne 'star' then printf,unit,format='("  GLAO WAVE-FRONT SENSOR FoV [asec]",d13.6)',GLAO_WFS.ANG
        if strlowcase(GLAO_WFS.TYPE) ne 'star' then printf,unit,format='("  GLAO WAVE-FRONT SENSOR FoV [amin]",d13.6)',GLAO_WFS.ANG/60.d
        if strlowcase(GLAO_WFS.TYPE) eq 'star' then printf,unit,format='("  GLAO WAVE-FRONT SENSOR FoV [asec]",d13.6)',2*mean(sqrt(GLAO_WFS.ANG(0,*)^2+GLAO_WFS.ANG(1,*)^2))
        if strlowcase(GLAO_WFS.TYPE) eq 'star' then printf,unit,format='("  GLAO WAVE-FRONT SENSOR FoV [amin]",d13.6)',2*mean(sqrt(GLAO_WFS.ANG(0,*)^2+GLAO_WFS.ANG(1,*)^2))/60.d
        if strlowcase(GLAO_WFS.TYPE) eq 'star' then printf,unit,format='("  NUMBER OF GLAO GUIDE STARS    [1]",i13)',n_elements(GLAO_WFS.ANG(0,*))
        if size(WFS_NEA,/type) ne 0 then printf,unit,format='("  WFS NOISE EQUIVALENT ANGLE [asec]",d13.6)',WFS_NEA
      endif
      printf,unit,''

      printf,unit,'DEFORMABLE MIRROR PARAMETERS'
      printf,unit,'----------------------------'
      if (size(DMTF))(0) eq 0 then printf,unit,'          DM TRANSFER FUNCTION TYPE      PERFECT'
      if (size(DMTF))(0) eq 2 then printf,unit,'          DM TRANSFER FUNCTION TYPE    REALISTIC'
      printf,unit,format='("           CONJUGATION ALTITUDE [m]",d13.6)',DM_HEIGHT
      printf,unit,format='("     ACTUATOR PITCH IN M1 PLANE [m]",d13.6)',ACTPITCH
      printf,unit,format='("   # OF ACTUATORS / M1 DIAMETER [1]",d13.6)',DM_NACT
      printf,unit,format='("    # OF ACTUATORS / M1 SURFACE [1]",d13.6)',TSC.SURF/ACTPITCH^2
      printf,unit,''
      printf,unit,'VARIANCE OF RESIDUAL PHASE COMPONENTS'
      printf,unit,'-------------------------------------'
      printf,unit,format='("   HIGH FREQUENCY (FITTING) [rad^2]",e13.3)',var(0)
      if var(1) ne 0 then printf,unit,format='("                   ALIASING [rad^2]",e13.3)',var(1)
      if var(2) ne 0 then printf,unit,format='("  ANISO-SERVO OR DISPERSION [rad^2]",e13.3)',var(2)
      if var(3) ne 0 then printf,unit,format='("                  WFS NOISE [rad^2]",e13.3)',var(3)
      printf,unit,format='("    TOTAL RESIDUAL VARIANCE [rad^2]",e13.3)',total(var)
      printf,unit,''
      printf,unit,'RMS OF RESIDUAL WAVEFRONT COMPONENTS'
      printf,unit,'------------------------------------'
      printf,unit,format='("      HIGH FREQUENCY (FITTING) [nm]",d13.3)',sqrt(var(0))*DIM.LAMBDA*1e3/2/!dpi
      if var(1) ne 0 then printf,unit,format='("                      ALIASING [nm]",d13.3)',sqrt(var(1))*DIM.LAMBDA*1e3/2/!dpi
      if var(2) ne 0 then printf,unit,format='("     ANISO-SERVO OR DISPERSION [nm]",d13.3)',sqrt(var(2))*DIM.LAMBDA*1e3/2/!dpi
      if var(3) ne 0 then printf,unit,format='("                     WFS NOISE [nm]",d13.3)',sqrt(var(3))*DIM.LAMBDA*1e3/2/!dpi
      printf,unit,format='("         TOTAL WAVEFRONT ERROR [nm]",d13.3)',sqrt(total(var))*DIM.LAMBDA*1e3/2/!dpi
      printf,unit,''
    endif

    if not keyword_set(ONLY_PSD) then begin
      printf,unit,'FINAL PSF CHARACTERISTICS'
      printf,unit,'-------------------------'
      printf,unit,format='("                STREHL RATIO    [1]",d13.6)',strehl
      printf,unit,format='("  RELATIVE CUTTING FREQUENCY    [1]",d13.6)',cutfreq
      if strlowcase(strtrim(MODE,2)) ne 'seli' then printf,unit,format='("     TRANSITION CORE -> HALO [asec]",d13.6)',co2halo
    endif
       
    printf,unit,'================================================'
    printf,unit,''
    if keyword_set(INFO) and keyword_set(LOGCODE) then begin
      if unit ne -1 then begin
        close,unit
        free_lun,unit
      endif
      INFO=0
      unit=-1
      goto,aff
    endif
    if unit ne -1 then free_lun,unit
  endif

end

;===============================================================================
;================= DEFINITION OF OPTICAL TURBULENCE PARAMETERS =================
;===============================================================================

;typical Paranal atmospheric turbulence conditions
;height=[48.0,162,324,649,1299,2598,5196,10392,20785]; 0 = ground level
height = [250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,3750,4000,4250,4500,4750,5000, $
          5250, 5500,  5750,  6000,  6250,  6500,  6750,  7000,  7250,  7500,7750,  8000,8250,  8500,  8750,  9000,  9250,  9500,  9750, 10000, $
          10250, 10500, 10750, 11000, 11250, 11500, 11750, 12000, 12250, 12500,12750, 13000, 13250, 13500, 13750, 14000, 14250, 14500, 14750, 15000,$
        15250, 15500, 15750, 16000, 16250, 16500, 16750, 17000, 17250, 17500,$
        17750, 18000, 18250, 18500, 18750, 19000, 19250, 24000, 24250, 24500,$
        24750,     0,     0,     0,     0,     0,     0,     0,     0,     0,$
         0,     0,     0,     0,     0,     0,     0,     0,     0,     0]
         
;Cn2 distribution

;dcn2=[53.28,1.45,3.5,9.57,10.83,4.37,6.58,3.71,6.71]
dcn2 =[5.781e-14, 2.049e-14, 9.398e-15, 1.493e-14, 7.930e-15, 1.440e-15, 3.704e-15,$
3.014e-15, 2.471e-15, 2.216e-15, 2.332e-15, 3.127e-15, 2.391e-15, 1.654e-15,$
2.205e-15, 2.345e-15, 1.744e-15, 1.169e-15, 7.218e-16, 2.748e-16, 9.485e-16,$
1.134e-15, 7.880e-16, 4.952e-16, 3.900e-16, 2.847e-16, 2.371e-16, 2.002e-16,$
1.452e-16, 6.931e-17, 9.618e-17, 1.208e-15, 2.320e-15, 3.553e-15, 4.835e-15,$
5.581e-15, 5.241e-15, 4.900e-15, 6.584e-15, 8.368e-15, 1.293e-14, 1.956e-14,$
2.618e-14, 3.276e-14, 3.368e-14, 3.460e-14, 2.645e-14, 1.296e-14, 2.826e-15,$
2.639e-15, 2.452e-15, 2.265e-15, 1.744e-15, 1.070e-15, 6.513e-16, 8.122e-16,$
9.730e-16, 1.134e-15, 8.285e-16, 3.591e-16, 6.884e-17, 9.437e-17, 1.199e-16,$
1.454e-16, 1.554e-16, 1.613e-16, 1.672e-16, 2.301e-16, 4.861e-16, 7.420e-16,$
6.457e-16, 4.873e-16, 3.289e-16, 1.983e-16, 1.379e-16, 7.751e-17, 1.713e-17,$
3.277e-17, 6.995e-17, 1.071e-16, 1.238e-16, 0.000e+00, 0.000e+00, 0.000e+00,$
0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,$
0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00, 0.000e+00,$
0.000e+00, 0.000e+00]
dcn2=dcn2/total(dcn2)

;wind velocity profile
;v=[15.0,13,13,9,9,15,25,40,21]
v = [ 0.0  ,   0.0  ,   0.0   ,  0.0   ,  0.0  ,   0.0  ,   0.0  ,   0.0  ,   0.0  ,   0.0,$
      0.0  ,   0.0  ,   0.0   ,  0.0   ,  0.0  ,   0.0   ,  0.0  ,   0.0  ,   0.0   ,  0.0,$
      0.0  ,   0.0  ,   0.0  ,  0.0   ,  0.0  ,   0.0   ,  0.0  ,   0.0   ,  0.0   ,  0.0,$
      0.0  ,   0.0  ,   0.0   ,  9.557 , 0.0   ,  0.0    , 0.0   ,  0.0    ,10.253 , 0.0,$
      0.0  ,   0.0  ,   0.0   , 16.225, 18.522, 19.789, 23.085 ,13.89 ,  0.0   ,  0.0,$
      0.0  ,   0.0  ,   0.0   ,  0.0  ,  16.663,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0,$
      0.0  ,   0.0  ,   0.0   ,  0.0  ,   0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0,$
      0.0  ,   0.0  ,   0.0   ,  0.0  ,   0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0,$
      0.0  ,   0.0  ,   0.0   ,  0.0  ,   0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0,$
      0.0  ,   0.0  ,   0.0   ,  0.0  ,   0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ,  0.0   ]
;windir=[0,10,-20,30,-60,-70,-80,90,0]*!dpi/180
windir = [  0.0,    0.0,    0.0,    0.0,    0.0 ,   0.0 ,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,$
            0.0,    0.0,    0.0,    0.0,    0.0 ,   0.0,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,$
            0.0,    0.0,    0.0,    0.0,    0.0 ,   0.0 ,   0.0,    0.0,    0.0,  160.,    0.0 ,   0.0,$
            0.0,    0.0 , 255.7,   0.0 ,   0.0  ,  0.0 ,   0.0,  209.9 ,211.3 ,209.3, 208.6 ,250.5,$
            0.0 ,   0.0 ,   0.0,    0.0 ,   0.0 ,   0.0,  235.9,   0.0,    0.0 ,   0.0,    0.0,    0.0,$
            0.0 ,   0.0 ,   0.0 ,   0.0 ,   0.0 ,   0.0 ,   0.0 ,   0.0 ,   0.0 ,  0.0  ,  0.0  ,  0.0,$
            0.0 ,   0.0 ,   0.0 ,   0.0 ,   0.0 ,   0.0  ,  0.0 ,   0.0 ,   0.0 ,   0.0  ,  0.0 ,   0.0,$
            0.0 ,   0.0 ,   0.0 ,   0.0 ,   0.0  ,  0.0  ,  0.0 ,   0.0  ,  0.0  ,  0.0  ,  0.0  ,  0.0,$
            0.0 ,   0.0 ,   0.0 ,   0.0 ]*!dpi/180
vx=v*cos(windir)
vy=v*sin(windir)
;wind=dblarr(9,2)
wind = dblarr(100,2)
wind(*,0)=vx & wind(*,1)=vy
w0=0.7d ; seeing angle @ 500 nm
L0=27.d  ; 27-m outer-scale
za=0 ; zenith angle 30 deg    30..90

;;5TH EXAMPLE====================================================================
;;================================= GLAO MODES ==================================
;;===============================================================================
;
;;we will use here a simple telescope architecture - just a monolithic mirror
mir=SEGPOS(8.0,1.2)    ;8.0..1.2
dxf=-1
n_psf=-1 ; default value = such that FoV = 8 times seeing limited PSF FWHM
lambda=1.25
wfspitch=-1 ; means that WFS pitch = r0 @ lambda
dmh=0 ; conjugation height of the DM, here pupil level
dm_params={dmtf:-1,actpitch:-1,dm_height:dmh}
wfs_params={wfs_pitch:wfspitch}
ang=20     ; science object off-axis angle [asec]
ori=60 ; science object position angle in deg/x-axis
wfs_int=10
lag=5 ; loop time lag (WFS reading + DM commands calculation) in msec
;Here we have 6 stars on a circle of radius 60 arcsec
FoVrad=360  ; [asec]       ;60..4200
;ang06 = transpose([[314*cos(2*!dpi/6*0.75),314*cos(2*!dpi/6*0.75),317*cos(2*!dpi/6*0.74),376*cos(2*!dpi/6*5.99),374*cos(2*!dpi/6*5.99),0],$
;                   [314*sin(2*!dpi/6*0.75),314*sin(2*!dpi/6*0.75),317*sin(2*!dpi/6*0.74),376*sin(2*!dpi/6*5.99),374*sin(2*!dpi/6*5.99),0]])

;ang06 = transpose([[0,0,0,0,0,0],[0,0,0,0,0,0]])
;ang06=transpose([[FoVrad*cos(2*!dpi/6*[0,1,2,3,4,5])],$
;  [FoVrad*sin(2*!dpi/6*[0,1,2,3,4,5])]])
ang06=transpose([[FoVrad*cos(2*!dpi/6*[0,5,2,4,4,5])],$
  [FoVrad*sin(2*!dpi/6*[0,5,2,4,4,5])]])
;ang06=transpose([[FoVrad*cos(2*!dpi/6*[0,0,0,0,0,0])],$
;  [FoVrad*sin(2*!dpi/6*[0,0,0,0,0,0])]])
;ang06=transpose([[FoVrad*cos(2*!dpi/6*[2,2,2,2,2,2])],$
;  [FoVrad*sin(2*!dpi/6*[2,2,2,2,2,2])]])
;ang06=transpose([[FoVrad*cos(2*!dpi/6*[1,1,1,1,1,1])],$
;    [FoVrad*sin(2*!dpi/6*[1,1,1,1,1,1])]])
print,"ang06=",(2*!dpi/6*[0,1,2,3,4,5])
glao_wfs={type:'star',ang:ang06} ; GLAO wavefront sensing structure variable
psg2=PIXMATSIZE(mir,dxf,n_psf,lambda,w0,L0,ZA,height,dcn2,wind,wfspitch,dmh,$
  ang,ori,wfs_int,lag,glao_wfs,/info)

;now we need the telescope OTF
tsc=PSFOTFTSC(mir,psg2)

;; GLAO modeling, giving a NEA for the WFS noise error
gs_weight=-1 ; all NGS of the constellation are given the same weight
wfs_nea=dblarr(6)+0.02 ; WFS Noise Equivalent Angle / NGS, in asec
print,"wfs_nea=",wfs_nea

;glao_star1=PAOLA1('glao',psg2,tsc,w0,L0,za,height,dcn2,wind,dm_params,wfs_params,$
;                  ang,ori,wfs_int,lag,'open',1,glao_wfs,gs_weight,wfs_nea,$
;                  /info,/otf,/psf,/sf,/only_psd,logcode='star1')

glao_star1=PAOLA1('glao',psg2,tsc,w0,L0,za,height,dcn2,wind,dm_params,wfs_params,$
                  ang,ori,wfs_int,lag,'open',1,glao_wfs,gs_weight,wfs_nea,$
                  /info,/otf,/psf,/sf,logcode='star1')
                  
end