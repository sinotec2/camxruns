C-- pei(執行em-cnty-p,需讀取emteds.inp)
C-- ndct=0:base case
C-- ndct=1~7:關閉各空品區排放
C-- ndct=8:關閉台中市排放(點線面全關)
C-- ndct=9:關閉台中港區排放(點線面全關)
C-- ndct=81:關閉台中市線源排放
C-- ndct=82:關閉台中市面源排放
C-- ndct=91:關閉台中港區線源排放
C-- ndct=92:關閉台中港區面源排放
C
C*** READ EMISSION PACKETS AND CREATE THE FILE for the fine/coarse grid
C**  the line/point sources are put into the the fine grid
C**  the area sources are in the coarse grid
C**  THE SOURCE ARE DIVIDED BY 5  REGION GROUPS
C
c-----注意:當網格設定有變時，記得要修改標籤175下5行左右的rx與ry
      PARAMETER (NGP=5)                                                 !OSAT
      CHARACTER*4   MEM(10), MSPEC(10,50) !CAMX4.01
      INTEGER MFID(60),MRUNID(60),MCON(10),MEND(10),
     +        MFIELD(10),MSIM(10),ISOP(10),MBLANK
      DIMENSION RAD(24),TRAF(2,0:24)
      DIMENSION EMC(40000),EMF(40000)                                   !OSAT
      DIMENSION EMLSOX (200,200), EMLNOX (200,200), EMLCO  (200,200)
      DIMENSION EMLPM  (200,200), EMLNMHC(200,200), EMLOLE (200,200)
      DIMENSION EMLPAR (200,200), EMLTOL (200,200), EMLXYL (200,200)
      DIMENSION EMLFORM(200,200), EMLALD2(200,200), EMLETH (200,200)
      DIMENSION EMLISOP(200,200), EMLNR(200,200)
      DIMENSION EMLETHA(200,200), EMLMEOH(200,200), EMLETOH(200,200)
      DIMENSION EMLIOLE(200,200), EMLTERP(200,200), EMLALDX(200,200)
      DIMENSION EMLPRPA(200,200), EMLBENZ(200,200), EMLETHY(200,200)
      DIMENSION EMLACET(200,200), EMLKET(200,200)

      DIMENSION EMPSOX (200,200), EMPNOX (200,200), EMPCO  (200,200)
      DIMENSION EMPPM  (200,200), EMPNMHC(200,200), EMPOLE (200,200)
      DIMENSION EMPPAR (200,200), EMPTOL (200,200), EMPXYL (200,200)
      DIMENSION EMPFORM(200,200), EMPALD2(200,200), EMPETH (200,200)
      DIMENSION EMPISOP(200,200), EMPNR(200,200)
      DIMENSION EMPETHA(200,200), EMPMEOH(200,200), EMPETOH(200,200)
      DIMENSION EMPIOLE(200,200), EMPTERP(200,200), EMPALDX(200,200)
      DIMENSION EMPPRPA(200,200), EMPBENZ(200,200), EMPETHY(200,200)
      DIMENSION EMPACET(200,200), EMPKET(200,200)

      DIMENSION EMASOX (200,200), EMANOX (200,200), EMACO  (200,200)
      DIMENSION EMAPM0 (200,200)
      DIMENSION EMAPM  (200,200), EMANMHC(200,200), EMAOLE (200,200)
      DIMENSION EMAPAR (200,200), EMATOL (200,200), EMAXYL (200,200)
      DIMENSION EMAFORM(200,200), EMAALD2(200,200), EMAETH (200,200)
      DIMENSION EMAISOP(200,200), EMANR(200,200)
      DIMENSION EMAETHA(200,200), EMAMEOH(200,200), EMAETOH(200,200)
      DIMENSION EMAIOLE(200,200), EMATERP(200,200), EMAALDX(200,200)
      DIMENSION EMAPRPA(200,200), EMABENZ(200,200), EMAETHY(200,200)
      DIMENSION EMAACET(200,200), EMAKET(200,200)

      dimension nx(300),ny(300)
      DIMENSION BIOOLE (200,200), BIOALD2(200,200), BIOPAR (200,200)
      DIMENSION BIOISOP(200,200)
      DIMENSION EMPNH3 (200,200),EMANH3 (200,200),EMLNH3 (200,200)
      DIMENSION LAND(200,200),RAT_L_NO(2)
      DIMENSION SPEC(20),FACT(4,7,0:24)

C**   THE LINE SOURCE VOC ARRAY
      PARAMETER (IX0=81,IY0=2421)
      PARAMETER (LX=351,LY=2802)
      PARAMETER (DX=3,DY=3)     ! 每三公里一格
      PARAMETER (NCV=20)         ! 碳鍵物種個數
      PARAMETER (NCP=4)         ! 其他物種個數
      PARAMETER (M=(LX-IX0)/DX)
      PARAMETER (N=(LY-IY0)/DY) ! THE X-Y MATRIX (FINE MESH)
      PARAMETER (LTYP=4)        ! 國/省/縣/鄉
      REAL      VOCL(M,N,LTYP,NCV) !g-mole/hr
      REAL      POLL(M,N,LTYP,NCP )
      dimension idct(IX0:LX,IY0:LY)
      common /dct/idct,ndct
C**   面源VOC矩陣
      parameter (Narea=32236)
      dimension IE(Narea),IN(Narea),S(20,Narea)

C**   THE GROUND LEVEL POINT SOURCES ARRAY
      REAL      QA(M,N,24)      !  (SNCP g/Hr; VOC 1-20 g-mole/Hr)
      REAL      QG(M,N,24,24) !  (SNCP g/Hr; VOC 1-20 g-mole/Hr)

c      COMMON /CORRFAC/ CORRH(35,56),CORRN(35,56) ! CORRECTION FACTOR
      COMMON /NOXYZGc/ NOXG,NOYG,NOZG,XORG, YORG, DELTAX, DELTAY
      COMMON /NOXYZGf/ NOX2,NOY2,NOZf,XORGf,YORGf,DELTAXf,DELTAYf
      COMMON /TIMEVAL/ NBD, NBT, NED, NET, TBEG, TEND
      REAL NOX,NMHC
      INTEGER NOXG,NOYG,NOZG
      CHARACTER*80 IPATH
      CHARACTER*7 NMGRP(NGP)
      COMMON /OSAT/NSG,MSPEC,L                                          !OSAT
C
      LOGICAL LREST, LSINK, LPTS, LRWY, LTEMP, LTERR, LCVAR
      integer imonth
      character*2 amm
      integer days(12)  !各月日數
      real obsOZ5(24) !zuoying obs/Z5
C
      character*60 pinp,linp_pol,linp_voc,ainp,fort013,fort113  !輸出入檔名
      real rpsnv(4),rpsnvi(4)!rpm,rnox,rvoc,rsox  !排放量調整係數
      real specm(40000)  !用於調整PAR, OLE
c
      DATA MEM  / 1HE, 1HM, 1HI, 1HS, 1HS, 1HI, 1HO, 1HN, 1HS, 1H /
      DATA MCON / 1HC, 1HO, 1HN, 1HT, 1HR, 1HO, 1HL, 1H , 1H , 1H /
      DATA MEND / 1HE, 1HN, 1HD, 1H , 1H , 1H , 1H , 1H , 1H , 1H /
      DATA ISOP / 1HI, 1HS, 1HO, 1HP, 1H , 1H , 1H , 1H , 1H , 1H /
      DATA MBLANK /1H  /
C
      DATA NIN  /2 /
      DATA NOU  /3 /
      DATA NUEM1/ 43/
      DATA NUEM0/ 44/
C
      DATA NMGRP/'kauscty', 'kauscnt','pintcnt','tainnan', 'otherse'/
      data days /31,28,31,30,31,30,31,31,30,31,30,31/

C***********************************************************************
C CBM(1)=OLE    C VOC(1)=ISOPRENE
C CBM(2)=PAR    C VOC(2)=A-PINENE
C CBM(3)=ALD    C VOC(3)=MONOTERPENE2
C CBM(4)=ISO    C VOC(4)=UNIDENTIFIEDP
C CBM(5)=NR
      REAL  CBM(4,5),mw(4)

      DO 11 MM=1,4
      DO 11 LL=1,5
11    CBM(MM,LL)=0.0
      CBM(1,4)=1.0
      CBM(2,1)=0.5
      CBM(2,2)=6.0
      CBM(2,3)=1.5
      CBM(3,1)=0.5
      CBM(3,2)=6.0
      CBM(3,3)=1.5
      CBM(4,1)=0.5
      CBM(4,2)=8.5
      CBM(4,5)=0.5

C***********************************************************************
C
C-- OPEN INPUT AND OUTPUT FILES
C
      OPEN(NIN,FILE='EMISSION.INP',STATUS='OLD')
      OPEN(NOU,FILE='EMISSION.CHK',STATUS='UNKNOWN')
      WRITE (NOU,5000)
c      call GETARG (1,Arg)
C**  讀取細網格之座標、網格數等資訊
      OPEN(1,FILE='emteds.inp',STATUS='OLD')
      read(1,*)i1
      read(1,*)i2
      read(1,*)j1
      read(1,*)j2
      read(1,*)nmesh
      read(1,*)imonth,ndct  !模擬月份，用於生物源,
      write(*,*)'模擬月份:',imonth
      read(1,*)pinp
      write(*,*)'點源輸入檔:',pinp
      read(1,*)linp_pol
      write(*,*)'線源輸入檔:',linp_pol
      read(1,*)linp_voc
      write(*,*)'線源輸入檔:',linp_voc
      read(1,*)ainp
      write(*,*)'面源輸入檔:',ainp
      read(1,*)fort013
      write(*,*)'輸出檔:',fort013
      read(1,*)fort113
      write(*,*)'輸出檔:',fort113
      read(1,*)rpsnvi(1) !rpm
      write(*,*)'PM調整係數:',rpsnvi(1)
      read(1,*)rpsnvi(3) !rnox
      write(*,*)'NOx調整係數:',rpsnvi(3)
      read(1,*)rpsnvi(4) !rvoc
      write(*,*)'VOC調整係數:',rpsnvi(4)
      read(1,*)rpsnvi(2) !rsox
      write(*,*)'SOx調整係數:',rpsnvi(2)
      read(1,*)rspecm
      write(*,*)'VOC物種調整係數:',rspecm
      close(1)

      if(imonth.lt.1.or.imonth.gt.12)then  !check
        write(*,*)'Month error.',imonth
        stop
      endif

      OPEN(NUEM0, FILE=fort013,FORM='UNFORMATTED'  
     +  ,convert='BIG_ENDIAN',STATUS='UNKNOWN')
      OPEN(NUEM1, FILE=fort113,FORM='UNFORMATTED'  
     +  ,convert='BIG_ENDIAN',STATUS='UNKNOWN')
C
C--READ CONTROL PACKET
C----PACKET HEADER
      READ  (NIN,5100) MFIELD
      WRITE (NOU,5200) MFIELD
      DO 20 I=1,10
         IF (MFIELD(I).EQ.MCON(I)) GOTO 20
         WRITE (NOU,6110) MCON
         GOTO 900
  20  CONTINUE
C
C----FILE NAME
      READ  (NIN,5100) MFIELD
      WRITE (NOU,5200) MFIELD
C      DO 30 I=1,10
C         IF (MFIELD(I).EQ.MEM(I)) GOTO 30
C         WRITE (NOU,5300) MEM
C        GOTO 900
C  30  CONTINUE
C
C----FILE IDENTIFIER
      READ  (NIN,5400) MFID
      WRITE (NOU,5500) MFID
      READ  (NIN,5600) NSG,NOSPEC
      READ  (NIN,5600) NBD, NBT, NED, NET
      WRITE (NOU,5700) NBD, NBT, NED, NET
      TBEG = FLOAT(NBT/100) + FLOAT(MOD(NBT,100))/60.
      TEND = FLOAT(NET/100) + FLOAT(MOD(NET,100))/60.
C----REGION
      READ(NIN,*) XUTM,   YUTM,   NZONE
      READ(NIN,*) XORG,   YORG
      READ(NIN,*) DELTAX, DELTAY
      READ(NIN,*) NOXG,   NOYG,   NOZG
      write(*,*)nox2,noy2,noxg,noyg
      noxc= noxg                                                        !coarse
      noyc= noyg                                                        !coarse
      READ(NIN,*) NVLOW,  NVUP,   DZSURF,  DZMINL,  DZMINU
      READ(NIN,*) ISEGX,  ISEGY,  ISEGX2,  ISEGY2
      DO 50 I = 1,NOSPEC
        READ  (NIN, 5800) (MSPEC(J,I),J=1,10)
        WRITE (NOU, 5100) (MSPEC(J,I),J=1,10)
  50  CONTINUE
C
C----END
      READ  (NIN,5100) MFIELD
      WRITE (NOU,5200) MFIELD
      DO 60 I=1,10
        IF (MFIELD(I).EQ.MEND(I)) GOTO 60
        WRITE (NOU,6110) MEND
        GOTO 900
  60  CONTINUE
c      READ(NIN,*) XORG,   YORG
c      READ(NIN,*) DELTAX, DELTAY
c      READ(NIN,*) NOXG,   NOYG,   NOZG
c      READ  (NIN,5600) NBD, NBT, NED, NET
      write(*,*)nox2,noy2,noxg,noyg
C
C======================================================================
C--READ THE TIME VARIATION FACTOR FOR EACH BIOGENIC VOC CHEM. SPEC.
      CALL VOCTIME(FACT)
C--READ THE LAND USE DATA
      CALL READLAND(LAND)

C--READ THE TRAFIC VOL. RATIO (VEL/DAY * TRAF/100.VEL/HR OR SUM(TRAF(I))=100)
      OPEN(1,FILE='TRAF2.VOL',STATUS='OLD')
100   READ(1,*,END=110)I,(TRAF(K,I),K=1,2)!1:local way/ 2:highway
      GOTO 100
110   CLOSE(1)

C
C*** NOW WRITE THE EMISSIONS FILE
C
      WRITE (NOU,5900)
C
C--FILE DESCRIPTION HEADER RECORD
C

      WRITE (NUEM0)    MEM, MFID, NSG, NOSPEC, NBD, TBEG, NED, TEND
      WRITE (NUEM1)    MEM, MFID, NSG, NOSPEC, NBD, TBEG, NED, TEND
C
C-- WRITE BOUNDARY REGION HEADER
C

      WRITE(*,*)'Fine',I1,I2,J1,J2,NMESH
      nox2=(I2-I1+1)*NMESH+2                                            !FINE
      noy2=(J2-J1+1)*NMESH+2                                            !FINE
      NOXY=NOX2*NOY2                                                    !FINE
C**i1,i2,j1,j2,nz,mesh 160       220       2460      2540      3
C**   IFINE/JFINE 為舊網格系統中，CHILD網格位置的起點(I,J)
      IFINE=I1
      JFINE=J1
      XORGf=(I1-1)*DELTAX+XORG-2*DELTAX/NMESH
      YORGf=(J1-1)*DELTAY+YORG-2*DELTAY/NMESH
      DELTAXC=DELTAX
      DELTAYC=DELTAY
      DELTAXF=DELTAX/NMESH
      DELTAYF=DELTAY/NMESH
      write(*,*)nox2,noy2,noxg,noyg
      noxc= noxg                                                        !coarse
      noyc= noyg                                                        !coarse
      write(*,*)noxC,noyC
      NOXYc=NOXc*NOYc                                                   !coarse

C**   寫入表頭
      WRITE (NUEM0)XUTM,YUTM,NZONE,XORG,YORG,DELTAX  ,DELTAY  ,        !COARSE
     +             NOXc,NOYc,NOZG, NVLOW,NVUP,DZSURF,DZMINL,DZMINU      !COARSE
      WRITE (NUEM1)XUTM,YUTM,NZONE,XORGF,YORGF,DELTAXF,DELTAYF    ,     !FINE
     +             NOX2,NOY2,NOZG, NVLOW,NVUP,DZSURF,DZMINL,DZMINU      !FINE
C**   TRADITION CONVENTION
C
C-- WRITE SEGMENT DEFN RECORD FROM BOUNDARY AND CALC MAXDIM
C
      WRITE (NUEM0)1,1,NOXc,NOYc                                        !COARSE
      WRITE (NUEM1)1,1,NOX2,NOY2                                        !FINE
C
C-- NOW WRITE THE SPECIES DEFN RECORD ON BOUNDARY
C
      WRITE (NUEM0)((MSPEC(I,J),I=1,10),J=1,NOSPEC)
      WRITE (NUEM1)((MSPEC(I,J),I=1,10),J=1,NOSPEC)
      WRITE (*, 5100) (MSPEC(J,5),J=1,10)
c      WRITE (*, '(A3)') Arg
C
C--READ THE AREA GROUND LEVEL SOURCE
C-----------80AEMISN.TXT----1997,7,16(v2.1)--
C      UTME,UTMN,QSOx,QNOx,QCO,QNMHC,QPM
C          NMHC DATA --> NEGLECT
C--------------------------------------------
C*************************************
C** SKIP AREA  SOURCE
C
!       GOTO 130
C*************************************
      
        WRITE(*,*)' READING AREA EMIS DATA'
        open(1,file='dict.txt',STATUS='OLD')
        NDICT=0
        idct(:,:)=99
        read(1,*)
101     read(1,'I3,5X,I4,4X,I4',end=102)i,j,K
        K=K/100
        if((K-1)*(K-11)*(K-31)*(K-32).eq.0) idct(i,j)=1
        if((K-12)*(K-33)*(K-35).eq.0) idct(i,j)=2
        if((K-17)*(K-36)*(K-37)*(K-38).eq.0) idct(i,j)=3
        if((K-21)*(K-22)*(K-39)*(K-40)*(K-41).eq.0) idct(i,j)=4
        if((K-2)*(K-42)*(K-43).eq.0) idct(i,j)=5
        if((K-34).eq.0) idct(i,j)=6
        if((K-45)*(K-46).eq.0) idct(i,j)=7
        if(ndct.eq.8.or.ndct.eq.81.or.ndct.eq.82.or.ndct.eq.9.or.
     +   ndct.eq.91.or.ndct.eq.92)then
         if((K-17)*(K-36).eq.0) idct(i,j)=8 !台中市
        endif

cc        print*,i,j,idct
cc        if(NDICT.eq.0) then
cc                NDICT=NDICT+1
cc                NX(NDICT)=I
cc                NY(NDICT)=J
cc                GOTO 101
cc        endif
cc        do k=1,NDICT
cc                if(NX(K).eq.I.and.NY(K).eq.J)goto 101
cc        enddo
cc                NDICT=NDICT+1
cc                NX(NDICT)=I
cc                NY(NDICT)=J
                GOTO 101
102     close(1)
130   continue
C*************************************
C--READ THE LINE GROUND LEVEL SOURCE
C**  2.順序為 PM,SOX,NOX,CO (kg/day)
C
C**  SKIP LINE  SOURCE
!      GOTO 151
      WRITE(*,*)' READING LINE EMIS DATA'
      OPEN(1,FILE='../line/'//linp_pol,FORM='UNFORMATTED',
     +   STATUS='old')
      READ(1,END=151)POLL
        KgpD2gph= 1000./24. !Kg/day -> g/Hr QQQ
        summ=0
      DO J =1,N
      DO I =1,M
      DO LT=1,LTYP
      DO L =1,NCP 
        POLL(I,J,LT,L)=POLL(I,J,LT,L) * KgpD2gph
      ENDDO
        summ=summ+POLL(I,J,LT,3)
      ENDDO
      ENDDO
      ENDDO
151   CLOSE(1)
        print*,summ/KgpD2gph/1000.*366.
c       stop

C*************************************

C-----------AVOCSPEC.TXT----1997,6,27---------
C     VOC speciate -- area Source emission
C           UTME,UTMN,(VOC(I),I=1,9)
C                      (g-mole/hr)
C--------------------------------------------
C*************************************
C** not SKIP AREA VOC  SOURCE
C
!      GOTO 170
C*************************************
      WRITE(*,*)' READING AREA EMIS DATA'
      OPEN(1,FILE='../area/'//ainp,
     +  FORM='UNFORMATTED',STATUS='old')
      READ(1,END=170)QA
      CLOSE(1)

      OPEN(1,FILE='areaCHK.DAT')
c      write(*,*)'H01',sum(qa(:,:,5:24))  !!

C**   THE UNIT IS g/HR-PNSC |gmole/HR other
      DO IX=1,M
      DO IY=1,N
        I=IX
        J=IY

C---pei(關閉面源)
        KmX=(I-1)*DX+IX0
        KmY=(J-1)*DY+IY0
        if(ndct.eq.82.and.idct(KmX,KmY).eq.8) qa(i,j,:)=0. !台中市面源        
        KmXR=(KmY-2161.6)/2.5861 !台中港區X邊界      
        if(ndct.eq.92)then !台中港區面源
          if(KmX.le.KmXR.and.KmY.ge.2677.and.KmY.le.2690) qa(i,j,:)=0.          
        endif        
C---pei

c-----加入新增排放源
c         if(i.eq.47.and.j.eq.117)then  !中科一二期
c           qa(i,j,5)=qa(i,j,5)+0.04
c           qa(i,j,6)=qa(i,j,6)+18691.1
c           qa(i,j,7)=qa(i,j,7)+115.24
c           qa(i,j,8)=qa(i,j,8)+106.35
c           qa(i,j,9)=qa(i,j,9)+40.74
c           qa(i,j,10)=qa(i,j,10)+1077.25
c           qa(i,j,11)=qa(i,j,11)+1.16
c           qa(i,j,12)=qa(i,j,12)+0.0
c           qa(i,j,13)=qa(i,j,13)+150.66
c         endif
c         if(i.eq.50.and.j.eq.121)then  !后里基地
c           qa(i,j,5)=qa(i,j,5)+0.0
c           qa(i,j,6)=qa(i,j,6)+4549.725
c           qa(i,j,7)=qa(i,j,7)+28.053
c           qa(i,j,8)=qa(i,j,8)+25.888
c           qa(i,j,9)=qa(i,j,9)+9.917
c           qa(i,j,10)=qa(i,j,10)+262.229
c           qa(i,j,11)=qa(i,j,11)+0.0
c           qa(i,j,12)=qa(i,j,12)+0.0
c           qa(i,j,13)=qa(i,j,13)+39.739
c         endif
c-----
         qa(i,j,5:24)=qa(i,j,5:24)*rpsnvi(4) !rvoc  !voc排放量乘調整係數

         EMAPM0 (I,J)=QA(I,J, 1)*rpsnvi(1) !rpm  !乘調整係數
         EMANOX (I,J)=QA(I,J, 2)*rpsnvi(3) !rnox !乘調整係數
         EMASOX (I,J)=QA(I,J, 3)*rpsnvi(2) !rsox !乘調整係數
         EMACO  (I,J)=QA(I,J, 4)
         EMAOLE (I,J)=QA(I,J, 5)
         EMAPAR (I,J)=QA(I,J, 6)
         EMATOL (I,J)=QA(I,J, 7)
         EMAXYL (I,J)=QA(I,J, 8)
         EMAFORM(I,J)=QA(I,J, 9)
         EMAALD2(I,J)=QA(I,J,10)
         EMAETH (I,J)=QA(I,J,11)
         EMAISOP(I,J)=QA(I,J,12)
         EMANR  (I,J)=QA(I,J,13)
         EMAETHA(I,J)=QA(I,J,14)
         EMAMEOH(I,J)=QA(I,J,15)
         EMAETOH(I,J)=QA(I,J,16)
         EMAIOLE(I,J)=QA(I,J,17)
         EMATERP(I,J)=QA(I,J,18)
         EMAALDX(I,J)=QA(I,J,19)
         EMAPRPA(I,J)=QA(I,J,20)
         EMABENZ(I,J)=QA(I,J,21)
         EMAETHY(I,J)=QA(I,J,22)
         EMAACET(I,J)=QA(I,J,23)
         EMAKET(I,J)=QA(I,J,24)

C**     CHECKING FOR NOT ZERO EMISSION
        SUMm=0
        DO I=1,24
             SUMm=SUMm+QA(IX,IY,I) !QG(IX,IY,I,IT)
        ENDDO
        IF(SUMm.GT.0) THEN
          XJX=(iX-.5)*DX+IX0
          YJY=(IY-.5)*DY+IY0
C**       THE UNIT IS Kg-mole/HR
          WRITE(1,51)XJX*1000,YJY*1000,(QA(IX,IY,K)/1000.,K=1,24)
51        FORMAT(2F10.0,24E15.4)
        ENDIF
      ENDDO
      ENDDO
170   CLOSE(1)
c      write(*,*)'H02',sum(qa(:,:,5:24))  !!


C*************************************
C     VOC speciate -- line Source emission
C                      (g-mole/hr)
C OLE  1/ PAR  2/ TOL  3/ XYL  4/ FORM 5/ ALD2 6/ ETH  7/ ISOP 8/ NR   9/
C** not SKIP LINE voc SOURCE
C
!     GOTO 171
C*************************************
      WRITE(*,*)' READING LINE VOC EMIS DATA'
      OPEN(1,FILE='../line/'//linp_voc,
     +  FORM='UNFORMATTED',STATUS='old')
      READ(1,END=171)VOCL
171   CLOSE(1)

C---pei(關閉線源)
      DO I=1,M
      DO J=1,N  
        KmX=(I-1)*DX+IX0
        KmY=(J-1)*DY+IY0
        if(ndct.eq.81.and.idct(KmX,KmY).eq.8)then !台中市線源
         POLL(i,j,:,:)=0.
         VOCL(i,j,:,:)=0.
        endif         
        KmXR=(KmY-2161.6)/2.5861 !台中港區X邊界      
        if(ndct.eq.91)then !台中港區線源
          if(KmX.le.KmXR.and.KmY.ge.2677.and.KmY.le.2690)then
           POLL(i,j,:,:)=0.
           VOCL(i,j,:,:)=0.
          endif           
        endif
      ENDDO
      ENDDO          
C---pei

c-----!乘調整係數 
!     poll(:,:,:,1)=poll(:,:,:,1)*rpm
!     poll(:,:,:,2)=poll(:,:,:,2)*rsox
!     poll(:,:,:,3)=poll(:,:,:,3)*rnox
!     vocl=vocl*rvoc
c-----!

c      TRF=1. ! qqq 000./365. g-mole/hr
c      DO LT=1,4
c        DO J =1,NOYG
c          DO I =1,NOXG
c            DO L =1,NCV
c              VOCL(I,J,LT,L)=VOCL(I,J,LT,L)*TRF
c            ENDDO
c          ENDDO
c        ENDDO
c      ENDDO

      write(amm,'(i2.2)')imonth

      OPEN(1,FILE='./bvoc/areabio'//amm           !*************************
     +  ,STATUS='OLD')
C*************************************
c1      UTM-E座標               UTME          Character       6               東向(公尺)
c2      UTM-N座標               UTMN          Character       7               北向(公尺)
c3      總VOC排放量             TBVOC         Numeric 12      5       計算之排放量公噸／年
c4      ISOPRENE排放量          ISOPRENE      Numeric 12      5       計算之排放量公噸／年
c5      MONOTERPEN排放量        MONOTERPEN    Numeric 12      5       計算之排放量公噸／年
c6      其他VOC排放量           OTVOC         Numeric 12      5       計算之排放量公噸／年
c       ref.    d:/teds4.0/point/cbm.dat
c                       k:/lin/camx/taipeif/bio.for
cname                          saroa     num cas             mw   ole   par   tol   xyl  form  ald2   eth  isop    nr
cTERPENES                       43123      19             136.23   0.5   6.0                    1.5
cISOPRENE                       43243      88     78795    68.12                                            1.0
cA-PINENE                       43256      96     80568   136.24   0.5   6.0                    1.5
cB-PINENE                       43257      97    127913   136.24
cunidentified                                                                                           140.0    0.5   8.5                    0.5
C
C** NOT SKIP BIO VOC   SOURCE
!     GOTO 175
C--------------------------------------------
      WRITE(*,*)' READING BIO VOC EMIS DATA'
      DO MM=1,4
        DO LL=1,5
          CBM(MM,LL)=0.0
        enddo
      enddo
      CBM(1,4)=1.0
      CBM(2,1)=0.5
      CBM(2,2)=6.0
      CBM(2,3)=1.5
      CBM(3,1)=0.5
      CBM(3,2)=6.0
      CBM(3,3)=1.5
      CBM(4,1)=0.5
      CBM(4,2)=8.5
      CBM(4,5)=0.5
      DATA MW   /68.0,136.24,136.23,140.0/

c-----
c  spec:1.ISOPRENE 2.PINENE 3.TERPENES 4.unidentified
c  CBM第二維: 1.OLE 2.PAR 3.ALD2 4.ISOP 5.NR
c-----
        
C     (g-mole/hr)
c      TpY2gpHr=1000.*1000./8760.
c      READ (1,36,END=175)
c162   READ (1,*,END=175)IUTM_E,IUTM_N,TBVOC,SPEC(1),SPEC(3),OTVOC  !teds5.0

      TpY2gpHr=1000./24./real(days(imonth))  !由kg/月 換算為 g/hr
      read(1,*)
162   READ (1,*,END=175)ii,ii,IUTM_E,IUTM_N,SPEC(1),SPEC(3),OTVOC,rr

      otvoc=otvoc+rr  !rr為甲基丁烯醇MBO,2-Methyl-3-Buten-2-ol，併入other vocs
      iutm_e=iutm_e
      iutm_n=iutm_n

 36   FORMAT(I7,I8,4F12.3)
      I=(IUTM_E-IX0)/DX+1
      J=(IUTM_N-IY0)/DY+1
      IF(I*(I-M).GE.0) GOTO 162
      IF(J*(J-N).GE.0) GOTO 162
C*************************************
      SPEC(2)=OTVOC/2.
      SPEC(4)=OTVOC/2.
      do isp=1,4
        SPEC(isp)=SPEC(isp)*TpY2gpHr/ mw(isp)
      enddo

      do icbm=1,5
        summ=0
        do isp=1,4
          summ=summ+SPEC(isp)*cbm(isp,icbm)
        enddo
        if(icbm.eq.1)BIOOLE (I,J) = BIOOLE (I,J) + summ
        if(icbm.eq.2)BIOPAR (I,J) = BIOPAR (I,J) + summ
        if(icbm.eq.3)BIOALD2(I,J) = BIOALD2(I,J) + summ
        if(icbm.eq.4)BIOISOP(I,J) = BIOISOP(I,J) + summ
      enddo

      GOTO 162
175   CLOSE(1)
      open(1,file='bio.dat')
      do i=1,200
        do j=1,200
          rx=68500.+(i-1)*3000.
          ry=2364500.+(j-1)*3000.
          ii=nint(rx)
          jj=nint(ry)
          summ=BIOOLE(I,J)+BIOPAR(I,J)+BIOALD2(I,J)+BIOISOP(I,J)
          if(summ.gt.0)then
            write(1,'(i6,1x,I7,4F10.3)')ii,jj,BIOOLE(I,J),BIOPAR(I,J),
     +       BIOALD2(I,J),BIOISOP(I,J)
          endif
        enddo
      enddo
      close(1)
C**   INPUT THE AMMONIA EMISSION ( IN UNIT OF TON/YR)
      OPEN(1,FILE='../nh3/nh3.dat',STATUS='OLD')
C** NOT SKIP Ammonia   SOURCE
C
!     GOTO  174
C--------------------------------------------
      WRITE(*,*)' READING Ammonia EMIS DATA'
        EMR=1. !0.2
173   READ(1,*,END=174)X,Y,EMSNH3
      X=X/1000.
      Y=Y/1000.
      I=(X-IX0)/DX+1
      J=(Y-IY0)/DY+1
      IF(I*(I-M).GE.0) GOTO 173
      IF(J*(J-N).GE.0) GOTO 173
      EMANH3(I,J)=EMSNH3 * 1000000./ 8760. + EMANH3(I,J) * EMR
        GOTO 173
174   CLOSE(1)
C
C*** WRITE TIME VARIANT PACKET
c     NINTERV=(NED-NBD)*24+(TEND-TBEG)
      NINTERV=24                ! 1 day (date ignored)
        JBD=NBD

        JED=JBD
        TB =TBEG-1
        TE =TEND

C--READ THE POINT GROUND LEVEL SOURCE
C*************************************
C** not SKIP POINT SOURCE
C
!     GOTO 190
C*************************************
      WRITE(*,*)' READING POINT EMIS DATA'
      OPEN(1,FILE='../ptsource/'//pinp,STATUS='OLD',form='unformatted')
C---------------from pt-camx.for
C      QSOx,QNOx,QCO,QPM,QVOC(1~20) (g-mole/hr)
C----------------------------------------------------
      READ (1)QG
190   CLOSE(1)

c-----注意:點源預設已於ptsource階段調整排放量，此處不再調整
c      do i=5,24
c        qg(:,:,i,:)=qg(:,:,i,:)*rvoc  !voc乘調整係數
c      enddo

      IF(M.NE.NOXG.OR.N.NE.NOYG)then    
        write(*,*)' GROUND LEVEL POINT SOURCE NG'
        write(*,*)' No. of grid not match.'
        write(*,*)m,n,noxg,noyg
c        stop
      endif
!zuoying
      read(fort013(9:10),'(A1,I1)')IPATH(1:1),iZ1
      if(IPATH(1:1).eq.'Z'.and.iZ1.eq.6) then
        open(1,file='obsOZ5.txt')
        read(1,*)
        DO IT=1,NINTERV
          read(1,*)i,avgg,obsOZ5(IT)
        enddo
        close(1)       
      endif

      DO 180 IT=1,NINTERV
         TB=TB+1
         IF (TB.GE.24) THEN
            JBD=JBD+1
            TB=TB-24
         ENDIF
         TE=TB+1
         JED=JBD
         IF (TE.GE.24) THEN
            JED=JED+1
            TE=TE-24
         ENDIF
C--(1) TIME INTERVAL RECORD
       WRITE (NUEM0)JBD, TB, JED, TE
       WRITE (NUEM1)JBD, TB, JED, TE
       WRITE (*,*)  JBD, TB, JED, TE
C qqq
      fa=0.   
      IF(TB.GE.8..AND.TB.LT.16.) fa=3.    !working 8 hr for NOx area sources
      RK=0.   
      IF(TB.GE.8..AND.TB.LT.16.) RK=3.  !time var. for VOCs area sources
      Rp=0.   
      IF(TB.GE.8..AND.TB.LT.16.) Rp=3.  !time var. for PM   area sources
        emR=1.0
      do 192 i=1,M
      do 192 j=1,N
         DO L=1,24
           IF(QG(I,J,L,IT).LT.0) STOP 'QG NG'
         ENDDO
         EMPPM  (I,J)=QG(I,J, 1,IT) *emR!*rpm   !->點源已在ptsource調整過，不在此調整
         EMPNOX (I,J)=QG(I,J, 2,IT)     !*rnox  !乘調整係數 
         EMPSOX (I,J)=QG(I,J, 3,IT)     !*rsox  !乘 桴舕Y數 
         EMPCO  (I,J)=QG(I,J, 4,IT)
         EMPOLE (I,J)=QG(I,J, 5,IT) * RK
         EMPPAR (I,J)=QG(I,J, 6,IT) * RK
         EMPTOL (I,J)=QG(I,J, 7,IT) * RK
         EMPXYL (I,J)=QG(I,J, 8,IT) * RK
         EMPFORM(I,J)=QG(I,J, 9,IT) * RK
         EMPALD2(I,J)=QG(I,J,10,IT) * RK
         EMPETH (I,J)=QG(I,J,11,IT) * RK
         EMPISOP(I,J)=QG(I,J,12,IT) * RK
         EMPNR  (I,J)=QG(I,J,13,IT) * RK
         EMPETHA(I,J)=QG(I,J,14,IT) * RK
         EMPMEOH(I,J)=QG(I,J,15,IT) * RK
         EMPETOH(I,J)=QG(I,J,16,IT) * RK
         EMPIOLE(I,J)=QG(I,J,17,IT) * RK
         EMPTERP(I,J)=QG(I,J,18,IT) * RK
         EMPALDX(I,J)=QG(I,J,19,IT) * RK
         EMPPRPA(I,J)=QG(I,J,20,IT) * RK
         EMPBENZ(I,J)=QG(I,J,21,IT) * RK
         EMPETHY(I,J)=QG(I,J,22,IT) * RK
         EMPACET(I,J)=QG(I,J,23,IT) * RK
         EMPKET(I,J)=QG(I,J,24,IT) * RK

         EMAPM  (I,J) =            +EMAPM0(I,J)*Rp
192   continue
!ZuoYing Station Addition
      X=176.679
      Y=2509.045 
      if(IPATH(1:1).eq.'Z') then
        I=(X-IX0)/DX+1
        J=(Y-IY0)/DY+1
        IF(I*(I-M).GE.0)stop 'X_ZuoYingSt. wrong' 
        IF(J*(J-N).GE.0)stop 'Y_ZuoYingSt. wrong'
        print*,'zuoying, PAL',TB,EMAOLE(I,J),EMLOLE(I,J),EMPOLE(I,J)
        rzuoy=1.
        if(iZ1.eq.6)rzuoy=obsOZ5(IT)
        EMPNOX (I,J)=EMPNOX(I,J)+1000.**(iZ1-1)*rzuoy
        EMPOLE (I,J)=EMPOLE(I,J)+1000.**(iZ1-1)*rzuoy
        print*,'zuoying, PAL',TB,EMAOLE(I,J),EMLOLE(I,J),EMPOLE(I,J)
!0-> 1.2t/y OLE, 2->12, 3->120, 4->1,200., 5->12,000 t/y
      endif  
!ZuoYing Station Addition
C
C***********************************************************************

C** 轉換VOCL/POLL
       DO J =1,N
       DO I =1,M
         EMLPM  (I,J) = 0.
         EMLSOX (I,J) = 0.
         EMLNOX (I,J) = 0.
         EMLCO  (I,J) = 0.
         EMLOLE (I,J) = 0.
         EMLPAR (I,J) = 0.
         EMLTOL (I,J) = 0.
         EMLXYL (I,J) = 0.
         EMLFORM(I,J) = 0.
         EMLALD2(I,J) = 0.
         EMLETH (I,J) = 0.
         EMLISOP(I,J) = 0.
         EMLNR  (I,J) = 0.
         EMLETHA(I,J)= 0.
         EMLMEOH(I,J)= 0.
         EMLETOH(I,J)= 0.
         EMLIOLE(I,J)= 0.
         EMLTERP(I,J)= 0.
         EMLALDX(I,J)= 0.
         EMLPRPA(I,J)= 0.
         EMLBENZ(I,J)= 0.
         EMLETHY(I,J)= 0.
         EMLACET(I,J)= 0.
         EMLKET(I,J)= 0.
       ENDDO
       ENDDO
       DO LT=1,4
         TRF=TRAF(1,INT(TB))/100.*24.
         IF(LT.EQ.1)TRF=TRAF(2,INT(TB))/100.*24.  !highway
         rpsnv=1.
         IF(TB.GE.8..AND.TB.LE.15.) rpsnv=rpsnvi!日間狀況
       DO J =1,N
       DO I =1,M
         EMLPM  (I,J) =EMLPM  (I,J)+POLL(I,J,LT,1)*TRF*rpsnv(1)
         EMLSOX (I,J) =EMLSOX (I,J)+POLL(I,J,LT,2)*TRF*rpsnv(2)
         EMLNOX (I,J) =EMLNOX (I,J)+POLL(I,J,LT,3)*TRF*rpsnv(3)
         EMLCO  (I,J) =EMLCO  (I,J)+POLL(I,J,LT,4)*TRF
         EMLOLE (I,J) =EMLOLE (I,J)+VOCL(I,J,LT,1)*TRF*rpsnv(4)
         EMLPAR (I,J) =EMLPAR (I,J)+VOCL(I,J,LT,2)*TRF*rpsnv(4)
         EMLTOL (I,J) =EMLTOL (I,J)+VOCL(I,J,LT,3)*TRF*rpsnv(4)
         EMLXYL (I,J) =EMLXYL (I,J)+VOCL(I,J,LT,4)*TRF*rpsnv(4)
         EMLFORM(I,J) =EMLFORM(I,J)+VOCL(I,J,LT,5)*TRF*rpsnv(4)
         EMLALD2(I,J) =EMLALD2(I,J)+VOCL(I,J,LT,6)*TRF*rpsnv(4)
         EMLETH (I,J) =EMLETH (I,J)+VOCL(I,J,LT,7)*TRF*rpsnv(4)
         EMLISOP(I,J) =EMLISOP(I,J)+VOCL(I,J,LT,8)*TRF*rpsnv(4)
         EMLNR  (I,J) =EMLNR  (I,J)+VOCL(I,J,LT,9)*TRF*rpsnv(4)
         EMLETHA  (I,J) =EMLETHA  (I,J)+VOCL(I,J,LT,10)*TRF*rpsnv(4)
         EMLMEOH  (I,J) =EMLMEOH  (I,J)+VOCL(I,J,LT,11)*TRF*rpsnv(4)
         EMLETOH  (I,J) =EMLETOH  (I,J)+VOCL(I,J,LT,12)*TRF*rpsnv(4)
         EMLIOLE  (I,J) =EMLIOLE  (I,J)+VOCL(I,J,LT,13)*TRF*rpsnv(4)
         EMLTERP  (I,J) =EMLTERP  (I,J)+VOCL(I,J,LT,14)*TRF*rpsnv(4)
         EMLALDX  (I,J) =EMLALDX  (I,J)+VOCL(I,J,LT,15)*TRF*rpsnv(4)
         EMLPRPA  (I,J) =EMLPRPA  (I,J)+VOCL(I,J,LT,16)*TRF*rpsnv(4)
         EMLBENZ  (I,J) =EMLBENZ  (I,J)+VOCL(I,J,LT,17)*TRF*rpsnv(4)
         EMLETHY  (I,J) =EMLETHY  (I,J)+VOCL(I,J,LT,18)*TRF*rpsnv(4)
         EMLACET  (I,J) =EMLACET  (I,J)+VOCL(I,J,LT,19)*TRF*rpsnv(4)
         EMLKET  (I,J) =EMLKET  (I,J)+VOCL(I,J,LT,20)*TRF*rpsnv(4)
       ENDDO
       ENDDO
       ENDDO
C
C***********************************************************************
C
C--(2) FOR EACH SPECEICES
C---- DEFINE THE NOx RATIO
       DATA RAT_P_NO/0.00/
       DATA RAT_A_NO/0.00/
       DATA RAT_L_NO/0.10,0.10/

       IF(TB.GE.2..AND.TB.LT.17.) then !日間狀況
        RAT_L_NO(1)=0.1  !線源no比例
        RAT_L_NO(2)=0.1  !此參數無用
        RAT_P_NO=0.10    !點源no比例
       else ! 夜間
        RAT_L_NO(1)=0.9  !線源no比例
        RAT_L_NO(2)=0.9  !此參數無用
        RAT_P_NO=0.9     !點源no比例
       endif

C----FOR THE CASE WHEN SPECIES NAME IS NO
       L=1
       IJ=1
       DO 1010 J=1,N
       DO 1010 I=1,M
          TMP=    EMPNOX(I,J)*RAT_P_NO/46. ! nox as no2
          TMP=TMP+EMLNOX(I,J)*RAT_L_NO(1)/46.
          EMc(IJ)=TMP+EMANOX(I,J)*RAT_A_NO/46.*fa
          EMF(IJ)=0.
          IJ=IJ+1
1010   CONTINUE
       call FINCORS(IFINE,JFINE,EMF,EMC)
Cs    OR THE CASE WHEN SPECIES NAME IS NO2
       L=2
       IJ=1
c       write(*,*)it,( 180000.-XORG)/DELTAX+1,(2500000.-YORG)/DELTAY+1
       DO 1020 J=1,N
       DO 1020 I=1,M
          TMP1=EMPNOX(I,J)*(1-RAT_P_NO)/46.
          TMP2=EMLNOX(I,J)*(1-RAT_L_NO(1))/46.
                TMP3=EMANOX(I,J)*(1-RAT_A_NO)/46.*fa 
          EMc(IJ)=TMP1+TMP2+TMP3
        if(i.ge.int( (176000.-XORG)/DELTAX+1).and.
     +   I.le.int(( 194000.-XORG)/DELTAX+1).and.
     +   j.ge.int((2487000.-YORG)/DELTAY+1).and.
     +   j.le.int((2505000.-YORG)/DELTAY+1))
c     + write(NOU,*)it,EMc(ij),TMP1,TMP2,TMP3
c       if(i.eq.int( (185000.-XORG)/DELTAX+1).and.
c     +   j.eq.int((2469000.-YORG)/DELTAY+1))
     +  write(NOU,*)it,i*100+J,TMP1,TMP2,TMP3
          EMF(IJ)=0.
          IJ=IJ+1
1020   CONTINUE
       call FINCORS(IFINE,JFINE,EMF,EMC)
C----FOR THE CASE WHEN SPECIES NAME IS O3
       L=3
       EMO3=0
       DO IJ=1,N*M
         EMF(IJ)=EMO3
         EMc(IJ)=EMO3
       ENDDO
       call FINCORS(IFINE,JFINE,EMF,EMC)

C-- FOR THE CASE WHEN SPECIES NAME IS SO2
       L=4
       IJ=1
       DO 1040 J=1,N
       DO 1040 I=1,M
          TMP=    EMPSOX(I,J)/64.
          TMP=TMP+EMLSOX(I,J)/64.
          EMc(IJ)=TMP+EMASOX(I,J)/64.
          EMF(IJ)=0.
          IJ=IJ+1
1040   CONTINUE
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS NH3
       L=5
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP=EMANH3(I,J)/17.
          TMP=TMP+EMPNH3(I,J)/17.
          TMP=TMP+EMLNH3(I,J)/17.
          EMc(IJ)=TMP   *   1.0         ! nh3-try
          EMf(IJ)=0
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS PNO3
       L=6
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP=(EMANOX(I,J)+EMLNOX(I,J)+EMPNOX(I,J))*0.01 ! 1% as NO3 in g/Hr
          EMc(IJ)=TMP 
          EMf(IJ)=0
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS PSO4
       L=7
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP=(EMASOX(I,J)+EMLSOX(I,J)+EMPSOX(I,J))*0.01 ! 1% as SO2 in g/Hr
          EMc(IJ)=TMP 
          EMf(IJ)=0
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS PNH4
       L=8
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP=0.
          EMc(IJ)=TMP 
          EMf(IJ)=0
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS POA
       L=9
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP=EMLPM(I,J)*0.01 ! 1% as LPM in g/Hr
          EMc(IJ)=TMP 
          EMf(IJ)=0
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS PEC
       L=10
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP=EMLPM(I,J)*0.01 ! 1% as LPM in g/Hr
          EMc(IJ)=TMP 
          EMf(IJ)=0
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS FPRM (FINE PRIMARY)
       L=11
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP=EMPPM(I,J)*0.2+EMLPM(I,J)      !燃燒
          EMc(IJ)=TMP
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS CPRM (COARSE PRIMARY)
       L=12
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP=EMPPM(I,J)*0.2     !黑煙
          EMc(IJ)=TMP
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS FCRS (FINE CRUSTAL)
       L=13
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP=EMAPM(I,J)*0.5+EMPPM(I,J)*0.3
          EMc(IJ)=TMP
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS CCRS (COARSE CRUSTAL)
       L=14
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP=EMAPM(I,J)*0.5+EMPPM(I,J)*0.3
          EMc(IJ)=TMP
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS SOA1~5
       EMS=0
       DO IJ=1,N*M
         EMF(IJ)=EMS
         EMC(IJ)=EMS
       ENDDO
       DO L=15,19
        call FINCORS(IFINE,JFINE,EMF,EMC)
       ENDDO
C-- FOR THE CASE WHEN SPECIES NAME IS PAR
       L=20
       IJ=1
       DO J=1,N
       DO I=1,M
          LT=LAND(I,J)
          TMP    =    EMPPAR(I,J)
          TMP    =TMP+EMLPAR(I,J)
          TMP    =TMP+EMAPAR(I,J) * RK
          EMc(IJ)=TMP+BIOPAR(I,J)*FACT(2,LT,int(TB))          
c-----修改PAR(調整成OLE)
          specm(ij)=emc(ij)*(1.-rspecm)
          emc(ij)=emc(ij)*rspecm
c          write(288,*)i,j,specm(ij)
          IJ=IJ+1
       ENDDO
       ENDDO
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS ETHA
       L=21
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPETHA(I,J)
          TMP    =TMP+EMLETHA(I,J)
          EMc(IJ)=TMP+EMAETHA(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS MEOH
       L=22
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPMEOH(I,J)
          TMP    =TMP+EMLMEOH(I,J)
          EMc(IJ)=TMP+EMAMEOH(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS ETOH
       L=23
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPETOH(I,J)
          TMP    =TMP+EMLETOH(I,J)
          EMc(IJ)=TMP+EMAETOH(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS ETH
       L=24
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPETH(I,J)
          TMP    =TMP+EMLETH(I,J)
          EMc(IJ)=TMP+EMAETH(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS OLE
       L=25
       IJ=1
       DO J=1,N
       DO I=1,M
          LT=LAND(I,J)
          TMP    =    EMPOLE(I,J)
          TMP    =TMP+EMLOLE(I,J)
          TMP    =TMP+EMAOLE(I,J) * RK
          EMc(IJ)=TMP+BIOOLE(I,J)*FACT(1,LT,int(TB))
c-----OLE調升(from PAR)
          emc(ij)=emc(ij)+specm(ij)
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS IOLE
       L=26
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPIOLE(I,J)
          TMP    =TMP+EMLIOLE(I,J)
          EMc(IJ)=TMP+EMAIOLE(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS ISOP
       L=27
       IJ=1
       DO J=1,N
       DO I=1,M
          LT=LAND(I,J)
          TMP    =    EMPISOP(I,J)
          TMP    =TMP+EMLISOP(I,J)
          TMP    =TMP+EMAISOP(I,J) * RK
          EMc(IJ)=TMP+BIOISOP(I,J)*FACT(4,LT,int(TB))
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS TERP
       L=28
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPTERP(I,J)
          TMP    =TMP+EMLTERP(I,J)
          EMc(IJ)=TMP+EMATERP(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS FORM
       L=29
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPFORM(I,J)
          TMP    =TMP+EMLFORM(I,J)
          EMc(IJ)=TMP+EMAFORM(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS ALD2
       L=30
       IJ=1
       DO J=1,N
       DO I=1,M
          LT=LAND(I,J)
          TMP    =    EMPALD2(I,J)
          TMP    =TMP+EMLALD2(I,J)
          TMP    =TMP+EMAALD2(I,J) * RK
          EMc(IJ)=TMP+BIOALD2(I,J)*FACT(3,LT,int(TB))
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS ALDX
       L=31
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPALDX(I,J)
          TMP    =TMP+EMLALDX(I,J)
          EMc(IJ)=TMP+EMAALDX(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS TOL
       L=32
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPTOL(I,J)
          TMP    =TMP+EMLTOL(I,J)
          EMc(IJ)=TMP+EMATOL(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS XYL
       L=33
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPXYL(I,J)
          TMP    =TMP+EMLXYL(I,J)
          EMc(IJ)=TMP+EMAXYL(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS PRPA
       L=34
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPPRPA(I,J)
          TMP    =TMP+EMLPRPA(I,J)
          EMc(IJ)=TMP+EMAPRPA(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS BENZ
       L=35
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPBENZ(I,J)
          TMP    =TMP+EMLBENZ(I,J)
          EMc(IJ)=TMP+EMABENZ(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS ETHY
       L=36
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPETHY(I,J)
          TMP    =TMP+EMLETHY(I,J)
          EMc(IJ)=TMP+EMAETHY(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS ACET
       L=37
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPACET(I,J)
          TMP    =TMP+EMLACET(I,J)
          EMc(IJ)=TMP+EMAACET(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
C-- FOR THE CASE WHEN SPECIES NAME IS KET
       L=38
       IJ=1
       DO J=1,N
       DO I=1,M
          TMP    =    EMPKET(I,J)
          TMP    =TMP+EMLKET(I,J)
          EMc(IJ)=TMP+EMAKET(I,J) * RK
          EMF(IJ)=0.
          IJ=IJ+1
       ENDDO	
       ENDDO	
       call FINCORS(IFINE,JFINE,EMF,EMC)
180   CONTINUE

C
C     *** THATS IT -- REWIND AND LEAVE
C
      GO TO 999
C
C     *** COME HERE IF THERE HAS BEEN AN ERROR
C
  900 WRITE (NOU,6100)
      GO TO 999
C
  999 CONTINUE
      STOP
C***********************************************************************
 5100 FORMAT (10A1)
 5800 FORMAT (10A1,F10.0)
 5400 FORMAT (60A1)
 5600 FORMAT (4I10)
 5000 FORMAT (1H1)
 5200 FORMAT (1X, 10A1)
 5300 FORMAT (43H **ERROR ON ABOVE CARD -- SHOULD HAVE BEEN , 10A1)
 5500 FORMAT (1X, 60A1)
 5700 FORMAT (1X, 4I10)
 5900 FORMAT (54H1*************** EMISSION FILE ***********************)
 6000 FORMAT (32H0FILE DESCRIPTION HEADER RECORD /
     + 5X, 17HFILE TYPE      = , 10A1 /
     + 5X, 17HFILE ID        = , 60A1 /
     + 5X, 17HNO OF SEGMENTS = , I2   /
     + 5X, 17HNO OF SPECIES  = , I2   /
     + 5X, 17HBEG DATE       = , I5   /
     + 5X, 17HBEG TIME       = , F3.0 /
     + 5X, 17HEND DATE       = , I5   /
     + 5X, 17HEND TIME       = , F3.0 )
 6100 FORMAT (43H0PROGRAM TERMINATING DUE TO ABOVE ERROR(S)/
     + 5X, 27H -FIX CARD(S) AND RESUBMIT )
 6110 FORMAT (43H **ERROR ON ABOVE CARD -- SHOULD HAVE BEEN , 10A1)
C
      END
      SUBROUTINE READLAND(LAND)
      PARAMETER(XMIN=140,XMAX=258,YMIN=2420,YMAX=2608)
      PARAMETER(M=(XMAX-XMIN)/2+1,N=(YMAX-YMIN)/2+1)
      CHARACTER CLAS(13)*1,A601*60,A602*60,A603*60,A1*1,A2*1,A3*1
      DIMENSION LAND(200,200)
      DO I=1,200
      DO J=1,200
       LAND(I,J)=5
      ENDDO
      ENDDO
      RETURN
      END
C**   FACT(SPEC,LAND,TIME)
C**   SPEC 1:OLE/2:PAR/3:ALD2/4:ISOP/
C**   W(ORG,CBM);ORG 1:ISOP/2:APIN/3:TERP/4:UNDF
      SUBROUTINE VOCTIME(FACT)
      DIMENSION  FORG(4,0:24)
      DIMENSION  FACT(4,7,0:24)
      REAL PERC(4,7),MW(4),W(4,4)
C--BIOGENIC EMISSION
C--PERC(SPEC,LAND-TYPE)
       DATA PERC/8*0.,20,25,25,30,
     +  20,25,25,30,   73,3,3,21,
     +  63,5,6,26,     24,21,23,32/MW/68.11,136.24,136.24,140.0/
      DATA W/0.,0.5,0.5,0.5,  0.,6.,6.,8.5,  0.,1.5,1.5,0.,
     +       1.,0.,0.,0./
C**   READ THE TIME VARIATION OF BIOGENIC EMISSION RATE
      OPEN(1,FILE='BIOFAC.DAT',STATUS='OLD')
      READ (1,*)
      DO 90 I=1,24
         READ(1,15) (FORG(L,I),L=1,4)
 15      FORMAT (6X,4F8.2)
 90   CONTINUE
      DO L=1,4
        FORG(L,0)=FORG(L,24)
      ENDDO
      CLOSE(1)
      WRITE(3,*)' W(2,3)=',W(2,3)
      WRITE(3,*)' W(2,5)=',PERC(2,5)
      DO IT=0,24
        DO LN=1,7
          TMP1=0
          DO IORG=1,4
            TMP1=TMP1+PERC(IORG,LN)
          ENDDO
          IF(TMP1.EQ.0) THEN
            DO ICBM=1,4
              FACT(ICBM,LN,IT)=0
            ENDDO
          ELSE
           DO ICBM=1,4
            TMP1=0
            TMP2=0
            DO IORG=1,4
             TMP1=TMP1+W(IORG,ICBM)*PERC(IORG,LN)*FORG(IORG,IT)/MW(IORG)
             TMP2=TMP2+W(IORG,ICBM)*PERC(IORG,LN)/MW(IORG)
            ENDDO
            FACT(ICBM,LN,IT)=TMP1/TMP2
           ENDDO
          ENDIF
        ENDDO
        WRITE(3,*)IT,(FACT(ICBM,5,IT),ICBM=1,4)
      ENDDO
      RETURN
      END
C***********************************************************************
      SUBROUTINE FINCORS(IFINE,JFINE,EMF,EMC)
C***********************************************************************
      PARAMETER (NGP=5)                                                 !OSAT
!     CHARACTER*4   MSPEC(10,50) !CAMX4.01
      dimension     MSPEC(10,50) !CAMX4.01
      DIMENSION tmp(200,200),EMF(40000),EMC(40000)
      DIMENSION EM (40000)                                              !OSAT
      COMMON /OSAT/NSG,MSPEC,L                                          !OSAT
      COMMON /NOXYZGc/ NOXG,NOYG,NOZG,XORG, YORG, DELTAX, DELTAY
      COMMON /NOXYZGf/ NOX2,NOY2,NOZf,XORGf,YORGf,DELTAXf,DELTAYf
c      COMMON /CORRFAC/ CORRH(35,56),CORRN(35,56) ! CORRECTION FACTOR
      PARAMETER (IX0=81,IY0=2421)
      PARAMETER (LX=351,LY=2802)
      PARAMETER (DX=3,DY=3)     ! 每三公里一格
      PARAMETER (M=(LX-IX0)/DX)
      PARAMETER (N=(LY-IY0)/DY) ! THE X-Y MATRIX (FINE MESH)
      dimension idct(IX0:LX,IY0:LY)
      common /dct/idct,ndct
      DATA NUEM1/ 43/
      DATA NUEM0/ 44/
       NOXYf=nox2*noy2
C**    STORE THE LINE/POINT/AREA SOURCES INTO THE TMP ARRAY
          TMP(:,:)=0.
       IJ=1
       DO 1280 J=1,N
       DO 1280 I=1,M
        KmX=(I-1)*DX+IX0
        KmY=(J-1)*DY+IY0
        II=(KmX*1000.-XORGf)/DELTAXf+1
        JJ=(KmY*1000.-YORGf)/DELTAYf+1
        if(ndct.ne.0.and.idct(KmX,KmY).eq.ndct) EMC(IJ)=0.
C-- pei(台中港區)
        KmXR=(KmY-2161.6)/2.5861 !台中港區X邊界      
        if(ndct.eq.9)then
          if(KmX.le.KmXR.and.KmY.ge.2677.and.KmY.le.2690) EMC(IJ)=0.          
        endif
C-- pei         
        if((II*(II-NOX2).LT.0).AND.(JJ*(JJ-NOY2).LT.0))
     +   TMP(II,JJ)=EMC(IJ)
          IJ=IJ+1
1280   CONTINUE
C**    CUT THE FINE GRID INTO emF(ij) TO WRITE OUT
       IJ=1
       DO 1281 J=1,NOY2
       DO 1281 I=1,NOX2
          EMF(IJ)=TMP(I,J)
          IJ=IJ+1
1281   CONTINUE
C**   SUM INTO THE COARSE GRID SYSTEM FOR THE WHOLE REGION
       IJ=1
       DO 1283 J=1,NOYG
       DO 1283 I=1,NOXG
        X=(I-1)*DELTAX+XORG
        Y=(J-1)*DELTAY+YORG
        KmX=X/1000.
        KmY=Y/1000.
        I1=(KmX-IX0)/DX+1
        J1=(KmY-IY0)/DY+1
        SUMm=0
         DO II=I1,I1+2
         DO JJ=J1,J1+2
         if((II*(II-200).LT.0).AND.(JJ*(JJ-200).LT.0)) 
     +   SUMm=SUMm+TMP(II,JJ)
         ENDDO
         ENDDO
         EMC(IJ)=SUMm
         IJ=IJ+1
1283   CONTINUE
       NOXY=NOXG*NOYG
       WRITE(NUEM0)NSG,(MSPEC(J,L),J=1,10),(EMc(IJ),IJ=1,NOXY)      !COARSE
       WRITE(NUEM1)NSG,(MSPEC(J,L),J=1,10),(EMf(IJ),IJ=1,NOXYf)     !FINE
      RETURN
      END

      subroutine readsurf(n,a)
      parameter (ni=35,nj=56)
      dimension a(ni,nj)
      do 41 L=1,5
        read (n,*)
 41   continue
      I1=NI/10+1
      I2=MOD(NI,10)
      IF(I2.EQ.0)I1=I1-1
      do 42 j=1,nj
        K1=1
        DO KK=1,I1
          K2=K1+9
          IF(KK.EQ.I1.AND.I2.NE.0) K2=K1+I2-1
          read(n,*) (a( k,j),K=K1,K2)
          K1=K2+1
        enddo
        read(n,*)
        do k=1,ni
          if (a(k,j).lt.0) a(k,j)=0.
        enddo
 42   continue
      return
      end

c      BLOCK DATA
c      COMMON /CORRFAC/ CORRH(35,56),CORRN(35,56) ! CORRECTION FACTOR
c      DATA CORRH/1960*1./  !***網格數改變時，此二行要改
c      DATA CORRN/1960*1./  !***
c      END
