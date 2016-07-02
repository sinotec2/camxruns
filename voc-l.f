C PROGRAM VOC-Line.for
C*INPUT the file out.?? (run write.for)
C*INPUT PROFIL.TXT
C 0068 43205 57  115071   100.00  F19890105
C ----       ---          ------
C (PRO_NO.,SPE_NO.,WT%)
C FORMAT(A4,7X,I3,10X,F6.2)
C
C*
C NAME(SOURCE,TYPE)
C SOURCE=1 > PIPE-TAIL
C SOURCE=2 > EVAPOR
C SOURCE=3 > R
C TYPE  =1 > PLDGV
C TYPE  =2 > PLDDV
C TYPE  =3 > BLDGV
C TYPE  =4 > BLDLPG
C TYPE  =5 > LDGT
C TYPE  =6 > LDDT
C TYPE  =7 > BUS
C TYPE  =8 > HDGV
C TYPE  =9   HDDT
C TYPE  =10  MC2
C TYPE  =11  MC4
C

      PARAMETER (NVTYP=11)
      PARAMETER (NETYP=3)
      PARAMETER (NEbyV= NETYP*NVTYP)
      CHARACTER DICT*4,ROAD_TYPE*1
      CHARACTER NAME(NETYP,NVTYP)*4 !三處排放口、六種車種之PROFILE NUMBER
      INTEGER   UTM_E,UTM_N
      INTEGER   IT(NEbyV)      !每車種排放物種之個數
      INTEGER   SPE(NEbyV,116) !物種的CAS NO.(1-734)
      REAL      W(NEbyV,116)
      CHARACTER PRO_NO*4
      CHARACTER VT*4,VTYP(NVTYP)*4 ! 車種名稱
      CHARACTER ET*4,ETYP(NETYP)*2 ! 排放種類名稱

      PARAMETER (NREC=17248)       ! num. of REC.
      PARAMETER (NVEH=NVTYP)       ! NUM. OF VEHICLE TYPES
      PARAMETER (NPOL=10)          ! NUM. OF POLLUTANTS
      DIMENSION EM(NVEH,NPOL,NREC) ! THE RESULTANT MATRIX

      INTEGER X(NREC) ! the utm (km)
      INTEGER Y(NREC) ! the utm (km)
      INTEGER R(NREC) ! the road type (1-4)

C**   CBM RELATIVE DATA
      PARAMETER (L=171,NO=734,NC=9)
      INTEGER   SPE_NO
      REAL      BASE(NO,NC),BAS(L,NC),MW(NO),MW1(L)

C**   THE RESULTANT ARRAY
      PARAMETER (IX0= 81,IY0=2421)
      PARAMETER (LX =351,LY =2802)
      PARAMETER (DX=3,DY=3)
      PARAMETER (M=(LX-IX0)/DX)
      PARAMETER (N=(LY-IY0)/DY) ! THE X-Y MATRIX (FINE MESH)
      PARAMETER (LTYP=4)
      REAL      VOCB(M,N,LTYP,NC)

      real sumt(nvtyp)  !for check

C      DATA VTYP/ 'PLDG', 'BLDG', 'LDGT', 'LDDT', 'BUS ',
C     +           'HDGV', 'HDDT', 'MC2 ', 'MC4 '/
      DATA VTYP/'PLDG','BLDG',   !此順序須與rdline.f中的vtypout相同
     +  'LDGT','LDDT','HDGV','HDDT',
     +  'BUS ','MC2 ','MC4 ','PLDD','BLDL'/
      DATA ETYP/ 'EX','E ','R '/

      OPEN(1,FILE='vocparm/ASSIGN-L.TXT',status='old')
      OPEN(4,FILE='vocparm/PROFIL-L.TXT',status='old')

C ASSING PROFILE NO.
      DO I=1,NETYP
      DO J=1,NVTYP
      NAME(I,J)='    '
      ENDDO
      ENDDO
103   READ(1,101,END=102)VT,ET,PRO_NO
101   FORMAT(A4,1X,A2,2X,A4)
      KET=0
      DO IET=1,NETYP
        IF(ET.EQ.ETYP(IET))KET=IET
      ENDDO
      IF(KET.EQ.0) STOP' WRONG EMISSION SOURCE TYPE'
      KVT=0
      DO IVT=1, NVTYP
        IF(VT.EQ.VTYP(IVT))KVT=IVT
      ENDDO
      IF(KVT.EQ.0) STOP' WRONG VEHICLE TYPE'
      NAME(KET,KVT)=PRO_NO
      GOTO 103
102   DO I=1,NETYP
      DO J=1,NVTYP
      IF(NAME(I,J).EQ.'    ') THEN
        WRITE(*,*)' NOT ASSIGNED',I,J
        STOP
      ENDIF
      ENDDO
      ENDDO
      CLOSE(1)

C READ PROFIL.TXT (PRO-NO.& WT%)
C W(PRO_NO,SPE_NO)
C TEST
C
c1115 43135 31             0.16  F19890105
C


C**   計算各PROFILE的物種個數、存入陣列中
      DO I=1,NETYP
      DO J=1,NVTYP
        K=(I-1)*NVTYP+J
        IT(K)=0
      ENDDO
      ENDDO
77    READ(4,300,END=777)PRO_NO,SPE_NO,WT ! CAS NO. / WT IN %
      WRITE(*,301)PRO_NO,SPE_NO,WT
300   FORMAT(A4,7X,I3,10X,T25,F6.2)
301   FORMAT(1X,A4,7X,I3,10X,F6.2)
      DO I=1,NETYP
      DO J=1,NVTYP
        IF(NAME(I,J).EQ.PRO_NO) THEN
          K=(I-1)*NVTYP+J
          IT(K)=IT(K)+1
          SPE(K,IT(K))=SPE_NO
          W(K,IT(K))=WT
        ENDIF
      ENDDO
      ENDDO
      GOTO 77
777   CLOSE(4)
      DO I=1,NEbyV
        WRITE(*,*)'I1=',IT(I)
      ENDDO

C**   READ THE CARBON BOND MECHANISM MAPPING ARRAY
      OPEN(3,FILE='vocparm/CBM.DAT',status='old')

c-----物種:1.ole 2.par 3.tol 4.xyl 5.form 6. ald2 7.eth 8.isop 9.nr
      WRITE(*,*)'READ CBM(SPE_NO,I),I=1,9 ! '
      READ(3,*)
10    READ(3,200,END=997)SPE_NO,MW(SPE_NO),(BASE(SPE_NO,I),I=1,NC)
200   FORMAT(41X,I3,13X,F6.2,9F6.1)
      GOTO 10
997   CLOSE(3)


C START TO READ INPUT INVENTORY FILE
      OPEN(1,FILE='cl96.bin',FORM='UNFORMATTED',STATUS='OLD')
      READ (1)EM
998   CLOSE(1)
C** EM(NVEH,NPOL,NREC) ! THE RESULTANT MATRIX
CTSP    1-9      1       CEXHC   1-9      6
CPM     1-9      2       CEHC    1-3 8-9  7
CSOx    1-9      3       CRHC    1-3 8-9  8
CNOx    1-9      4       CMHC    1-9      9
CCO     1-9      5       CPb     1-9     10

C**   READ THE SORTED SEQUENCE OF DB'S
C**   THE FILE IS READ BY L2-LOC.F AND SORTED IN HE6
C**   FORMAT CNT//TWN//X/Y/IRD(2A2,A3,A4,A1)
      OPEN(1,FILE='LE-loc12.KIN',status='old')
      READ (1,*)NUMREC
      IF(NUMREC.GT.NREC)STOP'WRONG IN LE-loc12.KIND FILE'
      DO J=1,numrec
        READ (1,'(4X,I3,I4,I1)')X(J),Y(J),R(J)
      ENDDO
      CLOSE(1)

C CACULTE VOC=HC*WT

C EM=ton/yr
C CMB(g-mole/hr)=VOC*ratio/MW_VOC*1,000,000./8760.

      UNIT=1000.*1000./(24.*365.)
      sumT=0
      DO Loc=1,NREC
      IF(MOD(Loc,500).EQ.0) WRITE(*,*)Loc
      ix=(x(Loc)-IX0)/DX+1
      iy=(y(Loc)-IY0)/DY+1
      IR=R(Loc)

      DO I=1,NETYP
      DO J=1,NVTYP
        HC=EM(J,5+I,Loc)  ! 6→EX/ 7→E /8→R
        sumT(j)=sumT(j)+hc   !for check
        K=(I-1)*NVTYP+J
        DO II=1,IT(K) !某路段車種排放口之所有可能物種
          VOCwt=HC*W(K,II)/100.
          IS=SPE(K,II)
          VOCmole=VOCwt/MW(IS)*UNIT
          DO LS=1,NC ! CARBON BOND MECHANISM SPECISE
            VOCB(IX,IY,IR,LS)=VOCB(IX,IY,IR,LS)+VOCmole*BASE(IS,LS)
          ENDDO
        ENDDO ! END OF THIS車種&排放口
      ENDDO ! VEH-TYPE
      ENDDO ! SOURCE TYPE
      ENDDO ! NEXT RECORD

C**   WRITE THE RESULTANT MATRIX
      OPEN (1,FILE='96LVOC.BIN',FORM='UNFORMATTED',STATUS='UNKNOWN')
      WRITE(1)VOCB
      CLOSE(1)

C**   CHECKING THE RESULT
      OPEN(1,FILE='VOC-L1.dat')
      OPEN(2,FILE='VOC-L2.dat')
      OPEN(3,FILE='VOC-L3.dat')
      OPEN(4,FILE='VOC-L4.dat')
      do i=1,M
        do j=1,N
          IX=(i-1)*DX+IX0
          IY=(j-1)*DY+IY0
          do ir=1,4
            rSUM=0
            DO K=1,NC
              rSUM=rSUM+VOCB(I,J,IR,K)
            ENDDO
            if(rSUM.gt.0.)then
              WRITE(IR,50)IX*1000,IY*1000,(VOCB(I,J,IR,K),K=1,9)
            endif
50          FORMAT(2I10,9F10.4)
          enddo
        enddo
      enddo
      write(*,*)'Input voc emission (t/y):'
      do i=1,nvtyp
        write(*,*)vtyp(i),':',sumt(i)
      enddo
      write(*,*)'Total:',sum(sumt)
      END
