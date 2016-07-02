C PROGRAM VOC-A
C*INPUT INPUT
C BASE-DATA = (ASC-80.TXT)
C NSC          Length = 6    Type = C
C NSC_SUB      Length = 1    Type = C
C DICT         Length = 4    Type = C
C SIC          Length = 4    Type = C
C BASICDATA    Length = 15   Type = N  Dec = 3
C UNIT         Length = 15   Type = C
C CITY_DIS     Length = 6    Type = C
C HSIANG_DIS   Length = 6    Type = C
C EEM_THC      Length = 15   Type = N  Dec = 3    > 0 ?
C EEM_MHC      Length = 15   Type = N  Dec = 3    > 0 ?
C*INPUT V_ASSIGN.TXT
C 10100101 1178 FE
C -------- ----
C (SCC,PRO_NO)
C FORMAT(A8,1X,A4)
C*INPUT V_PROFIL.TXT
C 0068 43205 57  115071   100.00  F19890105
C ----       ---          ------
C (PRO_NO.,SPE_NO.,WT%)
C FORMAT(A4,7X,I3,10X,F6.2)
      PARAMETER (NC=9)  !VOC 9物種
      REAL      VOCB(NC)
c      REAL      THC,VOC,WT,BASIC
      REAL      VOC,BASIC
      CHARACTER SCC*7
      CHARACTER NSC*6,NSC_S*1,SIC*4,CITY*6,HSIANG*6

      parameter (nrec=10000)
      character TYPE(nrec)*7
      character UTME(nrec)*6
      character UTMN(nrec)*7
      character DICT(nrec)*4
      real em(7,nrec) !SOX/TSP/PM/CO/NOX/THC/NMHC

C**   CBM RELATIVE PARAMETER
c      REAL      THC,MHC,HR
      INCLUDE 'SPECDATA.INC'

      PARAMETER (IX0=81,IY0=2421)
      PARAMETER (LX=351,LY=2802)
      PARAMETER (DX=3,DY=3)
      PARAMETER (M=(LX-IX0)/DX)
      PARAMETER (N=(LY-IY0)/DY) ! THE X-Y MATRIX (FINE MESH)
      REAL      QG(M,N,13)!in g-mole/HR PNSC
      DIMENSION ES(4) ! in g-mole/HR (PM g/hr)
      DIMENSION WM(4),MAP(4) !PNSC
      DATA WM/ 1.0, 46., 64.,28./
      DATA MAP/3,5,1,4/

      UNIT=1000*1000./365./24.        !T/Y->g/hr

      CALL SPECIATE(0,THC,SCC,VOCB,UNIT)  !0初始化，讀入speciate相關table

      open(2,file='96ag-x.bin.f95',status='OLD',form='unformatted')
      open(1,file='96AREA.f95',status='UNKNOWN',form='unformatted')

      open(22,file='warn.txt')

C *
C START TO READ SCC
C *
c10    READ(2,100,END=999) NSC,NSC_S,DICT,SIC,BASIC,CITY,HSIANG,THC
      DO 11 I=1,M
      DO 11 J=1,N
      DO 11 K=1,13
        QG(I,J,K)=0.
11    CONTINUE
      KREC=0
1     read(2,end=999)TYPE,UTME,UTMN,DICT,EM   !SOX/TSP/PM/CO/NOX/THC/NMHC
      print*,'pass1',KREC
      KREC=KREC+NREC
      DO 10 irec=1,NREC
        SCC=TYPE(irec)  !NSC//NSC_S
        IF(SCC.EQ.''.or.SCC(1:2).eq.'95')GOTO 10
        READ(UTME(irec),*)IUTNE
        READ(UTMN(irec),*)IUTNN
        IX=(IUTNE-IX0)/DX+1
        IY=(IUTNN-IY0)/DY+1
        IF(IX.LE.0.OR.IY.LE.0)then
c          PRINT*,IX,IY,UTME(irec),UTMN(irec)
          cycle
        endif
        IF(IX.GT.M.OR.IY.GT.N)then
c          PRINT*,IX,IY,UTME(irec),UTMN(irec)
          cycle
        endif
        DO IE=1,4       !PNSC
          QG(IX,IY,IE)=QG(IX,IY,IE)+EM(MAP(IE),Irec)*UNIT  !N/S/C在後續emission程序會除以分子量，故此處單位處理為g/hr
c          QG(IX,IY,IE)=QG(IX,IY,IE)+EM(MAP(IE),Irec)*UNIT/WM(IE)
        ENDDO
        THC=EM(6,Irec)
        IF(THC.GE.0.) THEN                    !有thc排放量的才做speciate
          CALL SPECIATE(1,THC,SCC,VOCB,UNIT)  !1表進行speciate作業,thc:排放量,scc:面源代碼,vocb:9種voc排放量(結果),unit:T/y->g/hr(換算為mole數在subroutine中處理)
          DO IE=1,NC  !VOC CB的9個物種
            QG(IX,IY,IE+4)=QG(IX,IY,IE+4)+VOCB(IE)
          ENDDO
        ENDIF
10    CONTINUE         !NEXT RECORD
      GOTO 1 !NEXT 10000 REC
999   CLOSE(2)
      WRITE(1)QG
      CLOSE(1)

      OPEN(1,FILE='areaCHK.DAT')

C**   THE UNIT IS g-mole/HR
      DO IX=1,M
        DO IY=1,N

C**     CHECKING FOR NOT ZERO EMISSION
          rSUM=0
          DO IE=1,13
            rSUM=rSUM+QG(IX,IY,IE)
          ENDDO
          IF(rSUM.GT.0) THEN
            XJX=(iX-.5)*DX+IX0
            YJY=(IY-.5)*DY+IY0
C**       THE UNIT IS Kg-mole/HR
            WRITE(1,50)XJX*1000,YJY*1000,(QG(IX,IY,K)/1000.,K=1,13)
50          FORMAT(2F10.0,13E15.4)
          ENDIF
        ENDDO
      ENDDO
      CLOSE(1)
      STOP
      END

      SUBROUTINE SPECIATE(IFLAG,THC,SCCin,VOCB,UNIT)
C**   CBM RELATIVE PARAMETER
      CHARACTER SCCin*7,PRO_NO*4,PRO_NO1*4
      PARAMETER (NC=9)  !9種CB物種
      DIMENSION VOCB(NC)
      INCLUDE 'SPECDATA.INC'
      COMMON /SPECDATA/BASE,MW,SCCsgn,PRO_NOsgn,PFtab,IBegProf,SPE_NO,WT

C-----此處使用之CB(Carbon Bond)物種與順序(ref. CBM.DAT)--------
c     1.OLE        6.ALD2
C     2.PAR        7.ETH
C     3.TOL        8.ISOP
C     4.XYL        9.NR
C     5.FORM
C
C *
C READ THE CBM.DAT
C *
      IF(IFLAG.EQ.0) THEN

      OPEN(5,FILE='CBM.DAT',status='old')
      read(5,*)
      WRITE(*,*)'READ CBM(SPE_NO,I),I=1,9 ! '
500   FORMAT(41X,I3,13X,F6.2,9F6.1)
11    READ(5,500,END=12)SPE,MW(SPE),(BASE(SPE,I),I=1,NC)
      GOTO 11
12    CLOSE(5)

C**   讀入SCC與PROFILE NO之對照表
      OPEN(3,FILE='ASSIGN-A.TXT',status='old')
      DO 20 I=1,Nsgn
        READ(3,'(A7,1X,A4)')SCCsgn(I),PRO_NOsgn(I)
20    continue
      CLOSE(3)

C CHECK PROFILE NO.
      OPEN(4,FILE='V_PROFIL.TAB',status='old')
      DO I=1,Kprof
        READ(4,'(A4,I10)')PFtab(I),IBegProf(I)  !起始的序號及對應的PROFILE NO
      ENDDO
      CLOSE(4)

      IBegProf(Kprof+1)=Nprof

C**
      OPEN(4,FILE='V_PROFIL.TXT',status='old')
      DO I=1,Kprof
        NumOfSp=IBegProf(I+1)-IBegProf(I)
        DO J=1,NumOfSp
          READ(4,300)PRO_NO1,SPE_NO(I,J),WT(I,J)
300       FORMAT(A4,7X,I3,10X,F6.2)
        ENDDO
      ENDDO
      CLOSE(4)

      ELSE ! doing the speciation

C**       MATCHING THE SCC's TO FIND PROFILE NUMBER
          PRO_NO1='    '
          DO ISC=1, Nsgn
           IF(SCCsgn(ISC).EQ.SCCin) THEN
            PRO_NO1 = PRO_NOsgn(ISC)
            GOTO 31
           ENDIF
          ENDDO
          IF(PRO_NO1.EQ.'    ')then
            write(22,*)' SCC NOT FOUND ',sccin
c            stop
            return  !找不到對應的會跳出，不處理
          endif

C**       LOOK-UP THE PROFILE NUMBER
31        JBegProf=0
          DO IPF=1,Kprof
           IF(PFtab(IPF).EQ.PRO_NO1) THEN
            JBegProf=IBegProf(IPF)
            GOTO 32
           ENDIF
          ENDDO
          IF(JBegProf.EQ.0) then
            print*,SCCin(1:7),'x',PRO_NO1
            STOP' PROF. NO. NOT FOUND'
          endif

C**       SPLIT THE SPECIATION
32        NumOfSp=IBegProf(IPF+1)-JBegProf
          DO ICBM=1,NC
            VOCB(ICBM)=0.
          ENDDO
          DO JS=1,NumOfSp
            IF(WT(IPF,JS).GT.0.) THEN !約有20%機會可能等於零，白作工
              VOCwt=THC*WT(IPF,JS)/100.
              IS=SPE_NO(IPF,JS)
              if(MW(IS).GT.0)VOCmole=VOCwt/MW(IS)*UNIT  !=>mole/hr，VOCmole與VOCwt為各單項voc物種
              DO ICBM=1,NC
                VOCB(ICBM)=VOCB(ICBM)+VOCmole*BASE(IS,ICBM)  !使用CBM.DAT資料lump VOC為9類物種
              ENDDO
            ENDIF
          ENDDO ! SUM ALL SPECIES IN THIS PROFILE NUMBER
       ENDIF   ! FLAG ALL
       RETURN
       END
