C*
C SOURCE=1 > PIPE-TAIL           CTSP    1-9      1
C SOURCE=2 > EVAPOR              CPM     1-9      2
C SOURCE=3 > R                   CSOx    1-9      3
C TYPE  =1 > PLDGV               CNOx    1-9      4
C TYPE  =2 > BLDGV               CCO     1-9      5
C TYPE  =3 > LDGT                CEXHC   1-9      6
C TYPE  =4 > LDDT                CEHC    1-3 8-9  7
C TYPE  =5 > BUS                 CRHC    1-3 8-9  8
C TYPE  =6 > HDGV                CMHC    1-9      9
C TYPE  =7   HDDT                CPb     1-9     10
C TYPE  =8   MC2
C TYPE  =9   MC4
C

      PARAMETER (NVTYP=11)
      CHARACTER DICT*4,ROAD_TYPE*1
c      CHARACTER VT*4,VTYP(NVTYP)*4 ! 車種名稱

      PARAMETER (NREC=17248)       ! num. of REC.
      PARAMETER (NVEH=NVTYP)       ! NUM. OF VEHICLE TYPES
      PARAMETER (NPOL=10)          ! NUM. OF POLLUTANTS
      DIMENSION EM(NVEH,NPOL,NREC) ! THE INPUT MATRIX

      INTEGER X(NREC) ! the utm (km)
      INTEGER Y(NREC) ! the utm (km)
      INTEGER R(NREC) ! the road type (1-4)

C**   THE RESULTANT ARRAY
      PARAMETER (IX0= 81,IY0=2421)
      PARAMETER (LX =351,LY =2802)
      PARAMETER (DX=3,DY=3)
      PARAMETER (M=(LX-IX0)/DX)
      PARAMETER (N=(LY-IY0)/DY) ! THE X-Y MATRIX (FINE MESH)
      PARAMETER (LTYP=4)
      PARAMETER (NC=4)! PM,SOX,NOX,CO
      REAL      POLB(M,N,LTYP,NC)

      INTEGER   MAP(NC)
      DATA MAP/2,3,4,5/

C START TO READ INPUT INVENTORY FILE
      OPEN(1,FILE='cl96.bin',ERR=998,FORM='UNFORMATTED'
     +  ,STATUS='OLD')
      READ (1,ERR=998)EM
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


C EM=ton/yr
C POL-L(Kg/day)=EM*1000./365.
      UNIT=1000./365.

      DO Loc=1,NREC
        IF(MOD(Loc,500).EQ.0) WRITE(*,*)Loc
        ix=(x(Loc)-IX0)/DX+1
        iy=(y(Loc)-IY0)/DY+1
        IR=R(Loc)
        DO IVTP=1,NVTYP
          DO LS=1,NC
            POLB(IX,IY,IR,LS)=
     +        POLB(IX,IY,IR,LS)+EM(IVTP,MAP(LS),Loc)*unit
          ENDDO ! END OF POLLUTANTS
        ENDDO ! VEH-TYPE
      ENDDO ! NEXT RECORD

C**   WRITE THE RESULTANT MATRIX
      OPEN (1,FILE='96LPOL.BIN',FORM='UNFORMATTED',STATUS='UNKNOWN')
      WRITE(1)POLB
c      write(*,*)'line sox:',sum(polb(:,:,:,2)),'kg/day'  !!
      CLOSE(1)

C**   CHECKING THE RESULT
      OPEN(1,FILE='POL-L1.dat')!國
      OPEN(2,FILE='POL-L2.dat')!省
      OPEN(3,FILE='POL-L3.dat')!縣
      OPEN(4,FILE='POL-L4.dat')!鄉
      do i=1,M
        do j=1,N
          IX=(i-1)*DX+IX0
          IY=(j-1)*DY+IY0
          do ir=1,4
            rSUM=0
            DO K=1,NC
              rSUM=rSUM+POLB(I,J,IR,K)
            ENDDO
            if(rSUM.gt.0)WRITE(IR,50)IX*1000,IY*1000,
     +        (POLB(I,J,IR,K),K=1,NC)  !PM,SOX,NOX,CO
50            FORMAT(2I10,9F10.4)
          enddo
        enddo
      enddo

      close(1)
      close(2)
      close(3)
      close(4)

c      OPEN(1,FILE='nox.dat')
c      do i=1,M
c        do j=1,N
c          IX=(i-1)*DX+IX0
c          IY=(j-1)*DY+IY0
c          do ir=1,4
c            rSUM=0
c            DO K=3,3 !PM,SOX,NOX,CO
c              rSUM=rSUM+POLB(I,J,IR,K)
c            ENDDO
c          enddo
c          K=3
c          if(rSUM.gt.0)WRITE(1,500)
c     +      IX*1000,IY*1000,(POLB(I,J,IR,K),ir=1,4),rsum
c500     FORMAT(2I10,9F10.4)
c        enddo
c      enddo
      END
