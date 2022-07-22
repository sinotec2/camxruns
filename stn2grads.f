c     this program use big_endian or little endian, must accord to epa.ctl
c-----------------------------------------------------------------------
      parameter(NPOL=11,LDATE=366)
      integer nst
      real,allocatable::x(:),y(:),C(:,:,:,:)
      real H(24)
      DIMENSION IDATE(LDATE)
      character date*4,YR*2,imo*2,namSt*2
      character*8 stn
      CHARACTER CH_nm*2,A2*2,A3*3
      dimension PM(4),ovm(2,0:11),novm(2,0:11)

      CHARACTER NAMPOL(0:NPOL)*3,myr*4
      character*4,allocatable::stnam(:)
      character*4 a4
      character cn
      character*500 aline
      logical lexist
      logical,allocatable::lsteff(:)  !▒▒O▒_▒b▒▒▒▒d▒▒

      DATA NAMPOL/'TMP','SO2','CMO','OZN','PMT',
     +  'NOX','P25','NO2','THC','NMH','WSP','WDR'/

      PARAMETER (NSP=6)  ! SP. NUM
      integer,allocatable::ISQ(:)     ! STATION COORDINATES (M)

      DATA NUAV/41/

      CHARACTER AD(12)*2,AM(12)*3,time*12
      DATA AD/'31','28','31','30','31','30','31','31'
     +  ,'30','31','30','31'/
      data AM/'jan','feb','mar','apr','may','jun',
     +        'jul','aug','sep','oct','nov','dec'/

      character*40 nam0(60),MRUNID

      open(111,file='abi_inp.txt',status='old')
      read(111,*)nam0(3)
      read(111,*)nam0(1),nam0(2)
      close(111)


      YR=NAM0(1)(1:2)
      READ(YR,'(I2)')iyr

      READ(NAM0(1)(1:2),*)IYB
      READ(NAM0(1)(3:4),*)IMB
      READ(NAM0(1)(5:6),*)IDB
      READ(NAM0(1)(7:8),*)IHB
      READ(NAM0(2)(1:2),*)IYE
      READ(NAM0(2)(3:4),*)IME
      READ(NAM0(2)(5:6),*)IDE
      READ(NAM0(2)(7:8),*)IHE
      IF(MOD(IYR,4).EQ.0)AD(2)='29'
      JdateB=IYB*100*100+IMB*100+IDB
      JdateE=IYE*100*100+IME*100+IDE
      IF(JdateB*100+IHB.GT.JdateE*100+IHE)then
        write(*,*)'Error, BEGIN TIME->END TIME',jdateB,jdateE
        stop
      ENDIF

C**   ▒M▒wnst
      open(21,file='cwb_epa.csv',
     +  status='old')
      read(21,*)
      nst=0
      do while(.true.)
        read(21,*,end=51)
        nst=nst+1
      enddo
100   format(i3,5x,a4,t16,i7,t25,i7,t33,i1)
51    continue
c-----Allocate
      allocate(x(nst))
      allocate(y(nst))
      allocate(C(nst,0:24,0:npol,ldate)) !0 for temp
      allocate(stnam(nst))
      allocate(isq(nst))
      allocate(lsteff(nst))

      rewind(21)
      read(21,*)
      do i=1,nst
        read(21,*,end=52)isq(i),x(i),y(i),StNam(i)
        if(isq(i).gt.460000)isq(i)=mod(isq(i)/10,1000)
      enddo
52    close(21)

      IID=0
      DO 98 JDATE=JdateB,JdateE
        ICAL=JDATE
        CALL juldate(ICAL)
        IF(ICAL.LE.0) GOTO 98
        IID=IID+1
        IF(IID.GT.LDATE) STOP 'IID.GT.LDATE'
        IDATE(IID)=JDATE
        print*,JDATE
98    CONTINUE
      MDATE=IID
        print*,MDATE

c-----Ū▒J▒ʴ▒▒▒
      DO jst=1,nst
      DO IMm=IYB*100+IMB,IYE*100+IME
        IM=mod(IMm,100)
        if(IM.le.0.or.IM.gt.12)cycle
        write(IMO,'(I2.2)')IM
          write(A3,'(i3.3)')isq(jst)
          write(myr,'(I2,I2.2)')20,IMm/100
          write(yr,'(I2.2)')IMm/100
          WRITE(*,*)'o:/users/4139/epa/'
     +      //myr//'/hs'//YR//IMO//AD(IM)//'.'//A3

          Open(1,file=
     +      '/st1/data/epa/'//myr//'/HS'//YR//IMO//AD(IM)//'.'
     +      //A3,STATUS='OLD',iostat=iio)
          if(iio.ne.0)then
            Open(1,file=
     +        '/st1/data/epa/'//myr//'/hs'//YR//IMO//AD(IM)//'.'
     +        //A3,STATUS='OLD',err=99)
           endif
 11       read(1,5,end=99)ist, ipollu,idat2,(H(I),i=1,24)
  5       format(I3,9x,i2,4x,i6,24(1x,f6.0))
          IF(IST.NE.isq(JST)) STOP 'Error01'
          IGG=1
          DO K=0,NPOL
            KK=K
            IF(K.eq.0)KK=14 ! temp_c
            IF(K.eq.6)KK=33 ! NO ->PM2.5
            IF(IPOLLU.EQ.KK) IGG=0
          ENDDO
c**     write(*,5)ist, ipollu,idat2,(H(I),i=1,24)

          IF(IGG.EQ.1) GOTO 11
          IF(IDAT2.LT.IDATE(1)    ) GOTO 11
          IF(IDAT2.GT.IDATE(MDATE)) GOTO 11
          DO I=1,MDATE
            IF(IDAT2.EQ.IDATE(I)) M=I
          ENDDO
          K=IPOLLU
          IF(K.eq.14)K=0 ! temp_c
          IF(K.eq.33)K=6 ! NO ->PM2.5
          DO I=1,24
            ih=i-1
            if(h(i).lt.0..and.K.ne.0)then
              C(JST,ih,K,M)=-99.
            else
              C(JST,ih,K,M)=H(I)
            endif
          ENDDO !hr
          GOTO 11
99        CLOSE(1)

        ENDDO !NEXT STATION
      ENDDO !NEXT MONTH

C**
C**   WRITE THE TIME SERIES FILE
C**
      open(21,file='stn2grads.dat',form='unformatted',
     &status='unknown' ,access='stream')
!    &status='unknown',access='SEQUENTIAL')
!,recordtype='STREAM',
!     &convert='big_endian' ,status='unknown')
!,access='stream')
!     &status='unknown' ,access='stream')
      nlev=1
      nflag=1
      itime=0
      Atime=0
      DO M=1,MDATE
        ICAL = IDATE(M)
        CALL juldate(ICAL)
        IY=ICAL/1000.
        ICAL=ICAL-IY*1000.
        DO I=0,23
          IF(M.EQ.1.AND.I.LT.IHB) GOTO 202        !NOT YET
          IF(M.EQ.MDATE.AND.I.GT.IHE) GOTO 202    !!END OF REQUIREMENT
!         Atime=Atime+0.001
          itime=itime+1
!       Atime=itime
          DO J=1,NST
            WS=C(J,I,10,M)
            WD=C(J,I,11,M)
            !print*,IDATE(M),isq(J),WD,WS
            IF(WD.LT.0)WD=WD+360
            if(wd.gt.360.)then
              wd=-99.
            endif
          write(stn,'(I3.3)')isq(J)
cthe time label is not used, all def by ctl
c         write(time,'(I2.2,A1,I2,A3,I4)')I,'Z',mod(IDATE(M),100),
c    +    AM(mod(IDATE(M)/100,100)),2000+IDATE(M)/10000
c 00z31dec2009
         u=-ws*sin(wd*3.14/180)
         v=-ws*cos(wd*3.14/180)
        if(sqrt(u*u+v*v).lt.0.1) then
        u=-99.99
        v=-99.99
        endif
          write(21)stn,y(J),x(J),Atime,nlev,nflag,u,v
!          write(21)ws
          ENDDO ! LOOP FOR STATION
          write(21)stn,0.0,0.0,0.0,0,0
202      CONTINUE
        ENDDO !NEXT HR
      ENDDO !LOOP FOR DAY
      print*,itime

      END
c-------------------------------------------------------------
      subroutine juldate(idate)
      dimension nday(12)
      data nday/31,28,31,30,31,30,31,31,30,31,30,31/
c
c-----Entry point
c
      iyear = idate/10000
      imonth = (idate - iyear*10000)/100
      IF(imonth.GT.12.OR.imonth.LE.0) THEN
        IDATE=-1
        RETURN
      ENDIF
      iday = idate - iyear*10000 - imonth*100
      IF(IDAY.GT.NDAY(IMONTH).OR.IDAY.LE.0) THEN
        IDATE=-1
        RETURN
      ENDIF
c
      nday(2) = 28
      if (mod(iyear,4).eq.0) nday(2) = 29

      mday = 0
      do 10 n = 1,imonth-1
        mday = mday + nday(n)
 10   continue
      jday = mday + iday
      idate = iyear*1000 + jday
c
      return
      end
c----------------------------------------------------------------
      function StVal(IX,IY,L,XRf,YRf)
      PARAMETER (NI=200,NJ=200,MXSP=90)
      COMMON /MTRX/A1(NI,NJ,MXSP)
      real MTX(0:11)
      DIMENSION A(4)
!            1   2    3    4   5    6     7   8
!           O3  NO2  SO2  VOC PM25 PM10   NO PM2ND (L=8,after shrink
!     DATA NAMPOL/'SO2','CMO','OZN','PMT','NOX',
!    +'P25','NO2','THC','NMH','WSP','WDR'/
      data MTX /   3,3 ,    8 ,   1 ,   6 ,   7 ,
     +  5 ,  2 ,   4,     4,    5,    6/ !(u &v)
      K=1
      ISP=MTX(L) !L=0, for temp_K in 3dMET file(location 3)
      if(ISP.eq.0) then
        StVal=0.
        return
      endif
      A=0
      DO I=IX,IX+1
        DO J=IY,IY+1
            A(K)=A1(I,J,ISP)
          K=K+1
        ENDDO
      ENDDO
      u11=a(1)*(1.-YRf)+a(2)*YRf
      u12=a(3)*(1.-YRf)+a(4)*YRf
      StVal=U11*(1.-XRf)+U12*XRf
      return
      end
c---------------------------------------------------------------------
      subroutine caldate(idate)
c
c-----CAMx v3.10 020410
c
c     CALDATE converts date from Julian (YYJJJ) format to calender
c     (YYMMDD) format
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        idate               julian date (YYJJJ)
c
c     Output arguments:
c        idate               calender date (YYMMDD)
c
c     Routines Called:
c        none
c
c     Called by:
c        AREAPREP
c        BNDPREP
c        CNCPREP
c        DRYDEP
c        PIGPREP
c        RDPTHDR
c        READZP
c        STARTUP
c
      dimension nday(12)
      data nday/31,28,31,30,31,30,31,31,30,31,30,31/
c
c-----Entry point
c
c-----If it is already in calender date, return
c
      if (idate.gt.100000) goto 9999
      iyear = idate/1000
      jday = idate - iyear*1000
c
      nday(2) = 28
      if (mod(iyear,4).eq.0) nday(2) = 29
      mday = 0
      do 10 imonth = 1,12
        mday = mday + nday(imonth)
        if (mday.ge.jday) go to 20
 10   continue
 20   iday = jday - (mday - nday(imonth))
      idate = iyear*10000 + imonth*100 + iday
c
 9999 return
      end
        subroutine UV2WsWd(wu,wv,ws,wd)
      pi=acos(-1.)
      ws=sqrt(wu**2+wv**2)
      if (ws.eq.0)then
          wd=0.
          return
      endif
      ang=wv/ws
      if (abs(ang).le.10e-5) then
         ang=0
      endif
      if (wu.ge.0 .and. wv.ge.0) then
         wd=asin(ang)*180./pi-90.
      endif
      if (wu.lt.0 .and. wv.ge.0) then
         wd=-asin(ang)*180./pi+90.
      endif
      if (wu.lt.0 .and. wv.lt.0) then
         wd=asin(-ang)*180./pi+90.
      endif
      if (wu.ge.0 .and. wv.lt.0) then
         wd=-asin(-ang)*180./pi-90.
      endif
          WD=-WD-180.
          if(WD.lt.0) WD=WD+360.
        return
        end
