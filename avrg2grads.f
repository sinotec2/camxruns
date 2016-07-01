c this routine calculate all time average of AVRG fields, and OP AS A ARG FILE
       INCLUDE  'PARAMS.CMD'
       INCLUDE  'CHPARM.CMD'
       INCLUDE  'CNTROL.CMD'
       INCLUDE  'FILCON.CMD'
       INCLUDE  'SEGTAB.CMD'
       INCLUDE  'NETDEP.CMD'
       INCLUDE  'BALANC.CMD'
       INCLUDE  'LOCPTR.CMD'
       INCLUDE  'MSCAL.CMD'
      integer,parameter::fmax=300
      integer,parameter::MXSPEC=40
      integer itmp(4),nti(2)
      character*4,allocatable:: SPNAME(:,:)
      CHARACTER*60 nam0(fmax),unit ! input/output file names
      character*4 fname(10)
      character note(60)*4,names*10,amonth(12)*3
      logical O3,GRD 
      integer,allocatable::ndate2(:),ndlast2(:)
      real,allocatable::ttime2(:),ttlast2(:)
      real,allocatable:: A1(:,:,:,:,:),tim(:,:)
      real,allocatable:: tm(:,:,:)
      integer,allocatable:: yyjul(:,:)
      real SCR(40000),r(4)
      data amonth/'jan' ,'feb' ,'mar' ,'apr' ,'may' ,'jun',
     &            'jul' ,'aug' ,'sep' ,'oct' ,'nov' ,'dec'/
      NUAV=41
      narg=iARGc ()
      if(narg.ne.3)stop  'avrg2grads_all File ispec, ihr'

      do i=1,narg
        call getarg(i,nam0(i))
      enddo
      read(nam0(2),*)ispec
      read(nam0(3),*)ihr
      narg=2
      i=1
          open(i+10,file=trim(nam0(i)),
     +      form='unformatted',convert='BIG_ENDIAN',STATUS='unknown')
        print*,trim(nam0(i))
      nfile=narg-1
      iout=narg+10

      allocate(ndate2(nfile))
      allocate(ndlast2(nfile))
      allocate(ttime2(nfile))
      allocate(ttlast2(nfile))
      DO IRD=1,nfile
        READ (IRD+10) fname, note, NOSEG, NOSPEC,
     +    NDATE2(ird), TTIME2(ird),
     $    NDLAST2(ird), TTLAST2(ird)
       enddo
      allocate(SPNAME(10,NOSPEC))
       ndate=minval(ndate2)
       ndlast=maxval(ndlast2)
       ttime=minval(ttime2)
       ttlast=maxval(ttlast2)
       numHr=(ndlast-ndate)*24+(ttlast-ttime)

       do ird=1,nfile
C
C--REGION DESCRIPTION HEADER
        READ (IRD+10) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $    NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
!       print*, XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY

!       if(ird.eq.1)then
!         write(*,*) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
!    $      NOXG, NOYG, NOZ, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
!       endif
C
C--SEGMENT DESCRIPTION HEADER
        READ (IRD+10) (Itmp(j), J=1,4)
!       print*, (Itmp(j), J=1,4)
C
C--SPECIES DESCRIPTION HEADER
        nti(IRD)=0
        READ (IRD+10) ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
!       print*, ((SPNAME(I,J), I=1,10), J=1,NOSPEC)
        do 
          READ (ird+10,END=30,err=30)i1,t1,i2,t2
!       print*,i1,t1,i2,t2
          DO  L=1,NOSPEC
            DO  K=1,NOZ
              READ (ird+10,END=30,err=30)
            enddo !k
          enddo !l
        nti(IRD)=nti(IRD)+1
        enddo !it
C
C         FIRST, WRITE TIME INTERVAL
C
30    rewind(IRD+10)
      do i=1,4!skip the header
              READ(10+ird)
      enddo
      ENDDO ! next IRD input file
      NXY=NOXG*NOYG
      NT=nti(1)!(NDLAST-NDATE)*24+(TTLAST-TTIME)
      allocate(A1(NOXG,NOYG,NOZ,NOSPEC,NT))
      allocate(tm(NXY,NOZ,NOSPEC))
      allocate(yyjul(2,NT),tim(2,NT))

      do  ird=1,nfile
        icount=0
        do it=1,NT
          READ (ird+10,END=33,err=33) yyjul(1,it),tim(1,it), yyjul(2,it),tim(2,it)
          DO  L=1,NOSPEC
            DO  K=1,NOZ
              READ(10+ird,err=33)ISEG,(SPNAME(I,L),I=1,10),((A1(i,j,k,l,it),i=1,NOXG),j=1,NOYG)
            enddo !k
          enddo !l
          icount=icount+1
        enddo !it
      enddo !ird
      goto 34
33    print*,'end of file'
34    if(NT.ne.icount) stop 'NT not right'
      L=ispec
        do i=1,10
        names(i:i)=SPNAME(I,L)(1:1)
        enddo
      nam0(2)=trim(nam0(1))//trim(names)//'.ga'
      OPEN(12,file=trim(nam0(2)),action='write',form='unformatted',access='direct',
     & recl=NOXG*NOYG,status='unknown')
!loop over timesteps
      irec =0 
      do it  = 1, nt
      do l  = 1, NOSPEC
      do k=1,NOZ
        irec = irec + 1
        WRITE(12,rec=irec) ((A1(i,j,k,l,it),i=1,NOXG),j=1,NOYG)
      enddo
       write(*,*)maxval(A1(:,:,:,L,it))
      enddo
      enddo
      close(12)
      nam0(3)=trim(nam0(1))//trim(names)//'.ctl'
      OPEN(13,file=trim(nam0(3)),status='unknown')
      data XLONC/120.99/PHIC/23.61/
      iutmzon=(180+XLONC)/6+1
      call utmgeo(0,iutmzon,rx4,ry4,XLONC,PHIC)
      Xcent=rx4*1000.           !center point UTM-m
      Ycent=ry4*1000.           !center point UTM-m
      r(1)=xorg+Xcent
      r(2)=xorg+real(noxg-1)*deltax+Xcent
      r(3)=yorg+Ycent
      r(4)=yorg+real(noyg-1)*deltay+Ycent
      do i=1,4
        r(i)=r(i)/1000.
      enddo
      call utmgeo(1,iutmzon,r(1),r(3),r1,r3)
      call utmgeo(1,iutmzon,r(2),r(4),r2,r4)
      dx=(r2-r1)/(noxg-1)
      dy=(r4-r3)/(noyg-1)
      ical=yyjul(1,1)
      call caldate(ical)
      IY=ical/100/100
      IM=mod(ical/100,100)
      ID=mod(ical,100)
      print*, int(tim(1,1))
      write(nam0(4),'(I2.2,A1,I2.2,A3,I4,I4,A2)')
     & int(tim(1,1)),'Z',ID,amonth(IM),2000+IY,ihr,'hr'
      write(13,*)'dset ',trim(nam0(2))
      write(13,*)'title ',(trim(fname(i)),i=1,10),' ',(trim(note(i)),i=1,60)
      !write(13,*)'options little_endian'
      write(13,*)'byteswapped'
      write(13,*)'undef -9999.'
      write(13,*)'xdef ',NOXG,' linear ',r1,dx
      write(13,*)'ydef ',NOYG,' linear ',r3,dy
      write(13,*)'zdef ',NOZ,'linear 1.0 1.0'
      write(13,*)'tdef ',nt,' linear ',trim(nam0(4))
      write(13,*)'vars ',ispec
      DO  L=ispec,ispec !1,NOSPEC
        unit=trim(names)//'  (ppb)'
        if(names(1:1).eq.'P')unit=trim(names)//'  (ug/M3)'
        if(names(1:4).eq.'PM10')unit='PM10  (ug/M3)'
        if(names(1:4).eq.'PM25')unit='PM2.5  (ug/M3)'
        do i=1,10
          if(names(i:i).eq.'_') then
            unit=names(i+1:10)
            do j=i+1,10
              names(j:j)=' '
            enddo
            exit
          endif
        write(13,*) trim(names),NOZ,' 99 ', trim(unit)
      enddo
      write(13,*)'endvars'
      end

