      integer,allocatable::ISQ(:)     ! STATION COORDINATES (M)
      real,allocatable::x(:),y(:),var(:,:,:,:),varout(:,:)
      character, allocatable:: A4(:)*12
      integer, allocatable::lstn(:)
      character stn*8, fmt*100
      parameter(NPOL=11,LDATE=366)
      CHARACTER NAMPOL(0:NPOL)*3,myr*4
      DATA NAMPOL/'TMP','SO2','CMO','OZN','PMT',
     +  'NOX','P25','NO2','THC','NMH','WSP','WDR'/
C** 
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
      allocate(isq(nst))
      allocate(A4(nst))
      allocate(lstn(nst))

      rewind(21)
      read(21,*)
      do i=1,nst
        read(21,*,end=52)isq(i),x(i),y(i)
        if(isq(i).gt.460000)isq(i)=mod(isq(i)/10,1000)
      enddo
52    close(21)
      open(1,file='ovm.dat',status='unknown')
      itime=0
      lstn=0
      read(1,300,end=98)
      DO 
        read(1,300,end=98)is
        js=0
        do i=1,nst
          if(is.eq.isq(i))then
            js=i
            exit
          endif
        enddo
        if(js.eq.0) then
        print*,is
        stop 'stn not found'
        endif
        if(lstn(js).eq.0)lstn(js)=1
       itime=itime+1
      enddo
98    rewind(1)
      
      ntime=(itime)/sum(lstn) !first line
      if(ntime*sum(lstn).ne.itime)stop 'line missing'
      allocate(var(2,9,nst,ntime)) 
      read(1,*)
      iss=0
      DO js=1,nst
        if(lstn(js).eq.0)cycle
      DO it=1,ntime 
        read(1,300,end=99)is,A4(js),Jul,ical,
     & (var(1,j,js,it),j=1,9),ws,wd,(var(2,j,js,it),j=1,8)
300     FORMAT(I3,A12,I8,1x,I6,7F7.1,2f8.0,9f7.1,f8.0)
          var(1,8,js,it)=amax1(0.,var(1,9,js,it))
        do j=1,8
          var(1,j,js,it)=amax1(0.,var(1,j,js,it))
        enddo
      ENDDO !LOOP FOR DAY
        iss=iss+1
        k=int(alog10(real(iss)))+1
        write(fmt,'(A,I1,A)')'(A4,I',k,',A1,I3.3,A4)'
!       write(*,trim(fmt))'nst.',iss,'=',is,A4(js)
      ENDDO ! LOOP FOR STATION
      do L=1,8
        k=0
      DO it=1,ntime
      DO JS=1,NST
      if(lstn(js).eq.0)cycle
        if(var(1,L,js,it).le.0.or.var(2,L,js,it).le.0) cycle
        if(isnan(var(1,L,js,it)).or.isnan(var(2,L,js,it))) cycle
        k=k+1
      ENDDO ! LOOP FOR STATION
      ENDDO !LOOP FOR DAY
      LL=L
      if(L.eq.8)LL=9
      print*,'total num=',NAMPOL(LL),k
      allocate (varout(2,k))
      open(21,file='ovm2grads_'//NAMPOL(LL)//'.dat',form='unformatted',
     &status='unknown' ,access='direct',recl=k)
          k=0
      DO it=1,ntime
      DO JS=1,NST
      if(lstn(js).eq.0)cycle
        if(var(1,L,js,it).le.0.or.var(2,L,js,it).le.0) cycle
        if(isnan(var(1,L,js,it)).or.isnan(var(2,L,js,it))) cycle
        k=k+1
             do i=1,2
               varout(i,k)=var(i,L,js,it)
             enddo
      ENDDO ! LOOP FOR STATION
      ENDDO !LOOP FOR DAY
        do i=1,2
          write(21,rec=i)varout(i,:)
        enddo
        close(21)
        deallocate(varout)
      open(21,file='scat_'//NAMPOL(LL)//'.ctl',status='unknown')
      write(21,*)'dset ^ovm2grads_'//NAMPOL(LL)//'.dat'
      write(21,*)'undef 99.99'
      write(21,*)'xdef ',k,' linear 1 1'
      write(21,*)'ydef 1   linear 1 1'
      write(21,*)'zdef 1   levels 1000'
      write(21,*)'tdef 1   linear 12Z31dec2009 1hr'
      write(21,*)'vars 2'
      write(21,*)'vo 0 99 SPEC (UNIT) obs'
      write(21,*)'vm 0 99 SPEC (UNIT) mdl'
      write(21,*)'endvars'
      close(21)
      enddo !L
      stop
99    print*,'ntime wrong'
      stop
      end


