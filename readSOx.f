	parameter(Xmax=149.5,Xmin=60.,Ymax=59.5,Ymin=-11.5)
        parameter(PHIC=23.60785)!   ; CENTRAL LATITUDE (minus for southern hemesphere)
        parameter(XLONC=120.9908)
	parameter(MX=200,MY=200)
	dimension TPY(MX,MY)
       INCLUDE  'CNTROL.CMD'
      PARAMETER (IX0=45,IY0=2421)
      PARAMETER (LX=351,LY=2802)
      PARAMETER (lengX=LX-IX0)
      PARAMETER (lengY=LY-IY0)
      dimension idct(IX0:LX,IY0:LY)
      common /dictxy/idct
	open(1,file='SO2_EM_TOTAL_2010_REF.0.5x0.5',status='unknown')
	xutm=-999
	XORG=+999
	yutm=-999
	YORG=+999
	do i=1,10
	read(1,*)
	enddo
	do
	read(1,*,end=1)alon,alat
	xutm=amax1(xutm,alon)
	XORG=amin1(XORG,alon)
	yutm=amax1(yutm,alat)
	YORG=amin1(YORG,alat)
	enddo
1	rewind(1)
	print*,xutm,XORG,yutm,YORG
	NOXG=(xutm-XORG)/0.5
	NOYG=(yutm-YORG)/0.5
	if(NOXG.gt.MX)stop 'MX not enough'
	if(NOYG.gt.MY)stop 'MY not enough'
	do i=1,10
	read(1,*)
	enddo
	do
	read(1,*,end=2)alon,alat,aa
	i=(alon-XORG)/0.5+1
	j=(alat-YORG)/0.5+1
	if(aa.eq.0) then
	  aa=-999
	else
	  aa=alog10(aa)
	endif
	TPY(i,j)=aa
	enddo
2	close(1)
!	goto 10
	iutmzon=(180+XLONC)/6+1
        call readDict
	i1=(118.-XORG)/0.5+1
	j1=(21.-YORG)/0.5+1
	i2=(123.-XORG)/0.5+1
	j2=(26.-YORG)/0.5+1
	tmin=9999999
	do i=i1,i2
	do j=j1,j2
	IsLand0=0
	do ii=i-1,i+1
	do jj=j-1,j+1
	  rlon4=(ii-1)*0.5+XORG
	  rlat4=(jj-1)*0.5+YORG
          call utmgeo(0,iutmzon,rx4,ry4,rlon4,rlat4)
          KmX=rx4
          KmY=ry4
	print*,rlon4,rlat4,KmX,KmY
          if(IsLand(KmX,KmY).ne.0)IsLand0=1
	enddo
	enddo
	if(IsLand0.eq.1) then
	  TPY(i,j)=-999
	else
	  tmin=amin1(tmin,TPY(i,j))
	endif
	enddo
	enddo
	open(1,file='taiwanBlank.dat',status='unknown')
	ib=1
	do i=i1,i2
	do j=j1,j2
	 if(TPY(i,j).eq.-999) then
	   TPY(i,j)=tmin
	   rlon4=(i-1)*0.5+XORG
	   rlat4=(j-1)*0.5+YORG
	   write(1,*)rlon4,rlat4,ib
	  ib=ib+1
	  endif
	enddo
	enddo
	close(1)
10	continue 
	open(1,file='sox-lt.grd',status='unknown')
	call surfer(1,TPY)
	close(1)
	stop
	end
      subroutine surfer(n,a)
       INCLUDE  'CNTROL.CMD'
      dimension a(200,200)
      DIMENSION SCR(40000)
	ni=NOXG
	nj=NOYG
        CMIN=99999
        CMAX=-99999
        J1=1
        J2=ni*nj
        JJ=1
        DO 221 J=1,nj
        DO 221 I=1,ni
          SCR(JJ)=a(i,j)
          CMIN=AMIN1(CMIN,SCR(JJ))
          CMAX=AMAX1(CMAX,SCR(JJ))
          JJ=JJ+1
221     CONTINUE
        WRITE(n,'(a4)')'DSAA'
        WRITE(n,*)NOXG,NOYG
        WRITE(n,*)XORG,XUTM !dot in grid
        WRITE(n,*)YORG,YUTM
        WRITE(n,*)CMIN,CMAX
        I1=NI/10+1
        I2=MOD(NI,10)
        IF(I2.EQ.0)I1=I1-1
C--MAKE THE BOUNDARY BLANK
        goto 102
        DO 30 I=1,NI
         SCR(I)=1.70141E+38
30      CONTINUE
        DO 31 I=J2-NI+1,J2
         SCR(I)=1.70141E+38
31      CONTINUE
        DO 32 I=J1,J2,NI
         SCR(I)=1.70141E+38
32      CONTINUE
        DO 33 I=J1+NI-1,J2,NI
         SCR(I)=1.70141E+38
33      CONTINUE
102     K1=J1
        DO 2 J=1,NJ
          DO 3 KK=1,I1
            K2=K1+9
            IF(KK.EQ.I1.AND.I2.NE.0) K2=K1+I2-1
            WRITE(N,101)(SCR(I),I=K1,K2)
4           K1=K2+1
3         CONTINUE
          WRITE(N,*)
2       CONTINUE
101   FORMAT(1x,10(1PG14.5E3))
      RETURN
      END

      subroutine readDict
      PARAMETER (IX0=45,IY0=2421)
      PARAMETER (LX=351,LY=2802)
      PARAMETER (lengX=LX-IX0)
      PARAMETER (lengY=LY-IY0)
      dimension idct(IX0:LX,IY0:LY)
      common /dictxy/idct
        WRITE(*,*)' READING AREA EMIS DATA'
        open(1,file='dict.txt',STATUS='unknown')
        idct=0
        read(1,*)
        do
101     read(1,*,end=102)i,j,k
        if((i-IX0)*(i-LX).gt.0)then
        print*,i
        stop
        endif
        if((j-IY0)*(j-LY).gt.0)then
        print*,j
        stop
        endif
        idct(i,j)=k
        enddo
102     close(1)
        do j=LY,IY0,-4
!        write(*,'(80I1)')(idct(i,j),i=IX0,LX,4)
        enddo
        return
        end

      integer function IsLand(i,j)
      PARAMETER (IX0=45,IY0=2421)
      PARAMETER (LX=351,LY=2802)
      dimension idct(IX0:LX,IY0:LY)
      common /dictxy/idct
      IsLand=0
        if((i-IX0)*(i-LX).gt.0)return
        if((j-IY0)*(j-LY).gt.0)return
      IsLand=idct(i,j)
        return
        end


