      CHARACTER*4 MEM(10), MSPEC(10,50) !CAMX4.01
      INTEGER MFID(60)
        parameter (NVARS3D=26)
        parameter (NOSPEC=50)
        CHARACTER SPNAM*4,A1*1
        COMMON /BIGV/SPNAM(10,NOSPEC),mw(NVARS3D),mtx(NOSPEC,NVARS3D)
        real mw,mtx
        parameter(MROWS=200,MCOLS=200)
      parameter(PHIC=23.60785)!   ; CENTRAL LATITUDE (minus for southern hemesphere)
      parameter(XLONC=120.9908)
      PARAMETER(IX0=81,IY0=2421)
      PARAMETER(LX=351,LY=2802)
      common /CentE/Em(-MCOLS:MCOLS,-MROWS:MROWS,NOSPEC,24)
      common /saveE/Ec(IX0:LX,IY0:LY,NOSPEC,24)
      dimension tmp(MCOLS,MROWS)
      DATA NUEM1/45/
      DATA NUEM3/46/
      call readDict             !input the Land/water MASK(dict)
      call ReadEmOld            !read previous data file
      do igrd=1,4
      write(A1,'(I1)')igrd      
      call cent(igrd,MX,MY)     !get Taiwan emiss. in igrd coord. sys.
      if(igrd.le.3) then        !read ems1-3 results
      print*,'olds/fortBE.'//A1//'13'  
      OPEN(NUEM1, FILE='olds/fortBE.'//A1//'13',FORM='UNFORMATTED'  
     +  ,convert='BIG_ENDIAN',STATUS='UNKNOWN')
      read(NUEM1) MEM, MFID, NSG, NSPEC, NBD, TBEG, NED, TEND
      read(NUEM1) XUTM,YUTM,NZONE,XORG,YORG,DELTAX,DELTAY,     !FINE
     +             NX,NY,NZ, NVLOW,NVUP,DZSURF,DZMINL,DZMINU    !FINE
      read(NUEM1)idum1,idum1,NX,NY
      read(NUEM1)((MSPEC(I,J),I=1,10),J=1,NSPEC)
      iMove=36000./DELTAX
      jMove=36000./DELTAY
      else                      !the Taiwan Island mesh #4
      call readIn(igrd,NX,NY,NZ,XORG,YORG,DELTAX,DELTAY)
      XORG=XORG*1000.
      YORG=YORG*1000.
      DELTAX=DELTAX*1000.
      DELTAY=DELTAY*1000.
      endif
        XUTM=XLONC
        YUTM=PHIC
        NZONE=(180+XLONC)/6+1
        NZ=1
      print*,'write output file',A1,' headers'
      OPEN(NUEM3, FILE='fortBE.'//A1//'13',FORM='UNFORMATTED'  
     +  ,convert='BIG_ENDIAN',STATUS='UNKNOWN')
      WRITE(NUEM3) MEM, MFID, NSG, NSPEC, NBD, TBEG, NED, TEND
      WRITE(NUEM3) XUTM,YUTM,NZONE,XORG,YORG,DELTAX,DELTAY,     !FINE
     +             NX,NY,NZ, NVLOW,NVUP,DZSURF,DZMINL,DZMINU    !FINE
      WRITE (NUEM3)1,1,NX,NY
      WRITE (NUEM3)((MSPEC(I,J),I=1,10),J=1,NSPEC)
      ic=(NX+1)/2               !define the center of domain
      jc=(NY+1)/2
      JBD=NBD
      JED=JBD
      DO it=1,24
      TB=IT-1
      TE=MOD(IT,24)
      if(igrd.le.3) read (NUEM1)JBD, TB, JED, TE
      WRITE(NUEM3)JBD, TB, JED, TE
      WRITE (*,*)JBD, TB, JED, TE
      do L=1,NSPEC
      if(igrd.le.3)then
        read(NUEM1)NSG,(MSPEC(J,L),J=1,10),((tmp(I,J),I=1,NX),J=1,NY)
        r=1.0
        if(L.eq.1.or.L.eq.2)r=0.5 !nox 
        if(L.eq.4)r=0.2 !sox 
        do i=NX,iMove+1,-1 !move north-eastern for 0.5 deg
        do j=NY,jMove+1,-1
          tmp(i,j)=tmp(i-iMove,j-jMove)*r   !cut for point source
        enddo
        enddo
      endif !igrd le 3
      do ix=-MX,MX              !center is at (0,0)
      do iy=-MY,MY
	if((ic+ix).le.0.or.(ic+ix).gt.NX) cycle
	if((jc+iy).le.0.or.(jc+iy).gt.NY) cycle
      if(igrd.le.2.and.Em(ix,iy,L,it).le.0)cycle !1/2 only land replace
      if(igrd.eq.3.and.(16/12*(ic+ix)-(jc+iy)+49).le.0) cycle
        tmp(ic+ix,jc+iy)=amax1(0.,Em(ix,iy,L,it))!3/4 all replace
      ENDDO
      ENDDO
      WRITE(NUEM3)NSG,(MSPEC(J,L),J=1,10),((tmp(I,J),I=1,NX),J=1,NY)
      ENDDO        !L
      ENDDO        !it
      do 1284 i=43,46
1284  close(i)
      ENDDO        !igrd
      stop
      end
      subroutine readIn(igrd,NX,NY,NZ,XORG,YORG,DELTAX,DELTAY)
      character fname*100,project*10 ,A1*1
      write(A1,'(I1)')igrd  
      open(1,file='/st0/Out08/d'//A1//'.in'
     +,status='unknown')
      do i=1,2
      read(1,*)
      enddo          
      read(1,'(20x,a10)') project
      write(*,'(a,1x,a10)') '                     Projection:',project
      do i=4,7
      read(1,*)
      enddo          
      read(1,'(20x,a)') fname !line 8th
      read(fname,*) NX,NY,NZ
      write(*,'(a,3i10)')'                 CAMx grid size:',NX,NY,NZ
      read(1,'(20x,a)') fname
      read(fname,*)DELTAX,DELTAY
       write(*,'(a,2f10.3)')
     &                  '              CAMx grid spacing:',DELTAX,DELTAY
       if (project.eq.'LCP       ') then
        iproj = 2
        read(1,'(20x,a)') fname
        read(fname,*) XORG,YORG,clonin,clatin,tlat1in,tlat2in
        write(*,'(a,2f10.3)')
     &                  '    CAMx LCP Origin (SW corner):',XORG,YORG
        write(*,'(a,4f10.3,/)')
     &  '    CAMx LCP Projection Params :',clonin,clatin,tlat1in,tlat2in
        endif 
        close(1)
        return
        end
