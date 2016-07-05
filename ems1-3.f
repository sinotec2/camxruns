      CHARACTER*4 MEM(10), MSPEC(10,50) !CAMX4.01
      INTEGER MFID(60)
        parameter (NVARS3D=26)
        parameter (NOSPEC=38)
        CHARACTER SPNAM*4,A1*1
        COMMON /BIGV/SPNAM(10,NOSPEC),mw(NVARS3D),mtx(NOSPEC,NVARS3D)
        real mw,mtx
        parameter(MX=200,MY=200,MZ=50)
        parameter(MROWS=200,MCOLS=200)
        dimension NCOLS(3),NROWS(3),i1(3),j1(3)
      parameter(PHIC=23.60785)!   ; CENTRAL LATITUDE (minus for southern hemesphere)
      parameter(XLONC=120.9908)
      common /saveE/EC(MX,MY,NOSPEC)
      DATA NUEM3/46/
      DATA MEM /1HE, 1HM, 1HI, 1HS, 1HS, 1HI, 1HO, 1HN, 1HS, 1H  /
      data i1/3,2,2/
      data j1/3,2,4/
	call readMTX
        call readTerr(NCOLS,NROWS)
      data NBD/10001/TBEG/0/NED/10365/TEND/2400./
      data XUTM/0/YUTM/0/NZONE/0/
      do igrd=1,3
      write(A1,'(I1)')igrd
      call readIn(igrd,NX,NY,NZ,XORG,YORG,DELTAX,DELTAY)
      ii=i1(igrd)
      jj=j1(igrd)
      call TXTbyMTX(igrd,NCOLS(igrd),NROWS(igrd),ii,jj,NX,NY)
      OPEN(NUEM3, FILE='olds/fortBE.'//A1//'13',FORM='UNFORMATTED'  
     +  ,convert='BIG_ENDIAN',STATUS='UNKNOWN')
        XUTM=XLONC
        YUTM=PHIC
        NZONE=(180+XLONC)/6+1
        NZ=1
        XORG=XORG*1000.
        YORG=YORG*1000.
        DELTAX=DELTAX*1000.
        DELTAY=DELTAY*1000.
      WRITE(NUEM3) MEM, MFID, NSG, NOSPEC, NBD, TBEG, NED, TEND
      WRITE(NUEM3) XUTM,YUTM,NZONE,XORG,YORG,DELTAX,DELTAY,     !FINE
     +             NX,NY,NZ, NVLOW,NVUP,DZSURF,DZMINL,DZMINU    !FINE
      WRITE (NUEM3)1,1,NX,NY
C
C-- NOW WRITE THE SPECIES DEFN RECORD ON BOUNDARY
C
  	do J=1,NOSPEC
	do I=1,10
	      MSPEC(I,J)=SPNAM(I,J)
	enddo
	enddo
      WRITE (NUEM3)((MSPEC(I,J),I=1,10),J=1,NOSPEC)
      JBD=NBD
      JED=NBD
      TB=0 
        sum=0
      do i=1,NX
        do j=1,NY
        do L=1,NOSPEC
        sum=sum+EC(i,j,l)
        enddo
        enddo
        enddo
      do it=0,23
      TE=TB+1
      if(TE.eq.24) then
        JED=JED+1
        TE=0
      endif  
      WRITE (*,*)JBD, TB, JED, TE,sum
      WRITE (NUEM3)JBD, TB, JED, TE
  	do L=1,NOSPEC
      WRITE(NUEM3)NSG,(MSPEC(J,L),J=1,10),((EC(I,J,L),I=1,NX),J=1,NY)
      ENDDO        !L
      JBD=JED
      TB=TE  
      ENDDO        !it
      do 1284 i=43,46
1284  close(i)
      ENDDO        !igrd
      stop
      end
	subroutine TXTbyMTX(igrd,NCOLS,NROWS,i1,j1,NX,NY)
        parameter (NVARS3D=26)
        parameter (NOSPEC=38)
        CHARACTER SPNAM*4,A1*1
        COMMON /BIGV/SPNAM(10,NOSPEC),mw(NVARS3D),mtx(NOSPEC,NVARS3D)
        real mw,mtx
        parameter(MX=200,MY=200,MZ=50)
        parameter(MROWS=200,MCOLS=200)
        common /saveE/EC(MX,MY,NOSPEC)
        dimension tpy(MCOLS,MROWS,NVARS3D)
        write(A1,'(I1)')igrd
        tpy=0.
        EC=0.
        open(1,file='/nfs/REAS/2010/2010REF'//A1//'.txt',status='old')
        read(1,*)
        do i=1,NCOLS
        do j=1,NROWS
          read(1,*)x,y,(tpy(i,j,l),l=1,NVARS3D)
        enddo
        enddo
        unit=1000.*1000./8760
        do i=1,NX
        do j=1,NY
          ii=i-1+i1
          jj=j-1+j1
          do k=1,NOSPEC
          sum=0
          do l=1,NVARS3D
            if(mtx(k,l).eq.0) cycle
            sum=sum+tpy(ii,jj,l)*mtx(k,l)
          enddo
          if(mw(k).ne.0) EC(i,j,k)=sum/mw(k)*unit      !gmole/hr
          enddo!k
        enddo
        enddo
        return
        end
	subroutine readMTX
        parameter (NVARS3D=26)
        parameter (NOSPEC=38)
        CHARACTER SPNAM*4
        COMMON /BIGV/SPNAM(10,NOSPEC),mw(NVARS3D),mtx(NOSPEC,NVARS3D)
        real mw,mtx
        character spnami(-1:NOSPEC)*20,emnam(NVARS3D)*20
        open(1,file='REAS-CB06.csv',status='unknown')
        read(1,*)(SPNAMi(i),i=-1,NOSPEC)
        do j=1,NVARS3D
        read(1,*)emnam(j),mw(j),(mtx(i,j),i=1,NOSPEC)
        do i=1,NOSPEC
          if(mtx(i,j).ne.0)then
                print*,emnam(j),spnami(i),mtx(i,j)
          endif
        enddo
        enddo
        close(1)
        do i=1,NOSPEC
	  do j=1,10
          SPNAM(j,i)=spnami(i)(j:j)
	  enddo
        enddo
        return
        end

        subroutine readTerr(NCOLS,NROWS)
        parameter (KW=3)
        character A130*130,KEYWORD(KW)*20,proj*10,A1*1,root*20
        dimension LENGTH(KW),jj(0:3)
        data KEYWORD/'NESTIX =','NESTJX =','DIS    ='/
        data LENGTH/8,8,8/
        dimension NCOLS(3),NROWS(3)

        open(1,file='t1-4.deck',status='old')
        do
          read(1,'(A)',end=1)A130
          do n=1,KW
            m=LENGTH(n)-1
          do i=1,130-m
            if(A130(i:i).eq.';')exit
            if(A130(i:i+m).eq.KEYWORD(n)(1:m+1))then
                jj(0)=i+m
                k=1
              do j=i+m+1,130-m
                if(A130(j:j).eq.',')then
                  jj(k)=j
                  k=k+1
                  if(k.eq.4)exit
                endif
              enddo
              do k=1,3
              if(n.eq.1)read(A130(jj(k-1)+1:jj(k)-1),*)NROWS(k)
              if(n.eq.2)read(A130(jj(k-1)+1:jj(k)-1),*)NCOLS(k)
              enddo
            endif !is keyword
          enddo !next phrase
          enddo !next keyword
        enddo !next record of input
1       close(1)
        NCOLS(1)=NCOLS(1)+4
        NROWS(1)=NROWS(1)+4
        print*,NCOLS,NROWS
        return
        end
      subroutine readIn(igrd,NX,NY,NZ,XORG,YORG,DELTAX,DELTAY)
      character fname*100,project*10 ,A1*1
      write(A1,'(I1)')igrd  
      open(1,file='/nfs/kuang/20100327CAMx/Domain1/mm05/d'//A1//'.in'
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
