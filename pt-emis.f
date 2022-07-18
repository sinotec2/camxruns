      program REAS_pt
      CHARACTER*4 MEM(10), MSPEC(10,50),ASPEC(10,50) !CAMX4.01
      INTEGER MFID(60),LL(50)
      PARAMETER (MP=80000)
      character*15 pid(MP)
      character*4 MPTS(10) !camx4.01
      character*15 plantid(80000)
      DIMENSION X(MP),Y(MP),IX(MP),IY(MP),D(MP),H(MP)
      DIMENSION T(MP),V(MP)
      DIMENSION ILOC(MP,24), IJPS(MP,24), KPTS(MP,24)
      DIMENSION FLOW(MP,24),EFPLH(MP,24)
      DIMENSION IPOPT(6),QPTS(50,MP)
      DIMENSION Q2(MP,24)
      DIMENSION HD(MP)    ! HOUR OF DAY
      DIMENSION QP(13,MP) ! EMISSION IN g-mole/HR
      COMMON /NT24/NTIME(24,24)
      parameter(MROWS=200,MCOLS=200)
      dimension tmp(MCOLS,MCOLS)
      dimension ems(MCOLS,MCOLS,3)
      ILOC=0
      IJPS=0
      KPTS=0
      FLOW=0.
      EFPLH=0.
      OPEN(1,FILE='PTSOURCE.BIN',FORM='UNFORMATTED',STATUS='UNKNOWN')
      read(1)pid,X,Y,H,D,T,V
      read(1)QP,HD,NOPTS
      read(1)NTIME
      close(1)
      data PHIC/23.61/!   ; CENTRAL LATITUDE (minus for southern hemesphere)
      data XLONC/120.99/
      iutmzon=(180+XLONC)/6+1
      call utmgeo(0,iutmzon,rx4,ry4,XLONC,PHIC)
      Xcent=rx4*1000.           !center point UTM-m
      Ycent=ry4*1000.
      print*,Xcent,Ycent
      do ip=1,NOPTS
!      print*,pid(ip),x(ip),y(ip),qp(ip,1),h(ip)
        X(ip)=X(ip)-Xcent       !new coord. vs center point
        Y(ip)=Y(ip)-Ycent
!      print*,pid(ip),x(ip),y(ip),qp(ip,1),h(ip)
      enddo
      call readIn(4,NX,NY,NZ,XORG,YORG,DELTAX,DELTAY)
      XORGD4=XORG*1000.
      YORGD4=YORG*1000.
      XLEND4=XORGD4+NX*DELTAX*1000.
      YLEND4=YORGD4+NY*DELTAY*1000.
      OPEN(NUEM1, FILE='../emission/olds/fortBE.113',FORM='UNFORMATTED'
     +  ,convert='BIG_ENDIAN',STATUS='UNKNOWN')
      read(NUEM1) MEM, MFID, NSG, NSPEC, NBD, TBEG, NED, TEND
      read(NUEM1) XUTM,YUTM,NZONE,XORG,YORG,DELTAX,DELTAY,
     +             NX,NY,NZ, NVLOW,NVUP,DZSURF,DZMINL,DZMINU
      print*, XUTM,YUTM,NZONE,XORG,YORG,DELTAX,DELTAY
      read(NUEM1)idum1,idum1,NX,NY
      read(NUEM1)((MSPEC(I,J),I=1,10),J=1,NSPEC)
        do i=1,10
        ASPEC(i,1)=MSPEC(i,1)
        ASPEC(i,2)=MSPEC(i,2)
        ASPEC(i,3)=MSPEC(i,4)
        enddo
!     DO it=1,24
!     goto 104
      read(NUEM1)JBD, TB, JED, TE
      do L=1,NSPEC
      read(NUEM1)NSG,(MSPEC(J,L),J=1,10),((tmp(I,J),I=1,NX),J=1,NY)
      K=0
      if(L.eq.1)K=1
      if(L.eq.2)K=2
      if(L.eq.4)K=3
      if(K.eq.0)cycle
        r=0.2 !for the SOx, 80% as point sources
        if(L.eq.1.or.L.eq.2)r=0.5 !for the NOx, half as point sources
      do i=1,NX
      do j=1,NY
        ems(i,j,K)=0.
        xx=(i-1)*deltax+xorg
        yy=(j-1)*deltay+yorg
!       print*,xx,yy,XORGD4,XLEND4,YORGD4,YLEND4
        if(((XX-XORGD4)*(XX-XLEND4).le.0).and.
     *   ((YY-YORGD4)*(YY-YLEND4).le.0))then
!       print*,xx,yy
        cycle
        endif
        ems(i,j,K)=tmp(i,j)*(1-r)
        enddo
        enddo
      ENDDO        !L
      close(NUEM1)
      u1=46*8760/1000/1000
      u2=64*8760/1000/1000
      k=NOPTS+1
      call random_seed()
      do i=1,NX
      do j=1,NY
      sum=sum+(ems(i,j,1)+ems(i,j,2))*u1+ems(i,j,3)*u2
        if(sum.eq.0) cycle
        h(k)=250.
        if(sum.le.1000) h(k)=50+sum*0.05 !max=100m
        if(sum.le.100) h(k)=10+sum*0.8 !max=50m
        call CORRECT(D(k),H(k),T(k),Q0,V(k))    !Q0=0
        Xsw=XORG+(i-1)*DELTAX
        Ysw=YORG+(j-1)*DELTAY
        call  RANDOM_NUMBER(x1) !0-1 noise
        call  RANDOM_NUMBER(y1)
        X(k)=Xsw+x1*DELTAX
        Y(k)=Ysw+y1*DELTAY      !in m
        HD(k)=8
        if(H(k).gt.50)HD(k)=24
        do l=1,3
          QP(l,k)=ems(i,j,l)
        enddo
        k=k+1
        enddo !next j
        enddo !next i
        NOPTS=k-1
104   continue !writing
      DATA MPTS / 1HP, 1HT, 1HS, 1HO, 1HU, 1HR, 1HC, 1HE, 1H , 1H /
      DATA NUPTS /14/
      OPEN(NUPTS, FILE='fortBE.14',FORM='UNFORMATTED'
     ,  ,convert='BIG_ENDIAN',STATUS='UNKNOWN')
      WRITE (NUPTS)    MPTS, MFID, NSG, NSPEC, NBD, TBEG, NED, TEND
      WRITE (NUPTS) XUTM, YUTM, NZONE, XORG, YORG, DELTAX, DELTAY,
     $         NOXG, NOYG, NOZG
     $, NVLOW, NVUP, DZSURF, DZMINL, DZMINU
      WRITE (NUPTS) 1,1,NOXG,NOYG
      WRITE (NUPTS) ((MSPEC(I,J),I=1,10),J=1,NSPEC)
      WRITE (NUPTS) NSEG, NOPTS
      WRITE (NUPTS) (X(K),Y(K),H(K),D(K),T(K),V(K),K=1,NOPTS)

      DO 680 IT=1,24  !▒p▒ɰj▒▒
        TB=IT-1
        TE=mod(IT,24)
        JED=NBD
        if(TE.eq.0)JED=JED+1
        WRITE ( *,* ) NBD, TB  , JED, TE
        WRITE (NUPTS) NBD, TB  , JED, TE
        WRITE (NUPTS) NSEG, NOPTS
        DO 410 IP=1,NOPTS
          if(NTIME(INT(HD(IP)),IT).EQ.0)then
            DO I=1,NSPEC
              QPTS(I,IP)=0.
            ENDDO
          ELSE
c            QPTS( 1,IP)=Q2*9./10.
c            QPTS( 2,IP)=Q2*1./10.
            QPTS( 1,IP)=QP(2 ,IP)* 8./10.
            QPTS( 2,IP)=QP(2 ,IP)* 2./10.
            QPTS( 4,IP)=QP(11,IP)              !ETH
            QPTS( 5,IP)=QP(5 ,IP)              !OLE
            QPTS( 6,IP)=QP(6 ,IP)              !PAR
            QPTS( 7,IP)=QP(7 ,IP)              !TOL
            QPTS( 8,IP)=QP(8 ,IP)              !XYL
            QPTS( 9,IP)=QP(9 ,IP)              !FORM
            QPTS(10,IP)=QP(10,IP)              !ALD2
            QPTS(19,IP)=QP(12,IP)              !ISOP
            QPTS(20,IP)=QP(4 ,IP)              !CO
            QPTS(24,IP)=QP(3 ,IP)              !SO2
            QPTS(27,IP)=QP(1 ,IP) *0.1 !FCRS
            QPTS(28,IP)=QP(1 ,IP) *0.0 !CCRS
            QPTS(29,IP)=QP(1 ,IP) *0.8 !FPRM
            QPTS(30,IP)=QP(1 ,IP) *0.1 !CPRM
          ENDIF
410     CONTINUE
        call  NewQPTS(IT,QPTS)
        WRITE (NUPTS)(ILOC(IP,IT),IJPS(IP,IT),KPTS(IP,IT),FLOW(IP,IT),
     $    EFPLH(IP,IT), IP=1,NOPTS)
        DO 500 I=1,NSPEC
          WRITE (NUPTS) NSEG,(MSPEC(J,I),J=1,10),
     $      (QPTS(I,IP), IP=1,NOPTS)
500     CONTINUE
680   CONTINUE
      STOP
      end
      subroutine NewQPTS(IT,QPTS)
      INTEGER LL(50)
      PARAMETER (MP=80000)
      DIMENSION QPTS(50,MP),tmp(38)
      common /LL31to8/LL
      if(IT.eq.1) then
      open(1,file='31to38.txt',status='unknown')
      do i=1,31
        read(1,*,end=1)LL(i)          !mapping of 38(new) to 31(old) system
      enddo
1     close(1)
      print*,LL
      endif
      do ip=1,MP
      tmp=0
      do l=1,31
       tmp(LL(l))=QPTS(l,ip)
      enddo
      do l=1,38
       QPTS(l,ip)=tmp(l)
      enddo
      enddo
      return
      end
C**********************************************************************
       SUBROUTINE CORRECT(DIA,HE,TMP,Q,VS)
C**▒▒▒Ƶ{▒▒▒N▒ɥ▒▒ϹD▒▒▒Ѽ▒
C**▒▒▒޿謰▒G 1.▒ϧw▒▒▒׬▒▒̰򥻪▒▒▒ơA▒▒▒M▒▒▒▒▒T
C**             2.▒L▒f▒|▒▒ƪ▒▒p▒▒w▒ʦ▒▒▒
C**             3.▒L▒ū׸▒▒(▒Τ▒▒▒▒Ф▒▒ū▒)
C**               ▒▒▒P▒ϧw▒W▒ұ▒▒@▒▒▒▒▒▒
C**             4.▒L▒y▒t▒▒▒(▒Τ▒▒▒▒Ф▒▒ū▒)
C**               ▒▒▒P▒ϧw▒W▒ұ▒▒@▒▒▒▒▒▒
                TEMP=TMP-273.
                if (dia.le.0) then
                    dia=he*0.3
                endif
                if (TEMP .le. 100) then
                    if (he.le.12.5) then
                        TEMP=444-273
                    elseif (he.le.50 ) then
                        TEMP=449-273
                    elseif (he.le.100) then
                        TEMP=475-273
                    elseif (he.le.150) then
                        TEMP=409-273
                    elseif (he.le.200) then
                        TEMP=403-273
                    elseif (he.gt.200) then
                        TEMP=363-273
                    endif
                endif
                tmp=real(TEMP+273)
                vs=q*4.*tmp/273./3.14159/dia/dia/60.
                if (he.le.12.5) then
                    if (vs.lt.4.5 .or. vs.gt.8.5 ) then
                        vs=6.5
                    endif
                elseif (he.le.50 ) then
                    if (vs.lt.5.3 .or. vs.gt.9.3 ) then
                        vs=7.3
                    endif
                elseif (he.le.100) then
                    if (vs.lt.7.1 .or. vs.gt.11.1) then
                        vs=9.1
                    endif
                elseif (he.le.150) then
                    if (vs.lt.14.1.or. vs.gt.18.1) then
                        vs=16.1
                    endif
                elseif (he.le.200) then
                    if (vs.lt.15.4.or. vs.gt.19.4) then
                        vs=17.4
                    endif
                elseif (he.gt.200) then
                    if (vs.lt.18.0.or. vs.gt.22.0) then
                        vs=20.0
                    endif
                endif
                RETURN
                END
      subroutine readIn(igrd,NX,NY,NZ,XORG,YORG,DELTAX,DELTAY)
      character fname*100,project*10 ,A1*1
      write(A1,'(I1)')igrd
      open(1,file='/nfs/kuang/20100327CAMx/Domain4/mm05/d'//A1//'.in'
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
