!netcdf includes
      INCLUDE 'PARMS3.EXT'     ! I/O parameters definitions
      INCLUDE 'FDESC3.EXT'     ! file header data structure
      INCLUDE 'IODECL3.EXT'    ! I/O definitions and declarations
      INCLUDE 'STATE3.EXT'     ! M3EXIT status codes
      INCLUDE 'ext/HGRD.EXT'       ! horizontal dimensioning parameters
      INCLUDE 'ext/VGRD.EXT'       ! vertical dimensioning parameters
      INCLUDE 'ext/COORD.EXT'      ! horizontal dimensioning parameters
      INCLUDE 'ext/GC_CONC.EXT'    ! gas chem conc file species and map table e
      CHARACTER*16 PNAME
      CHARACTER*16 FNAME1
      CHARACTER*16 FNAME2
      CHARACTER*16 FNAME3
      CHARACTER*999 XMSG
      DATA         PNAME / 'DRIVER' /
	real GC(NCOLS,NROWS,NLAYS,N_GC_CONC)
	real EM(NCOLS,NROWS,NLAYS)
        narg=iARGc ()
        if(narg.NE.1) then
      WRITE(*,*) ' What is the case name ?'
        write(*,*) ' A '
                stop
        ENDIF
                do i=1,narg
                      call getarg(i,FNAME1)
                enddo
	open(1,file=trim(FNAME1),status='unknown')
	read(1,*)
	do i=1,NCOLS
	do j=1,NROWS
	read(1,*)x,y,(GC(i,j,1,l),l=1,N_GC_CONC)
	enddo
	enddo
	close(1)
!netcdf header
      JDATE=2010001
      JTIME=0
      TSTEP=0	!hh:mm:ss
      call HEADER(JDATE,JTIME,TSTEP)
!netcdf opening file
      LOGDEV = INIT3( )
      print*, NVARS3D
      FNAME1=FNAME1(1:8)//'.nc'
      IF ( .NOT. OPEN3(FNAME1 , FSUNKN3,PNAME)) THEN
         XMSG = 'Could not open '// FNAME1     // ' file'
         CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
         END IF
         IF ( .NOT. DESC3( FNAME1 ) ) THEN
            XMSG = 'Could not get ' // FNAME1 // ' file description'
            CALL M3EXIT( PNAME, STDATE, STTIME, XMSG, XSTAT1 )
            END IF
	DO L=1,N_GC_CONC
!	print *,VNAME3D,GC_CONC( L )
	do j=1,NROWS
	do i=1,NCOLS
	EM(i,j,1)=GC(i,j,1,L)
	enddo
	enddo
C**     WRITE THE RESULTANT ARRAY
            IF ( .NOT. WRITE3( FNAME1, VNAME3D( L ),
     &         JDATE, JTIME, EM) ) THEN
               XMSG = 'Could not WRITE' // VNAME3D( L ) // ' to  '
     &               // FNAME1
               CALL M3EXIT( PNAME, JDATE, JTIME, XMSG, XSTAT1 )
               END IF
	ENDDO	!L

!netcdf closing
      IF ( SHUT3() ) THEN
         XMSG = ' ----  Program  DRIVER  completed successfully  ---- '
         WRITE ( *,'( //5X, A, // )' ) XMSG
         STOP
         ELSE
         XMSG = ' *** FATAL ERROR shutting down Models-3 I/O *** '
         WRITE ( *,'( //5X, A, // )' ) XMSG
         CALL EXIT (1)
         END IF
	stop
	END

      SUBROUTINE HEADER(JDATE,JTIME,TSTEP)
!-I/home/kuang/web31/ioapi3.1/ioapi/
      INCLUDE 'PARMS3.EXT'     ! I/O parameters definitions
      INCLUDE 'FDESC3.EXT'     ! file header data structure
      INCLUDE 'IODECL3.EXT'    ! I/O definitions and declarations
      INCLUDE 'STATE3.EXT'     ! M3EXIT status codes
      INCLUDE 'ext/HGRD.EXT'       ! horizontal dimensioning parameters
      INCLUDE 'ext/VGRD.EXT'       ! vertical dimensioning parameters
      INCLUDE 'ext/COORD.EXT'      ! horizontal dimensioning parameters
      INCLUDE 'ext/GC_CONC.EXT'    ! gas chem conc file species and map table e
      character*16 VGTPUN3D ! = VGTPUN_GD ! currently, not defined

      FTYPE3D = GRDDED3
      SDATE3D = JDATE
      STIME3D = JTIME
      TSTEP3D = TSTEP

      NVARS3D = N_GC_CONC
      NCOLS3D = NCOLS
      NROWS3D = NROWS
      NLAYS3D = NLAYS
      NTHIK3D = NTHIK
      GDTYP3D = GDTYP_GD
      P_ALP3D = P_ALP_GD
      P_BET3D = P_BET_GD 
      P_GAM3D = P_GAM_GD
      XORIG3D = XORIG_GD
      YORIG3D = YORIG_GD
      XCENT3D = XCENT_GD
      YCENT3D = YCENT_GD
      XCELL3D = XCELL_GD
      YCELL3D = YCELL_GD
      VGTYP3D = VGTYP_GD
      VGTOP3D = VGTOP_GD
      VGTPUN3D = VGTPUN_GD ! currently, not defined
      DO L = 1, NLAYS3D + 1
         VGLVS3D( L ) = VGLVS_GD( L )
         END DO
      GDNAM3D = GDNAME_GD

            VDESC3D( 1 ) = 'BC'
            VDESC3D( 2 ) = 'CO2'
            VDESC3D( 3 ) = 'NH3'
            VDESC3D( 4 ) = 'Ethane'
            VDESC3D( 5 ) = 'Propane'
            VDESC3D( 6 ) = 'Butanes'
            VDESC3D( 7 ) = 'Pentanes'
            VDESC3D( 8 ) = 'Other_Alkanes'
            VDESC3D( 9 ) = 'Ethene'
            VDESC3D(10 ) = 'Propene'
            VDESC3D(11 ) = 'Terminal_Alkenes'
            VDESC3D(12 ) = 'Internal_Alkenes'
            VDESC3D(13 ) = 'Acetylene'
            VDESC3D(14 ) = 'Benzene'
            VDESC3D(15 ) = 'Toluene'
            VDESC3D(16 ) = 'Xylenes'
            VDESC3D(17 ) = 'Other_Aromatics'
            VDESC3D(18 ) = 'HCHO'
            VDESC3D(19 ) = 'Oth_Ald_'
            VDESC3D(20 ) = 'Ketones'
            VDESC3D(21 ) = 'Halocarbons'
            VDESC3D(22 ) = 'Others'
            VDESC3D(23 ) = 'Total'
            VDESC3D(24 ) = 'NOX'
            VDESC3D(25 ) = 'OC'
            VDESC3D(26 ) = 'SO2'
         DO L=1,NVARS3D
	    GC_CONC_MAP(L)=L
            VTYPE3D( L ) = M3REAL
            VNAME3D( L ) = VDESC3D( L )        !long_name
            UNITS3D( L ) = 'Ton/Yr'            !units',
            ENDDO
      FDESC3D( 1 ) = 'Annual Emission rate  FIELD'
      DO L = 2, MXDESC3
         FDESC3D( L ) = ' '
         END DO
      RETURN
      END

