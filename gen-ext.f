	parameter (KW=3)
        character A130*130,KEYWORD(KW)*20,proj*10,A1*1,root*20
	dimension LENGTH(KW),x(200,200),y(200,200)
	integer GDTYP_GD,jj(0:3)
	data KEYWORD/'NESTIX =','NESTJX =','DIS    ='/
	data LENGTH/8,8,8/
C   1: LATGRD for lat-lon coordinates (unused)
C   2: LAMGRD for Lambert coordinates
C   3: MERGRD for Mercator coordinates
C   4: STEGRD for Stereographic coordinates
C   5: UTMGRD for UTM coordinates
	dimension NCOLS(3),NROWS(3),NLAYS(3),rlon(3),rlat(3),clon(3),clat(3)
	dimension XORIG_GD(3),YORIG_GD(3)
	dimension XCELL_GD(3),YCELL_GD(3)

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
	      if(n.eq.3)read(A130(jj(k-1)+1:jj(k)-1),*)XCELL_GD(k)
	      YCELL_GD(k)= XCELL_GD(k)
	      enddo !k
	    endif !is keyword
	  enddo !next phrase
	  enddo !next keyword
	enddo !next record of input
1	close(1)
	NCOLS(1)=NCOLS(1)+4
	NROWS(1)=NROWS(1)+4
	do k=1,3
	write(A1,'(I1)')k
	open(1,file='2010REF'//A1//'.txt',status='unknown')
	read(1,*)
	do i=1,NCOLS(k)
	do j=1,NROWS(k)
	read(1,*)X(i,j),Y(i,j)
	enddo
	enddo
4	close(1)
	ic=(NCOLS(K)+1)/2
	jc=(NROWS(K)+1)/2
	delLon=(x(ic+1,jc)-x(ic-1,jc))/2.
	delLat=(y(ic,jc+1)-y(ic,jc-1))/2.
	XORIG_GD(k)=x(ic,jc)-delLon*(ic-1)
	YORIG_GD(k)=y(ic,jc)-delLat*(ic-1)
	enddo
	print*,XORIG_GD,YORIG_GD
	print*,XCELL_GD,YCELL_GD,NCOLS,NROWS,NLAYS,root
	do i=1,3
	write(A1,'(I1)')i
	open(1,file='ext/COORD.bak',status='unknown')	
	open(2,file='ext/COORD.EXT',status='unknown')	
	do 
	  read(1,'(A)',end=2)A130
	  if(A130(7:27).eq.'DATA         P_GAM_GD') then
	    write(A130(31:40),'(F10.5)')XORIG_GD(i)
	  elseif(A130(7:27).eq.'DATA         XCENT_GD') then
	    write(A130(31:40),'(F10.5)')XORIG_GD(i)
	  elseif(A130(7:27).eq.'DATA         YCENT_GD') then
	    write(A130(31:40),'(F10.5)')YORIG_GD(i)
	  elseif(A130(7:27).eq.'DATA         XORIG_GD') then
	    write(A130(31:40),'(F10.5)')XORIG_GD(i)
	  elseif(A130(7:27).eq.'DATA         YORIG_GD') then
	    write(A130(31:40),'(F10.5)')YORIG_GD(i)
	  elseif(A130(7:27).eq.'DATA         XCELL_GD') then
	    write(A130(31:40),'(F10.1)')XCELL_GD(i)*1000.
	  elseif(A130(7:27).eq.'DATA         YCELL_GD') then
	    write(A130(31:40),'(F10.1)')YCELL_GD(i)*1000.
	  endif
	  write(2,'(A)')trim(A130)
	enddo
2	close(1)
	close(2)
	open(1,file='ext/HGRD.bak',status='unknown')	
	open(2,file='ext/HGRD.EXT',status='unknown')
	do 
	  read(1,'(A)',end=3)A130
	  if(A130(7:22).eq.'PARAMETER( NCOLS') then
	    write(A130(25:27),'(I3)')NCOLS(i)
	  elseif(A130(7:22).eq.'PARAMETER( NROWS') then
	    write(A130(25:27),'(I3)')NROWS(i)
	  endif
	  write(2,'(A)')trim(A130)
	enddo
3	close(1)
	close(2)
	result=system("rm -f txt-nc.o")
	result=system("make")
	result=system("mv txt-nc txt_"//A1)
	enddo !next grid
	stop
	end	
