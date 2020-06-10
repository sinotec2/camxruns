        parameter (nrec=10000)
        character TYPE(nrec)*7
        character type_tmp*7
        character UTME(nrec)*6
        character UTMN(nrec)*7
        character DICT(nrec)*4
        real em(10,nrec)        !TSP/PM(1-3)/SOX/NOX/THC/NMHC/CO/PB
        real em2(7,nrec)        !SOX/TSP/PM/CO/NOX/THC/NMHC
c        open(1,file='teds71_area_grid.txt',status='old')
        open(1,file='teds90_areagrid_twd97.sdf',status='old')
        open(2,file='102ag-x.bin',status='unknown',form='unformatted')
        irec=0
        sum=0. !pei
3       do i=1,nrec
                TYPE(i)=''
                UTME(i)=''
                UTMN(i)=''
                DICT(i)=''
                do k=1,10
                        em(k,i)=0
                enddo
        enddo

        do i=1,nrec
          read(1,100,end=2)TYPE(i),UTME(i),UTMN(i),DICT(i),
     +      (EM(k,i),k=1,10)
100       format(a7,a6,a7,a4,10f10.6)
          it=1
          type_tmp=''
          do j=1,7
            if(TYPE(i)(j:j).ne.' ') then
              type_tmp(it:it)=TYPE(i)(j:j)
              it=it+1
            endif
          enddo
          TYPE(i)=trim(type_tmp)
         
          j=i
          do ii=1,10
            if(em(ii,i).lt.0.)then
              write(*,*)irec+j,em(ii,i)
              em(ii,i)=0.
            endif
          enddo
        enddo
2       if(i.le.1) then
                close(1)
                close(2)
                stop
        endif
        do i=1,nrec
        em2(1,i)=em(5,i)  !sox
        em2(2,i)=em(1,i)  !tsp    
        em2(3,i)=em(2,i)  !pm10
        em2(4,i)=em(9,i)  !co
        em2(5,i)=em(6,i)  !nox
        em2(6,i)=em(7,i)  !thc
        em2(7,i)=em(8,i)  !nmhc
        sum=sum+em2(7,i) !pei
        enddo
        write(2)TYPE,UTME,UTMN,DICT,EM2
        irec=irec+j
        write(*,*)irec
        write(*,*)'nmhc sum=',sum !pei
        goto 3  
        end
