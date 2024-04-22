      program quadTree
      implicit none
      
      INTEGER :: capacity,i,j,pt1,pt2,lpt,lptmax,ip,ll,idiv,jdiv,lvlcellmax,limitlvl,io,nlog2_X,nlog2_Y
      INTEGER :: loop,l,nc,ncmax,npts,newnc,cmax,cellptsmax,di,i1,i2,i3,i4,ipts,strlen,node
      DOUBLE PRECISION :: int_pt,dec_pt,r,xi,yi,xe,ye,dx,dy,cx,cy,dum,lx,ly
      DOUBLE PRECISION,dimension(4) :: xtemp,ytemp
      INTEGER,pointer,dimension(:) :: cellpts,lvlcell,nlvlcell,Pcellnpts
      INTEGER,pointer,dimension(:,:) :: ptconnect,ptid
      DOUBLE PRECISION,pointer,dimension(:) :: xpts,ypts
      DOUBLE PRECISION,pointer,dimension(:,:) :: xp,yp
      DOUBLE PRECISION,pointer,dimension(:,:) :: nxp,nyp
      DOUBLE PRECISION,pointer,dimension(:,:) :: pt
      CHARACTER(200) :: fn,chb,chn,che,chsid,chsol,filename,outname,step_to_msh
      LOGICAL :: run,export_grid
      
      limitlvl=7
      export_grid = .TRUE.
      capacity=1
!      npts=100
      xi=0 ; xe=0.5
      yi=0 ; ye=0.5
      cx=(xe-xi)/2
      cy=(ye-yi)/2
      dx=0.01  ; dy=0.005

      filename='BulbDynamicBlades_Solid_Shell_10mm.STEP'

      ! GENERATE THE MSH IF THE IMPORT FILE IS STEP FILE:
      strlen = len(trim(adjustl(filename)))
      if(filename(strlen-3:) .eq. 'STEP')  filename = step_to_msh(dx,dy,filename)

      write(6,'(a)') trim(adjustl(filename))


      ! CREATE THE INITAL 2D FIELD, FIRST PARENT OF QUADTREE:
      call geometry_bound(filename,dx,dy,xi,xe,yi,ye,node,fn)    ! input filename

      idiv=int((xe-xi)/dx) 
      jdiv=int((ye-yi)/dy)
      write(6,'(a,2I10)') 'idiv,jdiv:',idiv,jdiv
      
      nlog2_X=int(log(idiv*1.d0)/log(2.d0))+1
      nlog2_Y=int(log(jdiv*1.d0)/log(2.d0))+1
      limitlvl = max(nlog2_X,nlog2_Y)
      write(6,'(a,3I10)') 'nlog2_X,nlog2_Y,limitlvl',nlog2_X,nlog2_Y,limitlvl
      
      idiv = 2**(limitlvl)
      jdiv = 2**(limitlvl)
      write(6,'(a,2I10)') 'idiv,jdiv:',idiv,jdiv
      
      
      lx = idiv*dx ; ly = jdiv*dy
      write(6,'(a,2F9.3)') 'lx,ly:',lx,ly

      idiv=1
      jdiv=1
      ncmax=(idiv)*(jdiv)

      write(6,'(a,3I10)') trim(adjustl(fn)), idiv,jdiv, ncmax

!      STOP
      
      allocate(xp(ncmax,2),yp(ncmax,2),cellpts(ncmax),ptconnect(ncmax,4))!,pt(npts,2))
      allocate(xpts(ncmax*4),ypts(ncmax*4),lvlcell(ncmax))
      
      ipts=0
      do j=1,jdiv ; do i=1,idiv
	ipts=ipts+1
	xp(ipts,1)=xi+(i-1)*lx ; xp(ipts,2)=xi+i*lx
	yp(ipts,1)=yi+(j-1)*ly ; yp(ipts,2)=yi+j*ly
	lvlcell(ipts)=1
	WRITE(6,'(I4,4F10.3)') ipts,xp(ipts,1),xp(ipts,2),yp(ipts,1),yp(ipts,2)
      end do ; end do
      if(ipts.ne.ncmax) WRITE(6,*) 'ERROR'
!      STOP

      ! IMPORT .DAT GEOMETRY:    
      write(6,'(a)') 'Improting '//trim(adjustl(fn))//' ...'
      open(20,file=fn)
        read(20,*) npts ; allocate(pt(npts,2))
        do l=1,npts
         read(20,*) pt(l,1),pt(l,2),dum
         pt(l,1)=pt(l,1)!+cx
         pt(l,2)=pt(l,2)!+cy
!	  write(6,'(I4,2F10.3)') l,(pt(l,di),di=1,2)
        end do
      close(20)
  
      run=.TRUE.
      newnc=ncmax
      loop=-1
      DO WHILE (run.eqv..TRUE.)
        loop=loop+1
	 cellptsmax=0
	 cellpts=0
        lvlcellmax=0
	
	 WRITE(6,*)
	 WRITE(6,'(a,I4)') 'loop:',loop

	 WRITE(6,'(a,3I4)') 'arrange data',loop,lptmax,ncmax
	 
        if(export_grid.eq..TRUE.) CALL grid_export(ncmax,loop,xp,yp)
!        STOP

        ! cell
        lvlcellmax=maxval(lvlcell)
	 WRITE(6,'(a,3I4)') 'loop,lvlmax,limit',loop,lvlcellmax,limitlvl
        if(lvlcellmax.lt.limitlvl) then
	   do nc=1,ncmax
	     do l=1,npts
!	       WRITE(6,'(I4,6F10.6)') l,pt(l,1),xp(nc,1),xp(nc,2),pt(l,2),yp(nc,1),yp(nc,2)
	       if(pt(l,1).gt.xp(nc,1) .and. pt(l,1).le.xp(nc,2) .and. pt(l,2).gt.yp(nc,1) .and. pt(l,2).le.yp(nc,2)) then
		  cellpts(nc)=cellpts(nc)+1
	       end if
	     end do

	     !  CALCULATES THE NUMBER OF NEW PARENTS:
	     if(cellpts(nc).gt.capacity)then
	       newnc=newnc+3
	     end if
	   end do
        else 
	   do nc=1,ncmax
	     do l=1,npts
!	       WRITE(6,'(I4,6F10.6)') l,pt(l,1),xp(nc,1),xp(nc,2),pt(l,2),yp(nc,1),yp(nc,2)
	       if(pt(l,1).gt.xp(nc,1) .and. pt(l,1).le.xp(nc,2) .and. pt(l,2).gt.yp(nc,1) .and. pt(l,2).le.yp(nc,2)) then
!		  Pcellnpts(nc)=Pcellnpts(nc)+1
!                ptid(nc,Pcellnpts(nc))=l   ! save the index of the point of interest
	       end if
	     end do
          end do
!          CALL IBPs_removal
	   WRITE(6,*) 'Maximum cell refinement was reached:',loop
	   WRITE(6,*) 'Program exit!'
	   STOP
        end if

	 allocate(nxp(newnc,2),nyp(newnc,2),xpts(newnc*4),ypts(newnc*4),nlvlcell(newnc))

	 ipts=1 ; lpt=0
	 do nc=1,ncmax 
	   ! CALCULATES THE BOUNDS OF THE NEW PARENTS:
	   if(cellpts(nc).le.capacity) then
	     nxp(ipts,1)=xp(nc,1)    ; nxp(ipts,2)=xp(nc,2)  ; nyp(ipts,1)=yp(nc,1)   ; nyp(ipts,2)=yp(nc,2)
	     nlvlcell(ipts)=lvlcell(nc)
	     ipts=ipts+1
	   
          else if(cellpts(nc).gt.capacity)then
	     WRITE(6,'(a,2I4,4F10.3)') 'divide,nc,cellnpts',nc,cellpts(nc),xp(nc,1),xp(nc,2),yp(nc,1),yp(nc,2)
	     
            i1=ipts ; i2=ipts+1 ; i3=ipts+2 ; i4=ipts+3
	    
	     nxp(i1,1)=xp(nc,1)			      ; nxp(i1,2)=xp(nc,1)+(xp(nc,2)-xp(nc,1))/2  
	     nyp(i1,1)=yp(nc,1)			      ; nyp(i1,2)=yp(nc,1)+(yp(nc,2)-yp(nc,1))/2
	     nlvlcell(i1)=lvlcell(nc)+1
	    
	     nxp(i2,1)=xp(nc,1)+(xp(nc,2)-xp(nc,1))/2   ; nxp(i2,2)=xp(nc,2)    
	     nyp(i2,1)=yp(nc,1)			      ; nyp(i2,2)=yp(nc,1)+(yp(nc,2)-yp(nc,1))/2
	     nlvlcell(i2)=lvlcell(nc)+1
	    
	     nxp(i3,1)=xp(nc,1)			      ; nxp(i3,2)=xp(nc,1)+(xp(nc,2)-xp(nc,1))/2  
	     nyp(i3,1)=yp(nc,1)+(yp(nc,2)-yp(nc,1))/2   ; nyp(i3,2)=yp(nc,2)
	     nlvlcell(i3)=lvlcell(nc)+1
	    
	     nxp(i4,1)=xp(nc,1)+(xp(nc,2)-xp(nc,1))/2   ; nxp(i4,2)=xp(nc,2)    
	     nyp(i4,1)=yp(nc,1)+(yp(nc,2)-yp(nc,1))/2   ; nyp(i4,2)=yp(nc,2)
	     nlvlcell(i4)=lvlcell(nc)+1
	    
	     ipts=ipts+4
	   end if
	 end do
      
!       do nc=1,newnc
!	  WRITE(6,'(a,I4,4F10.3)') 'nc,nxic,nxec,nyic,nyec:',nc,nxp(nc,1),nxp(nc,2),nyp(nc,1),nyp(nc,2)
!	end do

	 deallocate(xp,yp,cellpts,ptconnect,xpts,ypts,lvlcell)

	 allocate(xp(newnc,2),yp(newnc,2),cellpts(newnc),ptconnect(newnc,4),lvlcell(newnc))
	 allocate(xpts(newnc*4),ypts(newnc*4))
      
	 xp(:,:)=nxp(:,:) ; yp(:,:)=nyp(:,:) ; lvlcell(:)=nlvlcell(:)

	 do nc=1,newnc
	   WRITE(6,'(a,I4,4F10.3,I4)') 'nc,xic,xec,yic,nyec:',nc,xp(nc,1),xp(nc,2),yp(nc,1),yp(nc,2),lvlcell(nc)
	 end do

	 deallocate(nxp,nyp,nlvlcell)

	 ncmax=newnc

      END DO

      end program
!#######################################################################
      FUNCTION step_to_msh(dx,dy,filename) result(fn)

! Generate a .msh file from your CAD .step file using GMsh. 
! Run GMsh in bash using a slightly more refined mesh that your target
! resolution. Return the name of the file with the .msh extension.
!#######################################################################
      implicit none
      
      DOUBLE PRECISION,intent(in) :: dx,dy
      CHARACTER(200),intent(in) :: filename

      INTEGER :: sn
      DOUBLE PRECISION :: deltamin,delta
      CHARACTER(200) :: fn,chd,chd1
      CHARACTER(600) :: command
      
      write(6,*)     ! USER OUTPUT
      write(6,'(a)') '--------------------------------------'
      write(6,'(a)') ' >>>> step file to msh using GMsh >>>>'
      write(6,'(a)') '--------------------------------------'

      ! DETERMINE FINEST RESOLUTION: 
      deltamin=min(dx,dy)

      ! DETERMINE COEFFICIENT RESOLUTION FOR GMESH:
      delta = deltamin * 1000.d0 / 2.d0  ! m -> mm for gmsh
 
      write(chd,'(F10.3)') delta
      write(chd1,'(F10.3)') delta/1000.d0

      sn = len(trim(adjustl(filename)))
      fn = trim(adjustl(filename(:sn-5)))

      ! RUN THE GMSH COMMAND INTO TERMINAL TO CREATE .MSH FILE:
      command = 'gmsh -match  -clmin 0 -clmax '//trim(adjustl(chd))//' -1 -2 -3 ' &
                //trim(adjustl(fn))//'.STEP | tee '//trim(adjustl(fn))//'.log >> /dev/null'
      WRITE(6,'(a)') 'GMesh Resolution: '//trim(adjustl(adjustl(chd)))//' [mm] | '//trim(adjustl(chd1))//' [m]'
      WRITE(6,'(a)') 'GMesh Command: '//trim(adjustl(command))
      CALL SYSTEM(trim(adjustl(command)))
      
      ! RESULT NAME:
      fn = trim(adjustl(fn))//'.msh'

      END FUNCTION step_to_msh
!######################################################################
      SUBROUTINE geometry_bound(filename,dx,dy,xi,xe,yi,ye,node,fn)

!     Read GMsh File or dat file to evaluate the bound of the geometry
!######################################################################
      implicit none
      
      integer,intent(out) :: node
      double precision,intent(in) :: dx,dy
      double precision,intent(out) :: xi,xe,yi,ye
      character(200),intent(in) :: filename
      character(200),intent(out) :: fn

      integer :: ipts,ref1,ref2,ref3,filenodes ! FLAG RECENTER
      integer :: l,nodepackages,n,sn,dumlines
      integer :: entity1,entity2,entity3,entity4,entities,dum
      DOUBLE PRECISION :: xmin,xmax,ymin,ymax,zmin,zmax
      DOUBLE PRECISION :: zi,ze,Cxx,Cyy,Czz
      DOUBLE PRECISION,pointer,dimension(:) :: xf,yf,zf
      
      write(6,*)     ! USER OUTPUT
      write(6,'(a)') '--------------------------'
      write(6,'(a)') ' >>>> Geometry Bounds >>>>'
      write(6,'(a)') '--------------------------'

      sn = len(trim(adjustl(filename))) 
      
      write(6,'(a,I10)') trim(adjustl(filename(sn-3:))),sn

      IF(filename(sn-3:).eq.'.msh') THEN
      
        open(20,file=filename)
          read(20,*) ! $MeshFormat
          read(20,*) ! 4.1 0 8
          read(20,*) ! $EndMeshFormat
          read(20,*) ! $Entities

          ! READ THE ENTITIES SECTION:
          read(20,*) entity1,entity2,entity3,entity4
          entities=entity1+entity2+entity3+entity4
          do l=1,entities
           read(20,*) 
          end do
          read(20,*) !#EndEntities'
          
          ! READ THE NODES SECTION:
          read(20,*) !#Nodes'
          read(20,*) nodepackages,filenodes
          allocate(xf(filenodes),yf(filenodes),zf(filenodes))
          ipts=0
          do n=1,nodepackages
            read(20,*) dum,dum,dum,node
            do dumlines=1,node     ! read all node ID
              read(20,*) 
            end do
            do l=1,node	      
              ipts=ipts+1	  
              read(20,*) xf(ipts),yf(ipts),zf(ipts)
            end do 
          end do
        close(20)

        node = ipts
        
        Cxx=0.d0 ; Cyy=0.d0 ; Czz=0.d0

        ! CHANGE THE MM UNIT INTO THE M UNIT AND THE CENTROID LOCATION OF THE GEOMETRY
        xf(:)=xf(:)/1000+Cxx ; yf(:)=yf(:)/1000+Cyy ; zf(:)=zf(:)/1000+Czz
  !	 
        if(ipts.ne.filenodes) write(6,*) 'Error more nodes then expected'
        close(20) 
        
        fn = trim(adjustl(filename(:sn-4)))//'.xyz'
        write(6,'(a)') trim(adjustl(fn))
        open(20,file=fn) 
          write(20,*) node
          do l=1,node
            write(20,'(3F20.8)') xf(l),yf(l),zf(l)
          end do
        close(20)
      
      ELSE
        
        open(20,file=filename)   ! IF THE FILE IS 2D
          read(20,*) node
          do l=1,node
	     read(20,*) xf(l),yf(l),zf(l)
          end do
        close(20)

      END IF
      
      ! PERMUTE THE AXIS TO ENSURE THE CORRECT VALUES ARE ASSIGNED TO XF,YF,ZF
      call axis_permutation(node,xf,yf,zf)

      ! EVALUATE THE BOUNDARY BOX OF THE GEOMETRT OF THE .MSH:
      xmin=minval(xf) ; xmax=maxval(xf)
      ymin=minval(yf) ; ymax=maxval(yf)
      zmin=minval(zf) ; zmax=maxval(zf)
      
      ! CREATE A BOUNDARY THAT MULTIPLY OF THE TARGETED MESH:
      xi=int(xmin/dx)*dx-dx
      xe=int(xmax/dx)*dx+dx
      yi=int(ymin/dy)*dy-dy
      ye=int(ymax/dy)*dy+dy

!      zi=int(zmin/dz)*dz-dz
!      ze=int(zmax/dz)*dz+dz

      WRITE(6,'(a,6F10.5)') ' Bounds:',xi,xe,yi,ye!,zi,ze	
      
      RETURN
      END SUBROUTINE geometry_bound
! ######################################################################
      SUBROUTINE axis_permutation(node,xf,yf,zf)
! ######################################################################
      implicit none
      
      INTEGER,intent(in) :: node
      DOUBLE PRECISION,dimension(node),intent(inout) :: xf,yf,zf
      
      INTEGER :: ref1,ref2,ref3
      DOUBLE PRECISION,dimension(node) :: xtemp,ytemp,ztemp
      
      ref1 = 1 ; ref2 = 2 ; ref3 = 3
      xtemp(:) = xf(:) ; ytemp(:) = yf(:) ; ztemp(:) = zf(:)
      

      if (ref1.eq.2 .and. ref2.eq.1 .and. ref3.eq.3) then           !YXZ
        xf(:) = ytemp(:) ; yf(:) = xtemp(:) ; zf(:) = ztemp(:)
	 
      else if (ref1.eq.1 .and. ref2.eq.3 .and. ref3.eq.2) then	     !XZY
        xf(:) = xtemp(:) ; yf(:) = ztemp(:) ; zf(:) = ytemp(:)
	 
      else if (ref1.eq.3 .and. ref2.eq.1 .and. ref3.eq.2) then      !ZXY
        xf(:) = ztemp(:) ; yf(:) = xtemp(:) ; zf(:) = ytemp(:)

      else if (ref1.eq.2 .and. ref2.eq.3 .and. ref3.eq.1) then      !YZX
        xf(:) = ytemp(:) ; yf(:) = ztemp(:) ; zf(:) = xtemp(:)
	 
      else if (ref1.eq.3 .and. ref2.eq.2 .and. ref3.eq.1) then      !ZYX
        xf(:) = ztemp(:) ; yf(:) = ytemp(:) ; zf(:) = xtemp(:)
	 
      end if

      RETURN
      END SUBROUTINE axis_permutation
!#####################################################################
      SUBROUTINE grid_export(ncmax,loop,xp,yp)

!     Export the grid
!#####################################################################
      implicit none
      
      INTEGER,intent(in) :: ncmax,loop
      DOUBLE PRECISION,dimension(ncmax,2),intent(in) :: xp,yp
      
      INTEGER :: nc,ipts,lpt,lptmax,ip
      DOUBLE PRECISION :: area
      INTEGER,dimension(ncmax,4) :: ptconnect
      DOUBLE PRECISION,dimension(ncmax*8) :: xpts,ypts,area_norm
      CHARACTER(200) :: chb,chn,che,chsid,chsol,fn


      ! ARRANGE THE DATA FOR THE TECPLOT EXPORTATION:
      ipts=1
      do nc=1,ncmax 
	 xpts(ipts)=xp(nc,1)   ; ypts(ipts)=yp(nc,1)	; ptconnect(nc,1)=ipts
	 xpts(ipts+1)=xp(nc,1) ; ypts(ipts+1)=yp(nc,2)   ; ptconnect(nc,2)=ipts+1
	 xpts(ipts+2)=xp(nc,2) ; ypts(ipts+2)=yp(nc,2)   ; ptconnect(nc,3)=ipts+2
	 xpts(ipts+3)=xp(nc,2) ; ypts(ipts+3)=yp(nc,1)   ; ptconnect(nc,4)=ipts+3
        area = (xp(nc,2)-xp(nc,1)) * (yp(nc,2)-yp(nc,1))
        area_norm(ipts) = area / (0.01*0.01)
        area_norm(ipts+1) = area / (0.01*0.01)
        area_norm(ipts+2) = area / (0.01*0.01)
        area_norm(ipts+3) = area / (0.01*0.01)
	 ipts=ipts+4
      end do
      lptmax=ipts-1

      ! EXPORT GRID IN A TECPLOT FORMAT:
      write(chb,'(I4)') loop
      write(chn,'(I10)') lptmax
      write(che,'(I4)') ncmax
      write(chsid,'(I4)') loop+1
      write(chsol,'(I4)') loop

      fn='QuadMesh_'//trim(adjustl(chb))//'.dat'

      open(20,file=fn)
	 WRITE(20,'(a)') 'Variables=x,y,area_norm'
	 write(20,'(a)') 'ZONE T="Quad", NODES='//trim(adjustl(chn))// &
     &	 ', ELEMENTS='//trim(adjustl(che))//', STRANDID='//trim(adjustl(chsid))// &
     &	 ', SOLUTIONTIME='//trim(adjustl(chsol))//', DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'
	 do lpt=1,lptmax
	   write(20,'(3F10.3)') xpts(lpt),ypts(lpt),area_norm(lpt)
	 end do
	 do nc=1,ncmax
	   write(20,'(4I10)') (ptconnect(nc,ip),ip=1,4)
	 end do
      close(20)
	
      END SUBROUTINE
! !#####################################################################
!       subroutine IBPs_removal(ib,M,filename,filetype)

! ! Once the octree has located the cell that has IBPs,
! ! at the correct mesh refinement.
! ! This subroutine check if there is more than one IBP in each cell
! ! If IBP > 1 -> Remove the other IBPs and place the IBP in cell center
! !#####################################################################
!       implicit none
      
!       integer,intent(in) :: ib,M,filetype
!       character(200),intent(in) :: filename
!       integer :: nc,l,rpts,ipts,iptsrm,iptsmax,sn,n,c
!       double precision :: px,py,pz
!       character(200) :: fn,chb
!       logical :: centerCell
      
!       rpts=0
!       do nc=1,ncmax
        
! 	IF(npts(nc).eq.1) then         ! ONLY ONE POINT IN THE CELL
! 	  DO c=1,npts(nc)
! 	    if(centerCell) then
! 	      ! MOVE THE IBP IN THE CENTER OF THE CELL:
! 	      l=Pptid(nc,c) 
! 	      pt(l,1)=xp(nc,1)+(xp(nc,2)-xp(nc,1))/2d0
! 	      pt(l,2)=yp(nc,1)+(yp(nc,2)-yp(nc,1))/2d0
! 	    end if
! 	  END DO

! 	 ELSE IF(npts(nc).gt.1) then    ! MORE THAN ONE POINT IN A THE CELL
!           if(pt_average) then
!             do c=1,npts(nc)
! 	       l=Pptid(nc,c) 
!               avpt(l,1)=avpt(l,1)+pt(l,1)/npts(nc)
!               avpt(l,2)=avpt(l,1)+pt(l,2)/npts(nc)
!             end do
            
!             do c=1,npts(nc)
!              if(c.eq.1) then
!                pt(l,1)=avpt(l,1) ; pt(l,2)=avpt(l,2)
!              else
! 	        rpts=rpts+1
! 	        pt(l,:)=-1000d0  ! flag the number to not export this point
!              end if
            
!           else
!             do c=1,npts(nc)
! 	       if(c.eq.1) then
!                 if(cellcenter) then
! 		    pt(l,1)=Pxoct(nc,1)+(Pxoct(nc,2)-Pxoct(nc,1))/2d0
! 		    pt(l,2)=Pyoct(nc,1)+(Pyoct(nc,2)-Pyoct(nc,1))/2d0
!                 end if
! 	       else
! 	         rpts=rpts+1
! 	         pt(l,:)=-1000d0  ! flag the number to not export this point
! 	       end if
! 	  END DO !c
! 	END IF !cellpts
!       end do !nc
!       iptsmax=ibpsdom(ib)-rpts
!       iptsrm=nint(1.d0*rpts/ibpsdom(ib)*100)
! !      WRITE(6,'(2I8,2I15,2I17)') myrank,dom_id(ib),Pncells
! !     & ,iptsmax,ibpsdom(ib),iptsrm

!       ! EXPORT THE FINAL GEOMETRY POINT OF THE SUBDOMAINS:
! !      fn=filename
!       write(chb,'(I6)') dom_id(ib)
!       sn=len(trim(adjustl(chb)))
!       chb=repeat('0',(6-sn))//trim(adjustl(chb))
      
!       fn=trim(adjustl(filename))//'_'//trim(adjustl(chb))//'.pts'
!       open(20,file=fn)
! 	write(20,*) iptsmax,rpts
!        if(filetype.ne.3) then
! 	do l=ibs(ib),ibe(ib)
! 	  if(pt(l,1).ne.-1000.d0) then
! 	    write(20,'(3E27.8,2I10)') pt(l,1),pt(l,2),pt(l,3),myrank,dom_id(ib)
! 	  end if
! 	end do
!        else
! 	do l=ibs(ib),ibe(ib)
! 	  if(pt(l,1).ne.-1000.d0) then
! 	    write(20,'(6E27.8)') pt(l,1),pt(l,2),pt(l,3),pn(l,1),pn(l,2),pn(l,3)
! 	  end if
! 	end do
!        end if

!       close(20)

!       fn=trim(adjustl(filename))//'_'//trim(adjustl(chb))//'.dat'
!       WRITE(6,'(a)') trim(adjustl(fn))
!       open(20,file=fn)
!       if(filetype.ne.3) then
! 	write(20,*) 'Variables=x,y,z,proc,dom_id'
! 	do l=ibs(ib),ibe(ib)
! 	  if(pt(l,1).ne.-1000.d0) then
! 	    write(20,'(3E27.8,2I10)') pt(l,1),pt(l,2),pt(l,3),myrank,dom_id(ib)
! 	  end if
! 	end do
!       else 
! 	write(20,*) 'Variables=x,y,z,nx,ny,nz'
! 	do l=ibs(ib),ibe(ib)
! 	  if(pt(l,1).ne.-1000.d0) then
! 	    write(20,'(6E27.8)') pt(l,1),pt(l,2),pt(l,3),pn(l,1),pn(l,2),pn(l,3)
! 	  end if
! 	end do
!       end if
!       close(20)

!       END SUBROUTINE

