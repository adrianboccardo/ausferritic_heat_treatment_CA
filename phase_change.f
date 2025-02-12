      program phase_change
      
      implicit real*8 (a-h,o-z)
      integer :: row,column,page,nset
      integer :: icellgra
      character*15 :: filegraphset,yn,iplotf,fileinmicro
      integer, parameter :: n=5000,n1=100000000
      real*8, dimension(n,2) :: graphset
      
      real*8, dimension(:,:,:), allocatable :: grid,delgrid,cdiff,
     .     gridtamb             !grid for state to ambient temperature
      integer, dimension(:,:,:), allocatable :: capaux
      
      type :: capture_data
      integer, dimension(:,:,:), allocatable :: icap !capture index
      real*8, dimension(:,:,:), allocatable :: theta,alpha,ll,lt !angles and lengths
      real*8, dimension(:,:,:), allocatable :: cx,cy,cz !x, y, and z coordinate of cell centers
      end type capture_data
      type (capture_data) capture
      
      type :: chemical_comp_data
      real*8, dimension(:,:,:), allocatable :: carbon,carbon1,carbcap, !alloy elements
     .     Si,Mn,Cr,Ni,Cu,Mo,cgammao,cgammamax
      end type chemical_comp_data
      type (chemical_comp_data) chemicomp
      
      type :: nucleation_sites
      integer, dimension(:,:), allocatable :: graph,tips
      end type nucleation_sites
      type (nucleation_sites) nuclea
      
      type :: austenite_boundary
      integer :: ibound         !austenite boundary
      integer, dimension(6) :: icelflux !flux constraint in austenite bundary
      end type austenite_boundary
      
      type :: carbon_rejection
      real*8, dimension(:,:,:), allocatable :: cel2rej !number of cells that interface considers to refect carbon
      type (austenite_boundary), dimension(:,:,:), allocatable :: 
     .     iausboun
      end type carbon_rejection
      type (carbon_rejection) carbrejec
      
      type :: maximum_free_energy
      real*8, dimension(:,:,:), allocatable :: ako,ak1,ak2 !coefficient of free energy equation
      end type maximum_free_energy
      type (maximum_free_energy) freenergy

      
!     load input file
      call input_data(dtime,celld,igridd,inputd,timemax,dtimeplot,
     .     iplotf,ainis,aends,agraphs,tempgamma,temp,fgr,filegraphset,
     .     fileinmicro,csi,cmn,ccr,cni,ccu,cmo,ck,ck1,ck2,c2)
      
!     variable definition
      time=0.d0
      timeplot=0.d0
      timemart=0.d0
      istep=0
      
!     graphite set
      call graph_set(filegraphset,n,graphset,nset)
      
!     array dimensions
      call gridimension(celld,igridd,inputd,fgr,nset,n,graphset,
     .     dx,dy,dz,row,column,page,sup_inter,alen)

      
!     array allocation
      allocate(grid(row,column,page))
      allocate(gridtamb(row,column,page))
      allocate(delgrid(row,column,page))
      allocate(capaux(row,column,page))
      allocate(cdiff(row,column,page))
      allocate(capture%icap(row,column,page))
      allocate(capture%theta(row,column,page))
      allocate(capture%alpha(row,column,page))
      allocate(capture%ll(row,column,page))
      allocate(capture%lt(row,column,page))
      allocate(capture%cx(row,column,page))
      allocate(capture%cy(row,column,page))
      allocate(capture%cz(row,column,page))
      allocate(chemicomp%carbon(row,column,page))
      allocate(chemicomp%carbon1(row,column,page))
      allocate(chemicomp%carbcap(row,column,page))
      allocate(chemicomp%cgammao(row,column,page))
      allocate(chemicomp%cgammamax(row,column,page))
      allocate(chemicomp%Si(row,column,page))
      allocate(chemicomp%Mn(row,column,page))
      allocate(chemicomp%Cr(row,column,page))
      allocate(chemicomp%Ni(row,column,page))
      allocate(chemicomp%Cu(row,column,page))
      allocate(chemicomp%Mo(row,column,page))
      allocate(nuclea%graph(row*column*page,3))
      allocate(nuclea%tips(row*column*page,3))
      allocate(carbrejec%cel2rej(row,column,page))
      allocate(carbrejec%iausboun(row,column,page))
      allocate(freenergy%ako(row,column,page))
      allocate(freenergy%ak1(row,column,page))
      allocate(freenergy%ak2(row,column,page))   
      
!     define a seed for random numbers
      call seed
      
!     initialization of arrays
      capture%icap(:,:,:)=0
      capture%alpha(:,:,:)=0.d0
      capture%theta(:,:,:)=0.d0
      capture%ll(:,:,:)=0.D0
      capture%lt(:,:,:)=0.D0
      capture%cx(:,:,:)=0.D0
      capture%cy(:,:,:)=0.D0
      capture%cz(:,:,:)=0.D0
      carbrejec%cel2rej(:,:,:)=0.D0
      carbrejec%iausboun(:,:,:)%ibound=0
      grid(:,:,:)=ainis         !set initial state, then updated
      delgrid(:,:,:)=0.d0
      capaux(:,:,:)=0
      cdiff(:,:,:)=0.d0
      
      if(fileinmicro.eq.'y') then !from file
         write(*,*) 'input file from solidification model is used'
         open (unit=21,file='micro_results.out',status='old',err=5)
         
         do i=1,n1,1            !state and chemical composition
            read(21,*,end=6) irowaux,icolumnaux,ipageaux,cstateaux, !x y z state Si Mn Cr Ni Cu Mo
     .           csiaux,cmnaux,ccraux,cniaux,ccuaux,cmoaux
            grid(irowaux,icolumnaux,ipageaux)=cstateaux
            chemicomp%Si(irowaux,icolumnaux,ipageaux)=csiaux
            chemicomp%Mn(irowaux,icolumnaux,ipageaux)=cmnaux
            chemicomp%Cr(irowaux,icolumnaux,ipageaux)=ccraux
            chemicomp%Ni(irowaux,icolumnaux,ipageaux)=cniaux
            chemicomp%Cu(irowaux,icolumnaux,ipageaux)=ccuaux
            chemicomp%Mo(irowaux,icolumnaux,ipageaux)=cmoaux
         enddo
 6       continue
         do i=1,row,1           !carbon content
            do j=1,column,1
               do k=1,page,1
                  if(grid(i,j,k).ne.agraphs) then
c                     csi=chemicomp%Si(i,j,k)
c                     cmn=chemicomp%Mn(i,j,k)
c                     ccr=chemicomp%Cr(i,j,k)
c                     cni=chemicomp%Ni(i,j,k)
c                     ccu=chemicomp%Cu(i,j,k)
c                     cmo=chemicomp%Mo(i,j,k)
                     call carboncont(cgammamax,cgammao,tempgamma,temp,
     .                    csi,cmn,ccr,cni,ccu,cmo)
                     
                     chemicomp%cgammao(i,j,k)=cgammao
                     chemicomp%cgammamax(i,j,k)=cgammamax
                     chemicomp%carbon(i,j,k)=cgammao
                     chemicomp%carbon1(i,j,k)=cgammao
                     chemicomp%carbcap(i,j,k)=cgammao
                  endif
               enddo
            enddo
         enddo
      else                      !proposed values (uniform values)
         call microstructure(row,column,page,nset,n,dx,agraphs,graphset, !initial microstructure
     .        grid)
         gridtamb=grid          !set as initial microstructure, then updated
         
         call carboncont(cgammamax,cgammao,tempgamma,temp,csi,cmn,ccr, !austenite carbon concentration
     .     cni,ccu,cmo)
         chemicomp%Si(:,:,:)=csi
         chemicomp%Mn(:,:,:)=cmn
         chemicomp%Cr(:,:,:)=ccr
         chemicomp%Ni(:,:,:)=cni
         chemicomp%Cu(:,:,:)=ccu
         chemicomp%Mo(:,:,:)=cmo
         chemicomp%cgammao(:,:,:)=cgammao
         chemicomp%cgammamax(:,:,:)=cgammamax
         chemicomp%carbon(:,:,:)=cgammao
         chemicomp%carbon1(:,:,:)=cgammao
         chemicomp%carbcap(:,:,:)=cgammao
         
         write(*,*) cgammao,cgammamax
         
      endif
      
!     initial condition
      anc=0.d0                  !number of cell to nucleate
      
!     graphite interface cell of austenite
      call graph_inter(row,column,page,grid,nuclea,ainis,agraphs,
     .     icellgra,carbrejec)
      
!     constant incubation model ausferritic transformation
      do i=1,row,1
         do j=1,column,1
            do k=1,page,1
               if(grid(i,j,k).ne.agraphs) then
                  csi=chemicomp%Si(i,j,k)
                  cmn=chemicomp%Mn(i,j,k)
                  ccr=chemicomp%Cr(i,j,k)
                  cni=chemicomp%Ni(i,j,k)
                  ccu=chemicomp%Cu(i,j,k)
                  cmo=chemicomp%Mo(i,j,k)
                  call gmaxbha(temp,csi,cmn,ccr,cni,ccu,cmo,ako,ak1,ak2)
                  freenergy%ako(i,j,k)=ako
                  freenergy%ak1(i,j,k)=ak1
                  freenergy%ak2(i,j,k)=ak2
               endif
            enddo
         enddo
      enddo
      
!     plot at the beginning (output file at time=0)
      call pvplot(istep,row,column,page,grid,gridtamb,delgrid,cdiff,
     .     capture,chemicomp,carbrejec,time,freenergy,temp)
      
c!     run?
c      write(*,*) 'continue? (y/n)'
c      read(*,*) yn
c      if(yn.eq.'n') goto 20
      
!******************************
!     solve model
      write(*,*) 'start simulation'
      do while (time.lt.timemax)
         time=time+dtime
         timeplot=timeplot+dtime
         timemart=timemart+dtime
         istep=istep+1
         
!        nucleation at graphite interface
         call nucleation(row,column,page,dx,dy,dz,sup_inter,grid,
     .        capture,agraphs,nuclea,icellgra,chemicomp,dtime,anc,time,
     .        temp,ck,ck1,ck2,freenergy)
         
!        sheaves growth
         call growth(row,column,page,grid,delgrid,capture,dtime,dx,
     .        dy,dz,ainis,aends,chemicomp,time,temp,ck,ck1,ck2,
     .        freenergy)
         
!        fraction
         call fract(row,column,page,grid,capture,chemicomp,time,
     .        dx,temp,ainis,aends,agraphs,fracgraph,fracsh,carbrejec)
         
!        interface cells (possible subroutine position, it has a small
!        difference with the original version close to the end of the
!        transformation)
c         call interface(row,column,page,grid,capture,capaux,
c     .        chemicomp,dx,dy,dz,ainis,agraphs,freenergy,temp)
         
!        carbon diffusion
         call carbon_diff(row,column,page,dx,dtime,temp,ainis,aends,
     .        agraphs,fracgraph,fracsh,grid,delgrid,cdiff,
     .        carbrejec,chemicomp)
         
!        interface cells
         call interface(row,column,page,grid,capture,capaux,
     .        chemicomp,dx,dy,dz,ainis,agraphs,freenergy,temp)
         
!        output
         if (timeplot.ge.dtimeplot) then !output file at time t
            timeplot=0.d0
            call pvplot(istep,row,column,page,grid,gridtamb,delgrid,
     .           cdiff,capture,chemicomp,carbrejec,time,freenergy,
     .           temp)
         endif
         
!        martensite growth
c         dtimemart=5.d0
c         if (timemart.ge.dtimemart) then !output file at time t
c            timemart=0.d0
         call mart_growth(row,column,page,dx,temp,c2,grid,gridtamb,
     .        capture,chemicomp,time,ainis,aends,agraphs,fracgraph,
     .        fracsh)
c         endif
         
      enddo
!     model solved
!******************************
      
 20   continue
!     plot at the finish
      if(iplotf.eq.'y') then    !output file at the end
         call pvplot(istep,row,column,page,grid,gridtamb,delgrid,
     .        cdiff,capture,chemicomp,carbrejec,time,freenergy,
     .        temp)
      endif
      
!     plot results (state and alloy content)
      call results_out(row,column,page,alen,ainis,aends,grid,chemicomp)
      
!     end
      write(*,*) 'end of simulation'
      call exit
      
!     file didn't open
 5    continue
      write(*,*) 'error: input file has not been loaded'
      write(*,*) 'looking for: micro_results.out'
      
      end program
