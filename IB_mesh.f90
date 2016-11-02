subroutine  leggo_heart (nz,ny,nx,z,y,x)

integer nz,ny,nx
real, allocatable :: z(:), y(:), x(:)
integer fh, dummy, mi, Xi, mj, Xj, mk, Xk, Ns, err_alloc
character*256  buff
! double precision zD, yD, xD, dummy_d
real*8 zD, yD, xD, dummy_d

  fh = 20
  
  open (Unit = fh, file='LAW_032.rst', STATUS='OLD')
  do while (.TRUE.)
  read(Unit = fh, FMT='(a256)',END=11) buff
    select case  (buff(1:1))
      case ('Z')
        if (buff(2:3).eq.'  ') then
          read(buff(4:),*) inZONE
        elseif (buff(2:2).eq.'H') then
          continue
        else
          read(buff(2:3),*) B
          if(buff(4:5)=='I0'.or.buff(4:5)=='IN'.or.buff(4:5)=='J0'.or.buff(4:5)=='JN'.or.buff(4:5)=='K0'.or.buff(4:5)=='KN')then
            continue
          else
            if (B==0) then
             read(buff(4:),*)  dummy, mi, Xi, mj, Xj, mk, Xk
            endif
          end if
        end if
      case ('M')    ! Species: number, index of inert, names
        read(buff(4:),*) Ns
      case DEFAULT
    end select
  END DO
  11    close (fh)
! Ho letto dal .rst il numero di punti della mesh lungo le tre direzioni.  
  nz = Xi-mi+3;   ny = Xj-mj+3;   nx = Xk-mk+3

  allocate( z(nz), y(ny), x(nx), STAT=err_alloc)

  open (Unit = fh,file='LAW_032.ZONE0', FORM="UNFORMATTED", RECORDTYPE='STREAM', STATUS='OLD', CONVERT='BIG_ENDIAN')
  print*, "Reading grid from previous file"
  do k= 1, nx
     do j= 1, ny
        do i= 1, nz
          read (fh) zD, yD, xD ,dummy_d,dummy_d,dummy_d,dummy_d,dummy_d,dummy_d,(dummy_d,s=1,Ns),dummy_d ,dummy_d
!          z(i) = zD;   
          y(j) = yD;
!          x(k) = xD
        end do
     end do
  end do
  close(fh)
end subroutine leggo_heart

subroutine  leggo_Open_Smoke(nz,ny,nx,NsMax,zzz,yyy,xxx,Uread,Tread,Yis)
integer nz,ny,nx,n
real,allocatable :: zzz(:), yyy(:), xxx(:),Tread(:),Uread(:),Yis(:,:)
integer fh, mi, Xi, mj, Xj, mk, Xk, Ns, err_alloc
character*256  buff
! double precision zD, yD, xD, dummy_d
real*8 zD, yD, xD, dummy_d
CHARACTER*7 YiNAME(1:NsMAX)
CHARACTER*2000 prima_riga
REAL Ws(1:NsMAX), dummy
REAL*16 prova1,prova2




CHARACTER*7 YiFIND
INTEGER NsMAX, s, fh_chem
  fh = 20
  allocate( zzz(nz), yyy(ny), xxx(nx),Tread(nz),Uread(nz),Yis(NsMax,nz), STAT=err_alloc)
  open (Unit = fh, file='OPENSMOKE.dat', STATUS='OLD')
  read(fh,'(A2000)') prima_riga

  do i=1,nz
    read(fh,*) zzz(i),Uread(i),dummy,dummy,dummy,Tread(i), (dummy,n=1,(NsMax+15)),(Yis(n,i),n=1,NsMax)
  enddo


return
end subroutine  leggo_Open_Smoke



SUBROUTINE Read_chem(NsMAX,YiNAME,Ws)
IMPLICIT NONE
CHARACTER*32 YiNAME(1:NsMAX)
REAL Ws(1:NsMAX)

CHARACTER*32 YiFIND
INTEGER NsMAX, s, fh_chem
fh_chem=15
open (fh_chem,file='chem.dat',STATUS='OLD')
YiFIND = ""
do s=1,NsMAX
  read (fh_chem,53) YiFIND
  do 33 while (YiFIND.NE.YiNAME(s))
    read (fh_chem,53,END=1000) YiFIND
  33    continue
!     print*, s,YiFIND," ",YiNAME(s)
  read (fh_chem,*) Ws(s)                   !Molecular Weights [kg/mol]
  REWIND(fh_chem)
enddo
close(fh_chem)
RETURN
1000  WRITE (*,'(a,X,a,X, a)') "'Error on Input: '", Yiname(s), s,&
          " 'Invalid Name, Control input for species'"
close(fh_chem)
STOP
53 format(A)
END SUBROUTINE Read_chem

subroutine IB_create
  use file
  use IBmod

  integer tot_num_point,iunit,Read_zone
  integer IVM,Wfield,Compact,NZONE,izone
  real, allocatable :: zzz(:), yyy(:), xxx(:)
  character*20 file_name
  character*15 pre_file_name

  real :: B1,b2,up,Um,dl,l,W,SL

  real :: Strain_t(3,3), eigenvector(3,3), massimo, eigenvalues(3), modulo_rot, a
  integer :: cnt(3), rc, ind_massimo

  interface
     subroutine leggo_heart (nz,ny,nx,z,y,x)
          integer nz,ny,nx
          real, allocatable :: z(:), y(:), x(:)
     end subroutine
     subroutine leggo_Open_Smoke (nz,ny,nx,Nspecies,z,y,x,Uread, Tread, Yis)
          integer nz,ny,nx,Nspecies
          real, allocatable :: z(:), y(:), x(:), Uread(:), Tread(:), Yis(:,:)
     end subroutine
  end interface

!   Strain_t(:,:) = 0.
!   Strain_t(1,2) = 1.
!   Strain_t(2,1) = 1. 
! 
! 
!   call Jacobi(Strain_t, 3, eigenvalues, eigenvector, cnt)
!   print*, U_l(50.d0,0.70d0,0.d0,600.74d0,1.0184d0)
!   stop
!   print*, U_l(50.,0.7,0.,588.5,10.)
!   print*, U_l(50.,0.7,0.,588.5,40.)
!   stop!alfa,phi,Yu_res,T,P

  call readInt("Read_zone", Read_zone)
  if (Read_zone==1) then
  
!                    z   r  t  zz  rr  tt
  call leggo_heart (nx,ny,nz,zzz,yyy,xxx)
!  x = zzz;
  y = yyy; 
!  z = xxx;

  elseif (Read_zone==2) then
   call readInt("Nspecies", Nspecies)
   allocate(Uread(nx),Tread(nx),Yis(Nspecies,nx))
   call leggo_Open_Smoke(nx,ny,nz,Nspecies,zzz,yyy,xxx,Uread,Tread,Yis)
   x = zzz*0.01;
!  y = yyy; 
!  z = xxx;

  endif

  call readInt("NZONE", NZONE)
  print*,"NZONE", NZONE
 
  call readInt("Compact", Compact)
  call readInt("IVM", IVM)
  
  if (Compact==1) then
    iunit = iopen()
    ind_= index(Outfile,'_')
    file_name= Outfile(1:ind_-1)//"_COMPACT"
    open(iunit,file=file_name, form='unformatted',status='unknown',RECORDTYPE='STREAM',CONVERT='BIG_ENDIAN') 
    call compact_define(nz,ny,nx,z,y,x,1,iunit) 
    print*,"Writing compact x direction..."
    call compact_define(nz,ny,nx,z,y,x,2,iunit)
    print*,"Writing compact y direction..."
    call compact_define(nz,ny,nx,z,y,x,3,iunit)
    print*,"Writing compact z direction..."
    close(iunit)
   
        do izone = 1,NZONE 
         call write_multizone(nx,ny,nz,x,y,z,izone)
        enddo
!     call write_field(nx,ny,nz,x,y,z) 
  else
    if(IVM==1) then
     stl = STL_file
     call readgeo3ddim(nb,stl)
      
     allocate(xyzb(nb,9), barl(nb,3), nxyz(nb,3))
     
     call readgeo3d(nb,xyzb,barl,nb,stl)
     !deallocate(bar)
     
     call readnormals(nxyz,nb,stl)
     
     
     ! Create the vectors for the IB implementation
     tot_num_point = nz*ny*nx
      
     print*,'Total number of grid points',tot_num_point
     print*,'Number of stl triangles  ',nb
    endif     
    
    call readInt("Wfield", Wfield)
    print*,'Wfield',Wfield
    print*,'IVM',IVM
    
    if (IVM==1) then
     call ibgrid3d(nz,ny,nx,z,y,x,zm,ym,xm)
        do izone = 1,NZONE 
         call write_multizone(nx,ny,nz,x,y,z,izone)
        enddo
    elseif (IVM==0.and.Wfield==1) then

!       if(NZONE>1) then
        do izone = 1,NZONE 
         call write_multizone(nx,ny,nz,x,y,z,izone)
        enddo
!       else
!         call write_field(nx,ny,nz,x,y,z)
!       endif
     
    else
      continue
    endif
    
  endif
end subroutine






subroutine write_multizone(ni,nj,nk,zz1,rr1,tt1,izone)
  use file
  use IBmod
IMPLICIT NONE
      integer ni,nj,nk,ghost,NSMax,Ng,Nvar_S,i,j,k,s,n,ind, rec1,ip,jp,Npoints,KLEIN, ind_, izone, IUNIT, IERR, IOPEN, Read_zone
      parameter (ghost=1,NSMax= 130,Ng=1,Nvar_S=7)
      REAL  Ru,PI,RMAX(0:2),RMIN(0:2),Wmix1,in_Wmix2,p_in,T_in2,rho_in2,rhou2,A_Coax,inv_Wmix,Uinit,Wmix,dummy
      integer recl, reclen,reclen_S,fh,fh1,range_i,range_j, displ,stato,fh_BcIN,B, Nspecies,TVCZones,npb
      character*3 var
      character*20 file_name,file_name_S
      REAL  zz(1-ghost:ni-1),rr(1-ghost:nj-ghost), tt(1-ghost:nk-ghost)
      REAL  zz1(1:ni),rr1(1:nj), tt1(1:nk) !con i due ghost
      character*16 pre_file_name
      CHARACTER*1024 RIGA
      REAL , allocatable :: Uz_in(:,:), Ur_in(:,:),Ut_in(:,:), T_in(:,:), P_ing(:,:), Yi_in(:,:,:) , var_dum(:), Yc(:), Yj(:)
      REAL  buf(11+NSMax),rho,sump,Uj,Uc,cent,Thick,Thp,Tc,Tj,pc,pj, limitz(24)
      REAL , allocatable :: iBoundZ(:),xdir(:)
      INTEGER ix,fx,iy,fy,iz,fz,lowerInd,upperInd,findInd,NZONE
      CHARACTER*2048 string
      CHARACTER*32 , allocatable :: YiNAME(:)
      CHARACTER*7 Zlimit_,ZX

      REAL Ws(1:NsMAX),profile,senprofile
      CHARACTER*100 app_st
      

      INTEGER IMCWRK(20), ICKWRK(20),ind_NameSp, ind_N2
      REAL RMCWRK(20), CKWRK(20) 
      REAL DT(NSMax), COND(NSMax)
      INTEGER error,Dir,ndir
      LOGICAL trovato

!       REAL  ae(3,17),PE(3),GAM(3),BETA_UN,BETA_BU,ZeL(3),COEFF
!       INTEGER is,ie
! 
!       PE(1)= 12.011/1000.; PE(2)=1.00797/1000. ; PE(3)=15.994/1000. ;
!       !C
!       ae(1,1) =1.;ae(1,2) =1.;ae(1,3) =1.;ae(1,4) =1.;ae(1,5) =1.;ae(1,6) =1.;ae(1,7) =1.;ae(1,8) =1.;ae(1,9) =0.;ae(1,10) =0.;ae(1,11) =0.;ae(1,12) =0.;ae(1,13) =0.;ae(1,14) =0.;ae(1,15) =1.;ae(1,16) =0. ;ae(1,17) = 0.;  
!       !H
!       ae(2,1) =4.;ae(2,2) =3.;ae(2,3) =2.;ae(2,4) =2.;ae(2,5) =2.;ae(2,6) =1.;ae(2,7) =0.;ae(2,8) =0.;ae(2,9) =2.;ae(2,10) =1.;ae(2,11) =0.;ae(2,12) =0.;ae(2,13) =1.;ae(2,14) =1.;ae(2,15) =3.;ae(2,16) =2. ;ae(2,17) = 0.;
!       !O
!       ae(3,1) =0.;ae(3,2) =0.;ae(3,3) =0.;ae(3,4) =0.;ae(3,5) =1.;ae(3,6) =1.;ae(3,7) =2.;ae(3,8) =1.;ae(3,9) =0.;ae(3,10) =0.;ae(3,11) =2.;ae(3,12) =1.;ae(3,13) =1.;ae(3,14) =2.;ae(3,15) =1.;ae(3,16) =1. ;ae(3,17) = 0.;
!       !BILGER MIXTURE FRACTION
!       GAM(1) = 2./PE(1)
!       GAM(2) = 1./(2.*PE(2))
!       GAM(3) = -1./PE(3)


      call readInt("Read_zone", Read_zone)

     
      call readInt("Nspecies", Nspecies)
      call readInt("Dir", Dir)

      allocate (Yc(Nspecies),Yj(Nspecies),YiNAME(Nspecies))
      iunit = iopen()
      open (iunit, file=filename, iostat=ierr,status='unknown')
      trovato=.false.; error=0;
      do while (.not.trovato.or.error/=0)
       read(iunit,'(a)',IOSTAT=error) RIGA
       ind_NameSp=INDEX(RIGA,'NameSpecies')
       if(ind_NameSp>0) then
              read(iunit,'(a)') RIGA
              read(RIGA(2:),*) (YiNAME(i),i=1,Nspecies)
              trovato=.true.
       end if
      end do
      if(error/=0.and..not.trovato) then
       print *, 'IL RECORD NON ESISTE o ERRORE IN READ:',filename(1:len_trim(filename))
       STOP 
       end if 
      close(iunit)

      do s=1,Nspecies
      if("N2"==YiNAME(s)) then
      ind_N2 = s
      print*, "INERT INDEX: ",ind_N2
      endif
      enddo
     
      call archive_arraySize("pYcof", npb)
      if (Nspecies/=npb) then
          print*, "Nspecies cof /= npb",Nspecies ,npb
          stop
      endif
      call readArray("pYcof", Yc, npb)
      call archive_arraySize("pYjet", npb)
      if (Nspecies/=npb) then
          print*, "Nspecies jet /= npb",Nspecies ,npb
          stop
      endif
      call readArray("pYjet", Yj, npb)
      call readFloat("p_in", p_in)
      print*,"Pressione: ", p_in
     
      if(NsMAX<Nspecies)then
       print*, "NsMAX/=Nspecies",NsMAX,Nspecies
       stop
      endif
       
       Yi_in=0.
       print*, "Nspecies = ",Nspecies
       if (izone==1) then
        do s=1,Nspecies
         print*,  YiNAME(s)
        enddo
       endif
      Ru = 8.3144621d0 !8.314517000000000
      PI = dacos(-1.0)
      RMAX = 0.
      RMIN = 0.
      RMAX(2) = 0.0175
      RMIN(2) = 0.005
      Wmix1   = 0.028013400455
      in_Wmix2 = 1. /Wmix1
      T_in2 = 298.150
      rho_in2 = p_in / (Ru*in_Wmix2* T_in2)
      rhou2 = 0.42d0/60.d0
      A_Coax  = (RMAX(2)** 2.d0 - RMIN(2)** 2.d0) * PI
      rhou2   = rhou2 !/ rho_in2!A_Coax
      inv_Wmix = in_Wmix2 ![mol / kg]
      reclen=int((11+Nspecies)*8/4)

      reclen_S=int((Nvar_S+Ng)*8/4) !8/4 unita di misura in word (4 byte=1 word);
      pre_file_name = Outfile

      fh=10



      call readFloat("Uj", Uj)
      call readFloat("Uc", Uc)
      call readFloat("cent", cent)
      call readFloat("Thick", Thick)
      call readFloat("Thp", Thp)
      call readFloat("Tj", Tj)
      call readFloat("Tc", Tc)
      call readFloat("pj", pj)
      call readFloat("pc", pc)
      print*, pj,pc 

!      call readInt("NZONE", NZONE)

      call archive_arraySize("iBoundZ", NZONE)
      allocate(iBoundZ(NZONE))
      call readArray("iBoundZ", iBoundZ, NZONE)


      if(izone==1) then
        call readArray("Zlimit0", limitz, 6)
      elseif(izone==2) then
        call readArray("Zlimit1", limitz, 6)
      elseif(izone==3) then
        call readArray("Zlimit2", limitz, 6)
      elseif(izone==4) then
        call readArray("Zlimit3", limitz, 6)
      elseif(izone==5) then
        call readArray("Zlimit4", limitz, 6)
      elseif(izone==6) then
        call readArray("Zlimit5", limitz, 6)
      elseif(izone==7) then
        call readArray("Zlimit6", limitz, 6)
      elseif(izone==8) then
        call readArray("Zlimit7", limitz, 6)
      elseif(izone==9) then
        call readArray("Zlimit8", limitz, 6)
      elseif(izone==10) then
        call readArray("Zlimit9", limitz, 6)
      elseif(izone==11) then
        call readArray("Zlimit10", limitz, 6)
      elseif(izone==12) then
        call readArray("Zlimit11", limitz, 6)
      elseif(izone==13) then
        call readArray("Zlimit12", limitz, 6)
      elseif(izone==14) then
        call readArray("Zlimit13", limitz, 6)
      elseif(izone==15) then
        call readArray("Zlimit14", limitz, 6)
      elseif(izone==16) then
        call readArray("Zlimit15", limitz, 6)
      elseif(izone==17) then
        call readArray("Zlimit16", limitz, 6)
      elseif(izone==18) then
        call readArray("Zlimit17", limitz, 6)
      elseif(izone==19) then
        call readArray("Zlimit18", limitz, 6)
      elseif(izone==20) then
        call readArray("Zlimit19", limitz, 6)
      elseif(izone==21) then
        call readArray("Zlimit20", limitz, 6)
      elseif(izone==22) then
        call readArray("Zlimit21", limitz, 6)
      endif
    
      file_name = pre_file_name(:len_trim(pre_file_name))//".ZONE"//"X"
      B = izone-1
      write(file_name(index(file_name,"ZONE")+4:index(file_name,"ZONE")+5),"(i2.1)") izone-1
      print*, "Multizone: ",file_name
      if (Dir==1) then 
       allocate(xdir(0:nk-1))
      elseif(Dir==2) then 
       allocate(xdir(0:nj-1))
      elseif(Dir==3) then
       allocate(xdir(0:ni-1))
      endif

      do i=0,ni-1
        zz(i) = dble(zz1(i+1))
        if (Dir==3)  xdir(i) = zz(i)
      end do
      do j=0,nj-1
        rr(j) = dble(rr1(j+1))
        if (Dir==2)  xdir(j) = rr(j)
      end do
      do k=0,nk-1
        tt(k) = dble(tt1(k+1))
        if (Dir==1) xdir(k) = tt(k) 
      end do
      range_i= ni
      range_j= nj

      ix = findInd(limitz(1), tt(0:nk-1), nk) ! = Mink-1
      fx = findInd(limitz(2), tt(0:nk-1), nk) ! = Maxk
      iy = findInd(limitz(3), rr(0:nj-1), nj) ! = Minj-1
      fy = findInd(limitz(4), rr(0:nj-1), nj) ! = Maxj
      iz = findInd(limitz(5), zz(0:ni-1), ni) ! = Mini-1
      fz = findInd(limitz(6), zz(0:ni-1), ni) ! = Maxi
      !QUI FACCIO LA STAMPA DEI PUNTI INTERNI DA METTERE NEL FILE *.rst

      print*, "ZONE= ", B
      print*, "X direction ","ix= ",ix, "fx= ",fx+1
      write(*,21) "x(MINk-1)= ",tt(ix), "x(MAXk)= ",tt(fx)
! 
      print*, "Y direction ","iy= ",iy, "fy= ",fy+1
      write(*,21) "y(MINj-1)= ",rr(iy), "y(MAXj)= ",rr(fy)
! 
      print*, "Z direction ","iz= ",iz, "fz= ",fz+1
      write(*,21) "z(MINi-1)= ",zz(iz), "z(MAXi)= ",zz(fz)
      ZX="ZX"
      write(ZX(2:3),"(i2.2)") izone-1
      write(*,22) ZX, 1, iz+1,fz,iy+1,fy,ix+1,fx

21      FORMAT (1x,a,1x,e15.6,1x,a,1x,e15.6)
22      FORMAT (1x,a,1x,i1.1,6(1x,i4))
      open(UNIT=fh, file=file_name,FORM='UNFORMATTED',RECL=reclen, &
           STATUS='UNKNOWN',ACTION='WRITE',ACCESS='DIRECT' &
          ,CONVERT='BIG_ENDIAN',err=101,iostat=stato)

        ALLOCATE (Uz_in(0:nj-1,0:nk-1))
        ALLOCATE (Ur_in(0:nj-1,0:nk-1))
        ALLOCATE (Ut_in(0:nj-1,0:nk-1))
        ALLOCATE (T_in(0:nj-1,0:nk-1))
        ALLOCATE (P_ing(0:nj-1,0:nk-1))
        ALLOCATE (Yi_in(1:Nspecies,0:nj-1,0:nk-1)) 


19      FORMAT (<(Nspecies)>(1x,10A))

       call Read_chem(Nspecies,YiNAME,Ws)


! !BURNT
!       BETA_BU= 0.d0;
!       do ie=1,3
!        ZeL(ie)= 0.d0
!        do is=1,Nspecies
!         ZeL(ie) = ZeL(ie) + ae(ie,is)*PE(ie)*Yc(is)/Ws(is)
!        enddo
!        BETA_BU = BETA_BU +ZeL(ie)*GAM(ie) 
!       enddo 
!        print*, "BETA_BU",i,j,k,BETA_BU,ZeL(:)
! 
! !UNBURNT
!       BETA_UN= 0.d0;
!       do ie=1,3
!        ZeL(ie)= 0.d0
!        do is=1,Nspecies
!         ZeL(ie) = ZeL(ie) + ae(ie,is)*PE(ie)*Yj(is)/Ws(is)
!        enddo
!        BETA_BU = BETA_BU +ZeL(ie)*GAM(ie) 
!       enddo 
!        print*, "BETA_UN",i,j,k,BETA_BU,ZeL(:)
!        stop



        displ = 1
        do k= ix,fx+1 !0,nk-1 !X
          do j= iy,fy+1 !0,nj-1 !Y
            do i= iz,fz+1 !0,ni-1 !Z
              if (Dir==1) ip = k
              if (Dir==2) ip = j
              if (Dir==3) ip = i  
              sump = 0.q0
!             if(i>=100.and.i<=198) then
               Uz_in(j,k) = profile(Uj,Uc,xdir(ip),cent,Thick/2.,Thp*Thick) !& 
               !+ 35.0*senprofile(zz,i,1,ni)! + 0.01*((-1.)**(j))*(Uj+Uc)*0.5!senprofile(rr,j,0,nj-1)!*profile(1.,0.,zz(k),0.,-0.0116788820238/2.,0.001*2.*0.0116788820238)
!             else
!               Uz_in(j,k) = profile(Uj,Uc,xdir(ip),cent,Thick/2.,Thp*Thick) 
!             endif
!              if(zz(i)>-0.00)then
!               Uz_in(j,k) = 0.d0!Uc
!              endif
!               if(i<60.and.i>40) Uz_in(j,k)=profile(55.,Uc,xdir(ip),cent,Thick/2.,Thp*Thick) 
              if (Read_zone==2) then
                Uz_in(j,k) = Uread(i+1)*0.01
              endif
              if(Nspecies>1)then
                Yi_in(:,j,k) = 0.
                do n=1,Nspecies
                  Yi_in(n,j,k) = profile(Yj(n), Yc(n),xdir(ip),cent,Thick/2.,Thp*Thick)

!              if(zz(i)>-0.0)then
!                Yi_in(n,j,k) = Yc(n)
!              endif                   

                  if (Read_zone==2) then
                   Yi_in(n,j,k) = Yis(n,i+1)
                  endif 
                  sump = sump + Yi_in(n,j,k)
                enddo 
                if(sump.ne.1.q0) then
                 do s=1, Nspecies
!                   print*, sump,k,j,i
                 Yi_in(s,j,k) = (Yi_in(s,j,k) / sump) * (1.D0 - sump) + Yi_in(s,j,k)
                 end do
                end if 
              else
                Yi_in(1,j,k) = 1.
              endif
              T_in(j,k) =  profile(Tj,  Tc,xdir(ip),cent,Thick/2.,Thp*Thick)!senprofile(xdir,i,0,ni)

!              if(zz(i)>-0.0)then
!                 T_in(j,k) = 298.15q0
!                 Yi_in(1,j,k) = 1.q0
!                 Yi_in(2,j,k) = 0.q0
!              endif
 
              P_ing(j,k) = profile(pj,  pc,xdir(ip),cent,Thick/2.,Thp*Thick) !+ 100000.0*senprofile(zz,i,0,ni) 
              if (Read_zone==2) then
                T_in(j,k) = Tread(i+1)
              endif
              Ur_in(j,k) =  0.
              Ut_in(j,k) =  0.
              Wmix = 0.0
              do s=1, Nspecies
                Wmix = Wmix + Yi_in(s,j,k) / Ws(s)
!               print*,Yi_in(s,j,k), Ws(s), YiNAME(s)
              end do
              rho = P_ing(j,k)/(Ru * Wmix * T_in(j,k)) 
!               rho = rho + 0.01*((-1.)**(j))*rho

!               if (j<1)then
!                 print*,1./Wmix,rho,T_in(j,k),Yi_in(1:3,j,k)
!               endif
              WRITE(fh,rec=displ,err=100,iostat=stato) dbleq(zz(i)),dbleq(rr(j)),dbleq(tt(k)),dbleq(Uz_in(j,k))*dbleq(rho),dbleq(Ur_in(j,k))*dbleq(rho),dbleq(Ut_in(j,k))*dbleq(rho), &
                                                       dbleq(P_ing(j,k)),dbleq(T_in(j,k)),dbleq(rho),(dbleq(Yi_in(n,j,k)),n=1,Nspecies),dbleq(Wmix * T_in(j,k)),dbleq(1.)!buf
              displ = displ+1
            end do
          end do
        end do
        close(fh)
           
         fh_BcIN = 78787
         dummy = 0.d0
         IF( iBoundZ(izone) == 1.) THEN   
         ind_= index(pre_file_name,'_')
         file_name= pre_file_name(1:ind_-1)//".BCin_"//CHAR(B+48)
         print*,"File BC name: ",file_name
         recl = (7+Nspecies)*16 + 1
         print*, ""

! !    print '(I3,a,I3,a,I3,a,4(I5,x))',my_id, ' READ_INLET in ZONA ',B,' in DIR=', ind_dir, file_name,  mn_1,mx_1,mn_2,mx_2
!      open (UNIT=fh_BcIn ,file=file_name, STATUS='OLD', ACTION='READ', FORM='FORMATTED', RECL=rec_lenght, ACCESS='DIRECT')

         open (UNIT=fh_BcIN ,file=file_name, STATUS='UNKNOWN',ERR=100, FORM='FORMATTED',RECL=recl, ACCESS='DIRECT',iostat=stato)
         string(1:2048)=" "
!          write(string,*)
         write(string,13) p_in, fx-ix
         string(recl:recl)=CHAR(10)
         write(fh_BcIN,fmt='(a)',rec=1) string(1:recl)
         ind=2
         i = 0
         print*, "ix,fx+1",ix,fx+1
         print*, "iy,fy+1",iy,fy+1
         do_k: do k = ix,fx+1 !0,nk-1
         do_j: do j = iy,fy+1 !0,nj-1
              if (Dir==1) ip = k
              if (Dir==2) ip = j
              if (Dir==3) ip = i
              sump = 0.q0
              Uz_in(j,k) = profile(Uj,Uc,xdir(ip),cent,Thick/2.,Thp*Thick)
              if(Nspecies>1)then
                Yi_in(:,j,k) = 0.
                do n=1,Nspecies
                 Yi_in(n,j,k) =  profile(Yj(n), Yc(n),xdir(ip),cent,Thick/2.,Thp*Thick) 
                 sump = sump + Yi_in(n,j,k)
                enddo 
                if(sump.ne.1.q0) then
                 do s=1, Nspecies
                  Yi_in(s,j,k) = (Yi_in(s,j,k) / sump) * (1.q0 - sump) + Yi_in(s,j,k)
                 end do
                end if
              else
                Yi_in(1,j,k) = 1.
              endif
              T_in(j,k) =  profile(Tj,  Tc,xdir(ip),cent,Thick/2.,Thp*Thick)    
              Ur_in(j,k) =  0.0
              Ut_in(j,k) =  0.0
              string(1:2048)=""
              app_st=""
              write(app_st,'(a,i3,a)') "(e16.9,",6+Nspecies,"(1x,e15.6))"
              write(string,FMT=app_st) dbleq(dummy), dbleq(rr(j)), dbleq(tt(k)) , dbleq(Uz_in(j,k)), dbleq(Ur_in(j,k)), dbleq(Ut_in(j,k)), dbleq(T_in(j,k)), (dbleq(Yi_in(s,j,k)),s=1,Nspecies)
              string(recl:recl)=CHAR(10)
              write(unit=fh_BcIn, fmt='(a)',rec=ind) string(1:recl)
              ind = ind+1
             end do do_j
             end do do_k
        close(fh_BcIN)
        
        fh_BcIN = 78787
        dummy = 0.d0
        ind_= index(pre_file_name,'_')
        file_name= pre_file_name(1:ind_-1)//".BCin_"//CHAR(B+48)//"_OR"
        print*,"File BC_ORname: ",file_name
        rec1 = (7+Nspecies)*16 + 1
        open (UNIT=fh_BcIN ,file=file_name, STATUS='UNKNOWN', FORM='FORMATTED',RECL=rec1,err=100,iostat=stato, ACCESS='DIRECT')
         string(1:2048)=" "
!          write(string,*)
         write(string,13) p_in, fx-ix
         string(recl:recl)=CHAR(10)
         write(fh_BcIN,fmt='(a)',rec=1) string(1:recl)
         ind=2
         i = 0
        do_ki: do k = ix,fx+1 ! 0,nk-1
         do_ji: do j =  iy,fy+1 !0,nj-1
              if (Dir==1) ip = k
              if (Dir==2) ip = j
              if (Dir==3) ip = i
            sump = 0.q0
                  Uz_in(j,k) = profile(Uj,Uc,xdir(ip),cent,Thick/2.,Thp*Thick)
                  if(Nspecies>1)then
                    Yi_in(:,j,k) = 0.
                    do n=1,Nspecies
                     Yi_in(n,j,k) =  profile(Yj(n), Yc(n),xdir(ip),cent,Thick/2.,Thp*Thick) 
                     sump = sump + Yi_in(n,j,k)
                    enddo 
                    if(sump.ne.1.D0) then
                     do s=1, Nspecies
                      Yi_in(s,j,k) = (Yi_in(s,j,k) / sump) * (1.D0 - sump) + Yi_in(s,j,k)
                     end do
                    end if
                  else
                     Yi_in(1,j,k) = 1.
                  endif
                  T_in(j,k) = profile(Tj,  Tc,xdir(ip),cent,Thick/2.,Thp*Thick)   
                  Ur_in(j,k) =  0.0
                  Ut_in(j,k) =  0.0
!               print*, "ATTENZIONE NON STAI STAMPANDO NIENTE IN ",file_name 
!                   write (unit=fh_BcIn,fmt=11,rec=ind) dbleq(dummy), dbleq(rr(j)), &
!                   dbleq(tt(k)) ,dbleq(Uz_in(j,k)), dbleq(Ur_in(j,k)), dbleq(Ut_in(j,k)),dbleq(T_in(j,k)), (dbleq(Yi_in(s,j,k)),s=1,Nspecies)
              write(app_st,'(a,i3,a)') "(e16.9,",6+Nspecies,"(1x,e15.6))"
              write(string,FMT=app_st) dbleq(dummy), dbleq(rr(j)), dbleq(tt(k)) , dbleq(Uz_in(j,k)), dbleq(Ur_in(j,k)), dbleq(Ut_in(j,k)), dbleq(T_in(j,k)), (dbleq(Yi_in(s,j,k)),s=1,Nspecies)
              string(recl:recl)=CHAR(10)
              write(unit=fh_BcIn, fmt='(a)',rec=ind) string(1:recl)
           ind = ind+1
         end do do_ji
        end do do_ki

        ENDIF

11      FORMAT (e16.9,<(6+Nspecies)>(1x,e15.6))
12      FORMAT (<(7+Nspecies)*16>a)
13      FORMAT (f12.2,1x,I3,< (7+Nspecies)*16>a)
16      FORMAT (<(Nspecies)>(1x,e15.6))

        close(fh_BcIN)

      call readInt("KLEIN", KLEIN)

      if( (KLEIN == 1).and.(iBoundZ(izone) == 1.)) then
       call KLEIN_pre(B,ix+1,fx,iy+1,fy,rr,tt,nk,nj)
      endif
      return

100   print *," OPSSS write ", stato
101   print *," OPSSS open  ", stato

      deallocate(iBoundZ)

return
end subroutine write_multizone









SUBROUTINE KLEIN_pre(B,ix,fx,iy,fy,rr,tt,nk,nj)
use file
use IBmod

implicit none

      character*30 file_name
      REAL  rr(0:nj-1), tt(0:nk-1),Thick, Thp, Thickl,Thpl,Uf_j,Uf_c
      INTEGER rec, count_jj,j,k,ind_,MAXj,MAXjp1,MAXk,NsMAX,ind,nj,nk,ix,fx,iy,fy,B
      real Uz_max, ni, Re, INTENS
      real PI
      real r(0:4000), Rmax, r_m(0:4000)
      real Uz_COLL(0:4000), Uz_STAG(0:4000), &
                        dUzdr_w_COLL, U_star_COLL, Re_star_COLL, &
                        dUzdr_w_STAG, U_star_STAG, Re_star_STAG, &
                        y_plus_COLL(0:4000), Yplus_COLL, &
                        y_plus_STAG(0:4000), Yplus_STAG, &
                        UiUi, UiUi_IN, UiUi_OUT, TRASL, & 
                        R_uu_COLL(0:4000), R_uv_COLL(0:4000), &
                        R_vv_COLL(0:4000), R_ww_COLL(0:4000), &
                        R_uw_COLL(0:4000), R_vw_COLL(0:4000), &
                        R_uu_STAG(0:4000), R_uv_STAG(0:4000), &
                        R_vv_STAG(0:4000), R_ww_STAG(0:4000), &
                        R_uw_STAG(0:4000), R_vw_STAG(0:4000), &
                        AV_Ruu_COLL(0:4000), SG_Ruu_COLL(0:4000), &
                        AV_Rvv_COLL(0:4000), SG_Rvv_COLL(0:4000), &
                        AV_Rww_COLL(0:4000), SG_Rww_COLL(0:4000), &
                        AV_Ruu_STAG(0:4000), SG_Ruu_STAG(0:4000), &
                        AV_Rvv_STAG(0:4000), SG_Rvv_STAG(0:4000), &
                        AV_Rww_STAG(0:4000), SG_Rww_STAG(0:4000), &
                        L_zz_COLL(0:4000), L_rr_COLL(0:4000), L_tt_COLL(0:4000), &
                        Lzz_COLL, Lrr_COLL, Ltt_COLL, &
                        L_zz_STAG(0:4000), L_rr_STAG(0:4000), L_tt_STAG(0:4000), &
                        Lzz_STAG, Lrr_STAG, Ltt_STAG
      real U_j, U_c, L_j, L_c, prof, a, sigma, al, sigmal, appo

      call readFloat("Thick", Thick) !; print*,"Thick",Thick ;
      call readFloat("Thp", Thp)     !; print*,"Thp",Thp ;
      call readFloat("Thickl", Thickl) !; print*,"Thick",Thick ;
      call readFloat("Thpl", Thpl)     !; print*,"Thp",Thp ;

      call readFloat("Uf_j", Uf_j)   !; print*,"Uf_j",Uf_j ;
      call readFloat("Uf_c", Uf_c)   !; print*,"Uf_c",Uf_c ;
      call readFloat("L_j", L_j)     !; print*,"L_j",L_j ;
      call readFloat("L_c", L_c)     !; print*,"L_c",L_c ;


! NON FARE TURBOLENZA ISOTROPA
       MAXj  = nj-2 
       MAXjp1= MAXj + 1
       MAXk  = nk-2 
       NsMAX = 2  

 
 !*****************************************************
 !SON SANDIA: isotropic turbulence
 
       a   = Thick / 2.0
       sigma = Thp * Thick
       al  = Thickl / 2.0
       sigmal = Thpl * Thickl
       do j = iy,fy!1,MAXj
         prof = Uf_c + (Uf_j - Uf_c) / 2.0 * ( dtanh( (rr(j)+a) / sigma ) - dtanh( (rr(j)-a) / sigma ) )
         appo = prof * prof
         R_uu_STAG(j) = appo
         R_uv_COLL(j) = 0.00
         R_uw_STAG(j) = 0.00
         R_vv_COLL(j) = appo
         R_vw_COLL(j) = 0.00
         R_ww_STAG(j) = appo
         prof = L_c + (L_j - L_c) / 2.0 * ( dtanh( (rr(j)+al) / sigmal ) - dtanh( (rr(j)-al) / sigmal ) )
         L_zz_STAG(j) = prof
         L_rr_COLL(j) = prof * 0.50
         L_tt_STAG(j) = prof * 0.50
         !PRINT*, r(j), appo, prof
       enddo
 !*****************************************************
         print*, Outfile
         ind_= index(Outfile,'_')
         print*, "ind_",ind_
         file_name= Outfile(1:ind_-1)//".BCin_TURB_0_1"
         print*, "file_name ",file_name
         write(file_name(index(file_name,"0_1"):index(file_name,"0_1")),"(i1.1)") B
         print*,"File BC_TURB: ",file_name
         if(B==1)then
          continue
         endif
         
       rec   = (9)* 16 + 1
       open (UNIT=20 ,file=file_name,ACTION='WRITE', FORM='FORMATTED', RECL=rec, ACCESS='DIRECT')
         ind = 1
         do k=ix,fx
           do j=iy,fy
             write(unit=20,fmt=13,rec=ind) dbleq(R_uu_STAG(j)), dbleq(R_uv_COLL(j)), dbleq(R_uw_STAG(j)), &
                          dbleq(R_vv_COLL(j)), dbleq(R_vw_COLL(j)), dbleq(R_ww_STAG(j)), &
                          dbleq(L_zz_STAG(j)), dbleq(L_rr_COLL(j)), dbleq(L_tt_STAG(j)), CHAR(10)
             ind = ind + 1
           end do
         end do
       close(20)
 
 10    FORMAT(A100)
 11    FORMAT (e16.9,16(1x,e15.5))
 12    FORMAT (e16.9,10(1x,e15.5))
 13    FORMAT (9(1x,e15.5),a1)
 return
 end subroutine KLEIN_pre

integer function lowerInd(xvalue, x, nx)
    integer nx
    real x(0:nx-1),xvalue
   
    do i=0,nx-1
        if (xvalue <= x(i)) then
         continue
         exit
        endif
    end do
    lowerInd = i
return
end

integer function upperInd(xvalue, x, nx)
    integer nx
    real x(0:nx-1),dx(0:nx-1),xvalue
    do i=0,nx-1
     dx(i)=abs(x(i)-xvalue)
    enddo   
!     do i=0,nx-1
!         if (xvalue <= x(i)) then
!          continue
!          exit
!         endif
!     end do
    upperInd = minloc(dx,1)
return
end


integer function findInd(xvalue, x, nx)
    integer nx
    real x(0:nx-1),dx(0:nx-1),xvalue
    do i=0,nx-1
     dx(i)=abs(x(i)-xvalue)
    enddo  
    findInd = minloc(dx,1) - 1 !-1 perchè il vettore inizia da 0
return
end
!*************************************************************
real function profile(U_j,U_c,r,cent,a,sigma)
real U_j,U_c,r,a,sigma,sigmal,al,cent

sigmal= sigma*1.
al=a*1.

profile = U_c + (U_j-U_c)/2.*(tanh( (r+al-cent)/(sigmal) ) - tanh( (r-al-cent)/(sigmal) ))

end function


real function senprofile(r,i,a,b)
real U_j,U_c
integer a,b,i
real r(0:b)

senprofile= sin(3.14159265359d0/(r(b-1)- r(a))*8.*r(i))!+300.


end function
!*************************************************************
!*************************************************************
!*************************************************************
!*************************************************************
!*************************************************************


real function U_l(alfa,phi,Yu_res,T,p)
real alfa,phi,Yu_res,T,p
real T0,P0,alfa_T,alfa_p,alfa_res,f
U_l = U0_lam(alfa,phi) * ((T/300.)**alfa_T(alfa,phi,T)) * ((p/1.0)**alfa_p(alfa,phi,T,p)) * (1.-alfa_res(alfa,phi,T,p)*f(Yu_res))
end function

real function U0_lam(alfa,phi)
real alfa,phi,a0,a1,a2,a3,a4,a5,a6,gamma_1
a0l=2.72q-6
a1l=2.897
a2l=150.817
a3l=4.539
a4l=-2.448
a5l=-0.0017
a6l=-0.2248
U0_lam = ((1. + a0l*alfa**a1l)*a2l*(phi**a3l)*exp(a4l*(phi+a5l*alfa+a6l)**2.))*gamma_1(alfa,phi)

end function


real function gamma_1(alfa,phi)
real alfa,phi

gamma_1= (1. + g(alfa,70.,10.)*(A(alfa)-1.))*(1.+g(alfa,90.,10.)*(GPHI(phi)-1.))

end function
! 
real function A(alfa)
real alfa,lam0,lam1,lam2,lam3,lam4,lam5
lam0= 1.4
lam1= 3.390q-7
lam2=-1.17q-6
lam3= 1.17q-7
lam4=-3.75q-9 
lam5= 1.4q-10

A= lam0 + lam1*alfa     &
        + lam2*alfa**2. &
        + lam3*alfa**3. &
        + lam4*alfa**4. &
        + lam5*alfa**5.

end function

real function GPHI(phi)
real phi,phi0,phi1,phi2,phi3,phi4
phi0=1.750
phi1=-1.750
phi2=0.625
phi3=-0.038
phi4=0.138

GPHI= phi0 + phi1*phi     &
           + phi2*phi**2. &
           + phi3*phi**3. &
           + phi4*phi**4. 

end function

real function g(alfa,alfa0,Delta_alfa)
real alfa,alfa0,Delta_alfa

g= 0.5*(1. + qtanh((alfa-alfa0)/Delta_alfa))

end function

real function alfa_T(alfa,phi,T)
real alfa,phi,T,alfaT0,alfaT1,alfaT2,alfaT3,alfaT4,lambda_T0,lambda_T1
alfaT0= 3.2466
alfaT1=-1.0709
alfaT2= 0.1517
alfaT3=-3.201q-4
alfaT4=-1.0359
lambda_T0=0.5
lambda_T1=0.58
alfa_T= (alfaT0 + alfaT1*phi     &
                + alfaT2*phi**4. &
                + alfaT3*alfa    &
                + alfaT4*((300./T)**2.)*phi**2.)*(1.+g(alfa,90.,10.)*(exp(lambda_T1*(T/300.)**lambda_T0 - 1.) - 1.))
!  print*, "DENTRO", (alfaT0 + alfaT1*phi + alfaT2*phi**4. + alfaT3*alfa + alfaT4*(300.**2./T**2.)*phi**2.), &
!                    (1.+g(alfa,90.,10.)) , &
!                    (exp(lambda_T1*(T/300.)**lambda_T0 - 1.)-1.),&
!                    alfa_T

end function

real function alfa_p(alfa,phi,T,p)
real alfa,phi,T,p,alfap0,alfap1,alfap2,alfap3,alfap4,lambda_p0,lambda_p1,lambda_p2
alfap0=-0.5406
alfap1= 0.1347
alfap2=-0.0125
alfap3=-5.174q-4
alfap4= 2.289q-4
lambda_p0=-1.9026
lambda_p1=3.556q-2
lambda_p2=-1.630q-4
alfa_p= (alfap0 + alfap1*phi     &
              + alfap2*phi**4. &
              + alfap3*alfa*((phi+1./phi)**0.5)    &
              + alfap4*(T/300.)*(p/1.0)*phi**2.)*(1.+g(alfa,90.,10.)*(exp(lambda_p0+lambda_p1*(p/1.0)+lambda_p2*(p/1.0)**2.)-1.))

end function

real function alfa_res(alfa,phi,T,p)
real alfa,phi,T,p,alfa_res0,alfa_res1,alfa_res2,alfa_res3,alfa_res4

alfa_res0= 4.157
alfa_res1=-1.744
alfa_res2= 0.5124
alfa_res3=-0.0047
alfa_res4=-8.694q-4

alfa_res=alfa_res0+alfa_res1*phi &
                  +alfa_res2*phi**4. &
                  +alfa_res3*alfa&
                  +alfa_res4*(T/300.)*(p/1.0)*phi**2.
end function

real function f(Yu_res)
real f1,f2

f1=-1.115
f2= 1.323

f=Yu_res+f1*Yu_res**2.+f2*Yu_res**3.

end function
!*************************************************************
!* This subroutine computes all eigenvalues and eigenvectors *
!* of a real symmetric square matrix A(N,N). On output, ele- *
!* ments of A above the diagonal are destroyed. D(N) returns *
!* the eigenvalues of matrix A. V(N,N) contains, on output,  *
!* the eigenvectors of A by columns. THe normalization to    *
!* unity is made by main program before printing results.    *
!* NROT returns the number of Jacobi matrix rotations which  *
!* were required.                                            *
!* --------------------------------------------------------- *
!* Ref.:"NUMERICAL RECIPES, Cambridge University Press, 1986,*
!*       chap. 11, pages 346-348" [BIBLI 08].                *
!*************************************************************
Subroutine Jacobi(A,N,D,V,NROT)
integer N,NROT
real::  A(1:N,1:N),D(1:N),V(1:N,1:N)
real, pointer :: B(:), Z(:)
real  c,g,h,s,sm,t,tau,theta,tresh
integer ip,iq,i,j, ialloc

allocate(B(1:100),stat=ialloc)
allocate(Z(1:100),stat=ialloc)

  do ip=1, N    !initialize V to identity matrix
    do iq=1, N
      V(ip,iq)=0.d0
    end do
      V(ip,ip)=1.d0
  end do
  do ip=1, N
    B(ip)=A(ip,ip)
    D(ip)=B(ip)
    Z(ip)=0.d0
  end do
  NROT=0
  do i=1, 50
    sm=0.d0
    do ip=1, N-1     !sum off-diagonal elements
      do iq=ip+1, N
        sm=sm+DABS(A(ip,iq))
      end do
    end do
    if(sm==0.d0) then; deallocate(B,Z); return; end if  !normal return
    if(i.lt.4) then
      tresh=0.2d0*sm**2
    else
      tresh=0.d0
    end if
    do ip=1, N-1
      do iq=ip+1, N
        g=100.d0*DABS(A(ip,iq))
! after 4 sweeps, skip the rotation if the off-diagonal element is small
        if((i.gt.4).and.(DABS(D(ip))+g.eq.DABS(D(ip))) &
		.and.(DABS(D(iq))+g.eq.DABS(D(iq)))) then
		  A(ip,iq)=0.d0
        else if(DABS(A(ip,iq)).gt.tresh) then
	  h=D(iq)-D(ip)
	  if(DABS(h)+g.eq.DABS(h)) then
	    t=A(ip,iq)/h
          else
	    theta=0.5d0*h/A(ip,iq)
            t=1.d0/(DABS(theta)+DSQRT(1.d0+theta**2))
	    if(theta.lt.0.d0) t=-t
          end if
	  c=1.d0/DSQRT(1.d0+t**2)
	  s=t*c
          tau=s/(1.d0+c)
	  h=t*A(ip,iq)
	  Z(ip)=Z(ip)-h
	  Z(iq)=Z(iq)+h
	  D(ip)=D(ip)-h
	  D(iq)=D(iq)+h
	  A(ip,iq)=0.d0
	  do j=1, ip-1
	    g=A(j,ip)
	    h=A(j,iq)
	    A(j,ip)=g-s*(h+g*tau)
	    A(j,iq)=h+s*(g-h*tau)
          end do
	  do j=ip+1, iq-1
	    g=A(ip,j)
	    h=A(j,iq)
	    A(ip,j)=g-s*(h+g*tau)
	    A(j,iq)=h+s*(g-h*tau)
          end do
	  do j=iq+1, N
	    g=A(ip,j)
	    h=A(iq,j)
	    A(ip,j)=g-s*(h+g*tau)
	    A(iq,j)=h+s*(g-h*tau)
          end do
	  do j=1, N
	    g=V(j,ip)
	    h=V(j,iq)
	    V(j,ip)=g-s*(h+g*tau)
	    V(j,iq)=h+s*(g-h*tau)
          end do
          NROT=NROT+1
        end if !if ((i.gt.4)...
      end do !main iq loop
    end do !main ip loop
    do ip=1, N
      B(ip)=B(ip)+Z(ip)
      D(ip)=B(ip)
      Z(ip)=0.d0
    end do
  end do !main i loop
! pause ' 50 iterations !'
  deallocate(B,Z)
  return
END

! end of file ujacobi.f90

!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
!*********************************************
