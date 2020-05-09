c -----------------------------------------------------------------
c       This is the spectral code for three-D homogeneous shear 
c       turbulence using the message-passing in Origin 2000 machines
c       Oct. 1998 
c       nx = ny = nz = 2**n   ( n = 4,5,...)
c       nx = x-diretion mesh size in each node
c       ny = y-diretion mesh size in each node
c       nz = z-direction mesh size in each node
c       nproc = number of process (nodes)
c       lx=nx,ly=ny,lz=nz/nproc
c     
c       time = time - remap time transform included oct 22
c ----------------------------------------------------------------- 

       program three_spectral 

       include 'include.h'

       complex,allocatable,dimension(:,:,:)::vx,vy,vz
       complex,allocatable,dimension(:,:,:)::wx,wy,wz
       complex,allocatable,dimension(:,:,:)::wx_n,wy_n
       complex,allocatable,dimension(:,:,:)::wz_n,ww_n
       complex,allocatable,dimension(:,:,:)::ox,oy,oz
       complex,allocatable,dimension(:,:,:)::ww
       real,allocatable,dimension(:,:,:)::kx,ky,kz,k2,k_plus
       real,allocatable,dimension(:,:,:)::tmp,k2_e
       integer,allocatable,dimension(:)::iseed 
       complex, allocatable, dimension(:)::coeffy 
       real, allocatable, dimension(:)::coeffxz
       real, allocatable, dimension(:)::aaa1
       real force(2),time0, ss(500), v_square
       real vx_square, vy_square, vz_square, v_s

c  setup MPP environment

       call mpi_init(ierror)
       call mpi_comm_rank(mpi_comm_world,id,ierror)
       call mpi_comm_size(mpi_comm_world,nproc,ierror)
       nallgrp = mpi_comm_world

c   read parameters from node id = 0

       ak0 = 4.75683

       if(id.eq.0)then
         open(90,file='parameter.d',status='unknown')
         write(*,*)'enter nx,ny,nz'
         read(90,*)nx,ny,nz
         write(*,*)'enter number of processors'
         read(90,*)nprocs
         write(*,*) 'enter seed (321)'
         read(90,*) iseed_m0
         write(*,*)'enter nstep'
         read(90,*) nstep
         write(*,*)'enter tout'
         read(90,*) itout
         write(*,*)'enter eout'
         read(90,*) ieout
         write(*,*)'enter dt'
         read(90,*) dt
         write(*,*)'enter rnu'
         read(90,*) rnu
         write(*,*)'enter u0'
         read(90,*) u0
         write(*,*)'enter time'
         read(90,*) time
         write(*,*)'enter new'
         read(90,*) new
         write(*,*)'enter idp for rerun'
         read(90,*) idp
         write(*,*)'enter force'
         read(90,*) force
         write(*,*)'shear rate'
         read(90,*) srate
         write(*,*)nx,ny,nz,nprocs,iseed_m0,nstep,itout,ieout,dt,rnu,
     1     u0,time,new,idp,force,srate
       endif 

       call mpi_bcast(nx,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(ny,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(nz,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(nprocs,1,MPI_INTEGER,0,nallgrp,ierror)      
       call mpi_bcast(iseed_m0,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(nstep,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(itout,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(ieout,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(dt,1,MPI_REAL,0,nallgrp,ierror)
       call mpi_bcast(rnu,1,MPI_REAL,0,nallgrp,ierror)
       call mpi_bcast(u0,1,MPI_REAL,0,nallgrp,ierror)
       call mpi_bcast(time,1,MPI_REAL,0,nallgrp,ierror)
       call mpi_bcast(new,1,MPI_LOGICAL,0,nallgrp,ierror)
       call mpi_bcast(idp,1,MPI_INTEGER,0,nallgrp,ierror)
       call mpi_bcast(force,2,MPI_REAL,0,nallgrp,ierror)
       call mpi_bcast(srate,1,MPI_REAL,0,nallgrp,ierror)

       lx = nx
       ly = ny
       lz = nz/nproc
       lx1=lx+1
       lx2= lx*2
       lly=ly/nproc
       pi=4.*atan(1.)
       dx=1.0/pi/nx

       scale = 1.0/(nx*ny*nz*2.0)

       if(srate.gt.1.0e-10)then
         remap = 2.*nx/ny/srate
      else
         remap = 1.0e+10
      end if

c       remap = 1.0e+10

        time0 = 0.
     
c allocate memory 
      
       allocate( vx(lx1,ly,lz) )
       allocate( vy(lx1,ly,lz) )
       allocate( vz(lx1,ly,lz) )
       allocate( wx(lx1,ly,lz) )
       allocate( wy(lx1,ly,lz) )
       allocate( wz(lx1,ly,lz) )
       allocate( wx_n(lx1,ly,lz) )
       allocate( wy_n(lx1,ly,lz) )
       allocate( wz_n(lx1,ly,lz) )
       allocate( ww_n(lx1,ly,lz) )
       allocate( ox(lx1,ly,lz) )
       allocate( oy(lx1,ly,lz) )
       allocate( oz(lx1,ly,lz) )
       allocate( ww(lx1,ly,lz) )
       allocate( kx(lx1,ly,lz) )
       allocate( ky(lx1,ly,lz) )
       allocate( kz(lx1,ly,lz) )
       allocate( k2(lx1,ly,lz) )
       allocate( k2_e(lx1,ly,lz) )
       allocate( k_plus(lx1,ly,lz) )
       allocate( tmp(lx1,ly,lz) )
       allocate( iseed(nproc) )
       allocate( coeffy( ly+15 ) )
       allocate( coeffxz( (lx2+15)+2*(nz+15) ) )
       allocate( aaa1(nproc) )

       if(id.eq.0)print*,' memory allocated'

       is = 0
       iii = 1
       jjj = 1
       kkk = 1

       nek=0.9*0.5*sqrt(4.0*nx*nx+ny*ny+nz*nz)
        
       if(id.eq.0)then
         write(*,*)ieout,itout,nek
         open(11,file='energy.d',status='unknown',position='append')
         open(20,file='spectrum.d',status='unknown',position='append') 
         open(22,file='spectrum.x',status='unknown',position='append')
         open(24,file='spectrum.y',status='unknown',position='append')
         open(26,file='spectrum.z',status='unknown',position='append')
         open(36,file='ensthr.d',status='unknown',position='append') 
         open(70,file='initial.sp',status='unknown',position='append')
       end if
        
c---Generate initial condition -------

       if(new) then

c---Generate random number seed for each process ---

         if(id.eq.0)then
           nremesh=0
           call srand(iseed_m0)
           do i=1,nproc
             iseed(i) = irand()
           end do
         end if 
         call mpi_bcast(iseed,nproc,MPI_INTEGER,0,nallgrp,ierror)
         irr = iseed(id+1)
         call srand(irr)

         call gaussian(vx,pi)
         call gaussian(vy,pi)
         call gaussian(vz,pi)

         call wavenumber(kx,ky,kz,k2,k2_e,time,dt,rnu,id)

         where((k2.lt.256).or.(k2.gt.1024)) vx = (0.0,0.0) !modes k<16
         where((k2.lt.256).or.(k2.gt.1024)) vy = (0.0,0.0) !and
         where((k2.lt.256).or.(k2.gt.1024)) vz = (0.0,0.0) !k>32 rogal

         kx = vx * conjg(vx)            !   Measure tot. energy 
         vx_square = sum(kx) 
         kx = vy * conjg(vy)
         vy_square = sum(kx)
         kx = vz * conjg(vz)
         vz_square = sum(kx) 
         v_square = vx_square + vy_square + vz_square
      call mpi_reduce(v_square,v_s,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
         if (id.eq.0) write(*,100) istep, time, .5*v_s
         call mpi_bcast(v_s,1,MPI_REAL,0,nallgrp,ierror) 

         call wavenumber(kx,ky,kz,k2,k2_e,time,dt,rnu,id)

         tmp = (kx*real(vx) + ky*real(vy) + kz*real(vz))/k2
         vx = cmplx(real(vx) - kx*tmp, aimag(vx))
         vy = cmplx(real(vy) - ky*tmp, aimag(vy))
         vz = cmplx(real(vz) - kz*tmp, aimag(vz))
         tmp = (kx*aimag(vx) + ky*aimag(vy) + kz*aimag(vz))/k2
         vx = cmplx(real(vx), aimag(vx) - kx*tmp)
         vy = cmplx(real(vy), aimag(vy) - ky*tmp)
         vz = cmplx(real(vz), aimag(vz) - kz*tmp)
      
c         cc1 = u0 * sqrt(8.*sqrt(2./pi)/(3.*pi*ak0**5))
c         tmp = cc1 * sqrt(k2) * exp (-k2/(ak0*ak0))

c         vx = vx * tmp
c         vy = vy * tmp
c         vz = vz * tmp

         call symmetrize(vx,id)
         call symmetrize(vy,id)
         call symmetrize(vz,id)

c     calculating initial spectrum

         tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
         tmp(1,:,:) = 0.5*tmp(1,:,:)
         do i=1,nek
           ek = 0.0
           e_t = 0.0
           ek=sum(tmp(1:lx,:,:),mask=(abs(sqrt(k2)-i-0.49999).lt.0.5))
           call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
           call mpi_bcast(e_t,1,MPI_REAL,0,nallgrp,ierror)
           ss(i) = e_t
         end do
                  
         do i=1,nek
          if ((ss(i).gt.1.e-4).and.(i.ge.16).and.(i.le.32)) then
          where (abs(sqrt(k2)-i-0.49999).lt.0.5) vx=vx/sqrt(ss(i))
          where (abs(sqrt(k2)-i-0.49999).lt.0.5) vy=vy/sqrt(ss(i))
          where (abs(sqrt(k2)-i-0.49999).lt.0.5) vz=vz/sqrt(ss(i))
          elseif((i.ge.16).and.(i.le.32)) then
             if (id.eq.0) write(70,*) 'warning: e_t < 1.e-4'
          endif
         end do

         tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
         tmp(1,:,:) = 0.5*tmp(1,:,:)
         do i=1,nek
           ek = 0.0
           e_t = 0.0
           ek=sum(tmp(1:lx,:,:),mask=(abs(sqrt(k2)-i-0.49999).lt.0.5))
           call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
           if(id.eq.0)then
             write(70,*)i,e_t
           end if
         end do
         if(id.eq.0)close(70)
       else

c ------ Rreading from DISK ------ 

         call input(vx,vy,vz,idp,id,nallgrp,time,time0,nremesh)
         call symmetrize(vx,id)
         call symmetrize(vy,id)
         call symmetrize(vz,id)

         if(id.eq.0)then
           write(*,*)'after reading'
         end if
c
      endif

      call wavenumber(kx,ky,kz,k2,k2_e,time,dt,rnu,id)

      call cfft1di(ly,coeffy)
      call scfft2dui(lx2,nz,coeffxz)

      istep=-1
 
1     istep = istep + 1

c   write out the total energy in K space

       time1 = mpi_wtime()
       tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
       tmp(1,:,:)=0.5*tmp(1,:,:)
       ek = sum(tmp)
       call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
       if(id.eq.0)then
         write(*,*)istep,time,e_t
         write(11,*)time,e_t
       end if

       call wavenumber(kx,ky,kz,k2,k2_e,time,dt,rnu,id)

       wx = (0.,1.) * (ky*vz - kz*vy)
       wy = (0.,1.) * (kz*vx - kx*vz)
       wz = (0.,1.) * (kx*vy - ky*vx)
       tmp = wx*conjg(wx) + wy*conjg(wy) + wz*conjg(wz)
       tmp(1,:,:)=0.5*tmp(1,:,:)
       ek = sum(tmp)
       call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
       if(id.eq.0) write(36,*)istep,e_t
 
       if(mod(istep,itout).eq.0.and.istep.ne.0) then
        call output(vx,vy,vz,iii,id,nallgrp,time,time0,nremesh)
         iii = iii + 1
       endif

c     energy spectrum calculation 
     
      if(mod(istep,ieout).eq.0) then
        tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
        tmp(1,:,:)=0.5*tmp(1,:,:)
        if(id.eq.0)then   
          write(20,*)(kkk-1)
        end if
        do i=1,nek
          ek = 0.0
          e_t = 0.0
          ek=sum(tmp,mask=(abs(sqrt(k2)-i-0.49999).lt.0.5))
          call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
          if(id.eq.0)then
            write(20,*)i,e_t
          end if
        end do
       
         tmp = vx*conjg(vx)
         tmp(1,:,:) = .5 * tmp(1,:,:)
          if(id.eq.0)then
          write(22,*)(kkk-1)
        end if
          do i=1,nek
          ek = 0.0
          e_t = 0.0
          ek=sum(tmp,mask=(abs(sqrt(k2)-i-0.49999).lt.0.5))
          call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
          if(id.eq.0)then
            write(22,*)i,e_t
          end if
        end do

         tmp = vy*conjg(vy)
         tmp(1,:,:) = .5 * tmp(1,:,:)
          if(id.eq.0)then
          write(24,*)(kkk-1)
        end if
          do i=1,nek
          ek = 0.0
          e_t = 0.0
          ek=sum(tmp,mask=(abs(sqrt(k2)-i-0.49999).lt.0.5))
          call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
          if(id.eq.0)then
            write(24,*)i,e_t
          end if
        end do

         tmp = vz*conjg(vz)
         tmp(1,:,:) = .5 * tmp(1,:,:)
          if(id.eq.0)then
          write(26,*)(kkk-1)
        end if
          do i=1,nek
          ek = 0.0
          e_t = 0.0
          ek=sum(tmp,mask=(abs(sqrt(k2)-i-0.49999).lt.0.5))
          call mpi_reduce(ek,e_t,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
          if(id.eq.0)then
            write(26,*)i,e_t
          end if
        end do

        kkk = kkk + 1
 
      end if
        
      wx=vx
      wy=vy
      wz=vz

      ww_n = cmplx(cos(pi/ny*kplus), sin(pi/ny*kplus))

      wx_n = wx*ww_n
      wy_n = wy*ww_n
      wz_n = wz*ww_n

      ox=vx
      oy=vy
      oz=vz

c
      isign = 1
c    
      call mpifft(wx,isign,coeffy,coeffxz,id,nallgrp)
      call mpifft(wy,isign,coeffy,coeffxz,id,nallgrp)
      call mpifft(wz,isign,coeffy,coeffxz,id,nallgrp)      

      call mpifft(wx_n,isign,coeffy,coeffxz,id,nallgrp)
      call mpifft(wy_n,isign,coeffy,coeffxz,id,nallgrp)
      call mpifft(wz_n,isign,coeffy,coeffxz,id,nallgrp)

      aaa1=0.0
      do i=0,nproc-1
       if(id.eq.i) aaa1(i+1)=max(maxval(sqrt(real(wx)**2+real(wy)**2
     .         +real(wz)**2)),maxval(sqrt(aimag(wx)**2+aimag(wy)**2
     .         +aimag(wz)**2)))
      end do

      do i=0,nproc-1
        call mpi_bcast(aaa1(i+1),1,MPI_REAL,i,nallgrp,ierror)
      end do

      dt=dx/maxval(aaa1)*0.1
      dt=min(dt,remap-time0)
      if(id.eq.0) write(*,*)'dt=',dt

      call wavenumber(kx,ky,kz,k2,k2_e,time,dt,rnu,id)

      a1=(sum(real(wx)*real(wy))+sum(aimag(wx)*aimag(wy)))*scale
      call mpi_reduce(a1,a1_t,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror) 
      if(id.eq.0) write(35,*)istep,a1_t

c     v 

      tmp=(-k2+2.0*srate*time*kx*ky-srate**2*time*time*kx*kx)

      vx=-srate*(1.0+2.0*kx*kx/tmp)*vy
      vz=-2.0*kx*kz*srate/tmp*vy
      vy=-2.0*kx*srate*(ky-srate*time*kx)/tmp*vy

c     u*u

      ww=cmplx(real(wx)*real(wx),aimag(wx)*aimag(wx))
      ww_n=cmplx(real(wx_n)*real(wx_n),aimag(wx_n)*aimag(wx_n))

      isign = -1
      call mpifft(ww,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww,id)
      call mpifft(ww_n,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww_n,id)
      ww_n=ww_n*cmplx(cos(pi/ny*kplus),-sin(pi/ny*kplus))

      ww = 0.5*(ww + ww_n) 

      vx=vx-(0.,1.)*kx*(1.0+kx*kx/tmp)*ww
      vy=vy-(0.,1.)*(ky-srate*time*kx)*kx*kx/tmp*ww
      vz=vz-(0.,1.)*kz*kx*kx/tmp*ww

c     v*v

      ww=cmplx(real(wy)*real(wy),aimag(wy)*aimag(wy))
      ww_n=cmplx(real(wy_n)*real(wy_n),aimag(wy_n)*aimag(wy_n))

      isign = -1
      call mpifft(ww,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww,id)
      call mpifft(ww_n,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww_n,id)
      ww_n=ww_n*cmplx(cos(pi/ny*kplus),-sin(pi/ny*kplus))

      ww = 0.5*(ww + ww_n)

      vx=vx-(0.,1.)*kx*(ky-srate*time*kx)**2/tmp*ww
      vy=vy-(0.,1.)*((ky-srate*time*kx)**3/tmp+ky-kx*srate*time)*ww
      vz=vz-(0.,1.)*kz*(ky-srate*time*kx)**2/tmp*ww

c     w*w

      ww=cmplx(real(wz)*real(wz),aimag(wz)*aimag(wz))
      ww_n=cmplx(real(wz_n)*real(wz_n),aimag(wz_n)*aimag(wz_n))

      isign = -1
      call mpifft(ww,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww,id)
      call mpifft(ww_n,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww_n,id)
      ww_n=ww_n*cmplx(cos(pi/ny*kplus),-sin(pi/ny*kplus))

      ww = 0.5*(ww + ww_n)

      vx=vx-(0.,1.)*kx*kz*kz/tmp*ww
      vy=vy-(0.,1.)*(ky-srate*time*kx)*kz*kz/tmp*ww
      vz=vz-(0.,1.)*(kz*kz*kz/tmp+kz)*ww

c     u*v

      ww=cmplx(real(wx)*real(wy),aimag(wx)*aimag(wy))
      ww_n=cmplx(real(wx_n)*real(wy_n),aimag(wx_n)*aimag(wy_n))

      isign = -1
      call mpifft(ww,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww,id)
      call mpifft(ww_n,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww_n,id)
      ww_n=ww_n*cmplx(cos(pi/ny*kplus),-sin(pi/ny*kplus))

      ww = 0.5*(ww + ww_n)

      vx=vx-(0.,1.)*(ky-kx*srate*time)*(2.0*kx*kx/tmp+1.)*ww
      vy=vy-(0.,1.)*kx*(2.0*(ky-srate*time*kx)**2/tmp + 1.)*ww
      vz=vz-(0.,1.)*kz*2.0*kx*(ky-srate*time*kx)/tmp*ww

c     u*w

      ww=cmplx(real(wx)*real(wz),aimag(wx)*aimag(wz))
      ww_n=cmplx(real(wx_n)*real(wz_n),aimag(wx_n)*aimag(wz_n))

      isign = -1
      call mpifft(ww,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww,id)
      call mpifft(ww_n,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww_n,id)
      ww_n=ww_n*cmplx(cos(pi/ny*kplus),-sin(pi/ny*kplus))

      ww = 0.5*(ww + ww_n)

      vx=vx-(0.,1.)*kz*(2.0*kx*kx/tmp + 1.)*ww
      vy=vy-(0.,1.)*2.0*kx*kz*(ky-srate*time*kx)/tmp*ww
      vz=vz-(0.,1.)*kx*(2.0*kz*kz/tmp + 1.)*ww

c     v*w

      ww=cmplx(real(wy)*real(wz),aimag(wy)*aimag(wz))
      ww_n=cmplx(real(wy_n)*real(wz_n),aimag(wy_n)*aimag(wz_n))

      isign = -1
      call mpifft(ww,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww,id)
      call mpifft(ww_n,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww_n,id)
      ww_n=ww_n*cmplx(cos(pi/ny*kplus),-sin(pi/ny*kplus))

      ww = 0.5*(ww + ww_n)

      vx=vx-(0.,1.)*kx*2.0*kz*(ky-srate*time*kx)/tmp*ww
      vy=vy-(0.,1.)*kz*(2.0*(ky-srate*time*kx)**2/tmp + 1.)*ww
      vz=vz-(0.,1.)*(ky-srate*time*kx)*(2.0*kz*kz/tmp + 1.)*ww

c
      ww=ox
      ox=(ox+vx*dt/2.0)*k2_e
      vx=ww+vx*dt

      ww=oy
      oy=(oy+vy*dt/2.0)*k2_e
      vy=ww+vy*dt

      ww=oz
      oz=(oz+vz*dt/2.0)*k2_e
      vz=ww+vz*dt

c  2nd step of R-K

      time=time+dt
      time0 = time0 + dt

      wx=vx
      wy=vy
      wz=vz

      ww_n = cmplx(cos(pi/ny*kplus), sin(pi/ny*kplus))

      wx_n = wx*ww_n
      wy_n = wy*ww_n
      wz_n = wz*ww_n

c
      isign = 1
c
      call mpifft(wx,isign,coeffy,coeffxz,id,nallgrp)
      call mpifft(wy,isign,coeffy,coeffxz,id,nallgrp)
      call mpifft(wz,isign,coeffy,coeffxz,id,nallgrp)
      call mpifft(wx_n,isign,coeffy,coeffxz,id,nallgrp)
      call mpifft(wy_n,isign,coeffy,coeffxz,id,nallgrp)
      call mpifft(wz_n,isign,coeffy,coeffxz,id,nallgrp)

      tmp=(-k2+2.0*srate*time*kx*ky-srate**2*time*time*kx*kx)

      vx=-srate*(1.0+2.0*kx*kx/tmp)*vy
      vz=-2.0*kx*kz*srate/tmp*vy
      vy=-2.0*kx*srate*(ky-srate*time*kx)/tmp*vy

c     u*u

      ww=cmplx(real(wx)*real(wx),aimag(wx)*aimag(wx))
      ww_n=cmplx(real(wx_n)*real(wx_n),aimag(wx_n)*aimag(wx_n))

      isign = -1
      call mpifft(ww,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww,id)
      call mpifft(ww_n,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww_n,id)
      ww_n=ww_n*cmplx(cos(pi/ny*kplus),-sin(pi/ny*kplus))

      ww = 0.5*(ww + ww_n)

      vx=vx-(0.,1.)*kx*(1.0+kx*kx/tmp)*ww
      vy=vy-(0.,1.)*(ky-srate*time*kx)*kx*kx/tmp*ww
      vz=vz-(0.,1.)*kz*kx*kx/tmp*ww

c     v*v

      ww=cmplx(real(wy)*real(wy),aimag(wy)*aimag(wy))
      ww_n=cmplx(real(wy_n)*real(wy_n),aimag(wy_n)*aimag(wy_n))

      isign = -1
      call mpifft(ww,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww,id)
      call mpifft(ww_n,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww_n,id)
      ww_n=ww_n*cmplx(cos(pi/ny*kplus),-sin(pi/ny*kplus))

      ww = 0.5*(ww + ww_n)

      vx=vx-(0.,1.)*kx*(ky-srate*time*kx)*(ky-srate*time*kx)/tmp*ww
      vy=vy-(0.,1.)*((ky-srate*time*kx)**3/tmp+ky-kx*srate*time)*ww
      vz=vz-(0.,1.)*kz*(ky-srate*time*kx)*(ky-srate*time*kx)/tmp*ww

c     w*w

      ww=cmplx(real(wz)*real(wz),aimag(wz)*aimag(wz))
      ww_n=cmplx(real(wz_n)*real(wz_n),aimag(wz_n)*aimag(wz_n))

      isign = -1
      call mpifft(ww,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww,id)
      call mpifft(ww_n,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww_n,id)
      ww_n=ww_n*cmplx(cos(pi/ny*kplus),-sin(pi/ny*kplus))

      ww = 0.5*(ww + ww_n)

      vx=vx-(0.,1.)*kx*kz*kz/tmp*ww
      vy=vy-(0.,1.)*(ky-srate*time*kx)*kz*kz/tmp*ww
      vz=vz-(0.,1.)*(kz*kz*kz/tmp+kz)*ww

c     u*v

      ww=cmplx(real(wx)*real(wy),aimag(wx)*aimag(wy))
      ww_n=cmplx(real(wx_n)*real(wy_n),aimag(wx_n)*aimag(wy_n))

      isign = -1
      call mpifft(ww,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww,id)
      call mpifft(ww_n,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww_n,id)
      ww_n=ww_n*cmplx(cos(pi/ny*kplus),-sin(pi/ny*kplus))

      ww = 0.5*(ww + ww_n)

      vx=vx-(0.,1.)*(ky-kx*srate*time)*(2.0*kx*kx/tmp+1.)*ww
      vy=vy-(0.,1.)*kx*(2.0*(ky-srate*time*kx)**2/tmp + 1.)*ww
      vz=vz-(0.,1.)*kz*2.0*kx*(ky-srate*time*kx)/tmp*ww

c     u*w

      ww=cmplx(real(wx)*real(wz),aimag(wx)*aimag(wz))
      ww_n=cmplx(real(wx_n)*real(wz_n),aimag(wx_n)*aimag(wz_n))

      isign = -1
      call mpifft(ww,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww,id)
      call mpifft(ww_n,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww_n,id)
      ww_n=ww_n*cmplx(cos(pi/ny*kplus),-sin(pi/ny*kplus))

      ww = 0.5*(ww + ww_n)

      vx=vx-(0.,1.)*kz*(2.0*kx*kx/tmp + 1.)*ww
      vy=vy-(0.,1.)*2.0*kx*kz*(ky-srate*time*kx)/tmp*ww
      vz=vz-(0.,1.)*kx*(2.0*kz*kz/tmp + 1.)*ww

c     v*w

      ww=cmplx(real(wy)*real(wz),aimag(wy)*aimag(wz))
      ww_n=cmplx(real(wy_n)*real(wz_n),aimag(wy_n)*aimag(wz_n))

      isign = -1
      call mpifft(ww,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww,id)
      call mpifft(ww_n,isign,coeffy,coeffxz,id,nallgrp)
      call symmetrize(ww_n,id)
      ww_n=ww_n*cmplx(cos(pi/ny*kplus),-sin(pi/ny*kplus))

      ww = 0.5*(ww + ww_n)

      vx=vx-(0.,1.)*kx*2.0*kz*(ky-srate*time*kx)/tmp*ww
      vy=vy-(0.,1.)*kz*(2.0*(ky-srate*time*kx)**2/tmp + 1.)*ww
      vz=vz-(0.,1.)*(ky-srate*time*kx)*(2.0*kz*kz/tmp + 1.)*ww

c
      vx=ox+vx*dt/2.0
      vy=oy+vy*dt/2.0
      vz=oz+vz*dt/2.0

c Remesh

      if(nremesh.eq.0) then
        if(time.le.(remap).and.(time0-0.5*remap).ge.0)then
          time0 = time0-0.5*remap
          time = time - remap
          nremesh=nremesh+1
          if(id.eq.0)write(*,*)'nremesh=',nremesh
          ox = vx
          oy = vy
          oz = vz

          vx=0.0
          vy=0.0
          vz=0.0

          do i=1,nx
            do j=1,ny/2
              kk=j+1-i
              if(kk.le.0) kk=kk+ny
              vx(i,kk,:)=ox(i,j,:)
              vy(i,kk,:)=oy(i,j,:)
              vz(i,kk,:)=oz(i,j,:)
            end do

            do j=ny/2+1,ny
              kk=j+1-i
              if(kk.gt.ny/2)then
                vx(i,kk,:)=ox(i,j,:)
                vy(i,kk,:)=oy(i,j,:)
                vz(i,kk,:)=oz(i,j,:)
              end if
            end do
          end do

        end if
      end if

      if((time0-remap).ge.0)then

        time0 = time0-remap
        time = time - remap
        nremesh=nremesh+1
        if(id.eq.0)write(*,*)'nremesh=',nremesh
        ox = vx
        oy = vy
        oz = vz

        vx=0.0
        vy=0.0
        vz=0.0

          do i=1,nx
            do j=1,ny/2
              kk=j+1-i
              if(kk.le.0) kk=kk+ny
              vx(i,kk,:)=ox(i,j,:)
              vy(i,kk,:)=oy(i,j,:)
              vz(i,kk,:)=oz(i,j,:)
            end do

            do j=ny/2+1,ny
              kk=j+1-i
              if(kk.gt.ny/2)then
                vx(i,kk,:)=ox(i,j,:)
                vy(i,kk,:)=oy(i,j,:)
                vz(i,kk,:)=oz(i,j,:)
              end if
            end do
          end do

      end if

c     call supfor(vx,vy,vz,tmp,k2,force,id,nallgrp)

       time2=mpi_wtime()
        if(id.eq.0) then
          write(*,fmt='(5X,A,F20.4,A)')
     &   'elapsed time =',1.0D3*(time2-time1),' msec'
        end if

      if(istep.le.nstep) goto 1

      if(id.eq.0)write(*,*)'nremesh=',nremesh

      write(*,*)id,'finished '

      if(id.eq.0)then
        close(11)
        close(20)
        close(22)
        close(24)
        close(26)
        close(36)
      end if

100   format(' step',i4,':  t = ',e14.6,' E = ',e14.6,' W = ',e14.6)
200   format(8e13.5)
102   format(2(1x,E16.7))

      call MPI_FINALIZE(ierror)
      stop
      end

c --------------------------------------------------------------
      subroutine supfor(vx,vy,vz,tmp,k2,force,id,nallgrp)

      include 'include.h'

      complex,dimension(lx1,ly,lz)::vx,vy,vz
      real,dimension(lx1,ly,lz)::k2,tmp
      real force(2),ek,ff

      tmp = vx*conjg(vx) + vy*conjg(vy) + vz*conjg(vz)
      tmp(1,:,:)=0.5*tmp(1,:,:)
      do is = 1, 2
        ek = 0.0
        ff = 0.0
        ek=sum(tmp(1:lx,:,:),mask=(abs(sqrt(k2)-is-0.49999).lt.0.5))
        call mpi_reduce(ek,ff,1,MPI_REAL,MPI_SUM,0,nallgrp,ierror)
        call mpi_bcast(ff,1,MPI_REAL,0,nallgrp,ierror)
        ff = sqrt(force(is)/ff)
        where(abs(sqrt(k2)-is-0.49999).lt.0.5)
          vx=vx*ff
          vy=vy*ff
          vz=vz*ff
        end where
      end do

      return
      end

c-----------------------------------------------------------------
       subroutine wavenumber(kx,ky,kz,k2,k2_e,time,dt,rnu,id)

       include 'include.h'

       real,dimension(lx1,ly,lz)::kx,ky,kz,k2,k2_e

       do i=1,lx1
         kx(i,:,:)= i-1 
       end do

       do j=1,ly
         ky(:,j,:) = mod(j-1+ly/2,ly)-ly/2
       end do
   
       do k=1,lz
         k1=id*lz+k
         kz(:,:,k) = mod(k1-1+nz/2,nz)-nz/2
       end do
 
       k2=kx*kx+ky*ky+kz*kz
       k2_e = exp(-rnu*dt*(k2-kx*ky*srate*(2.0*time+dt)
     .     +srate*srate*kx*kx*(3.0*time*time+3.0*time*dt+dt*dt)/3.0))
 
       if(id.eq.0)then
         k2(1,1,1)=1.e-5
         k2_e(1,1,1) = 1.0
       end if

       return
       end

c --------------------------------------------------------------------
      subroutine gaussian(u,pi)

      include 'include.h'

      complex,dimension(lx1,ly,lz) :: u

      u = 0.0

      do i=1,lx1
       do j=1,ly
        do k=1,lz
          t1 = rand()
          t2 = rand()
          if(t1.le.1.e-10) t1 = 1.e-10
          if(t2.le.1.e-10) t2 = 1.e-10
          t2 = 2*pi*t2
          u(i,j,k)=sqrt(-2.0*log(t1))*cmplx(cos(t2),sin(t2))
        end do
       end do
      end do
  
      return
      end

c ---------------------------------------------------------------------
      subroutine symmetrize(c,id)

      include 'include.h'
 
      complex,dimension(lx1,ly,lz)::c
      complex,dimension(ly,lz)::cc

c   zero k2=ny/2 modes
      c(:,ly/2+1,:) = (0.0,0.0)
c   zero k3=nz/2 modes
      nid=(nz/2+1)/lz
      nip=nz/2+1-nid*lz
      if(id.eq.nid) c(:,:,nip) = (0.0,0.0)
c   zero k1=nx modes
      c(lx1,:,:)=(0.0,0.0)
c  zero k1=k2=k3=0 mode
      if(id.eq.0) c(1,1,1) = (0.0,0.0)

      cc=c(1,:,:)

c for k1=0 plane, zero all modes with k2<0

c     do j = ly/2+2,ly 
c       c(1,j,:) = (0.0,0.0)
c     end do

c for the line k1=0 and k2=0,  zero all modes with k3<0
c     do k = 1,lz
c       k1=k+id*lz
c       if(k1.gt.nz/2.and.k1.le.nz) c(1,1,k)=(0.0,0.0)
c     end do

      return
      end

c ----------------------------------------------------------------------
      subroutine output(ux,uy,uz,idump,id,nallgrp,time,time0,nremesh)

      include "include.h"

      complex,dimension(lx1,ly,lz)::ux,uy,uz
      integer i1d,i2d,i3d,i4d,i5d,nremesh

      character*60 fnm1

      i1d=int(id/100)
      i2d=int((id-i1d*100)/10)
      i3d=id-i1d*100-i2d*10

      if(idump.lt.10) then
        fnm1='/scratch/nko/shear/vel'//char(idump+48)//'.'
     .        //char(i1d+48)//char(i2d+48)//char(i3d+48)
      else
        i4d=int(idump/10)
        i5d=mod(idump,10)
        fnm1='/scratch/nko/shear/vel'//char(i4d+48)
     .        //char(i5d+48)//'.'//char(i1d+48)
     .        //char(i2d+48)//char(i3d+48)
      endif

      open(10,file=fnm1,status='unknown',form='unformatted')
      write(10)time,time0,nremesh
      write(10)ux
      write(10)uy
      write(10)uz
      close(10)

      return
      end

c ---------------------------------------------------------------------
      subroutine input(ux,uy,uz,idp,id,nallgrp,time,time0,nremesh)

      include "include.h"

      complex,dimension(lx1,ly,lz)::ux,uy,uz
      integer i1d,i2d,i3d,i4d,i5d,nremesh

      character*60 fnm1

      i1d=int(id/100)
      i2d=int((id-i1d*100)/10)
      i3d=id-i1d*100-i2d*10
      if(idp.lt.10) then
        fnm1='/scratch/nko/shear/vel'//char(idp+48)//'.'
     .        //char(i1d+48)//char(i2d+48)//char(i3d+48)
      else
        i4d=int(idp/10)
        i5d=mod(idp,10)
        fnm1='/scratch/nko/shear/vel'//char(i4d+48)
     .        //char(i5d+48)//'.'//char(i1d+48)
     .        //char(i2d+48)//char(i3d+48)
      endif

      open(10,file=fnm1,status='unknown',form='unformatted')
      read(10)time,time0,nremesh
      read(10)ux
      read(10)uy
      read(10)uz
      close(10)

      return
      end

c -------------------------------------------------------------------
        subroutine transpose(ux,id,nallgrp)

        include 'include.h'

        complex,dimension(lx+1,ly,lz)::ux
        complex,dimension(lx+1,ly/nproc,lz):: tmp1,tmp2
        complex,dimension(lx+1,ly/nproc,lz*nproc)::tmp
        integer isize,nzm,nzp,status(MPI_STATUS_SIZE,2),req(2)

        isize=(lx+1)*ly*lz/nproc
        lly=ly/nproc
        do i=1,nproc-1
          nzp=mod(id+i,nproc)
          nzm=mod(id-i+nproc,nproc)
          js=nzp*lly
          do j=1,lly
            j1=js+j
            tmp1(:,j,:)=ux(:,j1,:)
          end do
          call mpi_isend(tmp1,isize,MPI_COMPLEX,nzp,i,nallgrp,req(1),
     .                   ierror)
          call mpi_irecv(tmp2,isize,MPI_COMPLEX,nzm,i,
     .               nallgrp,req(2),ierror)
          call MPi_Waitall(2,req,status,ierror)
          ks=nzm*lz
          do k=1,lz
            k1=ks+k
            tmp(:,:,k1)=tmp2(:,:,k)
          end do
        end do

         ks=id*lz
         js=id*lly
        do k=1,lz
          k1=ks+k
          do j=1,ly/nproc
            j1=js+j
            tmp(:,j,k1)=ux(:,j1,k)
          end do
        end do

        do k=1,lz
          do j=1,ly
             ux(:,j,k)=tmp(:,k,j)
           end do
        end do

        return
        end

c ------------------------------------------------------------------------
       subroutine mpifft(ux,isign,coeffz,coeffxy,id,nallgrp)

       include "include.h"

       complex,dimension(lx+1,ly,lz)::ux
       complex,dimension(lx+1,ly/nproc,lz):: tmp1,tmp2
       complex,dimension(lx+1,ly/nproc,lz*nproc)::tmp
       complex, dimension(nz+15)::coeffz
       real, dimension((lx2+15)+2*(ly+15))::coeffxy
       real alpha1,alpha2

       if(isign.eq.-1) then

C      Step 1, 2D fft, real-to-complex in x, complex-to-complex in y

        alpha2=1.0/lx2/ly
c       time1=mpi_wtime()
        do k=1,lz
          call scfft2du(isign,lx2,ly,ux(:,:,k),lx2+2,coeffxy)
        end do
        ux=ux*alpha2
c       time=mpi_wtime()-time1
c       if(id.eq.0)write(*,*)'time for 2-d fft',time

c      Step 2, Transpose between z and y

c       time2=mpi_wtime()
        call transpose(ux,id,nallgrp)
c       time=mpi_wtime()-time2
c       if(id.eq.0)write(*,*)'time for transpose',time

C      Step 3, 1d complex-to-complex in z

c       time3=mpi_wtime()
        do j=1,ly/nproc 
          do i=1,lx
            call cfft1d(isign,nz,ux(i,:,j),1,coeffz)
          end do
        end do
        ux=ux/nz
c       time=mpi_wtime()-time3
c       if(id.eq.0)write(*,*)'time for 1-d fft',time

       else if(isign.eq.1)then

C      Step 1, 1d complex-to-complex in z

c       time1=mpi_wtime()
        do j=1,ly/nproc
        do i=1,lx
            call cfft1d(isign,nz,ux(i,:,j),1,coeffz)
c           ux(i,:,j)=ux(i,:,j)/nz
          end do
        end do
c       time=mpi_wtime()-time1
c       if(id.eq.0)write(*,*)'time for 1-d fft',time

c      Step 2, Reverse transpose between z and y

c       time2=mpi_wtime()
        call transpose(ux,id,nallgrp)
c       time=mpi_wtime()-time2
c       if(id.eq.0)write(*,*)'time for transpose',time

C      Step 3, 2D fft, complex-to-complex in y, complex-to-real in x

        alpha2=1.0/lx2/ly
c       time3=mpi_wtime()
        do k=1,lz
          call csfft2du(isign,lx2,ly,ux(:,:,k),lx2+2,coeffxy)
c         ux(:,:,k)=ux(:,:,k)*alpha2
c         call sscal2d(lx,ly,alpha2,ux(:,:,k),lx+2)
        end do
c       time=mpi_wtime()-time3
c       if(id.eq.0)write(*,*)'time for 2-d fft',time

       end if

       return
       end




