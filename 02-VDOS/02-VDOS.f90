{\rtf1\ansi\ansicpg936\cocoartf1345\cocoasubrtf380
{\fonttbl\f0\froman\fcharset0 TimesNewRomanPSMT;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural

\f0\fs24 \cf0 \
!======================================\
! To calculate VDOS via FFT\
! By Yanlei Wang@THU,Mar 2013\
! Reference: Goncalves S, Bonado H, Vibrational densities of states from molecular-dynamics calculations, Physical Review B, 1992, 46:12019\
!======================================\
\
program vdos_w\
  use, intrinsic :: iso_c_binding\
  include 'fftw3.f03'\
\
!  IMPLICIT NONE\
  INTEGER :: nsteps,i\
  INTEGER, DIMENSION(:), ALLOCATABLE :: tag \
  DOUBLE PRECISION :: timestep,dt,sumcorr,distance\
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: time,u,vdos\
    DOUBLE COMPLEX,DIMENSION(:),ALLOCATABLE :: c\
  TYPE(C_PTR) :: plan_t,datain1,dataout1\
  REAL(C_DOUBLE),POINTER :: in1(:)\
  COMPLEX(C_DOUBLE_COMPLEX),POINTER :: out1(:)\
  DOUBLE COMPLEX :: sum\
  DOUBLE PRECISION :: sum2,x,y\
  DOUBLE PRECISION,DIMENSION(:,:), ALLOCATABLE :: energy\
  DOUBLE PRECISION :: pi,mb,time_p,pol_p\
  \
  OPEN (11,FILE='1nor_facf.data',STATUS='OLD\'92) !input, the output file of VACF.f90\
  OPEN (12,FILE='vdos.data',STATUS='UNKNOWN') !output, VDOS\
\
!======================================\
! Read The information of VACF\
!======================================\
 read(11,*) nsteps\
\
 ALLOCATE(time(nsteps),u(nsteps),c(nsteps),vdos(nsteps))\
 \
 do i=1,nsteps\
 read(11,*) time(i),u(i)\
 end do\
 \
 datain1=fftw_alloc_real(INT(nsteps,C_SIZE_T))\
 dataout1=fftw_alloc_complex(INT((nsteps/2+1),C_SIZE_T))\
 CALL c_f_pointer(datain1,in1,[nsteps])\
 CALL c_f_pointer(dataout1,out1,[nsteps/2+1])\
\
!======================================\
! FFT r2c \
!======================================\
 plan_t = fftw_plan_dft_r2c_1d(nsteps,in1,out1, FFTW_ESTIMATE)\
    do i=1,nsteps\
	in1(i)= u(i)\
	end do\
	call fftw_execute_dft_r2c(plan_t,in1,out1)\
    do i=1,(nsteps/2+1)\
	c(i)=out1(i)\
	end do\
  \
 do i=1,nsteps/2+1\
    x=real(c(i))\
    y=aimag(c(i))\
    vdos(i) = x*x+y*y\
 end do\
\
 time_p=1.d0/(nsteps*0.002)\
!======================================\
! Output\
!======================================\
 write (12,*) "frequency vdos"\
 do i=1,(nsteps/2+1)\
  write (12,'(F20.5,F20.5)') i*time_p,vdos(i)\
 end do\
\
end program vdos_w\
}