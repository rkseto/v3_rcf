c...parameters and common block for RQMD/S
      parameter (mx = 1500)
      parameter (rho0 = 0.168d0)
      parameter (pi=3.141592653589793d0)
      parameter (hc=0.19732705d0)
c     parameter (pmu1=2.35d0*0.19732705d0, pmu2=0.4d0*0.19732705d0)
cc    parameter (vex1=-0.25854d0, vex2=0.3756d0) ! original Marur. par.
c     parameter (vex1=-0.277d0,   vex2=0.663d0 ) ! isse tuned
      parameter (clam=2.d0) ! cut-off par. used in jamrqmd.f
c     parameter(g=1.66667d0,el=2.05d0, a=-21.21d-3,  bt=97.93d-3)  ! Maru Hard
c     parameter(g=1.33333d0,el=2.1d0,  a=-92.86d-3,  bt=169.28d-3) ! Maru Medium
c     parameter(g=1.16667d0,el=2.1d0,  a=-223.56d-3, bt=298.78d-3) ! Maru Soft
c     parameter(g=1.333d0,  el=1.75d0, a=-127.86d-3, bt=204.28d-3) ! Hirata
c     parameter(g=1.66667d0,el=2.05d0, a=-33.d-3,    bt=110.d-3)   ! New Hard
c     parameter(g=1.33333d0,el=2.1d0,  a=-116.d-3,   bt=193.d-3)   ! New Medium
c     parameter(g=1.16667d0,el=2.1d0,  a=-268.d-3,   bt=345.d-3)   ! New Soft
c     parameter(g=2.000d0,  el=1.08d0, a=-124.00d-3, bt= 70.50d-3) ! Aich Hard
c     parameter(g=1.166d0,  el=1.08d0, a=-356.00d-3, bt=303.00d-3) ! Aich Soft

      common /potpar/t1,t3,t3f,pmu1,pmu2,vex1,vex2
      save /potpar/
c     common /eul/r0(3,mx),p0(3,mx),fp(3,mx),fr(3,mx),
c    & fr0(3,mx),fp0(3,mx),fr1(3,mx),fp1(3,mx)
      common /eul/fp(3,mx),fr(3,mx)
      common /mat/rb(mx,mx),r22(mx,mx)
      common /ener/vsky(mx),vmom(mx)
      common /mom/pm(mx,mx),p22(mx,mx),pmom1(mx,mx)
      common /cut/rhocut(mx)
      save /mat/,/ener/,/mom/
      save /eul/
      save /cut/
