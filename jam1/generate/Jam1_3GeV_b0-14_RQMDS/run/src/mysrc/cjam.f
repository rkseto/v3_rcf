
      subroutine c_jammain
      include 'jam1.inc'
      include 'jam2.inc'
      include 'cjam.inc'
      character frame*8,proj*8,targ*8,cwin*15
      save frame, proj, targ, cwin, mevent, nstep, bmin, bmax, dt
      integer i, idx_a
      double precision :: nonpart, npart, fpx, fpy, fpz, fpt, fpt2,
     $ fpcm, fycm ,feta, fecm, fmass, fpcm2, fchge,
     $ frx, fry, frz, frt,
     $ frfx, frfy, frfz, frft
      integer, intent(in) :: c_i, iev, inv
      integer, intent(out) :: iout
      double precision, intent(out) :: dout
      double precision, intent(in) :: din
      integer, intent(inout) :: cmevent, cnstep, cnpart
      double precision, intent(inout) :: cbmin, cbmax, cdt, cmentry,
     $ cmxcoll
      character, dimension(*),intent(inout)::ccwin, cframe, cproj, ctarg
      double precision, dimension(*), intent(inout) :: cd_ary
      real, dimension(*), intent(inout) ::
     $ px_ary, py_ary, pz_ary, e_ary, m_ary, chge_ary, pt_ary,
     $ ycm_ary, eta_ary,
     $ x_ary, y_ary, z_ary, t_ary,
     $ fx_ary, fy_ary, fz_ary, ft_ary
      integer, dimension(*), intent(inout) :: ks_ary, pid_ary
      return

      entry c_t1(din)
      return

      entry c_para(cmevent, cbmin, cbmax, ccwin, cframe, cproj, ctarg,
     $ cdt, cnstep, cnv, cnbary, cnmeson, cmentry, cmxcoll)
      cmevent=mevent
      cbmin=bmin
      cbmax=bmax
      do i = 1, 8
        ccwin(i)=cwin(i:i)
        cframe(i)=frame(i:i)
        cproj(i)=proj(i:i)
        ctarg(i)=targ(i:i)
      end do
      do i = 9, 15
        ccwin(i)=cwin(i:i)
      end do
      cdt=dt
      cnstep=nstep
      cnv=nv
      cnbary=nbary
      cnmeson=nmeson
      cmentry=mentry
      cmxcoll=mxcoll
      return

      entry c_s_para(cmevent, cbmin, cbmax, ccwin, cframe, cproj, ctarg,
     $ cdt, cnstep)
      mevent=cmevent
      bmin=cbmin
      bmax=cbmax
      do i = 1, 8
        cwin(i:i)=ccwin(i)
        frame(i:i)=cframe(i)
        proj(i:i)=cproj(i)
        targ(i:i)=ctarg(i)
      end do
      do i = 9, 15
        cwin(i:i)=ccwin(i)
      end do
      dt=cdt
      nstep=cnstep
      return

      entry c_s_mstc(c_i, din)
      mstc(c_i) = din
      return

      entry c_s_parc(c_i, din)
      parc(c_i) = din
      return

      entry c_s_mstd(c_i, din)
      mstd(c_i) = din
      return

      entry c_s_pard(c_i, din)
      pard(c_i) = din
      return

      entry c_s_pare(c_i, din)
      pare(c_i) = din
      return

      entry c_s_mste(c_i, din)
      mste(c_i) = din
      return

      entry c_s_mdcy(c_i, inv, din)
      mdcy(c_i, inv) = din
      return

      entry c_s_mdme(c_i, inv, din)
      mdme(c_i, inv) = din
      return

      entry c_mstc(c_i, dout)
      dout = mstc(c_i)
      return

      entry c_parc(c_i, dout)
      dout = parc(c_i)
      return

      entry c_mstd(c_i, dout)
      dout = mstd(c_i)
      return

      entry c_pard(c_i, dout)
      dout = pard(c_i)
      return

      entry c_pare(c_i, dout)
      dout = pare(c_i)
      return

      entry c_mste(c_i, dout)
      dout = mste(c_i)
      return

      entry c_jaminit
      call jaminit(mevent,bmin,bmax,dt,nstep,frame,proj,targ,cwin)
      return

      entry c_mevent(iout)
      iout=mevent
      return

      entry c_jamevt(iev)
      !call cjamevt(iev)
      call jamevt(iev)
      return

      entry c_jamfin
      call jamfin
      return

      entry c_jam_test(c_i, iout, dout)
      iout=mstc(c_i)
      dout=win
      return

      entry c_p(c_i, inv, dout)
      dout = p(c_i, inv)
      return

      entry c_k(c_i, inv, dout)
      dout = k(c_i, inv)
      return

      entry c_r(c_i, inv, dout)
      dout = r(c_i, inv)
      return

      entry c_v(c_i, inv, dout)
      dout = v(c_i, inv)
      return

      entry c_wp(c_i, cd_ary)
      do i = 1, nv
        cd_ary(i) = p(c_i, i)
      end do
      return

      entry c_wk(c_i, cd_ary)
      do i = 1, nv
        cd_ary(i) = k(c_i, i)
      end do
      return

      entry c_wr(c_i, cd_ary)
      do i = 1, nv
        cd_ary(i) = r(c_i, i)
      end do
      return

      entry c_wv(c_i, cd_ary)
      do i = 1, nv
        cd_ary(i) = v(c_i, i)
      end do
      return

      entry c_wt(c_i, cd_ary)
      do i = 1, nv
        cd_ary(i) = v(c_i, i)
      end do
      return

      entry c_wff(c_i, cd_ary)
      do i = 1, nv
        cd_ary(i) = v(c_i, i)
      end do
      return

      entry c_cdata(cnpart, pid_ary, px_ary, py_ary, 
     $ pz_ary)
      nonpart=0
      do i = 1, nv
        pid_ary(i) = int(k(2, i))
        px_ary(i) = real(p(1, i))
        py_ary(i) = real(p(2, i))
        pz_ary(i) = real(p(3, i))
        if(abs(k(7,i)) .eq. 1) then
                nonpart=nonpart+1
        end if
      end do
      npart=mstd(11)-nonpart
      cnpart=int(npart)
      return 

      entry c_wdata(cnpart, pid_ary, px_ary, py_ary, 
     $ pz_ary, e_ary, m_ary, chge_ary)
      nonpart=0
      do i = 1, nv
        pid_ary(i) = k(2, i)
        px_ary(i) = real(p(1, i))
        py_ary(i) = real(p(2, i))
        pz_ary(i) = real(p(3, i))
        e_ary(i) = real(p(4, i))
        m_ary(i) = real(p(5, i))
        chge_ary(i) = real(jamchge(k(2, i)) / 3)
        if(abs(k(7,i)) .eq. 1) then
                nonpart=nonpart+1
        end if
      end do
      npart=mstd(11)-nonpart
      cnpart=int(npart)
      return 

      entry c_wedata(cnpart,
     $ pid_ary, px_ary, py_ary,
     $ pz_ary, e_ary, m_ary, chge_ary, pt_ary,
     $ ycm_ary, eta_ary)
      nonpart=0
      do i = 1, nv
        pid_ary(i) = k(2, i)
        fpx = p(1, i)
        px_ary(i) = real(fpx)
        fpy = p(2, i)
        py_ary(i) = real(fpy)
        fpz = p(3, i)
        pz_ary(i) = real(fpz)
        fecm = p(4, i)
        e_ary(i) = real(fecm)
        fmass = p(5, i)
        m_ary(i) = real(fmass)
        fchge = jamchge(k(2, i)) / 3
        chge_ary(i) = real(fchge)
        fpt2 = fpx*fpx + fpy*fpy
        fpt = sqrt(fpt2)
        fpcm2 = fpt2 + fpz*fpz
        fpcm = sqrt(fpcm2)
        fecm = sqrt(fpcm2 + fmass*fmass)
        fycm = 0.5 * log((fecm + fpz) / (fecm - fpz))
        if (fycm .ne. fycm) then 
          fycm = -999.0
        end if
        feta = 0.5 * log((fpcm + fpz) / (fpcm - fpz))
        if (feta .ne. feta) then 
          feta = -999.0
        end if
        pt_ary(i) = real(fpt)
        ycm_ary(i) = real(fycm)
        eta_ary(i) = real(feta)
        if(abs(k(7,i)) .eq. 1) then
                nonpart=nonpart+1
        end if
      end do
      npart=mstd(11)-nonpart
      cnpart=int(npart)
      return

      entry c_wxdata(cnpart,
     $ ks_ary, pid_ary, x_ary, y_ary, z_ary, t_ary, ft_ary,
     $ px_ary, py_ary, pz_ary,
     $ e_ary, m_ary, chge_ary, pt_ary,
     $ ycm_ary, eta_ary)
      nonpart=0
      do i = 1, nv
        ks_ary(i) = k(1, i)
        pid_ary(i) = k(2, i)
c       if (abs(k(7,i)) .eq. 1) print *, pid_ary(i), k(2, i)
c       print *, pid_ary(i)
        frx = r(1, i)
        x_ary(i) = real(frx)
        fry = r(2, i)
        y_ary(i) = real(fry)
        frz = r(3, i)
        z_ary(i) = real(frz)
        frt = r(4, i)
        t_ary(i) = real(frt)
        frft = r(5, i)
        ft_ary(i) = real(frft)
        fpx = p(1, i)
        px_ary(i) = real(fpx)
        fpy = p(2, i)
        py_ary(i) = real(fpy)
        fpz = p(3, i)
        pz_ary(i) = real(fpz)
        fecm = p(4, i)
        e_ary(i) = real(fecm)
        fmass = p(5, i)
        m_ary(i) = real(fmass)
        fchge = jamchge(k(2, i)) / 3
        chge_ary(i) = real(fchge)
        fpt2 = fpx*fpx + fpy*fpy
        fpt = sqrt(fpt2)
        fpcm2 = fpt2 + fpz*fpz
        fpcm = sqrt(fpcm2)
        fecm = sqrt(fpcm2 + fmass*fmass)
        fycm = 0.5 * log((fecm + fpz) / (fecm - fpz))
        if (fycm .ne. fycm) then
          fycm = -999.0
        end if
        feta = 0.5 * log((fpcm + fpz) / (fpcm - fpz))
        if (feta .ne. feta) then
          feta = -999.0
        end if
        pt_ary(i) = real(fpt)
        ycm_ary(i) = real(fycm)
        eta_ary(i) = real(feta)
        if(abs(k(7,i)) .eq. 1) then
                nonpart=nonpart+1
        end if
      end do
      npart=mstd(11)-nonpart
      cnpart=int(npart)
      return

      entry c_cldata(cnpart,
     $ ks_ary, pid_ary, x_ary, y_ary, z_ary, t_ary, 
     $ fx_ary, fy_ary, fz_ary, ft_ary, 
     $ px_ary, py_ary, pz_ary)
      nonpart=0
      do i = 1, nv
        ks_ary(i) = k(1, i)
        pid_ary(i) = k(2, i)
c       if (abs(k(7,i)) .eq. 1) print *, pid_ary(i), k(2, i)
c       print *, pid_ary(i)
        frx = r(1, i)
        x_ary(i) = real(frx)
        fry = r(2, i)
        y_ary(i) = real(fry)
        frz = r(3, i)
        z_ary(i) = real(frz)
        frt = r(4, i)
        t_ary(i) = real(frt)
        frfx = v(1, i)
        fx_ary(i) = real(frfx)
        frfy = v(2, i)
        fy_ary(i) = real(frfy)
        frfz = v(3, i)
        fz_ary(i) = real(frfz)
        frft = v(4, i)
        ft_ary(i) = real(frft)
        fpx = p(1, i)
        px_ary(i) = real(fpx)
        fpy = p(2, i)
        py_ary(i) = real(fpy)
        fpz = p(3, i)
        pz_ary(i) = real(fpz)
        if(abs(k(7,i)) .eq. 1) then
                nonpart=nonpart+1
        end if
      end do
      npart=mstd(11)-nonpart
      cnpart=int(npart)
      return

      entry c_kchg(c_i, inv, dout)
      dout = kchg(c_i, inv)
      return

      entry c_pmas(c_i, inv, dout)
      dout = pmas(c_i, inv)
      return

      entry c_parf(c_i, dout)
      dout = parf(c_i)
      return

      entry c_vckm(c_i, inv, dout)
      dout = vckm(c_i, inv)
      return

      entry c_mdcy(c_i, inv, dout)
      dout = mdcy(c_i, inv)
      return

      entry c_mdme(c_i, inv, dout)
      dout = mdme(c_i, inv)
      return

      entry c_brat(c_i, dout)
      dout = brat(c_i)
      return

      entry c_kfdp(c_i, inv, dout)
      dout = kfdp(c_i, inv)
      return

      entry c_kq(c_i, inv, dout)
      dout = kq(c_i, inv)
      return

      entry c_vq(c_i, inv, dout)
      dout = vq(c_i, inv)
      return

      entry c_nv(dout)
      dout = nv
      return

      entry c_nbary(dout)
      dout = nbary
      return

      entry c_nmeson(dout)
      dout = nmeson
      return

      entry c_jamcomp(c_i, iout)
      iout = jamcomp(c_i)
      return

      entry c_mxcoll(iout)
      iout = mxcoll
      return

      entry c_mentry(dout)
      dout = mentry
      return

      entry c_icoll(c_i, inv, dout)
      dout = icoll(c_i, inv)
      return

      entry c_coll(c_i, inv, dout)
      dout = coll(c_i, inv)
      return

      entry c_on_timeana
        cjam_record_time=1
      return

      entry c_off_timeana
        cjam_record_time=0
      return

      entry c_on_skiptime
        cjam_specified_time = 1
      return

      entry c_off_skiptime
        cjam_specified_time = 0
      return

      entry c_specified_time(iout)
        iout = cjam_specified_time
      return

      entry c_specify_time_series(ks_ary, inv)
        cjam_specified_time = 1
        cjam_timepoint_size = inv
        idx_a = inv
        do i = 1, idx_a
          cjam_tseries(i) = ks_ary(i)
        end do
      return

      entry c_quiet_file_out(inv)
        fname(inv)='/dev/null'
      return

      end ! subroutine cjam

      subroutine c_tphase(catime, cnpart, cmul, ks_ary, pid_ary,
     $rx_ary, ry_ary, rz_ary, rt_ary, rft_ary,
     $px_ary, py_ary, pz_ary, e_ary, m_ary, chge_ary,
     $pt_ary, ycm_ary, eta_ary)
c     copied from jamana.f --> jamout3

c...Purpose: to out put phase space data.  2011/12/9
      include 'jam1.inc'
      include 'jam2.inc'
      integer :: i, ri
      integer, intent(inout) :: cnpart, cmul
      double precision :: atime k1, kf, x, y, z, frt, frft
      double precision :: nonpart, npart, fpx, fpy, fpz, fpt, fpt2,
     $ fpcm, fycm ,feta, fecm, fmass, fpcm2, fchge, atime
      real, intent(inout) :: catime
      real, dimension(*), intent(inout) ::
     $ px_ary, py_ary, pz_ary, e_ary, m_ary, chge_ary, pt_ary,
     $ ycm_ary, eta_ary
      integer, dimension(*), intent(inout) :: ks_ary, pid_ary
      real, dimension(*), intent(inout) ::
     $ rx_ary, ry_ary, rz_ary, rt_ary, rft_ary
      data id/0/
c...ioot=0:
      data iopt/1/
      save id,iopt

      atime = pard(1)
      catime = real(atime)

      !write(iunit,'(''# time '',f10.4,'' fm/c'')')atime

c     etot=0
c     px=0
c     py=0
c     pz=0
c     icharge=0
c....Loop over all particles.
      nonpart=0
      ri=0 ! Real Mul to record
      do i=1,nv
        k1=k(1,i)
        kf=k(2,i)
        if(k1.gt.10) continue    ! dead particle
        if(iopt.eq.1.and. k1.le.0) continue   ! within a formation time(newly)

        dt=atime-r(4,i)
c... (some const. quarks) within a formation time
        if(iopt.eq.1.and.dt.lt.0.0d0) continue

c        etot = etot+p(4,i)
c        px = px+p(1,i)
c        py = py+p(2,i)
c        pz = pz+p(3,i)
c        icharge=icharge+jamchge(kf)

        ri=ri+1
        x=r(1,i)+dt*p(1,i)/p(4,i)
        y=r(2,i)+dt*p(2,i)/p(4,i)
        z=r(3,i)+dt*p(3,i)/p(4,i)
        frt=r(4,i)+dt
        frft=r(5,i)
        rx_ary(ri) = real(x)
        ry_ary(ri) = real(y)
        rz_ary(ri) = real(z)
        rt_ary(ri) = real(frt)
        rft_ary(ri) = real(frft)
c       write(iunit,820)kf,jamchge(kf),x,y,z,(p(j,i),j=1,5)
        !write(iunit,820)kf,k(7,i),x,y,z,(p(j,i),j=1,5)
        ks_ary(ri) = k(1, i)
        pid_ary(ri) = k(2, i)
c       if (abs(k(7,i)) .eq. 1) print *, pid_ary(ri), k(2, i)
        fpx = p(1, i)
        px_ary(ri) = real(fpx)
        fpy = p(2, i)
        py_ary(ri) = real(fpy)
        fpz = p(3, i)
        pz_ary(ri) = real(fpz)
        fecm = p(4, i)
        e_ary(ri) = real(fecm)
        fmass = p(5, i)
        m_ary(ri) = real(fmass)
        fchge = jamchge(k(2, i)) / 3
        chge_ary(ri) = real(fchge)
        fpt2 = fpx*fpx + fpy*fpy
        fpt = sqrt(fpt2)
        fpcm2 = fpt2 + fpz*fpz
        fpcm = sqrt(fpcm2)
        fecm = sqrt(fpcm2 + fmass*fmass)
        fycm = 0.5 * log((fecm + fpz) / (fecm - fpz))
        if (fycm .ne. fycm) then
          fycm = -999.0
        end if
        feta = 0.5 * log((fpcm + fpz) / (fpcm - fpz))
        if (feta .ne. feta) then
          feta = -999.0
        end if
        pt_ary(ri) = real(fpt)
        ycm_ary(ri) = real(fycm)
        eta_ary(ri) = real(feta)
        if(abs(k(7,i)) .eq. 1) then
                nonpart=nonpart+1
        end if
        npart=mstd(11)-nonpart
        cnpart=int(npart)
      end do
      cmul = ri

c     etot=etot/mstd(11)
c     px=px/mstd(11)
c     py=py/mstd(11)
c     pz=pz/mstd(11)
c     print *,'px=',px,'py=',py,'pz=',pz,'e=',etot,'ich=',icharge/3

      end subroutine ! tphase
