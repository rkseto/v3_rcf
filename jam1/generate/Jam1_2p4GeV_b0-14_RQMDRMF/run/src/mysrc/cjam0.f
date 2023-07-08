
c***********************************************************************
      subroutine cjamevt(iev)

c...Purpose: to simulate complete one event.

      include 'jam1.inc'
      include 'jam2.inc'
      include 'cjam.inc'
c...Dirted by a student
      integer idx_a, kdt_assigned
c...For particle multiplicities.
      common/jamxml/mult(-500:500)
      save /jamxml/
c...For random seed
      common/jydatr/mrpy(6),rrpy(100)
      save /jydatr/
      common/rseed/iseed
      save /rseed/

c...Save current event number.
      mstd(21)=iev

      imev=0
 3000 continue
      imev=imev+1
      if(imev.gt.30) call jamerrm(30,0,'(jamevt:) infinit loop???')

c...Initialize JAM parameters, etc for each event.
      call jaminie(mrun)

cTABARA
c     call ttinitchk
      if(mrun.eq.1) return

c...Save current random seed.
      mstd(22)=mrpy(1)

c...Global time.
      pard(1)=0.0d0

c...Some initialization for summary.
      if(mstc(161).ge.1) call jamana(1)

c...Rest collision counters.
      mstd(29)=0   ! collision counter
      mstd(30)=0   ! number of dead particle
      mstd(41)=0   ! elastic
      mstd(42)=0   ! inelastic
      mstd(43)=0   ! absorption
      mstd(44)=0   ! BB collision
      mstd(45)=0   ! MB collision
      mstd(46)=0   ! MM collision
      mstd(47)=0   ! antiBB collision
      mstd(48)=0   ! parton-hadron collision
      mstd(49)=0   ! parton-parton collision
      mstd(50)=0   ! decay
      mstd(51)=0   ! Pauli block
      mstd(52)=0   ! low energy cut
      mstd(53)=0   ! decay after simul.
      mstd(54)=0   ! Econ block
      mstd(55)=0   ! Number of hard scatt.

      mste(40)=0


c...Option: output phase space data.
      if(mstc(164).eq.1) call jamout1(pard(1))
      if(mstc(164).eq.2) call jamout2(pard(1))

c...Reset time dependent analysis.
      if(mstc(3).eq.1) then
c...Option: output particles on the terminal.
c      if(mstc(8).ne.0) call jamdisp(6,1)
        npri=1
        nann=1
        if(mstc(161).ge.1) call jamanat(0)
      else
        npri=max(1,int(parc(8)/parc(2)))
        nann=max(1,int(parc(7)/parc(2)))
      endif

c...Start time evolution.
c======================================================================*
      kdt=1
      do 10000 while (kdt.le.mstc(3))
c======================================================================*

        mstd(23)=kdt

c...Cascade.
        if(mstc(6).ge.0) then

          if(mstc(3).ge.2.and.mod(kdt-1,npri).eq.0) then
            if(mstc(8).ge.1) then
              call jamdisp(6,kdt)
c             write(*,*) 'jam(2)'
            endif
c           if(mstc(161).ge.1) call jamana(2)
          endif
crqmd
          if(mstc(3).ge.2) then
            if(mod(kdt-1,nann).eq.0) then
              if(mstc(161).ge.1) call jamana(2)
            endif
c....I am sorry to dirty this code for ugly using
            if(cjam_record_time .eq. 1) call c_jam_time_record(pard(1))
          endif

          call jamcoll

c...LPC:Collision according to mean free path.
        else if(mstc(6).le.-100) then
          call jamfpath

c...Multipule AA collision by Glauber geometry.
        else
          call jamglaub(igl)
          if(igl.eq.1) goto 3000

c...Option for final state hadron cascade.
          if(mod(abs(mstc(6))/10,10).eq.2) then
            mstc55=mstc(55)
            mstc(55)=1        ! resonance decay forbidden
            call jamfdec      ! fragment strings
            call jamcoll      ! hadron cascade
            mstc(55)=mstc55
          endif
        endif

        if (cjam_specified_time .eq. 1) then
          kdt_assigned = 0
          do idx_a = 1, cjam_timepoint_size
            if (kdt .lt. cjam_tseries(idx_a)) then
              kdt = cjam_tseries(idx_a)
              kdt_assigned = 1
              exit
            end if
          end do
          if (kdt_assigned .eq. 0) then
            if (kdt .lt. mstc(3)) then
              kdt = mstc(3)
            else
              kdt = mstc(3) + 1
            end if
          end if
        else !cjam_specified_time .ne. 1
          kdt = kdt + 1
        end if

c======================================================================*
10000 continue  !******************** End time evolution
c======================================================================*

c...1+1 simulation: only inelastic event is allowed.
      if(mstc(17).eq.1.and.mstd(42).eq.0) then
        call jamerrm(1,0,' Elastic event in 1+1 simulation')
        goto 3000
      endif

c...Transport particles that are still within a formation time.
      do ip=1,nv
        call jamtrspt(ip,r(5,ip))
      end do

c...Final decay of resonances delta, n*, rho ....
      if(mstc(41).eq.1) call jamfdec
      if(mstc(8).ne.0) then
        call jamdisp(6,100)
        write(*,*) 'jam(1)'
      endif

c...2014/1/27 add mstc(18)
c     if(mstc(18).eq.0.and.mstd(30).gt.0) call jamedit
      if(mstc(18).eq.0.and.mstd(30).gt.0) then
      call jamedit
      write(mstc(38),*)'jamedit? mstc18 mstd30',mstc(18),mstd(30)
      endif

c----------------------------------------------------------------------*
c.....Summary and final output for each run
c----------------------------------------------------------------------*
      pard(1)=pard(1)+parc(2)
      if(mstc(8).ge.2.and.mstc(6).ge.-100)
     $  call jamcheck('After Simul. i.e. Final Check')

c...Output analysis results.
      if(mstc(161).ge.1) call jamana(3)

c...Count hadron multiplicities.
      do i=1,nv
        if(k(1,i).gt.10) cycle
        kf=k(2,i)
        kc=jamcomp(kf)
        if(kc.le.0) cycle
        kch=jamchge(kf)
        mult(0)=mult(0)+1
        if(kch.ne.0) mult(41)=mult(41)+1
        if(kch.lt.0) mult(42)=mult(42)+1
        if(kch.gt.0) mult(43)=mult(43)+1
        mult(kc*isign(1,kf))=mult(kc*isign(1,kf))+1
      end do

c...Option for deuteron coalescence.
      if(mstc(45).ge.1) call jamdeut

c...Update collision counters.
      pard(71)=pard(71)+mstd(41)   ! elastic
      pard(72)=pard(72)+mstd(42)   ! inelastic
      pard(73)=pard(73)+mstd(43)   ! absorption
      pard(78)=pard(78)+mstd(50)   ! decay
      pard(79)=pard(79)+mstd(51)   ! Pauli block
      pard(80)=pard(80)+mstd(52)   ! low energy cut
      pard(81)=pard(81)+mstd(53)   ! decay after simul.
      pard(82)=pard(82)+nv
      pard(83)=pard(83)+nmeson
      pard(84)=pard(84)+mstd(54)   ! Econ block
      pard(85)=pard(85)+mstd(48)   ! h-p coll.
      pard(86)=pard(86)+mstd(49)   ! p-p coll.
      pard(87)=pard(87)+mstd(55)   ! hard scatt.
      pard(74)=pard(74)+mstd(44)   ! BB collision
      pard(75)=pard(75)+mstd(45)   ! MB collision
      pard(76)=pard(76)+mstd(46)   ! MM collision
      pard(77)=pard(77)+mstd(47)   ! antiBB collision

      if(mstc(8).ge.3) then
        ih=mstc(38)
        write(ih,*)'Event=           ',mstd(21)
        write(ih,*)'Total collisions ',(mstd(41)+mstd(42))/mstc(5)
        write(ih,*)'inelastic        ',mstd(42)/mstc(5)
        write(ih,*)'elastic          ',mstd(41)/mstc(5)
        write(ih,*)'absorption       ',mstd(43)/mstc(5)
        write(ih,*)'decays           ',mstd(50)/mstc(5)
        write(ih,*)'h-p coll.        ',mstd(48)/mstc(5)
        write(ih,*)'p-p coll.        ',mstd(49)/mstc(5)
        write(ih,*)'Pauli blocks     ',mstd(51)/mstc(5)
        write(ih,*)'low energy cuts  ',mstd(52)/mstc(5)
        write(ih,*)'Econ block       ',mstd(54)/mstc(5)
        write(ih,*)'decay after simul.',mstd(53)/mstc(5)
        write(ih,*)'Number of hard scatt.',mstd(55)/mstc(5)
        write(ih,*)'Total particle   ',nv/mstc(5)
        write(ih,*)'Total mesons     ',nmeson/mstc(5)
      endif

      if(nv.ne.nbary+nmeson) then
        write(check(1),'(''nv='',i9)')nv
        write(check(2),'(''nbary='',i9)')nbary
        write(check(3),'(''nmeson='',i9)')nmeson
        call jamerrm(11,3,'nv .ne.nbary+nmeson')
      endif
 
      end ! cjamevt
