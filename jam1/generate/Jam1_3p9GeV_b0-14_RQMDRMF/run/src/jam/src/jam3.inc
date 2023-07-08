c...common block for buu 

c...Number of spatial grid.
c     parameter(maxx=20,maxz=24)
      parameter(maxx=1,maxz=1)

c...Force acting on hadrons.
      common/jamevf/force(3,mxv),forcer(3,mxv),vpot(mxv),dtfree(mxv)
     $ ,tforce

      parameter (mxw = 1500)
      common/jamevf2/rho(mxv),rhom(mxw,mxw),pmom2(mxw,mxw)

      ! for Coulomb potential
      common/jamevf3/rhoc(mxv),rhc(mxw,mxw)

      ! for symmetry energy
      common/jamevf4/rhos(mxv),rhs(mxw,mxw)

c...Baryon density.
      common/jambuu1/rhob(-maxx:maxx,-maxx:maxx,-maxz:maxz)

c...Charge density.
      common/jambuu2/rhoz(-maxx:maxx,-maxx:maxx,-maxz:maxz)

c....Coulomb potential.
      common/jambuu3/cop(-maxx:maxx,-maxx:maxx,-maxz:maxz)

c...4-momentum of the grid.
      common/jambuu4/avp(4,-maxx:maxx,-maxx:maxx,-maxz:maxz)


c....end of jam3.inc

