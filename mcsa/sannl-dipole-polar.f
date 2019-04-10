      program simulann
c use saparam3.prm 
c  revised: 2014-OCT-06 by mingjun
c  PEML based on the algorithm of Corona et al. ACM Trans Math 13(3). 262 (1987) 
      implicit none

      integer nmax, neps, ncons, ngrp
c  change the value of ncons, weight()!!!!
c  INPUT:
c  ncons, dipoleqm.dat, dipolemm.dat, saparams.prm (only for
c  charge parameters: the no.of grp type)
c  OUTPUT:
c  para.prm, optimized values
      parameter (nmax = 90, neps = 4, ncons = 1, ngrp=1)
      double precision weight(8), chggrp(nmax)

      double precision lowbnd(nmax), upbnd(nmax), x(nmax), xopt(nmax),
     $                 cvec(nmax), stplen(nmax),fstar(neps), xtry(nmax),
     $                 temp, eps, rt, fopt, x0(nmax)

      double precision value(nmax)
c !!!!!!!!!!!!!!!!!!!!!
c dxqm need to change
c dxmm need to change
c !!!!!!!!!!!!!!!!!!!!!
      double precision dxqm(ncons), dyqm(ncons), dzqm(ncons),dqm(ncons),
     $ alphqm1(ncons),alphqm2(ncons),alphqm3(ncons),alphqm4(ncons)
      double precision dxmm(ncons), dymm(ncons), dzmm(ncons),dmm(ncons),
     $ alphmm1(ncons),alphmm2(ncons),alphmm3(ncons),alphmm4(ncons)
     

      integer i, j, k, n, n1, n2, nacp(nmax), ns, nt, nfcnev, idum,
     $        maxevl, nacc, nobds
      integer grpinfo(nmax,nmax,2), grpinfo2(nmax,nmax,2),
     $        nparagrp(nmax), nparagrp2(nmax)
      integer  nchg, npolthol,
     $        igrpx(nmax), natmgrp(nmax), nn, nk

      character name(nmax)*6

      external fctn 

c  Open input file with number of variables, initial values, limits and step length
c  x: parameters for fitting; lowbnd & upbnd: low and up bound of the parameters; 
c  stplen: step length for the SA algorithm
c  nmax: maximum no. of parameters

      open (10, file='saparams-dipole-polar.prm')
c   README:
c   1. USE the same parameter for equivalent atoms!
c   2. write the charge parameters followed by alpha and Thole
c   3. two notices: (1) the no. of parameters for a group >= 2
c   4. limitation: the 1st parameter in each group MUST belong to
c      this group only
c   Input:
c   npara nchg npolthol
c   parameters(for CHARGE: only n-1 parameters are needed to change;
c              since the summation of n parameters is unchanged;)
c   lowbounds
c   upbounds
c   stepsize
c   reference parameters
c   temperature
c   weights: dx dy dz dtot
c   grp1  nparagrp1 
c   x1 n1 x2 n2 .. ..
c   grp2  nparagrp2
c   x1 n1 x2 n2 .. ..
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  grp1: group no.
c  nparagrp1: no. of parameters in this group
c  x1: the no. of parameter; n1: equivalent atoms for x1 in this group
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      read(10,*) n, nchg, npolthol
      read(10,*) (x(i),i=1,n) 
      read(10,*) (lowbnd(i),i=1,n)
      read(10,*) (upbnd(i),i=1,n)
      read(10,*) (stplen(i),i=1,n)
c      read(10,*) (x0(i),i=1,n)
      read(10,*) temp
      read(10,*) (weight(i), i=1, 8)
c      do i=1, nchg
c         do j=1, nmax
c             grpinfo(i, j, 1) = 1
c             grpinfo(i, j, 2) = 1
c         end do
c         nparagrp(i) = 1
c      end do
      do i=1, ngrp
          read(10,*) n1, nparagrp(i)
          read(10,*) (grpinfo(n1, j,1), grpinfo(n1, j,2), 
     $                j=1, nparagrp(i))
      end do
      close(10)
c     print out the charge information
      write(6,*) 'group info for each group:'
      do i=1, ngrp
         write(6,*) i, nparagrp(i)
         write(6,*) (grpinfo(i, j,1), grpinfo(i, j,2), j=1, nparagrp(i))
      end do
c     compute group charge; 
      do i=1, ngrp
         chggrp(i)=0.0 
         do j=1, nparagrp(i)
            n1=grpinfo(i,j,1)
            n2=grpinfo(i,j,2)
            chggrp(i)=chggrp(i)+ x(n1)*dble(n2)
         end do
         write(6,*) i, 'group charge:', chggrp(i)
      end do
c     assign grpinfo2(i,j,k)
      do i=1, nchg
         nparagrp2(i)=0
         do j=1, ngrp
            do k=1, nparagrp(j)
               if(grpinfo(j,k,1) .eq. i) then
                       nparagrp2(i) = nparagrp2(i)+1
                       n1=nparagrp2(i)
                       grpinfo2(i,n1,1)=j
                       grpinfo2(i,n1,2)=grpinfo(j,k,2)
               end if
            end do
         end do
      end do
c     print out the informaiton in grpinfo2
      write(6,*) 'group info for each atom:'
      do i=1, nchg
         write(6,*) i, nparagrp2(i), (grpinfo2(i, j,1),grpinfo2(i, j,2),
     $       j=1, nparagrp2(i))
      end do
c     read the qm dipoles
      open (10, file='dipoleqm.dat')
      do i=1, ncons

c !!!!!!!!!!!!!!!!!!!!!
c   read the qm dipoles
c   dxqm(i)....
c !!!!!!!!!!!!!!!!!!!!!
      read(10,*) dxqm(i), dyqm(i), dzqm(i), dqm(i), 
     $ alphqm1(i),alphqm2(i),alphqm3(i),alphqm4(i)
      end do
      close(10)

C  Define input parameters.
c  maxevl: maximum no. of evaluations
c  eps: convergence criterion
c  idum: random seed
c  rt: temp=temp*rt
c  adjust the step length after ns*n MC steps
c  nt: decrease temperature after nt*ns*n MC steps
c  neps: used for convergence check besides the stplen(i) < 1.D-4
      maxevl = 999999
      eps = 1.0d-3
      rt = .75
      idum = 1234
      ns = 100
      nt =  8
c      ns = 2
c      nt = 2

      do i = 1, n
         cvec(i) = 2.0d0
      end do

C  Set initial temperature                         

c      temp = 500.0d0

C  Print initial conditions
      write(6,*) 'no. of group types:', ngrp
      write(6,*) 'no. of conformations:', ncons
      write(6,*) 'no. of charge parameters:', nchg
      write(6,*) 'no. of alpha+Thole:', npolthol
      write(6,*) 'Weights:', weight(1), weight(2), weight(3), weight(4),
     $ weight(5), weight(6), weight(7),weight(8)
c      write(6,*) 'No. of equivalent atoms in one group for each para'
c      do i = 1, n
c         write(6,*) 'parameter: ', i, natmgrp(i) 
c      end do

      write(6,11) n, temp, rt, eps, ns, nt, neps, maxevl, 
     $              idum

  11  format(/,' number of parameters: ',i3,
     $       /,' temperature: ', f8.2, '   rt: ',f8.2, '   eps: ',e8.2,
     $       /,' ns: ',i3, '   nt: ',i2, '   neps: ',i2,
     $       /,' maxevl: ',i8, '   idum: ',i4)
  
      call wrtvec(x,n,'initial trial vector')
c      call wrtvec(x0,n,'reference vector')
c      call wrtvec(stplen,n,'initial step length')
c      call wrtvec(lowbnd,n,'lower bound')
c      call wrtvec(upbnd,n,'upper bound')
c      call wrtvec(cvec,n,'c vector')
c      stop              

C  call simulated annealing routine

      call work(n,x,rt,eps,ns,nt,neps,maxevl,lowbnd,upbnd,cvec,
     $        idum,temp,stplen,xopt,fopt,nacc,nfcnev,nobds,
     $        fstar,xtry,nacp,ncons,dxqm,dyqm,dzqm,dqm,
     $        alphqm1,alphqm2,alphqm3,alphqm4, weight,
     $        ngrp,nchg,npolthol,x0,grpinfo, grpinfo2, nparagrp,
     $        nparagrp2, chggrp)
     
      write(6,*) 'Final results after sa'
     
      call wrtvec(xopt,n,'solution')
      call wrtvec(stplen,n,'final step length')
     
      write(6,12) fopt, nfcnev, nacc, nobds, temp 
     

  12  format(/,' optimal function value: ',f20.6
     1 /,' number of function evaluations: ',i10,
     2 /,' number of accepted evaluations: ',i10,
     3 /,' number of out of bound evaluations: ',i10,
     4 /,' final temp: ', f20.6)
    
      end

C  main routine
      subroutine work(n,x,rt,eps,ns,nt,neps,maxevl,lowbnd,upbnd,cvec,
     $           idum,temp,stplen,xopt,fopt,nacc,nfcnev,nobds,
     $          fstar,xtry,nacp,ncons,dxqm,dyqm,dzqm,dqm,
     $           alphqm1,alphqm2,alphqm3,alphqm4, weight,
     $           ngrp,nchg,npolthol,x0,grpinfo, grpinfo2, nparagrp,
     $           nparagrp2, chggrp)
      implicit none
      integer nmax
      parameter (nmax = 90)

      double precision x(n), lowbnd(n), upbnd(n), cvec(n), stplen(n), 
     $                 fstar(n), xopt(n), xtry(n), temp, eps, rt, fopt
      integer nacp(n), n, ns, nt, neps, nacc, maxevl, 
     $        nobds, nfcnev, idum
      integer ngrp, nchg, npolthol
      double precision f, ftry, metrop, random, ratio
      integer nup, ndown, nrjct, nnew, lnobds, totalmov, 
     $        k, i, j, m, m0, m1, m2, ieps, nopt
      logical quit
      double precision expn, sx, na1, nb1
      integer ncons, n0, n1, n2, k0, k1, k2, kn, natmgrp(n), n3, n4
      integer kk, nkk, kk1, nkk1, nkk0, nkk2, kk0, kk2
c     double precision dxqm(ncons), dyqm(ncons), dzqm(ncons), dqm(ncons)
      double precision dxqm(ncons), dyqm(ncons), dzqm(ncons),dqm(ncons),
     $ alphqm1(ncons),alphqm2(ncons),alphqm3(ncons),alphqm4(ncons)

      double precision weight(8), x0(n), chggrp(ngrp)
      integer grpinfo(nmax,nmax,2), grpinfo2(nmax,nmax,2),
     $        nparagrp(n), nparagrp2(n)
      real ran1
      character*80 line
      character*40 filenam
      integer nbr

C initialize counters
c nopt: no. of minimized values found
c nacc: number of accepted evaluations
c nobds: number of out of bound evaluations
c nfcnev: number of function evaluations
      nopt = 0
      nacc = 0
      nobds = 0
      nfcnev = 0
      ieps = 1     ! inicialize IEPS: counter of neps 

      do i = 1, n
         xopt(i) = x(i)
         nacp(i) = 0
      end do  

      do i = 1, neps
         fstar(i) = 1.0d+30
      end do  

C  check temperatue and limits

      if (temp .le. 0.0) then
         write(6,*) 'initial temperature is negative'       
         stop
      end if

c      do i = 1, n
c         if ((x0(i) .gt. upbnd(i)) .or. (x0(i) .lt. lowbnd(i))) then
c            write(6,*)'reference values out of bounds ',
c     $                 i,x0(i),upbnd(i),lowbnd(i)
c            stop
c         end if
c      end do
      
      do i = 1, n
         if ((x(i) .gt. upbnd(i)) .or. (x(i) .lt. lowbnd(i))) then
            write(6,*)'initial values out of bounds ',
     $                 i,x(i),upbnd(i),lowbnd(i)
            stop
         end if
      end do

C  calculte the initial function

      call fctn(ncons, dxqm, dyqm, dzqm, dqm, 
     $     alphqm1,alphqm2,alphqm3,alphqm4, 
     $     weight, n,x,f)

      nfcnev = nfcnev + 1
      fopt = f
      fstar(1) = f

      write(6,'(//)')
      call wrtvec(x,n,'initial vector')
      write(6,'(a10,f14.6)') 'initial f: ',f

c      stop ! stop after function
100   continue

      nnew = 0     ! nnew: number of new minima
      nup = 0      ! ndown: number of higher function evaluations
      ndown = 0    ! nup: number of lower function evaluations
      nrjct = 0     ! nrjct: number of rejected evaluations
      lnobds = 0    ! lnobds: local counter of out-of-bound moves

      do m = 1, nt
         do j = 1, ns
            do k = 1, n

               do i = 1, n                ! copy curent vector to trial vector
                  xtry(i) = x(i)
               end do
c               write(6,*) 'ran1(idum) =', ran1(idum)
               kn=0
               do n0=1, ngrp
                  if(k .eq. grpinfo(n0, 1, 1)) kn=1
               end do
               if(kn .ne. 1) then !1st charge if
               xtry(k) = x(k) + (ran1(idum)*2.- 1.) * stplen(k)   ! update trial vector
c              treat the charge paramters
               if(k .le. nchg) then
                       do m0=1, nparagrp2(k)
                          m1=grpinfo2(k,m0,1) !which grp this parameter belongs to
                          n1=nparagrp(m1) !no. of parameters for this grp
c                          n0=1
c                         nkk: no. of equivalent atoms in this group;
c                          kk: position of this parameter in this group;
c                          do k0=1, nparagrp(m1)
c                             if(grpinfo(m1,k0,1) .eq. k) then
c                                     kk=k0
c                                     nkk=grpinfo(m1,k0,2)
c                             end if
c                          end do
c                          kk1=kk
c                          do while (kk1 .eq. kk)
c                               kk1=(n0+n1*ran1(idum))/1.
c                               write(6,*) 'CHECK1:', kk, kk1, n0, n1
c                          end do
c                          nkk1=grpinfo(m1,kk1,2)
c                          k1=grpinfo(m1,kk1,1)
c                          sx=dble(nkk)/dble(nkk1)
c
c                          xtry(k1)=x(k1)- (xtry(k)-x(k))*sx
c                          if ((xtry(k) .lt. lowbnd(k)) .or.
c     $                     (xtry(k) .gt. upbnd(k)) .or.
c     $                     (xtry(k1) .lt. lowbnd(k1)) .or. 
c     $                     (xtry(k1) .gt. upbnd(k1))) then
c                               lnobds = lnobds + 1
c                               nobds = nobds + 1
c                               write(6,*)'reset charge parameters',m1
c                               do k2=n0, n0+n1-1
c                                   xtry(k2)=x0(k2)
c                               end do
c                           end if
c  place restraint of total charge for this group
c                     round the charges to 3 digits
                           do k2=1, nchg
                            if(xtry(k2) .gt. 0) then
                                    n3=(xtry(k2)*1000+0.5)/1.
                            else if(xtry(k2) .lt. 0) then
                                    n3=(xtry(k2)*1000-0.5)/1.
                            else
                                    n3=xtry(k2)*1000
                            end if
                            xtry(k2)=dble(n3)/1000.
c                            write(6,*) 'TEST1:', k2, n3, xtry(k2) 
                           end do
c restraint the 1st charge parameter
                           na1=0.0
                           n0=grpinfo(m1,1,1)
                           nkk0=grpinfo(m1,1,2)
                           do k2=2, nparagrp(m1)
                              nkk2=grpinfo(m1,k2,2)
                              kk2=grpinfo(m1,k2,1)
                              na1=na1+xtry(kk2)*dble(nkk2)
                           end do
                           xtry(n0)=(chggrp(m1)-na1)/dble(nkk0)
                               write(6,*)'xtry(n0) charges', xtry(n0)
c check the boundary
                          if ((xtry(k) .lt. lowbnd(k)) .or.
     $                     (xtry(k) .gt. upbnd(k)) .or.
     $                     (xtry(n0) .lt. lowbnd(n0)) .or. 
     $                     (xtry(n0) .gt. upbnd(n0))) then
                               lnobds = lnobds + 1
                               nobds = nobds + 1
                               write(6,*)'reset charges', m1, k, n0
                               do k1=1, nchg
                                  xtry(k1)=x(k1)
                               end do
                               goto  102
                           end if
                       end do !m0

               end if

C  check new vector boundaries               
               if ((k .gt. nchg) .and. ((xtry(k) .lt. lowbnd(k)) .or. 
     $             (xtry(k) .gt. upbnd(k)))) then
                  lnobds = lnobds + 1
                  nobds = nobds + 1
                  write(6,'(//)')
c                  call wrtvec(x,n,'current vec')
c                  write(6,'(2x,a17,f14.6)') 'current function:',f 
c                 call wrtvec(xtry,n,'trial vec')
                  write(6,*)'Point out of bounds. Generate new point'
                  xtry(k) = lowbnd(k) + 
     $                         (upbnd(k) - lowbnd(k))*ran1(idum)
               end if

C  calculate the function 
               call fctn(ncons, dxqm, dyqm, dzqm, dqm,
     $         alphqm1,alphqm2,alphqm3,alphqm4, weight, 
     $                   n,xtry,ftry)

               nfcnev = nfcnev + 1
         
               write(6,'(//)')       
c               call wrtvec(x,n,'current vec')
c               write(6,'(2x,a17,f14.6)') 'current function:',f
               call wrtvec(xtry,N,'trial vec')
               write(6,'(2x,a13,f20.5)') 'new function:',ftry


C  check maximum number of function evaluations                    
               if (nfcnev .ge. maxevl) then
                  write(6,*) ' '
                  write(6,*) 'MAXEVL exceeded'
                  stop  
               end if

C  accept vector if function decreases

               if (ftry .le. f) then         ! accept if function value decreases: minimization
                  write(6,'('' point accepted'')')
                  do i = 1, n
                     x(i) = xtry(i)
                  end do
                  f = ftry
                  nacc = nacc + 1          ! increase nacc: number of accepted tries
                  nacp(k) = nacp(k) + 1    ! increase nacp(k): number of accepted tries for each variable. 
                  ndown = ndown + 1

C  If smaller than any other point, record as new optimum
                  if (ftry .lt. fopt) then

c                     call wrtfil(n,xtry,'optimum.str')

                     do i = 1, n
                        xopt(i) = xtry(i)
c                        x0(i) = xopt(i)
                     end do
                     fopt = ftry
                     nnew = nnew + 1
                     nopt = nopt + 1

                 write(6,'(a13,i5,1x,f12.3)') 'new OPTIMUM #',nnew,fopt 
                     open(11,file='num')
                     write(11,*) nopt
                     rewind(11)
                     read(11,'(a40)') line
c                     write(6,*) 'LINE =',line
                     close (11)
                     nbr = len_trim(line)
c                     write(6,*) 'nbr =', nbr
                     if (nopt .le. 9) then
                        filenam = 
     $                  'OPTIMA/optimum_'//line(nbr:nbr)//'.str'
                     else if ((nbr .ge. 10) .and. (nbr .le. 99)) then
                        filenam = 
     $                  'OPTIMA/optimum_'//line(nbr-1:nbr)//'.str'
                     else if ((nbr .ge. 100) .and. (nbr .le. 999)) then
                        filenam = 
     $                  'OPTIMA/optimum_'//line(nbr-2:nbr)//'.str'
                     else
                        filenam = 
     $                  'OPTIMA/optimum_'//line(nbr-3:nbr)//'.str'
                     end if
c                     write(6,*) filenam
                     call wrtfil(n,xopt,filenam)

                  end if

C  If the point is higher use the Metropolis criterium to decide on acceptance or rejection

               else

                  metrop = expn((f - ftry)/temp)
                  random = ran1(idum)
                  if (random .lt. metrop) then
                     write(6,*) '  point accepted despite higher'
                     do i = 1, n
                        x(i) = xtry(i)
                     end do
                     f = ftry
                     nacc = nacc + 1
                     nacp(k) = nacp(k) + 1
                     nup = nup + 1              ! increment accepted up moves
                  else
                     nrjct = nrjct + 1
                     write(6,*)'  rejected hgher point'
                  end if
               end if
               end if !1st charge if
102            continue
            end do   ! close h loop 
         end do      ! close m loop

C  adjust stplen to have half of the evaluations accepted

         do i = 1, n

            ratio = dfloat(nacp(i)) /dfloat(ns)
            if (ratio .gt. 0.6d0) then
               stplen(i) = 
     $         stplen(i)*(1.0d0 + cvec(i)*(ratio - 0.6d0)/0.4d0)
            else if (ratio .lt. 0.4d0) then
               stplen(i) = 
     $         stplen(i)/(1.0d0 + cvec(i)*((0.4d0 - ratio)/0.4d0))
            end if
            if (stplen(i) .gt. (upbnd(i)-lowbnd(i))) then     ! if interval is bigger than limits set interval to the limits
               stplen(i) = upbnd(i) - lowbnd(i)
            end if
         end do  

C  check convergence based on step length magnitude
         do i = 1, n
            if (stplen(i) .le. 5.0d-4) then
               write(6,*)'Convergence achieved based on step length'
               return
            end if
         end do

         write(6,*) 'results after adjustment of step length'  
         call wrtvec(stplen,n,'new step length')
         call wrtvec(xopt,n,'current optimum')
         call wrtvec(x,n,'current vec')
         write(6,*)' '

         do i = 1, n
            nacp(i) = 0
         end do  

      end do  !end do m = 1, nt

      totalmov = nup + ndown + nrjct
      write(6,*) 'intermediate results before temperature reduction'  
      write(6,'(a30,f12.5)') '  current temperature:        ', temp
      write(6,'(a30,f12.5)') '  current minimum:            ', fopt
      write(6,'(a30,i8)') '  total tries:                ', totalmov
      write(6,'(a30,i8)') '     downhill:                ', ndown
      write(6,'(a30,i8)') '     accepted uphill:         ', nup
      write(6,'(a30,i8)') '     rejected uphill:         ', nrjct
      write(6,'(a30,i8)') '  trials out of bounds:       ', lnobds
      write(6,'(a30,i8)') '  new minima at current temp: ', nnew
      call wrtvec(xopt,n,'current optimal vector')
      call wrtvec(stplen,n,'current step length')
      write(6,*) ' '      

C  Check termination criteria
      quit = .false.
      fstar(ieps) = f
      if ((fopt - fstar(1)) .le. eps) quit = .true.
      do i = 1, neps-1
         if (abs(fstar(i) - fstar(i+1)) .gt. eps) quit = .false.
      end do   

c Terminate SA if appropriate.
      if (quit) then
         do i = 1, n
            x(i) = xopt(i)
         end do   
         write(6,*) 'Converged!!!'
         return
      end if

C  If termination criteria is not met, prepare for another loop
      temp = rt*temp
      
      if (ieps .lt. neps) then   ! increment ieps if less than neps or reset to 1
         ieps = ieps + 1
      else
         ieps = 1
      end if

      f = fopt
      do i = 1, n
         x(i) = xopt(i)
      end do  

C  Loop again

      go to 100

      end

      function expn(f)
      implicit none

      double precision f, expn

      if (f .gt. 174.) then
         expn = 3.69d+75
      else if (f .lt. -180.) then
         expn = 0.0d0
      else
         expn = exp(f)
      end if

      return
      end

      function ran1(idum)
c numerical recipes
      implicit none
      integer idum, ia, im,iq, ir, ntab, ndiv
      real ran1, am, eps, rnmx
      parameter (ia=16807,im=2147483647, am=1.0/im, iq=127773, ir=2836,
     $           ntab=32, ndiv=1+(im-1)/ntab, eps=1.2e-7, rnmx=1.0-eps)
      integer j, k, iv(ntab),iy
      save iv,iy
      data iv /ntab*0/, iy /0/
      if (idum .le. 0 .or. iy .eq. 0) then
         idum = max(-idum,1)
         do j = ntab+8,1,-1
            k = idum/iq
            idum = ia*(idum-k*iq)-ir*k
            if (idum .lt. 0) idum = idum + im
            if (j .le. ntab) iv(j) = idum
         end do
      end if
      k = idum/iq
      idum = ia*(idum-k*iq)-ir*k
      j = 1 + iy/ndiv
      iy = iv(j)
      iv(j) = idum 
      ran1 = min(am*iy,rnmx)
      return
      end

      subroutine wrtvec(x,n,title)
      implicit none
      integer nelem
      parameter  (nelem=12)
      integer i, j, n, nline, lbegin, lend
      double precision x(n)
      character *(*) title

      write(6,'(2x,a)') title

      nline = int(n/nelem)     

      if (nline .ge. 1) then
         do i = 1, nline
            lbegin = nelem*(i - 1) + 1
            lend = lbegin + nelem - 1 
            write(6,100) (x(j), j = lbegin, lend)      
         end do  
         lbegin = lend + 1
         if (lbegin .le. n) then
            write(6,100) (x(j), j = lbegin, n)
         end if
      else
         lbegin = 1
         lend = n
         write(6,100) (x(j),j = lbegin, n)
      end if

  100 format( 12(f8.3,1X))

      return
      end

      subroutine wrtfil(n,x,name)
      implicit none
      integer n
      double precision x(n)
      character*(*) name

      include 'stream.txt'

      return
      end

C  Subroutine to calculate the current function

      subroutine fctn(ncons, dxqm, dyqm, dzqm, dqm, 
     $           alphqm1,alphqm2,alphqm3,alphqm4,
     $           weight,n,x,f)
      implicit none

      integer n, i, ncons

      double precision x(n), f,dx,dy,dz,d,s,w1, w2, w3, w4, w5, w6, w7, 
     $                 di1, di2, alph1, alph2, alph3, alph4
      double precision dxqm(ncons), dyqm(ncons), dzqm(ncons),dqm(ncons),
     $ alphqm1(ncons),alphqm2(ncons),alphqm3(ncons),alphqm4(ncons)
      double precision weight(10)

      character(len=80) :: task*80


      call wrtfil(n,x,'para.str')

c  generate the file dipolemm.dat
      task = './script-charmm-dip-alpha'
      call system(task)

c  compute RMSD btw QM and MM dipoles
      open(17,file='dipolemm.dat',status='old')
      f = 0.0
      do i = 1, ncons
         read(17,*) dx, dy, dz, d, alph1, alph2, alph3, alph4
         s=weight(1)*(dx-dxqm(i))**2 + 
     $     weight(2)*(dy-dyqm(i))**2 + 
     $     weight(3)*(dz-dzqm(i))**2 + 
     $     weight(4)*(d-dqm(i))**2  +
     $     weight(5)*(alph1-alphqm1(i))**2 +  
     $     weight(6)*(alph2-alphqm2(i))**2 + 
     $     weight(7)*(alph3-alphqm3(i))**2 + 
     $     weight(8)*(alph4-alphqm4(i))**2  
         f = f + s
      end do

      f = sqrt(f/ncons)
      write(6,*) 'RMSD:', f
      close(17)

c      stop

      return
      end

