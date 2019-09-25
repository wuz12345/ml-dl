      MODULE setstuff_mod
      USE precision
      USE commons2
      USE clips_library
      USE angkin2_mod
      USE angles_mod
      USE getvj_mod
      USE gridx_mod
      USE potent_mod
      USE jacobi_mod
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: setstuff

      CONTAINS

      SUBROUTINE setstuff(rbc,psibc,vbc,proc,rmin,rmax,srmin,srmax,&
     &                     jbstart,jbend)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reads in various details and sets up problem
! Some further explanations :
! ------------------------------------------------------------
!  brab,srab,cbrab,csrab   :
! ------------------------------------------------------------
! absorption parameters such that absorption occurs for
! R >= brab with strength exp(-cbrab*(R-brab)**2), and
! a similar absorption for r paramterized by srab,csrab
! ------------------------------------------------------------
! myjbig,myip, and jsym pertain to angular momentum and symmetry: 
! ------------------------------------------------------------
! myjbig= value of total a.m. quantum number J on this processor
! myip=  Parity of (-1)**J. If ip=0,then parity is
!          such that K=0 is the lowest possible projection
!          quantum number of the total angular momentum;  if
!          ip=1 then parity is such that K=1 is the lowest
!          projection quantum number.  (For J=0, ip can
!          ONLY be 0)
!
! jsym= 0,1 or 2 depending on whether or not the BC rotational
!          states in ABC are  mixed even and odd (jsym=0), only
!          odd (jsym=1) or only even (jsym=2).  jsym=1 and jsym=2
!          are appropriate only if B=C, i.e. if there is exchange
!          symmetry 
!
! kqm(i) is the k quantum number of the i'th basis function
!        kqm(1) is sometimes 0 and sometimes 1 depending on the parity.
! 
! inj(k) is the index for the first allowed j state for a
! given k (j>=k is the rule).  
! 
! kstart is the initial value of the K or helicity quantum number. It 
! must be between 0 and min(jbig,jstart). For some parities it starts
! at 1 not zero.  i.e. if (-1)**jbig is even even parity has a minimum 
! value of K=0 and odd parity has a minimum value of K=1.
! ------------------------------------------------------------
! nj, nr and nv determine the vibrational basis set : 
! ------------------------------------------------------------
!    nj=no. of j states associated with BC 
!    nr=number of evenly spaced grid pts associated with R
!    nv=number of evenly spaced grid pts associated with r
!    (v is old notation from when a vib. basis set was used)
! ------------------------------------------------------------
! rmin,rmax ; srmin,srmax : evenly spaced grid limits, with
! grids being defined by nx points in (xmin,xmax) -- neither
! xmin nor xmax is an explicit grid point, 
! x=R or r 
! ------------------------------------------------------------
!  vcut and cencut : potential and centrifugal cut-offs in au 
!  see subroutine potgrid for their use. 
! ------------------------------------------------------------
!  ncost=number of quadrature points for cos(gamma) quadratures 
! ---------------------------------------------------------------
!  rma,rmb,rmc=masses in amu 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      USE jacobi_mod 
      IMPLICIT NONE
      INTEGER, PARAMETER :: lft_FORLEN=16
      INTEGER, PARAMETER :: jmx=2*njmx
      REAL(REAL8), DIMENSION(:), INTENT(OUT) :: rbc,psibc,vbc
      INTEGER, INTENT(OUT) :: proc
      INTEGER :: k,j,kk,kkk,jj,nsize,ir,isr,iv,mjabmax,nk
      INTEGER :: jbstart,jbend,i,j1,j2,ipend,ic
      INTEGER :: mvabmax,mdimt,iab,ibc,ier,myk
      REAL(REAL8) :: amu,vcut
      REAL(REAL8), INTENT(OUT) :: rmin,rmax,srmin,srmax
      REAL(REAL8) :: brab,cbrab,srab,csrab,rbcmin,rbcmax,vmin
      REAL(REAL8) :: vv,ediss,rj2,xx,xx2,vsmall,rmbc,vxadd
      REAL(REAL8) :: vyadd,tkin,tmax,dele,rmuabc
      REAL(REAL8) :: blen,bexp,slen,sexp
      REAL(REAL8), DIMENSION(3) :: rij
      REAL(REAL8), DIMENSION(nxmx) :: evec
      DATA amu/1822.88734d0/
      REAL(REAL8), ALLOCATABLE, DIMENSION(:) :: v
      REAL(REAL8), ALLOCATABLE, DIMENSION(:) :: veff
      REAL(REAL8), ALLOCATABLE, DIMENSION(:,:) :: psimat
      REAL(REAL8), ALLOCATABLE, DIMENSION(:,:) :: vic
      REAL(REAL8), ALLOCATABLE, DIMENSION(:,:) :: vk
      REAL(REAL8), ALLOCATABLE, DIMENSION(:,:,:) :: pjk
      REAL(REAL8), ALLOCATABLE, DIMENSION(:) :: fleg
!      lft_FOR_7777=
!        dimension wcost(ncosmx)
! grid representation of centrifugal repulsion:
! amu --> au : 
! Basic parameters that define the problem.   See explanations
! that follow. 
!
!.....set default values for range of "ip" or "parity"
	ipend=1
!
      WRITE(6,*) ' Jacobi product coordinates propagation '
! ***** ALL reads for program : *******************************
      READ(5,*) rma,rmb,rmc
      READ(5,*) rdagger
      READ(5,*) emin,emax
      READ(5,*) nj,ncost,jsym
      READ(5,*) rmin,rmax,nr
      READ(5,*) srmin,srmax,nv
! read in a parameter to decide whether explicit R/r dvr matrices
! are used for ke operator (0) or whether sin fft is used (1). 
      READ(5,*) ikin
      READ(5,*) vcut
      READ(5,*) brab,cbrab
!      IF(brab>=rmax)CALL clips_error(1,'SETSTUFF',' brab.ge.rmax')
      READ(5,*) srab,csrab
!      IF(srab>=srmax)CALL clips_error(1,'SETSTUFF',' srab.ge.srmax')
      READ(5,*) nstep
      READ(5,*) etran
      READ(5,*) br0,alpha,bbeta
      READ(5,*) ivstart,jstart,kstartbc
	if(jstart.eq.0)ipend=0
!.....test kstartbc
!	IF(kstartbc.gt.jstart)CALL clips_error(1,'SETSTUFF',' kstartbc.gt.jstart')
      READ(5,*) nvab,njab
!      IF(njab>nj)CALL clips_error(1,'SETSTUFF',' njab gt nj')
!      IF(nvab>nv)CALL clips_error(1,'SETSTUFF',' nvab gt nv')
      READ(5,*)jbstart,jbend
      write(6,*)' jbstart,jbend=',jbstart,jbend
      call flush(6)
!.....test kstartbc
!	IF(kstartbc.gt.jbstart)CALL clips_error(1,'SETSTUFF',' kstartbc.gt.jbstart')
! write out the basic information : 
      WRITE(6,*) ' rdagger : ',rdagger
      WRITE(6,*) ' emin, emax : ',emin,emax
      WRITE(6,*) ' no. j states ',nj,' ncost=',ncost
      WRITE(6,*) ' R grid ',rmin,rmax,nr
      WRITE(6,*) ' r grid ',srmin,srmax,nv
      WRITE(6,*) ' ikin=',ikin
      WRITE(6,*) ' pot/cent cutoff ',vcut
      WRITE(6,*) ' abs in R : ',brab,cbrab
      WRITE(6,*) ' abs in r : ',srab,csrab
      WRITE(6,*) ' no. Cheb. steps ',nstep
      WRITE(6,*) ' initial tran E/eV ',etran
      WRITE(6,*) ' initial Gaussian R0, width, smoothing ',br0, &
     &              alpha,bbeta
! adjust jstart for biggest K
!      jstart=jbend
      WRITE(6,*) ' v0,j0,k0: ',ivstart,jstart,kstartbc
      WRITE(6,*) ' no. ab vib states examined ',nvab
      WRITE(6,*) ' no. ab rot states examined ',njab
! emin/emax are for Ch scaling: 
! which has coefficient data for prod analysis written to it
! Other parameters that are of relevance and could be changed
! depending on system and information desired:
      ip=0
!ggbk        jsym=0
!
! Product channel c :            AB + C 
! Jacobi coords AB,C system <--> OD + H  products so
!          rma=15.9949d0
!          rmb=2.0141d0
!          rmc=1.00783d0
! 
! masses now read in -- see second read(5,*) above
!
! need to also set the minimum/maximum values of rbc, reactant distance,
! and the number of points, nbc, to be used in defining the initial
! state.  Nb grid is defined later on in this subroutine
      rbcmin=0.02d0
      rbcmax=8.0d0
!ggbk      nbc=200
!
      IF(nbc>nxmx)THEN
         WRITE(6,*) 'increase nxmx to :',nbc
!         CALL clips_error(1,'SETSTUFF','nxmx too small')
      END IF
! also need set range/no. pts. associated with reactant channel R
! variable, Ra. 
      ramin=0.0d0
      ramax=20.0d0
      nra=500
! N.B. neither the rbc or ra type variables above are important
! for size/computation time.  They are involved in initial condition
! routine only. 
      WRITE(6,*) ' '
      WRITE(6,*) ' mases (amu) : '
      WRITE(6,*) ' rma ',rma
      WRITE(6,*) ' rmb ',rmb
      WRITE(6,*) ' rmc ',rmc
      WRITE(6,*) ' '
      WRITE(6,*) ' potential, centrifugal cut-offs (a.u.) : '
      WRITE(6,*) vcut
      WRITE(6,*) ' '
! some conversions to au
      rma=rma*amu; rmb=rmb*amu; rmc=rmc*amu
! Jacobi coordinate masses for the product channel c:
      rm=rma*rmb/(rma+rmb); rmu=rmc*(rma+rmb)/(rmc+rma+rmb)
! some mass ratios used in conversions between Jacobi coords.
      aoab=rma/(rma+rmb); boab=rmb/(rma+rmb)
      bobc=rmb/(rmb+rmc); cobc=rmc/(rmb+rmc)
      rmabc=rma+rmb+rmc
! **************  end of problem details specification ***************
! check to be sure dimensions are correct:
! (need change parameters.inc if not correct)
!      IF(nr>nrmx)CALL clips_error(1,'SETEXP','increase nrmx')
!      IF(nv>nvmx)CALL clips_error(1,'SETEXP','increase nvmx')
!      IF(nj>njmx)CALL clips_error(1,'SETEXP','increase njmx')
!      IF(nvmx>nxmx)CALL clips_error(1,'SETEXP','increase nxmx')
!      IF(jsym<0.OR.jsym>2)CALL clips_error(2,'SETEXP','jsym is wrong')
      CALL init_parallel(jbstart,jbend,proc,ipend)
      CALL alloc_all()
! K values truncated because of a.m. coupling
      IF(myip==0)THEN
! if myip=0, then p=+ (even parity) :
        DO k=1,MIN(myjbig+1,njmx)
          kqm(k)=k-1
        END DO
      ELSE
! odd parity
        DO k=1,MIN(myjbig,njmx-1)
          kqm(k)=k
        END DO
      END IF
! BC j states are determined :
      IF(jsym==0)THEN
        DO j=1,nj
          jqm(j)=j-1
        END DO
      ELSE IF(jsym==1)THEN
        DO j=1,nj
          jqm(j)=2*(j-1)+1
        END DO
      ELSE IF(jsym==2)THEN
        DO j=1,nj
          jqm(j)=2*(j-1)
        END DO
      END IF
! Output the input information :
      WRITE(6,*) ' calculation for total angular momentum ',myjbig
      WRITE(6,*) ' and parity such that k = ',kqm(1),' is lowest value '
      WRITE(6,*) ' '
      WRITE(6,*) ' number of k states ',mynk,' and the k states :'
      WRITE(6,'(1x,10i4)') kqm(:)
      WRITE(6,*) ' '
      WRITE(6,*) ' number of j states ',nj,' and the j states: '
      WRITE(6,'(1x,10i4)') (jqm(k),k=1,nj)
      WRITE(6,*) ' '
      WRITE(6,*) ' number of quadrature points for gamma',ncost
      WRITE(6,*) ' '
      WRITE(6,*) ' R dvr basis size ',nr
      WRITE(6,*) ' rmin, rmax ',rmin,rmax
      WRITE(6,*) ' '
      WRITE(6,*) ' r dvr basis size ',nv
      WRITE(6,*) ' srmin, srmax ',srmin,srmax
      WRITE(6,*) ' '
      WRITE(6,*) '  nstep ',nstep
      WRITE(6,*) ' '
! set up quantum number index arrays for vectors. 
! A vector is Sum_k Sum_j=inj(k) Sum_v Sum_i C_{jkvi} |jk>|v>|i>
! where inj(k) is the index for the first allowed j state for a
! given k (j>=k is the rule). (v,i denote BC vibrations and R grid) 
      nk=myjbig+1-myip
      DO k=1,nk
        kk=kqm(k)
        DO j=1,nj
          jj=jqm(j)
          IF(jj>=kk)EXIT
        END DO
        inj(k)=j
      END DO
! actual size of the problem :
      nsize=mynk*nv*nr*ncost
      WRITE(6,*) ' '
      WRITE(6,*) ' ************************************** '
      WRITE(6,*) ' size of problem = ',nsize
      WRITE(6,*) ' ************************************** '
      WRITE(6,*) ' '
! get grids and the potential matrix V(R,r,j,jp,K) 
      ALLOCATE(pjk(0:jmx,ncosmx,0:nkmx+1), STAT=ier)
      write(6,*)''
      write(6,*)rmin,rmax,srmin,srmax,vcut
      write(6,*)'going into potmat'
      CALL potmat(rmin,rmax,srmin,srmax,vcut,vmin,pjk)
! determine absorption details passed to common/absetc/ :
! determine absorption parameters:
! should really be in SETEXP ...
      iabbr=nr
      blen=wkbr(nr)-brab
      bexp=2.0d0*blen
      DO ir=nr,1,-1
        fabbr(ir)=1.0d0
! Valentina noticed a possible divide by zero
        IF(wkbr(ir)>brab)THEN
!quad            fabbr(ir)=EXP(-cbrab*(wkbr(ir)-brab)**2)
          fabbr(ir)=EXP(-cbrab*EXP(-bexp/(wkbr(ir)-brab)))
          iabbr=ir
        END IF
      END DO
      iabsr=nv
      slen=wksr(nv)-srab
      sexp=2.0d0*slen
      DO isr=nv,1,-1
        fabsr(isr)=1.0d0
! Valentina noticed a possible divide by zero
        IF(wksr(isr)>srab)THEN
!quad            fabsr(isr)=EXP(-csrab*(wksr(isr)-srab)**2)
          fabsr(isr)=EXP(-csrab*EXP(-sexp/(wksr(isr)-srab)))
          iabsr=isr	 
        END IF
      END DO
      WRITE(6,*) ' start grid for R abs ',iabbr
      WRITE(6,*) ' start grid for r abs ',iabsr
      WRITE(6,*) ' at end of R grid, abs fac is ',fabbr(nr)
      WRITE(6,*) ' at end of r grid, abs fac is ',fabsr(nv)
      WRITE(6,*) ' '
      ALLOCATE(v(nv), STAT=ier)
! Set up single reactant vib. state and the various
! asymptotic product states for determination of 
! reaction probabilities:
! First, let's form all the relevant psiab (product)
! vib-rot states and energies:
!cmh calculate the wavefunctions far out and at the analysis line
!cmh take the energies from the far out wavefunctions
!mh calculate the shift in the analysis program
!cmh....calculate psiab far out

! form bare v=vab :
      DO iv=1,nv
! rab =
         rij(1)=wksr(iv)
! rbc =
         rij(2)=50.0d0
! rac =
         rij(3)=rij(1)+rij(2)
! 
         CALL surfac(vv,rij)
         v(iv)=vv
      END DO
! save dissociation energy
      ediss=v(nv)
      WRITE(6,*) ' dissociation energy of product AB=',ediss
! get the vib states for the given j state:
      WRITE(6,*) ' '
!      WRITE(6,*)'do ab far out '
      CALL flush(6)
      mjabmax=0; mvabmax=0
      ALLOCATE(psimat(nv,nv), STAT=ier)
      ALLOCATE(veff(nv), STAT=ier)
! now allocate evjab - gets deallocated in main.f90
      CALL alloc_evjab(nvab,njab)
      CALL alloc_devj(nvab,njab)
      DO j=1,njab
! add centrifugal term to v
         jj=jqm(j)
         rj2=jj*(jj+1)
         DO iv=1,nv
            xx=srmin+iv*dsr; xx2=xx*xx
            veff(iv)=v(iv)+rj2/(2.0d0*rm*xx2)
         END DO
         mdimt=ncosmx*ncosmx*nvmx*nkmx
         CALL getvj(veff,rm,srmin,srmax,nv,psimat,evec)
!cmh save the energies
!cmh 05.09.2005
! save them as difference to be able to do the analysis properly 
!         evjab(1:nvab,j)=evec(1:nvab)
         devj(1:nvab,j)=evec(1:nvab)
!mh 18.09.08 do this here now to make sure it is using far out energies
         do k=1,nvab
          IF(((devj(k,j))<ediss).AND.((devj(k,j))<vcut))THEN
           IF(j>mjabmax)mjabmax=j
           IF(k>mvabmax)mvabmax=k
          ENDIF 
         enddo
      END DO
!      WRITE(6,*)'done ab energies'
      CALL flush(6)
      WRITE(6,*) ' Evj (AB) product energies far out: '
      DO j=1,njab
         WRITE(6,*) ' jqm=',jqm(j)
         WRITE(6,'('' Evj '',5d13.5)')(devj(iab,j),iab=1,nvab)
         CALL flush(6)
      END DO
!cmh find now rdagger grid point for product analysis
       idagger=0
       DO ir=1,nr
        IF(wkbr(ir).LE.rdagger)idagger=ir
       END DO
!cmh now get the wavefunctions at Rdagger
!cmh need to allocate vic
       ALLOCATE(vic(nv,ncost), STAT=ier)
       DO iv=1,nv
        DO ic=1,ncost
         CALL bondc(wkbr(idagger),wksr(iv),cost(ic),rij)
         CALL surfac(vv,rij)
         vic(iv,ic)=vv
        END DO
       END DO
!this stuff here still is not quite right!!!!
!      jm=jqm(njab)
      ALLOCATE(vk(nv,mynk), STAT=ier)
      DO kk=1,mynk
       myk=mykstart+kk-myip
!       kkk=kqm(myk)
       DO j=inj(myk),njab
        jj=jqm(j)
        DO iv=1,nv
         vk(iv,kk)=0.0d0
         DO ic=1,ncost
!cmh have to get k in here as well (?) and also need pjk
          vk(iv,kk)=vk(iv,kk)+vic(iv,ic)*(pjk(jj,ic,kk)**2)*wcost(ic)
         END DO
        END DO
        rj2=jj*(jj+1)
        DO iv=1,nv
         xx=srmin+iv*dsr
         xx2=xx*xx
         veff(iv)=vk(iv,kk)+rj2/(2.d0*rm*xx2)
        END DO
        CALL getvj(veff,rm,srmin,srmax,nv,psimat,evec)
!mh do this here not in the anal program
        DO k=1,nvab
         evjab(k,j)=evec(k)
!mh 10.07.2008         devj(k,j)=devj(k,j)-evec(k)
!         WRITE(20)devj(k,j)
! have it in commons and write it out with evjab in main later
        END DO
         DO k=1,nvab
!mh 10.07.2008            IF(((evec(k))<ediss).AND.((evec(k))<vcut))THEN
!mh devj contains energies far out without j average
            IF((j.le.mjabmax).AND.(k.le.mvabmax))THEN
!               IF(j>mjabmax)mjabmax=j
!               IF(k>mvabmax)mvabmax=k
               psiab(1:nv,j,k,kk)=psimat(1:nv,k)
            ELSE
              psiab(1:nv,j,k,kk)=0.0d0
            END IF
!mh 10.07.2008
            devj(k,j)=devj(k,j)-evec(k)
         END DO
! nb indices for psiab:  iv <--> r, j <--> jqm(j), k <--> ab vib
      END DO
! end loop over mynk
      END DO
      WRITE(6,*)'done the ab wavefunctions'
      CALL flush(6)
      call flush(6)
      call flush(6)
      DEALLOCATE(psimat,v,veff,vic,vk,pjk)
! finished generation of product diatomic wavefunctions.
      IF(njab>mjabmax)THEN
         WRITE(6,*) ' should reducing njab from',njab,' to ',mjabmax
!         njab=mjabmax
      END IF
      IF(nvab>mvabmax)THEN
         WRITE(6,*) ' should reducing nvab from',nvab,' to ',mvabmax
!         nvab=mvabmax
      END IF
! now find outer "turning point".
! examinine only highest J state for highest vib state.
      vsmall=1.0d-4/SQRT(REAL(nv))
      ivtp=0
      Do kk=1,mynk
!      write(6,*)' finding outer turning point'
      DO j=1,njab
         DO iv=nv,1,(-1)
!      write(6,*)iv,psiab(iv,j,nvab,kk),vsmall
            IF(abs(psiab(iv,j,nvab,kk))>vsmall)THEN
               IF(iv>ivtp)ivtp=iv
               GOTO 50
            END IF
         END DO
50       CONTINUE 
      END DO
      WRITE(6,*) ' outer turning point for vib analysis =',ivtp
!      IF(ivtp==0)CALL clips_error(1,'SETEXP','outer turning point')
      END DO
      WRITE(6,*) ' Evj (AB) product energies at Rdagger: '
      DO j=1,njab
         WRITE(6,*) ' jqm=',jqm(j)
         WRITE(6,'('' Evj '',5d13.5)')(evjab(iab,j),iab=1,nvab)
         CALL flush(6)
         call flush(6)
      END DO
      call flush(6)
! make it so psiab has normed st dsr * Sum psi**2=1 :
      psiab=psiab/SQRT(dsr)
! Now get psibc(nbc) (which is used in subroutine initial) 
! form bare v=vab and also rbc vector which is needed in
! interpolations for initial condition :
!cmh do this for large A----BC distance and for R=R_0
      drbc=(rbcmax-rbcmin)/(nbc+1)
      ALLOCATE(v(nbc),STAT=ier)
      DO ibc=1,nbc
         rbc(ibc)=rbcmin+ibc*drbc
! rab =
         rij(1)=60.0d0
! rbc =
         rij(2)=rbcmin+ibc*drbc
! rac =
         rij(3)=rij(1)+rij(2)
! get electronic potential part
         CALL surfac(vv,rij)
         v(ibc)=vv
      END DO
! get the vib states for the given j state:
      WRITE(6,*) ' '
! add centrifugal term to v :
      rmbc=rmb*rmc/(rmb+rmc)
      rj2=jstart*(jstart+1)
!      write(6,*)' veff--'
      ALLOCATE(veff(nbc), STAT=ier)
      DO ibc=1,nbc
         xx=rbcmin+ibc*drbc; xx2=xx*xx
         veff(ibc)=v(ibc)+rj2/(2.0d0*rmbc*xx2)
!	write(6,*)ibc,xx,veff(ibc)
      END DO
      ALLOCATE(psimat(nbc,nbc), STAT=ier)
      CALL getvj(veff,rmbc,rbcmin,rbcmax,nbc,psimat,evec)
!mh this evj0 is the shifted one used later in the analysis program
      evj0large=evec(1+ivstart)
      WRITE(6,*) ' initial state energy ',evj0large
! cmh      psibc(:)=psimat(:,(1+ivstart))/SQRT(drbc)
!need psibc at R=R_0 for initital condition
      DEALLOCATE(v)
!cmh still need these
!
!cmh do now initial state at R=R_0
!
      ALLOCATE(fleg(0:jmx), STAT=ier)
      DO ibc=1,nbc
       vbc(ibc)=0.0d0
       DO ic=1,ncost
        CALL bonda(br0,rbc(ibc),cost(ic),rij)
        CALL surfac(vv,rij)
        CALL alegen(jstart,kstartbc,cost(ic),fleg)
        vbc(ibc)=vbc(ibc)+vv*(fleg(jstart)**2)*wcost(ic)
       END DO
      END DO
      DEALLOCATE(fleg)
      WRITE(6,*)' B--C pot at R_0'
      rj2=jstart*(jstart+1)
      DO ibc=1,nbc
       xx=rbcmin+ibc*drbc
       xx2=xx*xx
       veff(ibc)=vbc(ibc)+rj2/(2.d0*rmbc*xx2)
      END DO
      CALL getvj(veff,rmbc,rbcmin,rbcmax,nbc,psimat,evec)
      evj0=evec(1+ivstart)
!      WRITE(6,*)'initial state energy ',evec(1+ivstart)
      devj0=evj0large-evj0
!mh Also have this in commons and write this out later in main
      WRITE(6,*)'diff of intital state energy at infinity and at R_0= '
      WRITE(6,*)devj0
      DO ibc=1,nbc
       psibc(ibc)=psimat(ibc,(1+ivstart))/dsqrt(drbc)
      END DO       
!
      ivtpbc=0
      WRITE(6,*)'finding outer turning point'
       DO iv=nbc,1,-1
        IF(abs(psibc(iv)).gt.vsmall)then
         IF(iv.gt.ivtpbc)ivtpbc=iv
         GO TO 51
        END IF
       END DO
51    CONTINUE
!      IF(ivtpbc==0)CALL clips_error(1,'SETEXP','outer turning point')
      WRITE(6,*)'outer turning point initial wavefn= ',ivtpbc
!
!cmh have to think about this WRITE 
!get ivtpbc, rbc, v and psinrm to main and write it out with the rest
!like evjab
!mh      WRITE(20)ivtpbc
!mh      DO ibc=1,ivtpbc
!mh       WRITE(20)rbc(ibc)
!mh       WRITE(20)v(ibc)
!mh       psinrm=psibc(ibc)*dsqrt(drbc)
!mh       WRITE(20)psinrm
!mh      END DO
      DEALLOCATE(psimat,veff)
!mh still need v, is needed for analysis so has to be written out in main
!
!cmh done idagger already to calculate AB wavefunctions
! find rdagger grid point for product analysis:
!      idagger=0
!      DO ir=1,nr
!         IF(wkbr(ir)<=rdagger)idagger=ir
!      END DO
      WRITE(6,*) ' '
      WRITE(6,*) ' amp analysis on about rdagger ',rdagger
      WRITE(6,*) ' actually, idagger=',idagger
      WRITE(6,*) ' or Rdagger=',wkbr(idagger)
! Now get tropr and tropsr matrices.
! These are the grid representations of the centrifugal repulsions.
      CALL angkin2(vcut,vxadd,vyadd)
! yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
! make sure emax is not too small :
! write out kinetic energies
      tkin=(pi**2)*(1.0d0/(rmu*dbr**2)+1.0d0/(rm*dsr**2))/2.0d0
      WRITE(6,*)                                               &
     &     ' sum of grid kinetic energies in big and small r =',tkin
!9.11.99      tmax=tkin+vxadd + vyadd
      tmax=tkin+vxadd + vyadd + 2.0d0*vcut
!      tmax=tkin+vxadd + vyadd + cencut
! form scalings for Ch recursion:
! test emax
      IF(emax<((tmax+vcut)+1.000000D-1))THEN
         WRITE(6,*) ' CHANGING emax from ',emax,' to',&
     &       (tmax+vcut)+1.000000D-1
         emax=tmax+vcut+0.1d0
      END IF
! test emin
      IF(emin>vmin)THEN
         WRITE(6,*) ' emin GT vmin!!'
         WRITE(6,*) ' emin =',emin,' vmin=',vmin
         WRITE(6,*) ' reducing emin too:=',vmin
         emin=vmin
      END IF
! form scalings for Ch recursion:
      dele=emax-emin
! note: as here forces Hmin --> -1 and Hmax --> 1
      as=2.0d0/dele; bs=1.0d0-as*emax
      WRITE(6,*) ' **************************** '
      WRITE(6,*) ' as,bs ',as,bs
      WRITE(6,*) ' **************************** '
!
! yyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyyy
! scale potential and kinetic energy :
! s.t. new H=as H + bs
      vgrid=as*vgrid+bs
! scale tropr() and tropsr()
!      write(6,*)' tropr -- after scaling'
      tropr=as*tropr
!      write(6,*)' tropsr -- after scaling'
      tropsr=as*tropsr
! need to also scale the centrifugal and Coriolis terms:
      DO k=1,mynk
        myk=mykstart+k-myip
!old        fplus(inj(myk):nj,k)=as*fplus(inj(myk):nj,k)
!old        fminus(inj(myk):nj,k)=as*fminus(inj(myk):nj,k)
!        j=nj+1
	DO i=1,nr
          DO j1=1,ncost
            DO j2=1,ncost
              fplus(j1,j2,i,k)=as*fplus(j1,j2,i,k)
              fminus(j1,j2,i,k)=as*fminus(j1,j2,i,k)
            END DO
          END DO
        END DO
      END DO
!test mh .............................
!	DO i=1,nr
!          DO j1=1,ncost
!            DO j2=1,ncost
!              fplus(j1,j2,i,k)=0.0d0
!              fminus(j1,j2,i,k)=0.0d0
!            END DO
!          END DO
!        END DO
!      END DO
! mh .................................
      IF(ikin==0)THEN
         WRITE(6,*) ' '
         WRITE(6,*) ' ke done via matrix-vector method'
         WRITE(6,*) ' '
!         CALL clips_error(1,'SETEXP','bail out...')
!ggbk        do ir=1,nr
!ggbk         do irp=1,nr
!ggbk         hbrmat(ir,irp)=hbrmat(ir,irp)*as
!ggbk         enddo
!ggbk        enddo
!
!ggbk       do iv=1,nv
!ggbk        do ivp=1,nv
!ggbk        hsrmat(iv,ivp)=hsrmat(iv,ivp)*as
!ggbk       enddo
!ggbk       enddo
      ELSE
         WRITE(6,*) ' '
         WRITE(6,*) ' ***************************** '
         WRITE(6,*) ' ke done via sin fft from vfftpk '
         WRITE(6,*) ' ***************************** '
         WRITE(6,*) ' '
! multiply k**2/2mu factors appropriately
         rkinbr(1:nr)=rkinbr(1:nr) * as
         rkinsr(1:nv)=rkinsr(1:nv) * as
      END IF
! determine rk0 (passed in common to mainch2c where it is
! written to unit 20) 
! a-channel mu is that for A,BC:
      rmuabc=rma*(rmb+rmc)/(rma+rmb+rmc)
      rk0=SQRT(2.0d0*rmuabc*etran/autev)
      END SUBROUTINE setstuff

      SUBROUTINE potmat(rmin,rmax,srmin,srmax,vcut,vmin,pjk)
! **********************************************************************
!
! input arguments: rmin,rmax,srmin,srmax=R,r quadrature ranges,
!                  nrq,nsrq=no. R,r quadrature points
!                  vcut=potential cut-off
!                  ncost=no. gamma quad. pts.
!
! (Other variables, e.g., nr, nv, nj,nk, etc.  are supplied in common)
!
! must change nrmx, etc. if basis set is appropriately incremented
! in subroutine initial 
!        parameter (mxnrnv=max(nrmx,nvmx))
! for reactant channel vib state interpolations:
! for product channel analysis:
!
! Centrifugal terms and related arrays:
!ggbk	dimension bsmlr(nvmx)
!
! some quantities required for angular functions:
!
! grid information for angular quadratures:
!        dimension wcost(ncosmx)
!
! for LSTH potential used since it needs vector of bond coords
! **************************************************************
      IMPLICIT NONE
      INTEGER, PARAMETER :: jmx=2*njmx
      REAL(REAL8), INTENT(IN) :: rmin,rmax,srmin,srmax,vcut
      REAL(REAL8), INTENT(OUT) :: vmin
      INTEGER :: iv,ir,jm,k,kk,ic,j,jj,is,myk
      REAL(REAL8) :: del,x,vmax,sr,br,vv
      REAL(REAL8), DIMENSION(3) :: rij
! automatic work array
      REAL(REAL8), DIMENSION(0:jmx,ncosmx,0:nkmx+1), INTENT(OUT) :: pjk
      REAL(REAL8), DIMENSION(0:jmx) :: fleg
! Determines dvrs for R,r and  V(R,r,j,j',k) matrix elements
! n.b. V is defined here to include the diagonal Centrifugal contributio
! j(j+1)*[B(R)+b(r)] +[J(J+1) - 2 K^2] B(R)
! get Cartesian grid for r and associated kinetic energy matrix hs
! (srmin is assumed not necess.=0 for this variable
      write(6,*)srmin,srmax
      CALL gridx(rm,srmin,srmax,nv,nvmx,wksr)
      dsr=wksr(2)-wksr(1)
      bsmlr(1:nv)=0.5d0/(rm*wksr(1:nv)**2)
! determine dvrs for big r, br:
      write(6,*)rmin,rmax
      CALL gridx(rmu,rmin,rmax,nr,nrmx,wkbr)
      dbr=wkbr(2)-wkbr(1)
      bbigr(1:nr)=0.5d0/(rmu*wkbr(1:nr)**2)
! if ikin=1,then hbrmat and hsrmat are not used and
! a (real) sin fft algorithm is used in hpsi.  Here
! is the determination of the required terms in that
! case:
      IF(ikin/=0)THEN
! specialized to sin fft ke's.  The ir=0 (rkinbr=0) case not 
! evaluated.
         del=rmax-rmin
         DO ir=1,nr
            rkinbr(ir)=(0.5d0/rmu)*(REAL(ir)*pi/del)**2
         END DO
         del=srmax-srmin
         DO iv=1,nv
            rkinsr(iv)=(0.5d0/rm)*(REAL(iv)*pi/del)**2
         END DO
! initializions for vfftpk sin fft :
         CALL vsinti(nr,wsavbr)
         CALL vsinti(nv,wsavsr)
      END IF
!        write(6,*)' '
! Gauss-Legendre quadrature gamma (bending angle) integrals
! so find the special gauss-legendre points:
      CALL gauleg(-1.0d0,1.0d0,cost,wcost,ncost)
! now determine all the associated legendre polynomials, and
      jm=jqm(nj)
looppjk:      DO k=0,mynk+1
        myk=mykstart+k-myip
      if((myk.lt.1).or.(myk.gt.(myjbig+1)))cycle looppjk
        kk=kqm(myk)
! determine all Legendre polynomials from 0,..,jm
! (actually, only j>=K in array count)
        DO ic=1,ncost
          x=cost(ic)
          CALL alegen(jm,kk,x,fleg)
          pjk(kk:jm,ic,k)=fleg(kk:jm)
        END DO
      END DO  looppjk
! pjk(j,ic,k) is such that j is the actual rot quantum number,
! ic is the cos(gamma) grid and k is the index for the k quantum
! number (i.e, kqm(k) is the k quantum number corresponding to k).
! now form tmat's :
      tmat=0.0d0
!      write(6,*)' tmat: '
loop1:      DO k=0,mynk+1
        myk=mykstart+k-myip
      if((myk.lt.1).or.(myk.gt.(myjbig+1)))cycle loop1
        DO j=inj(myk),nj
          jj=jqm(j)
          tmat(1:ncost,j,k)=pjk(jj,1:ncost,k)*SQRT(wcost(1:ncost))
!	write(6,*)(ic,j,tmat(ic,j,k),ic=1,ncost)
        END DO
      END DO loop1
! for some analysis stuff -- tailored to J=0 :
! temporary
! array needed to get average value of angle
!ggbk       do j=1,nj
!ggbk       jj=jqm(j)
!ggbk       do jp=j,nj
!ggbk       jjp=jqm(jp)
!ggbk       gjj(j,jp)=0.0d0
!ggbk       g2jj(j,jp)=0.0d0
!ggbk        do ic=1,ncost
!ggbk        gamma=acos(cost(ic))
!ggbk        gjj(j,jp)=gjj(j,jp)+
!ggbk     >     wcost(ic)*gamma*pjk(jj,ic,1)*pjk(jjp,ic,1)
!ggbk        g2jj(j,jp)=g2jj(j,jp)+
!ggbk     >     wcost(ic)*(gamma**2)*pjk(jj,ic,1)*pjk(jjp,ic,1)
!ggbk        enddo
!ggbk       gjj(jp,j)=gjj(j,jp)
!ggbk       g2jj(jp,j)=g2jj(j,jp)
!ggbk       enddo
!ggbk       enddo
!
! First, some useful quantities for diagonal Centrifugal terms.
! these are kept in commons2 for use in other places
! now form grid of potential values, used in matrix element
! quadratures or wp propagation:
!      call prepot
      write(6,*)'do potential matrix'
      vmin=999999.0d0; vmax=-999999.0d0
      DO ic=1,ncost
         DO is=1,nv
            sr=wksr(is)
            DO ir=1,nr
               br=wkbr(ir)
! convert from R,r,cos(gamma) to rab,rbc,rac:
! here we are converting from channel c (product) Jacobi's --
! AB + C to bond distances
               CALL bondc(br,sr,cost(ic),rij)
               CALL surfac(vv,rij)
! vcut is applied to true full potential:
               IF(vv<vmin)vmin=vv
               IF(vv>vmax)vmax=vv
               IF(vv>vcut)vv=vcut
               vgrid(ir,is,ic)=vv
!               write(66+myk,*)ir,is,ic,vv
            END DO
         END DO
      END DO
      WRITE(6,*) ' min/max values of potential ',vmin,vmax
      WRITE(6,*) ' before imposing vcut=',vcut,' cut-off'
!
! -----------------------------------------------------------------
! determine for each bigr the corresponding rot-vibational matrix
! elements:   A potential grid is used and so no explicit quadratures
! are used here.  (Contained implicitly in hpsi now) 
! -----------------------------------------------------------------
!
! Form the  Centrifugal contribution (diagonal in j and k) as
! (note inclusion of centrifugal cut-offs)
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ggbk        do k=1,mynk
!ggbk        rk=real(kqm(k))
!ggbk        rk2=rk*rk
!ggbk        do j=inj(k),nj
!
!ggbk        rj=real(jqm(j))
!ggbk        rj2=rj*(rj+1.0d0)
! perhaps more memory efficient ways to save this --
! however this way makes cut-off easiest
!ggbk        do iv=1,nv
!ggbk        do ir=1,nr
!
!ggbk        cent=rj2* bsmlr(iv)+(rj2+bigj2-2.0d0*rk2)*bbigr(ir)
!ggbk        if(cent.gt.cencut)cent=cencut
!ggbk        centri(ir,iv,j,k)=cent
!
!ggbk        enddo
!ggbk        enddo
!ggbk        enddo
!ggbk        enddo
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      WRITE(6,*) ' done potmat '
      END SUBROUTINE potmat
!
!      
      SUBROUTINE init_parallel(jbstart,jbend,proc,ipend)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PARALLELISATION STRATEGY
!
! all loops in the program (except some i/o) have the structure:
!         DO k=1,nk
!ggbk       do j=inj(k),nj
!            DO j=1,ncost
!               DO iv=1,nv
!                  DO ir=1,nr
!                  END DO
!               END DO
!            END DO
!         END DO
!
! Program currently calculates for one jbig over range of projections k
! we can do a different jbig,k,ip on each processor.
! Need to store jbig,kstart,kend,ip for each processor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: jbstart,jbend,ipend
      INTEGER, INTENT(OUT) :: proc
      INTEGER :: ibig,i,nk,iip,ier
! now count and assign jig,ip,nk values to each processor
      proc=0
      WRITE(6,'('' following J,k,p values will be processed '')')
      write(6,*)' ipend=',ipend
      DO iip=0,ipend
        DO ibig=jbstart,jbend
! k values truncated because of a.m. coupling
! need to check convergence of this :
          nk=MIN(ibig+1,njmx)-iip
! always make this the inner loop so neighbouring processors have
! |k> and |k+-1> values - simplifies message passing.
          DO i=1,nk
            ip(proc)=iip; jbig(proc)=ibig
            kstart(proc)=i-1+iip; kend(proc)=i-1+iip
            WRITE(6,'(1x,3i8,'' on proc '',i4)')ibig,kstart(proc),iip,proc
            proc=proc+1
          END DO
        END DO
      END DO
      WRITE(6,'('' need '',i4,'' processors '')')proc
!      IF(proc>nnode)CALL clips_error(1,'INIT_PARALLEL', &
!     &     ' more jbig,k values than nodes ')
      IF(inode+1>proc)THEN
        WRITE(6,'('' no work for this processor '')')
        CALL mpi_finalize(ier); STOP
      END IF
! to test, make them all do the same thing
! also need to put IF(.FALSE.)THEN in hpsi2.f to prevent
! communication
!                this test worked OK on SP2 8/3/99
!      myjbig=0
!      mykstart=0
!      mykend=0
!      myip=0
      myjbig=jbig(inode); myip=ip(inode)
      mykstart=kstart(inode); mykend=kend(inode)
      mynk=mykend-mykstart+1
      IF(mynk>1)WRITE(6,'('' mynk > 1 on proc '')')inode
      WRITE(6,'('' inode,myjbig,mykstart,mykend '',4i4)')inode,myjbig,mykstart,mykend
      END SUBROUTINE init_parallel
      
      END MODULE setstuff_mod
