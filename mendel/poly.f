C  October 2010
C  Program for computing breast cancer penetrance using family data.
C  TP53 data for  Clare Turnbull. 
C  This version assumes  fixed background incidences independent of birth cohort. 
C  Breast Cancer incidences for England and Wales


c 25/10/2010 : This code contstraints the overall incidence and estimates 
c frequency and RR simultaneously while allowing for a polygenic background

C  Outstanding 30/8/2010: Cohort specific incidences? Are years of bith available?

C Recode the genotype data to appear at the locus field rather than sens?
C                        



C     MODIFIED VERSION OF MENDEL TO ALLOW A POLYGENIC COMPONENT - 
C     HYPERGEOMETRIC POLYGENIC MODEL 

C     ARRAY AND VARIABLE DECLARATIONS.
C
      PARAMETER(LENC=100,LENI=100000,LENL=100,LENR=9000000,MXLOCI=10
     :,NEXTRA=1)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION EXTRA(NEXTRA),RARRAY(LENR)
      INTEGER IARRAY(LENI),COND
      CHARACTER*8 CARRAY(LENC),LNAME(MXLOCI),BASE,BATFIL*40,LOCFIL*40
     :,MUTLOC,OUTFIL*40,PEDFIL*40,TITLE*40,TRAVEL,XXSIGN,XYSIGN
      LOGICAL LARRAY(LENL),ASYCV,ECHO,LUMP,ORDERD,PMODE,STAND

C
C     DEFAULT VALUES FOR SOME VARIABLES.  TO RUN THE PROGRAM IN
C     BATCH MODE FILL IN THE NAME OF THE BATCH FILE BATFIL.  SET
C     PMODE TO TRUE TO DO CALCULATIONS IN PRODUCT MODE AND LUMP
C     TO TRUE TO AMALGAMATE ALLELES AT EACH LOCUS. 
C

      DATA MXTWIN,ABSENT/10,-1.0D20/
      DATA CONV,NCONV,MXSTEP,DP/1.0D-4,4,3,1.0D-7/
      DATA BATFIL,LUMP,ORDERD,PMODE/'',
     &.FALSE.,.FALSE.,.FALSE./
      
C     DEFAULT VALUES FOR THE PROBLEM MENU.
C
      DATA TITLE,LOCFIL/'TP53_data', 'locus3.dat'/
      DATA PEDFIL,OUTFIL/'pedigree_poly.txt','poly_out.dat'/
      DATA ECHO,XXSIGN,XYSIGN,NVAR/.FALSE.,'F','M',7/
      DATA MUTLOC,XXRATE,XYRATE,COND/'PGENE',0.0D0,0.0D0,2/
      DATA EXTRA/0.0D0/
      DATA LNAME/'MAJORG','PGENE',8*' '/
      DATA BASE,STAND,TRAVEL,NPOINT/'E',.FALSE.,'SEARCH',1/
      DATA NPAR,NCNSTR,ASYCV,MXITER/1,0,.TRUE.,150/
C
C     THE ARRAYS AND VARIABLES BEGIN A LONG DESCENT INTO THE
C     PROGRAM.  SAY GOODBYE TO THEM AND WISH THEM LUCK.
C
      CALL MENDEL(EXTRA,RARRAY,IARRAY,CARRAY,LNAME,LARRAY,ABSENT
     :,CONV,DP,XXRATE,XYRATE,COND,LENC,LENI,LENL,LENR,MXITER,MXLOCI
     :,MXSTEP,MXTWIN,NCNSTR,NCONV,NEXTRA,NPAR,NPOINT,NVAR,BASE,BATFIL
     :,LOCFIL,MUTLOC,OUTFIL,PEDFIL,TITLE,TRAVEL,XXSIGN,XYSIGN,ASYCV
     :,ECHO,LUMP,ORDERD,PMODE,STAND)
      END

      SUBROUTINE INITAL(ALLFRQ,CNSTR,CVALUE,EXTRA,GRID,PAR,PARMAX
     1,PARMIN,PNAME,XLINK,XXRATE,XYRATE,MAXALL,MUTATE,NCNSTR,NEXTRA
     2,NLOCI,NPAR,NPOINT,NVAR,PROBLM,UNIT3,TRAVEL)



C
C     IN THIS SUBROUTINE THE USER SHOULD DEFINE THE INITIAL
C     PARAMETER VALUES, THE PARAMETER BOUNDS, AND THE LINEAR
C     EQUALITY CONSTRAINTS FOR A LIKELIHOOD SEARCH.  WHEN A
C     GRID OF LIKELIHOOD VALUES IS DESIRED, THEN ONLY DEFINE 
C     THE ARRAY GRID.  PARAMETER NAMES CAN BE OPTIONALLY INPUT.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION ALLFRQ(NLOCI,MAXALL),CNSTR(NCNSTR,NPAR)
     1,CVALUE(NCNSTR),EXTRA(NEXTRA),GRID(NPOINT,NPAR),PAR(NPAR)
     2,PARMAX(NPAR),PARMIN(NPAR)
      INTEGER PROBLM,UNIT3
      CHARACTER*8 PNAME(NPAR),TRAVEL
      LOGICAL XLINK(NLOCI)

      integer pgene, mgene, nall2
      double precision inc(0:79), popbr(0:79), surv(0:80)
      double precision fdens(0:79),sumsurv

      common/breast/nall2,pgene, mgene
      common/rates/inc,fdens
C
C     Define the polygenic locus as the locus MUTATE

      pgene=mutate
      if(nloci.eq.2) then
       if(mutate.eq.1) then
        mgene=2
        pgene=1
       else
        pgene=2
        mgene=1
       endif
      else
       if(mutate.eq.1) then
       pgene=1
       mgene=0
       else
       pgene=0
       mgene=1
       endif
      endif
      write(*,*) 'mutate', mutate
      write(*,*) 'pgene=', pgene
      write(*,*) 'mgene=', mgene


C NUMBER OF ALLELES IN POLYGENIC COMPONENT
C NALL2 MUST BE IN A COMMON BLOCK 

      if(pgene.eq.0) then
       nall2=1
      else
      do 1 i=1,maxall
      if(allfrq(mutate,i).lt.1D-7) goto2
1     continue
      i=maxall+1
2     nall2=i-1
      endif
c      write(*,*) 'maxall', maxall
c      write(*,*) 'nall2=',nall2

      PNAME(1)='logRR'

      PAR(1)=2.0d0	
      PARMIN(1)=0.0d0
      PARMAX(1)=3.0d0



C----------------------------------------------------------------
C  Read in the population incidence rates for England and Wales 1993-97
c  per 100000 population.
C----------------------------------------------------------------

      data  popbr/0.050204,0.050204,0.050204,0.050204,
     :0.050204,0.0124487,0.0124487,0.0124487,0.0124487,
     :0.0124487,0,0,0,0,0,0.1671877,0.1671877,0.1671877,
     :0.1671877,0.1671877,1.191661,1.191661,1.191661,
     :1.191661,1.191661,8.00212,8.00212,8.00212,8.00212,
     :8.00212,26.28934,26.28934,26.28934,26.28934,26.28934,
     :57.95005,57.95005,57.95005,57.95005,57.95005,106.2854,
     :106.2854,106.2854,106.2854,106.2854,173.896,173.896,173.896,
     :173.896,173.896,244.4749,244.4749,244.4749,244.4749,244.4749,
     :245.3651,245.3651,245.3651,245.3651,245.3651,269.8981,269.8981,
     :269.8981,269.8981,269.8981,233.0805,233.0805,233.0805,
     :233.0805,233.0805,278.9083,278.9083,278.9083,278.9083,
     :278.9083,292.8277,292.8277,292.8277,292.8277,292.8277 /


CC Population breast cancer inc. rates for England and Wales, 1983-1987


      sumsurv=0.0d0
      surv(0)=1.0d0
      do  iage=0,79
          inc(iage)=popbr(iage)/100000.0d0
          fdens(iage)=surv(iage)*inc(iage)
          sumsurv=sumsurv+inc(iage)
          surv(iage+1)=exp(-sumsurv)
      end do






      END



      SUBROUTINE OUTPUT(EXTRA,PAR,SCORE,PNAME,LOGLIK,FINAL,ITER,MAXPAR
     1,NEXTRA,NPAR,NSTEP,UNIT3,BASE,TRAVEL,STAND,UMOVE)
C
C     THIS SUBROUTINE OUTPUTS THE LOGLIKELIHOOD AND PARAMETERS
C     AT EACH ITERATION.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION EXTRA(NEXTRA),PAR(MAXPAR),SCORE(MAXPAR),LOGLIK
      INTEGER FINAL,UNIT3
      CHARACTER*8 PNAME(MAXPAR),BASE,TRAVEL
      LOGICAL STAND,UMOVE
      SAVE BEST,IBEST,START

      double precision inc(0:79), fdens(0:79)
      double precision  cumbrrisk
      double precision  rrbr(0:69), survbr(0:70)


      common/rates/inc,fdens

 
C
      IF (ITER.EQ.1) THEN
      BEST=LOGLIK
      START=LOGLIK
      WRITE(UNIT3,10) (PNAME(I),I=1,NPAR)
      write(*,10) (PNAME(I),I=1,NPAR)
 10   FORMAT(/,' ITER  NSTEP  LOGLIKELIHOOD',(T28,4(4X,A8),:))
      END IF
      IF (LOGLIK.GE.BEST) THEN
      BEST=LOGLIK
      IBEST=ITER
      END IF
      IF (STAND) LOGLIK=LOGLIK-START
      IF (BASE.EQ.'10') LOGLIK=LOG10(EXP(1.0D0))*LOGLIK
      WRITE(UNIT3,20) ITER,NSTEP,LOGLIK,(PAR(I),I=1,NPAR)
      write(*,20) ITER,NSTEP,LOGLIK,(PAR(I),I=1,NPAR) 
 20   FORMAT(/,I4,3X,I3,3X,D14.7,(T28,4(1X,D11.4),:))





c          if(iter.eq.final) then
c            write(*,*) '                                         ' 
c            write(UNIT3,*) '                                         ' 
c            write(*,*) 'AGE        Cum BC Risk        '
c            write(UNIT3,*) 'AGE        Cum BC Risk        '
c            do i=1,79
c               if(mod(i,10).eq.0 .and. i.ge.20) then
c                write(*,23) i, 1.0d0-survbr(i)
c                write(UNIT3,23) i, 1.0d0-survbr(i)
c23              format(/, I3,8x,F11.4)
c               endif  
c            end do    
c           endif 




     
      IF (ITER.EQ.FINAL) WRITE(UNIT3,30) IBEST
 30   FORMAT(/,' THE MAXIMUM LOGLIKELIHOOD OCCURS AT ITERATION',I4,'.')
      END



      SUBROUTINE NEWLIK(EXTRA,PAR,LOGLIK,COND,ITER,NEXTRA,NPAR,NPED,PED)
C
C     THIS SUBROUTINE ALLOWS THE USER TO MODIFY THE LOGLIKELIHOOD
C     OF A PARTICULAR PEDIGREE.  FOR INSTANCE, IN GENETIC COUNSELING
C     PROBLEMS IT IS MAY BE NECESSARY TO FORM CONDITIONAL PROBABILITIES.
C     COND TELLS WHICH PEDIGREE, IF ANY, TO CONDITION ON.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION EXTRA(NEXTRA),PAR(NPAR),LOGLIK
      INTEGER COND,PED
C
      if(ped.gt.nped/2) loglik=-loglik
      END


      SUBROUTINE APEN(EXTRA,PAR,PEN,VAR,GENES,XLINK,ABSENT,XXRATE
     1,XYRATE,FIRST,LAST,MUTATE,NEXTRA,NGTYPE,NLOCI,NPAR,NVAR,PED
     2,PER,MALE,like)
C



C     APEN SUPPLIES THE PENETRANCE PROBABILITIES FOR EACH PERSON
C     IN A PEDIGREE.  THE CURRENT VERSION USES COX PROPORTIONAL
C     HAZARDS MODEL
C
      IMPLICIT NONE

      double precision absent
      double precision xxrate, xyrate
      integer          first,last
      integer          nextra
      integer          ngtype
      integer          nloci
      integer          npar
      integer          nvar
      integer          ped, per
      logical          male
      integer          mutate
      integer          like
      integer          maxlik

      DOUBLE PRECISION EXTRA(NEXTRA),PAR(NPAR),PEN(NGTYPE),VAR(NVAR)
      INTEGER GENES(FIRST:LAST,2,NGTYPE)
      LOGICAL XLINK(NLOCI)

      parameter(maxlik=12)




      double precision inc(0:79),ff(maxlik,0:80,0:2),rr(0:79,1:2)
      double precision lambda(maxlik,0:79,0:2), lampalb2(maxlik,0:79)
      double precision p0,p1,p2
      double precision sum0, sum1, sum2, sumpalb2
      double precision sumweight, sumfdens
      double precision fdens(0:79),frr(0:79), ffpalb2(maxlik,0:80)
      double precision avfrr,resvar,resstd
      double precision f0,f1,f2,g0,g1,g2
      double precision sp, xp, pp
      double precision bin
      double precision poprr(maxlik,0:79)


      integer i,k, kk, jj
      integer idis
      integer agelfu, agebc, age, agebc2, ageoc, agedeath
      integer ageother, ageother2
      integer isex
      integer kage, igen
      integer unit3

      integer pgene, mgene, nall2, mm     

 
      common/breast/nall2,pgene, mgene
      common/rates/inc,fdens

C  Polygnic component. Total variance is predifined, based on previous 
C  segregation analyses.

      mm=nall2-1
      sp=1.29d0

      

c---------------------------------------------------------------------
c The risks will correspond to the risk of a first BC
c Censoring occurs at the first of bc dx, oc dx, mastectomy or agelfu
c Note that this is the minimum of the age at last follow-up
c All censored at age 80. 
c Ignore oophorectomy [CAN BE CHANGED LATER]

c Define the disease status, idis(=0 if unaffected, =1 if BC)

c Define the censoring variable in the data. If icens=1 then the individual
c is censored at birth.
c--------------------------------------------------------------------- 

       agebc  = var(1)
       agebc2 = var(2)
       ageoc  = var(3)
       ageother=var(4)
       ageother2=var(5)
       agelfu = var(6)
       agedeath= var(7)
  

c By default everyone censored at 0, unless age information is available:

       age = 0
       idis=0


	if (ageother.eq. 999) then
	 	if(agelfu.gt.0 .and. agedeath.gt.0) then	 
			ageother=min(agelfu,agedeath)
		elseif(agelfu.gt.0 .and. agedeath.eq.0) then
			ageother=agelfu
		elseif(agelfu.eq.0 .and. agedeath.gt.0) then
			ageother=agedeath
		else
			ageother=0
		endif
	endif


	if (agebc.eq. 999) then
	 	if(agelfu.gt.0 .and. agedeath.gt.0) then	 
			agebc=min(agelfu,agedeath)
		elseif(agelfu.gt.0 .and. agedeath.eq.0) then
			agebc=agelfu
		elseif(agelfu.eq.0 .and. agedeath.gt.0) then
			agebc=agedeath
		else
			agebc=0
		endif
	endif

	if (ageoc.eq. 999) then
	 	if(agelfu.gt.0 .and. agedeath.gt.0) then	 
			ageoc=min(agelfu,agedeath)
		elseif(agelfu.gt.0 .and. agedeath.eq.0) then
			ageoc=agelfu
		elseif(agelfu.eq.0 .and. agedeath.gt.0) then
			ageoc=agedeath
		else
			ageoc=0
		endif
	endif



C Censor at the first cancer. 

	if(agebc.gt.0 .and. ageoc.gt.0 .and. ageother.gt.0) then
	 age=min(agebc, ageoc, ageother)
	elseif(agebc.gt.0 .and. ageoc.gt.0 .and. ageother.eq.0) then
	 age=min(agebc, ageoc)
	elseif(agebc.gt.0 .and. ageoc.eq.0 .and. ageother.gt.0) then
	 age=min(agebc, ageother)
	elseif(agebc.eq.0 .and. ageoc.gt.0 .and. ageother.gt.0) then
	 age=min(ageoc, ageother)
	elseif(agebc.gt.0 .and. ageoc.eq.0 .and. ageother.eq.0) then
	 age=agebc
	elseif(agebc.eq.0 .and. ageoc.gt.0 .and. ageother.eq.0) then
	 age=ageoc
	elseif(agebc.eq.0 .and. ageoc.eq.0 .and. ageother.gt.0) then
	 age=ageother
	elseif(agelfu.gt.0 .and. agedeath.eq.0) then
		age=agelfu
	elseif(agelfu.eq.0 .and. agedeath.gt.0) then
	        age=agedeath
	elseif(agelfu.gt.0 .and. agedeath.gt.0) then
		age=min(agelfu, agedeath)
	endif

	if (age.eq.agebc .and. agebc.gt.0) idis=1

c	write(*,*) agebc,ageoc,ageother,agelfu,agedeath,age,idis
  

c---------------------------------------------------------------------
c Censor at age 80
c---------------------------------------------------------------------
       
       if(age.ge.80) then
          idis=0
          age=80
       endif   
       
c--------------
C GENERAL MODEL
c-------------

      p0=(1.0d0-0.0004d0)**2
      p1=2*0.0004d0*(1.0d0-0.0004d0)
      p2=0.0004d0**2




c---------------------------------------------------------------------
c Define the relative risk parameters and residual variance:
c---------------------------------------------------------------------

      sumweight=0.0d0
      sumfdens=0.0d0


       do i=0,79

c Assign relative risks
          if(i.lt.20) then 
            rr(i,1) = 1.0d0
            rr(i,2) = 1.0d0
          elseif(i.lt.25) then
            rr(i,1) = exp(par(1))
            rr(i,2) = exp(par(1))
          elseif(i.lt.30) then
            rr(i,1) = exp(par(1))
            rr(i,2) = exp(par(1))
          elseif(i.lt.35) then
            rr(i,1) = exp(par(1))
            rr(i,2) = exp(par(1))
          elseif(i.lt.40) then
            rr(i,1) = exp(par(1))
            rr(i,2) = exp(par(1))
          elseif(i.lt.45) then
            rr(i,1) = exp(par(1))
            rr(i,2) = exp(par(1))
          elseif(i.lt.50) then
            rr(i,1) = exp(par(1))
            rr(i,2) = exp(par(1))
          elseif(i.lt.55) then
            rr(i,1) = exp(par(1))
            rr(i,2) = exp(par(1))
          elseif(i.lt.60) then
            rr(i,1) = exp(par(1))
            rr(i,2) = exp(par(1))
          elseif(i.lt.65) then
            rr(i,1) = exp(par(1))
            rr(i,2) = exp(par(1))
          elseif(i.lt.70) then
            rr(i,1) = exp(par(1))
            rr(i,2) = exp(par(1))
           elseif(i.lt.80) then
            rr(i,1) = exp(par(1))
            rr(i,2) = exp(par(1))
          endif

ccc Compute the residual variance based on the current parameter estimates.

          frr(i)=(p0+p1*rr(i,1)**2+p2*rr(i,2)**2)/
     :              (p0+p1*rr(i,1)+p2*rr(i,2))**2


c Weighted by density distribution

          if (i.ge.20) then
             sumweight=sumweight+frr(i)*fdens(i)
             sumfdens=sumfdens+fdens(i)
          end if   
 
       end do    
  
c Residual variance

       avfrr=sumweight/sumfdens
       
       resvar=sp**2-log(avfrr)

       if(ped.eq.1 .and. per.eq.1) then
          write(*,*) "sp^2:",sp**2,"   resvar:",resvar
       end if   
       resstd=sqrt(resvar)

c       resstd=sp

c-------------------------
          if(like.gt.maxlik) write(*,*) 'Like too large',like


          if(ped.eq.1) then
C-------------------------




C------------------
C SURVIVAL PROBABILITIES BY AGE AND GENOTYPE (0=dd,1=Dd,2=DD)
C ASSUMING POLYGENIC COMPONENT =0
C-----------------

      ff(like,0,0)=1.0d0
      ff(like,0,1)=1.0d0
      ff(like,0,2)=1.0d0

      sum0=0.0d0
      sum1=0.0d0
      sum2=0.0d0
      sumpalb2=0.0d0


C Here the baseline hazard is computed iteratively. The polygenic component 
C is approximated by the binomial distribution.


ccc Compute the residual variance based on the current parameter estimates.


      do kage=0,79

        f0=0.0d0
        f1=0.0d0
        f2=0.0d0
        g0=0.0d0
        g1=0.0d0
        g2=0.0d0

        do  k=0,2*nall2-2
         xp=exp(resstd*(dfloat(k)-mm)/sqrt(mm/2.0)) 
         pp=bin(2*nall2-2,k)

         f0=f0+pp*FF(like,kage,0)**xp
         f1=f1+pp*FF(like,kage,1)**xp
         f2=f2+pp*FF(like,kage,2)**xp
         g0=g0+pp*FF(like,kage,0)**xp*xp
         g1=g1+pp*FF(like,kage,1)**xp*xp
         g2=g2+pp*FF(like,kage,2)**xp*xp

        end do 
      

      lambda(like,kage,0)=inc(kage)*(f0*p0+f1*p1+f2*p2)/
     :             (g0*p0+g1*p1*rr(kage,1)+g2*p2*rr(kage,2))

      lambda(like,kage,1)= lambda(like,kage,0)*rr(kage,1)
      lambda(like,kage,2)= lambda(like,kage,0)*rr(kage,2)


      lampalb2(like,kage)=lambda(like,kage,1)*g1/f1	

      poprr(like,kage)=lampalb2(like,kage)/inc(kage)	

      sum0=sum0+lambda(like,kage,0)
      sum1=sum1+lambda(like,kage,1)
      sum2=sum2+lambda(like,kage,2)

      sumpalb2=sumpalb2+lampalb2(like,kage)

      ff(like,kage+1,0)=exp(-sum0) 
      ff(like,kage+1,1)=exp(-sum1)
      ff(like,kage+1,2)=exp(-sum2) 

      ffpalb2(like,kage+1)=exp(-sumpalb2)

      end do

	if(per.eq.1) write(*,*) 20, 1-ffpalb2(like,20)
	if(per.eq.1) write(*,*) 30, 1-ffpalb2(like,30)
	if(per.eq.1) write(*,*) 40, 1-ffpalb2(like,40)
	if(per.eq.1) write(*,*) 50, 1-ffpalb2(like,50)
	if(per.eq.1) write(*,*) 60, 1-ffpalb2(like,60)
	if(per.eq.1) write(*,*) 70, 1-ffpalb2(like,70)
	if(per.eq.1) write(*,*) 80, 1-ffpalb2(like,80)
 
 	if(per.eq.1) write(*,*) 25, poprr(like,25)
	if(per.eq.1) write(*,*) 35, poprr(like,35)
	if(per.eq.1) write(*,*) 45, poprr(like,45)
	if(per.eq.1) write(*,*) 55, poprr(like,55)
	if(per.eq.1) write(*,*) 65, poprr(like,65)
	if(per.eq.1) write(*,*) 75, poprr(like,75)


      endif



       if(male) then
        isex=1
       else
        isex=2
       endif





      do  igen=1,NGTYPE

        if(isex.eq.1) then
            pen(igen)=1.0
        else



c     polygenic component

         kk=genes(pgene,1,igen)+genes(pgene,2,igen)-2
         xp=exp(resstd*(dfloat(kk)-mm)/sqrt(mm/2.0))

c Major genotypes

      if(genes(mgene,1,igen).eq.1.and.genes(mgene,2,igen).eq.1) then
        k=0
      elseif(genes(mgene,1,igen).eq.2.and.genes(mgene,2,igen).eq.2) then
        k=2
      else
        k=1
      endif



        pen(igen)=ff(like,age,k)**xp


        if(idis.eq.1) then

          pen(igen)=pen(igen)*xp*lambda(like,age,k)
     
        endif


      endif

      end do


      END




      SUBROUTINE APRIOR(ALLFRQ,EXTRA,PAR,PRIOR,VAR,GENES,XLINK
     1,ABSENT,XXRATE,XYRATE,FIRST,LAST,MAXALL,MUTATE,NEXTRA,NGTYPE,
     2NLOCI,NPAR,NVAR,PED,PER,MALE)

C
C     PRIOR SUPPLIES THE PRIOR PROBABILITIES FOR EACH FOUNDER IN
C     A PEDIGREE.  THE CURRENT VERSION IS VALID IF HARDY-WEINBERG
C     AND LINKAGE EQUILIBRIUM HOLD AT ALL LOCI.  LOCI MAY BE AUTOSOMAL
C     OR X-LINKED.
C
c the prior probabilities for each founder for the poligenic component 
c are etermined by the binomial distribution

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION ALLFRQ(NLOCI,MAXALL),EXTRA(NEXTRA),PAR(NPAR)
     1,PRIOR(NGTYPE),VAR(NVAR)
      INTEGER FIRST,GENES(NLOCI,2,NGTYPE),PED,PER
      LOGICAL XLINK(NLOCI),MALE,TWICE
      integer pgene, mgene, nall2
      common/breast/nall2,pgene, mgene


      allfrq(mgene,2)=0.0004d0
      allfrq(mgene,1)=1-0.0004d0



       MM=NALL2-1
       DO 2 I=1,nall2
2      allfrq(pgene,i)=bin(mm,i-1)



      DO 10 I=1,NGTYPE
      P=1.0D0
      TWICE=.FALSE.
      DO 20 LOCUS=1,NLOCI
      IG1=GENES(LOCUS,1,I)
      IF (XLINK(LOCUS).AND.MALE) THEN
      P=P*ALLFRQ(LOCUS,IG1)
      ELSE
      IG2=GENES(LOCUS,2,I)
      P=P*ALLFRQ(LOCUS,IG1)*ALLFRQ(LOCUS,IG2)
        END IF
 20   CONTINUE
 10   PRIOR(I)=P

      END


      SUBROUTINE ATRANS(EXTRA,PAR,TRANS,VARI,VARJ,GAMETE,GENES,XLINK
     1,ABSENT,XXRATE,XYRATE,FIRST,LAST,MUTATE,NEXTRA,NGTYPE,NLOCI
     2,NPAR,NVAR,PED,PERI,PERJ,MALEI,MALEJ)
C
C     TRANS SUPPLIES THE TRANSMISSION PROBABILITIES FOR EACH
C     PARENT-OFFSPRING PAIR IN A PEDIGREE.  THE SUFFIX I INDICATES
C     THE PARENT AND THE SUFFIX J THE CHILD.  SEVERAL MULTIPLE
C     LOCUS GENOTYPES OF THE PARENT ARE PASSED VIA THE ARRAY GENES.
C     THE ARRAY GAMETE REPRESENTS ONE OF THE TWO GAMETES MAKING UP
C     A MULTIPLE LOCUS GENOTYPE OF THE CHILD.  

C	THIS VERSION IMPLEMENTS AN APPROXIMATION TO THE MIXED MODEL,
C	WHERE THE POLYGENIC COMPONENT (GIVEN BY THE SECOND "LOCUS") 
C	IS APPROXIMATED BY A BINOMIAL DISTRIBUTION WITH PARAMETERS (2N,1/2)
C	TRANSMISSION PROBABILITIES ARE GIVEN BY A HYPERGEOMETRIC DISTRIBUTION
C	
C	THE FIRST LOCUS IS A SIMPLE M ALLELE SYSTEM (USUALLY M=2)
C	WITH STANDARD MENDELIAN TRANSMISSION


      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION EXTRA(NEXTRA),PAR(NPAR),TRANS(NGTYPE)
     1,VARI(NVAR),VARJ(NVAR)
      INTEGER FIRST,GAMETE(FIRST:LAST),GENES(FIRST:LAST,2,NGTYPE)
     1,PED,PERI,PERJ
      LOGICAL XLINK(NLOCI),MALEI,MALEJ
      integer pgene, mgene, nall2      
      common/breast/nall2,pgene, mgene

      

      MM=NALL2-1
      do  I=1,NGTYPE

         TRANS(I)=1.0
         
         T=0.0D0
         IF (GENES(mgene,1,I).EQ.GAMETE(1)) T=T+0.5D0
         IF (GENES(mgene,2,I).EQ.GAMETE(1)) T=T+0.5D0
         TRANS(I)=TRANS(I)*T
 
       K=GENES(pgene,1,I)+GENES(pgene,2,I)-2
       K1=GAMETE(pgene)-1

      IF(K1.GT.K.or.K-K1.GT.MM) THEN
       TRANS(I)=0.0D0
      ELSE
       TRANS(I)=TRANS(I)*BINC(K,K1)*BINC(2*MM-K,MM-K1)
     :     /BINC(2*MM,MM)  
      end if




      end do   

      END






      FUNCTION BIN(N,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,R,M

      IF (R.GT.N/2) THEN
       M=N-R
      ELSE
       M=R
      ENDIF
      BIN=0.5**N
      DO 10 I=1,M
      BIN=BIN*DFLOAT(N-I+1)/DFLOAT(I)
10    CONTINUE
      RETURN
      END

      FUNCTION BINC(N,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,R,M

      IF (R.GT.N/2) THEN
       M=N-R
      ELSE
       M=R
      ENDIF
      BINC=1.0
      DO 10 I=1,M
      BINC=BINC*DFLOAT(N-I+1)/DFLOAT(I)
10    CONTINUE
      RETURN
      END
