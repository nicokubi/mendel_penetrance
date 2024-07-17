C
C  {{ label }}
C
C  Program for computing breast cancer penetrance using family data.
C  TP53 data for  Clare Turnbull.
C  This version assumes  fixed background incidences independent of birth cohort.
C  Breast Cancer incidences for England and Wales


      PARAMETER(LENC=100,LENI=100000,LENL=100,LENR=1000000,MXLOCI=10
     :,NEXTRA=100)
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
      DATA CONV,NCONV,MXSTEP,DP/1.0D-7,4,3,1.0D-7/
      DATA BATFIL,LUMP,ORDERD,PMODE/' ',.FALSE.,.FALSE.,.FALSE./
C
C     DEFAULT VALUES FOR THE PROBLEM MENU.
C
      DATA TITLE,LOCFIL/'TP53_data','locusbr.dat'/
      DATA PEDFIL,OUTFIL/'pedigree.txt','single_out.dat'/
      DATA ECHO,XXSIGN,XYSIGN,NVAR/.FALSE.,'F','M',7/
      DATA EXTRA/100*0.0D0/
      DATA LNAME/'MAJOR',9*' '/
      DATA MUTLOC,XXRATE,XYRATE,COND/' ',0.0d0,0.0d0,2/
      DATA BASE,STAND,TRAVEL,NPOINT/'E',.FALSE.,'{{ travel|default("SEARCH") }}',1/
      DATA NPAR,NCNSTR,ASYCV,MXITER/{{ npar }},0,.TRUE.,{{ mxiter|default(599) }}/
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

      double precision   popbr(0:{{ maxage - 1 }})

      common /popincid/popbr



C-----------------------------------------------------------------------
C The log-likelihood is parameterised in terms of the log relative risks
C Assume all individuals are censored at 70
C-----------------------------------------------------------------------



      pname(1)='RR'


      DO I=1,NPAR
        par(I) = {{ '%0.1f'| format(parinit|default(0.5)) }}d0
        parmin(I) =  {{ '%0.1f'| format(parmin|default(0)) }}d0
        parmax(I) =  {{ '%0.1f'| format(parmax|default(5)) }}d0
      END DO



C----------------------------------------------------------------
C  Read in the population incidence rates for England and Wales 1993-97
c  per 100000 population.
C----------------------------------------------------------------
      {% if cancer_type == "colorectal" %}
      data  popbr/0.0010, 0.0050, 0.0090, 0.0220, 0.0430,
     :0.0650, 0.0860, 0.1160, 0.1990, 0.2910,
     :0.3830, 0.4750, 0.5730, 0.7090, 0.8510,
     :0.9930, 1.1350, 1.2730, 1.3820, 1.4870,
     :1.5920, 1.6960, 1.8240, 2.0860, 2.3720,
     :2.6570, 2.9430, 3.2620, 3.7800, 4.3320,
     :4.8840, 5.4360, 6.0360, 6.9270, 7.8660,
     :8.8060, 9.7450, 10.7850, 12.4340, 14.1830,
     :15.9330, 17.6820, 19.5540, 22.1620, 24.8920,
     :27.6220, 30.3510, 33.5010, 39.1800, 45.2780,
     :51.3740, 57.4690, 62.9320, 64.6190, 65.6750,
     :66.7300, 67.7830, 69.2050, 72.8420, 76.8460,
     :80.8460, 84.8440, 89.0500, 94.5240, 100.2040,
     :105.8790, 111.5490, 117.1350, 122.2510, 127.2830,
     :132.3070, 137.3240, 142.7000, 150.2750, 158.2030,
     :166.1160, 174.0110, 181.8510, 189.4490, 196.9850,
     :204.4950, 211.9770, 219.2250, 225.2180, 230.9740,
     :236.6940, 242.3780, 246.8160, 243.9830, 239.9320 /
      {% elif cancer_type == "brain" %}
      data  popbr/3.7450, 3.9430, 4.1700, 4.1870, 3.9950,
     :3.8020, 3.6100, 3.4240, 3.2740, 3.1300,
     :2.9860, 2.8420, 2.7030, 2.5980, 2.4990,
     :2.4000, 2.3010, 2.2110, 2.1790, 2.1570,
     :2.1340, 2.1120, 2.1100, 2.2330, 2.3770,
     :2.5200, 2.6640, 2.8020, 2.9140, 3.0200,
     :3.1270, 3.2330, 3.3330, 3.3910, 3.4430,
     :3.4940, 3.5460, 3.6050, 3.7070, 3.8170,
     :3.9260, 4.0360, 4.1550, 4.3280, 4.5120,
     :4.6950, 4.8780, 5.0940, 5.5140, 5.9670,
     :6.4190, 6.8720, 7.3290, 7.8080, 8.2920,
     :8.7750, 9.2580, 9.7450, 10.2540, 10.7670,
     :11.2790, 11.7920, 12.3220, 12.9590, 13.6130,
     :14.2670, 14.9210, 15.5690, 16.1840, 16.7940,
     :17.4030, 18.0120, 18.5990, 19.0550, 19.4890,
     :19.9220, 20.3550, 20.7620, 21.0140, 21.2390,
     :21.4640, 21.6890, 21.8390, 21.5450, 21.1770,
     :20.8100, 20.4420, 19.9830, 18.9740, 17.8740 /
      {% elif cancer_type == "breast" %}
      data  popbr/0.0220, 0.0140, 0.0050, 0.0000, 0.0010,
     :0.0020, 0.0030, 0.0040, 0.0040, 0.0030,
     :0.0020, 0.0010, 0.0020, 0.0180, 0.0360,
     :0.0540, 0.0720, 0.1070, 0.2420, 0.3950,
     :0.5470, 0.6990, 0.9420, 1.7310, 2.6120,
     :3.4920, 4.3720, 5.4060, 7.3600, 9.4670,
     :11.5750, 13.6820, 16.0160, 19.7130, 23.6360,
     :27.5580, 31.4810, 36.1100, 44.9770, 54.5500,
     :64.1220, 73.6920, 83.0390, 91.0500, 98.8380,
     :106.6240, 114.4080, 121.7820, 126.7060, 131.2210,
     :135.7340, 140.2440, 144.5710, 147.8110, 150.8680,
     :153.9230, 156.9750, 160.7430, 168.8230, 177.6160,
     :186.4040, 195.1860, 204.0540, 213.4700, 222.9710,
     :232.4630, 241.9470, 250.8900, 256.6350, 261.8400,
     :267.0370, 272.2260, 276.4580, 274.9920, 272.5720,
     :270.1490, 267.7210, 264.9980, 260.5260, 255.7610,
     :250.9930, 246.2230, 241.4770, 236.8880, 232.3200,
     :227.7470, 223.1710, 218.1060, 210.1350, 201.6820 /
     {% elif cancer_type == "leukemia" %}
      data  popbr/4.5370, 6.1280, 7.9460, 8.3480, 7.3330,
     :6.3180, 5.3030, 4.3910, 4.1020, 3.9160,
     :3.7300, 3.5440, 3.3780, 3.3340, 3.3100,
     :3.2860, 3.2620, 3.2320, 3.1650, 3.0920,
     :3.0190, 2.9460, 2.8860, 2.9050, 2.9370,
     :2.9690, 3.0010, 3.0420, 3.1340, 3.2360,
     :3.3370, 3.4380, 3.5460, 3.6920, 3.8450,
     :3.9970, 4.1490, 4.3230, 4.6220, 4.9430,
     :5.2640, 5.5840, 5.9250, 6.3860, 6.8660,
     :7.3470, 7.8280, 8.3480, 9.1060, 9.9040,
     :10.7020, 11.4990, 12.3370, 13.4160, 14.5360,
     :15.6550, 16.7740, 17.9650, 19.5830, 21.2720,
     :22.9610, 24.6500, 26.4400, 28.8440, 31.3500,
     :33.8540, 36.3580, 38.8850, 41.5580, 44.2550,
     :46.9490, 49.6420, 52.3860, 55.4460, 58.5560,
     :61.6630, 64.7660, 67.7850, 70.3220, 72.7750,
     :75.2230, 77.6660, 80.0070, 81.7670, 83.4270,
     :85.0820, 86.7320, 88.0530, 87.4330, 86.4890 /
     {% elif cancer_type == "lfs" %}
      data  popbr/9.1,
     :13.4,
     :13.4,
     :13.4,
     :13.4,
     :8.3,
     :8.3,
     :8.3,
     :8.3,
     :8.3,
     :7.700000000000001,
     :7.700000000000001,
     :7.700000000000001,
     :7.700000000000001,
     :7.700000000000001,
     :7.9,
     :7.9,
     :7.9,
     :7.9,
     :7.9,
     :7.800000000000001,
     :7.800000000000001,
     :7.800000000000001,
     :7.800000000000001,
     :7.800000000000001,
     :12.899999999999999,
     :12.899999999999999,
     :12.899999999999999,
     :12.899999999999999,
     :12.899999999999999,
     :24.4,
     :24.4,
     :24.4,
     :24.4,
     :24.4,
     :47.599999999999994,
     :47.599999999999994,
     :47.599999999999994,
     :47.599999999999994,
     :47.599999999999994,
     :88.1,
     :88.1,
     :88.1,
     :88.1,
     :88.1,
     :144.9,
     :144.9,
     :144.9,
     :144.9,
     :144.9,
     :202.40000000000003,
     :202.40000000000003,
     :202.40000000000003,
     :202.40000000000003,
     :202.40000000000003,
     :282.2,
     :282.2,
     :282.2,
     :282.2,
     :282.2,
     :397.79999999999995,
     :397.79999999999995,
     :397.79999999999995,
     :397.79999999999995,
     :397.79999999999995,
     :534.5,
     :534.5,
     :534.5,
     :534.5,
     :534.5,
     :658.5,
     :658.5,
     :658.5,
     :658.5,
     :658.5,
     :743.5,
     :743.5,
     :743.5,
     :743.5,
     :743.5 /
     {% endif %}



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
C

      SAVE BEST,IBEST,START
C
      IF (ITER.EQ.1) THEN
      BEST=LOGLIK
      START=LOGLIK
      write(*,10) (pname(i),i=1,npar)
      WRITE(UNIT3,10) (PNAME(I),I=1,NPAR)
 10   FORMAT(/,' ITER  NSTEP  LOGLIKELIHOOD',(T28,4(4X,A8),:))
      END IF
      IF (LOGLIK.GE.BEST) THEN
      BEST=LOGLIK
      IBEST=ITER
      END IF
      IF (STAND) LOGLIK=LOGLIK-START
      IF (BASE.EQ.'10') LOGLIK=LOG10(EXP(1.0D0))*LOGLIK
      write(*,20) iter,nstep,loglik,(par(i),i=1,npar)
      WRITE(UNIT3,20) ITER,NSTEP,LOGLIK,(PAR(I),I=1,NPAR)
 20   FORMAT(/,I4,3X,I3,3X,D14.7,(T28,4(1X,D11.4),:))


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
C     THE LOGLIKELIHOOD FOR EACH PEDIGREE GIVEN FREE RECOMBINATION
C     IS STORED IN THE ARRAY EXTRA.  THE LIKELIHOOD IS THEN ADJUSTED
C     TO BE A CONVEX COMBINATION OF THIS LIKELIHOOD AND THE CURRENT
C     LIKELIHOOD WITH THE RECOMBINATION PARAMETERS TAKING VALUES LESS
C     THAN .5.



c Ascertainment part should follow as the second half of the predigrees:

C     if(ped.gt.nped/2) loglik=-loglik


      END






      SUBROUTINE APEN(EXTRA,PAR,PEN,VAR,GENES,XLINK,ABSENT,XXRATE
     1,XYRATE,FIRST,LAST,MUTATE,NEXTRA,NGTYPE,NLOCI,NPAR,NVAR,PED
     2,PER,MALE,like)
C
C     PEN SUPPLIES THE PENETRANCE PROBABILITIES FOR EACH PERSON
C     IN A PEDIGREE.  THE CURRENT VERSION IS VALID ONLY FOR SIMPLE
C     EITHER/OR TRAITS.
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
      integer          mutate, like

      DOUBLE PRECISION EXTRA(NEXTRA),PAR(NPAR),PEN(NGTYPE),VAR(NVAR)
      INTEGER GENES(FIRST:LAST,2,NGTYPE)
      LOGICAL XLINK(NLOCI)

      double precision  rrbr(0:{{ maxage - 1}})
      double precision  cumbrrisk
      double precision  popbr(0:{{ maxage - 1}})
      double precision  ffbr(0:{{ maxage }})
      double precision  ffncbr(0:{{ maxage }})
      double precision  cumncbr
      double precision  p1, p2
      double precision  lambda(0:{{ maxage - 1}},0:1)


      integer i
      integer idis
      integer agelfu, agebc, age, agebc2, ageoc, agedeath
      integer ageother, ageother2
      integer isex
      integer imut
      integer is
      integer iage





      common /popincid/popbr
      common /risk/ffbr



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
       idis = 1
	elseif(agebc.gt.0 .and. ageoc.gt.0 .and. ageother.eq.0) then
	 age=min(agebc, ageoc)
       idis = 1
	elseif(agebc.gt.0 .and. ageoc.eq.0 .and. ageother.gt.0) then
	 age=min(agebc, ageother)
       idis = 1
	elseif(agebc.eq.0 .and. ageoc.gt.0 .and. ageother.gt.0) then
	 age=min(ageoc, ageother)
       idis = 1
	elseif(agebc.gt.0 .and. ageoc.eq.0 .and. ageother.eq.0) then
	 age=agebc
       idis = 1
	elseif(agebc.eq.0 .and. ageoc.gt.0 .and. ageother.eq.0) then
	 age=ageoc
       idis = 1
	elseif(agebc.eq.0 .and. ageoc.eq.0 .and. ageother.gt.0) then
	 age=ageother
       idis = 1
	elseif(agelfu.gt.0 .and. agedeath.eq.0) then
		age=agelfu
	elseif(agelfu.eq.0 .and. agedeath.gt.0) then
	        age=agedeath
	elseif(agelfu.gt.0 .and. agedeath.gt.0) then
		age=min(agelfu, agedeath)
	endif

C	if (age.eq.agebc .and. agebc.gt.0) idis=1

c	write(*,*) agebc,ageoc,ageother,agelfu,agedeath,age,idis


c---------------------------------------------------------------------
c Censor at age {{ maxage }}
c---------------------------------------------------------------------

       if(age.ge.{{ maxage }}) then
          idis=0
          age={{ maxage }}
       endif


c---------------------------------------------------------------------
c For carriers we estimate the incidence rates using the
c relative risks. In the fixed incidence version, no Rel Risk applies
c to non-carriers.
c---------------------------------------------------------------------

      {% if rr_model == "linear" %}

       do i=0,{{ maxage - 1}}
          if(i.lt.{{ age_cutoffs[1] }}) then
            rrbr(i) = exp( (par(2) - par(1))/{{ age_cutoffs[1] }} * (i - 0) + par(1) )
          {% for N in range(3, npar) -%}
          elseif(i.lt.{{ age_cutoffs[N-1] }}) then
            rrbr(i) = exp( (par({{ N }}) - par({{ N-1 }}))/{{ age_cutoffs[N-1] - age_cutoffs[N-2] }} * (i - {{ age_cutoffs[N-2] }}) + par({{ N-1}}) )
          {% endfor -%}
          elseif(i.lt.{{ maxage }}) then
            rrbr(i) = exp( (par({{ npar }}) - par({{ npar - 1 }}))/{{ maxage - age_cutoffs[npar - 2] }} * (i - {{ age_cutoffs[npar - 2] }}) + par({{ npar - 1 }}) )
          endif
       end do
      {% elif rr_model == "piecewise" %}
       do i=0,{{ maxage - 1}}
          if(i.lt.{{ age_cutoffs[1] }}) then
            rrbr(i) = exp(par(1))
          {% for N in range(2, npar) -%}
          elseif(i.lt.{{ age_cutoffs[N] }}) then
            rrbr(i) = exp(par({{ N }}))
          {% endfor -%}
          else
            rrbr(i) = exp(par({{ npar }}))
          endif
       end do
      {% endif %}


c---------------------------------------------------------------------
c This code constraints overall incidence to the population rates to obtain
c the female non carriers have incidence rates.
c---------------------------------------------------------------------

c non carrier frequency
       p1=(1-0.0004)**2

c carrier frequency:
       p2=1-p1

c initialise survival probabilities:

c non-carriers:

       ffncbr(0)=1.0d0

c carriers:

       ffbr(0)=1.0d0

c initialise sums:

       cumncbr=0.0d0
       cumbrrisk=0.0d0


       do iage=0,{{ maxage - 1 }}

          lambda(iage,0)=(popbr(iage)/100000.0d0)*
     :             (p1*ffncbr(iage)+p2*ffbr(iage))/
     :              (p1*ffncbr(iage)+p2*rrbr(iage)*ffbr(iage))

          lambda(iage,1)=lambda(iage,0)*rrbr(iage)



c---------------------------------------------------------------------
c Cpmpute the survivor function for non carriers
c incidence rates
c---------------------------------------------------------------------

        cumncbr = cumncbr + lambda(iage,0)

        ffncbr(iage+1) = exp(-cumncbr)


c---------------------------------------------------------------------
c Compute the survivor function for carriers,
c---------------------------------------------------------------------


          cumbrrisk = cumbrrisk+lambda(iage,1)

          ffbr(iage+1) = exp(-cumbrrisk)


       end do



c---------------------------------------------------------------------
c Code the sex of each member of the pedigree
c---------------------------------------------------------------------

      if(male) then
       isex=1
      else
       isex=2
      end if




c---------------------------------------------------------------------
c Assume a  dominant model. So carriers denoted
c by is=2 and non carriers by  is=1
c---------------------------------------------------------------------

      DO 10 I=1,NGTYPE

       if(genes(1,1,I).eq.1.and.genes(1,2,I).eq.1) then
        is=1
       else
        is=2
       endif




c---------------------------------------------------------------------
c Assume that males do not develop either cancer and specify the penetrance
c for each male as 1*fact
c---------------------------------------------------------------------

{% if cancer_type != "breast" %}C{% endif %}      if(isex.eq.1) then
{% if cancer_type != "breast" %}C{% endif %}         pen(i)=1.0d0
{% if cancer_type != "breast" %}C{% endif %}      else

c---------------------------------------------------------------------
c Non-Carriers develop the cancers according to population specific rates.
c---------------------------------------------------------------------

c=====================================================================
         if(is.eq.1) then
c=====================================================================


c---------------------------------------------------------------------
c Compute the penetrance for females
c---------------------------------------------------------------------


c---------------------------------------------------------------------
c For disease free females:
c---------------------------------------------------------------------

            if(idis.eq.0) then
              pen(i) = ffncbr(age)


c---------------------------------------------------------------------
c If she develops  breast cancer
c---------------------------------------------------------------------
            elseif(idis.eq.1) then
              pen(i) = ffncbr(age)*lambda(age,0)


            endif


c---------------------------------------------------------------------
c Carriers develop the cancers according to the fixed background
c incidence times the relative risks.
c---------------------------------------------------------------------

c=====================================================================
         else
c=====================================================================


c---------------------------------------------------------------------
c For disease free carriers:
c---------------------------------------------------------------------

            if(idis.eq.0) then
               pen(i) = ffbr(age)

c---------------------------------------------------------------------
c If carriers  develops  breast cancer
c---------------------------------------------------------------------
            elseif(idis.eq.1) then
               pen(i) = ffbr(age)*lambda(age,1)

            endif

c=====================================================================
         endif
c=====================================================================



{% if cancer_type != "breast" %}C{% endif %}      endif

c	write(*,*) ped, per, age, pen(i)

10    continue
      return
      END




      SUBROUTINE APRIOR(ALLFRQ,EXTRA,PAR,PRIOR,VAR,GENES,XLINK
     1,ABSENT,XXRATE,XYRATE,FIRST,LAST,MAXALL,MUTATE,NEXTRA,NGTYPE
     2,NLOCI,NPAR,NVAR,PED,PER,MALE)
C
C     PRIOR SUPPLIES THE PRIOR PROBABILITIES FOR EACH FOUNDER IN
C     A PEDIGREE.  THE CURRENT VERSION IS VALID IF HARDY-WEINBERG
C     AND LINKAGE EQUILIBRIUM HOLD AT ALL LOCI.  LOCI MAY BE AUTOSOMAL
C     OR X-LINKED.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION ALLFRQ(NLOCI,MAXALL),EXTRA(NEXTRA),PAR(NPAR)
     1,PRIOR(NGTYPE),VAR(NVAR)
      INTEGER FIRST,GENES(FIRST:LAST,2,NGTYPE),PED,PER
      LOGICAL XLINK(NLOCI),MALE
C

      allfrq(1,2)=0.0004d0
      allfrq(1,1)=1.0d0-allfrq(1,2)



      DO 10 I=1,NGTYPE
      P=1.0D0
      DO 20 LOCUS=FIRST,LAST
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


      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION EXTRA(NEXTRA),PAR(NPAR),TRANS(NGTYPE)
     1,VARI(NVAR),VARJ(NVAR)
      INTEGER FIRST,GAMETE(FIRST:LAST),GENES(FIRST:LAST,2,NGTYPE)
     1,PED,PERI,PERJ
      LOGICAL XLINK(NLOCI),MALEI,MALEJ

      DO 10 I=1,NGTYPE

         TRANS(I)=1.0

         T=0.0D0
         IF (GENES(1,1,I).EQ.GAMETE(1)) T=T+0.5D0
         IF (GENES(1,2,I).EQ.GAMETE(1)) T=T+0.5D0
         TRANS(I)=TRANS(I)*T



10    CONTINUE

      END
