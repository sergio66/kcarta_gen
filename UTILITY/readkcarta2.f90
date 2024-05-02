! this program takes in the UNFORMATTED output from atmos.x and allows the
! user to read the final radiance and save to a FORMATTED or UNFORMATTED file

! either do Gaussian SRF,AIRS SRF or NOAA SRF -- 
! see paramSRF.f for AIRS details, paramNOAA for NOAA details

! NOTE : Absoft behaves very wierdly with fseek,ftell
! so for Linux machines, use g77 or PDF Fortran
! this file is for g77 (which uses subroutines fseek instead of functions)
!          iDummy=fseek(iIOUN,offset,1)       !! if using absoft
!          CALL fseek(iIOUN,OffSet,1)         !! if using g77
!
! example : ../BIN/readkcarta.x fin=../WORK/junk.dat

      PROGRAM readkcarta
!      use ifport             ! for getenv, fseek, ftell

      IMPLICIT NONE

      include 'convolve.param'

! the definitions of all these variables are found in the subroutines
      CHARACTER*80 caInName
      CHARACTER*130 caaMixedPathInfo(kMixFilRows)
      INTEGER iIOUN,iNumPathsOut
      INTEGER iNpmix,iNumMixPathsOut,iNatm
      INTEGER iaNumLayersOut(kMaxAtm),iaNumLayersInAtm(kMaxAtm)

      REAL rFr1,rFr2

      INTEGER iLorS,iNOAA,iInstrType
      INTEGER iInstr,rInPointSp,rInStart,rFWHM
      INTEGER iBinOrNC

      INTEGER iMainType,iSubMainType,iNumberOut,ikMaxPts,iWhichStore
      REAL rFrLow,rFrHigh,rDelta
      INTEGER iOutNumber,iaNumOut(kMaxPrint),iDummy,iDatapoints
      INTEGER iTotal,iStart,iStore,iFreqPts
      REAL rTemp,raFreq(kMaxEntire)
      REAL raaEntire(kMaxEntire,100),raConvolve(kMaxEntire)
      INTEGER iI,iJ,iK,iLoop,iStartOffset,iJump,iDataSize,iFileWhere

!      INTEGER ftell

! for vararg
      CHARACTER*80 FIN, FOUT, mode
      INTEGER IARGC,IOERR
      INTEGER J
      INTEGER NARGS
 
      CHARACTER*80 BUF
      CHARACTER*80 VAL
      CHARACTER*80 VAR

!************************************************************************
! see /asl/packages/rtpV201/utils/rdinfo_subset.f
!      ------------
!      Set defaults
!      ------------
       IOERR = 6
       FIN =  'kcarta.dat'           ! input filename
       FOUT = 'kcarta.txt'           ! output filename

!      -----------------------------------------------------------------
!      Loop on program parameters
!      --------------------------
!      Determine the number of command-line arguments
       NARGS=IARGC()

       IF (NARGS .EQ. 0) THEN
          WRITE(IOERR,1010)
 1010     FORMAT('readkcarta.f90 must be run with at least 1 argument')
          WRITE(IOERR,1011)
 1011     FORMAT('   fin  = <filename> txt input file  {mandatory}','//'    &
                 '   fout = <filename> rtp output file {optional, set to kcarta.txt}')
          STOP
       ENDIF

!      Loop over the command-line arguments
       DO iI = 1, NARGS

!         Pull out the ith argument
          CALL GETARG(iI, BUF)

!         Find the "=" character in the command-line argument string
          J=INDEX(BUF, '=')

          IF (J .NE. 0) THEN

!            Name of variable
             VAR = BUF(1:J-1)
             CALL UPCASE(VAR)

!            Specified value
             VAL = BUF(J+1:LEN(BUF))

!            Big "IF" to set parameters
!            ----------------------------
             IF (VAR(1:3) .EQ. 'FIN') THEN
                FIN=VAL
             ELSEIF (VAR(1:4) .EQ. 'FOUT') THEN
                FOUT=VAL
             ELSE
                WRITE(IOERR,1020) VAR
 1020           FORMAT('Unknown command-line argument: ',A6)
                STOP

             ENDIF

          ENDIF
       ENDDO  ! end of loop over command-line arguments
       write(*,'(A,A)') 'input  file name = ',fin
       write(*,'(A,A)') 'output file name = ',fout

!************************************************************************

 90   CONTINUE
      iInstrtype=0

      print *,'Do you want a (+1) binary output (0) text output (-1) netcdf output? '
      read *,iBinOrNC

       caInName = fin

! first read header to figure out iStore ******************
      iIOUN=10
!      OPEN(UNIT=iIOUN,FILE=caInName,STATUS='OLD',FORM='UNFORMATTED',RECL=4)
      OPEN(UNIT=iIOUN,FILE=caInName,STATUS='OLD',FORM='UNFORMATTED')

      CALL readmainheader(iIOUN,rFr1,rFr2,iLorS)

      CALL readgaspaths(iIOUN,iNumPathsOut,iLorS)

      CALL readmixedpaths(iIOUN,iNpmix,iNumMixPathsOut,iLorS)

      IF (iNpmix .EQ. 0) THEN
        iNumMixPathsOut=0
      END IF

      CALL readatmospheres(iIOUN,iNpmix,iNatm,iaNumLayersInAtm,iaNumLayersOut)

      IF (iNatm .EQ. 0) THEN
        DO ii=1,iNatm
          iaNumLayersOut(ii)=0
        END DO
      END IF

      IF (iNatm .GT.  0) THEN
        DO ii=1,iNatm
          IF (iaNumLayersOut(ii) .EQ. -1) THEN
            iaNumLayersOut(ii)=iaNumLayersInAtm(ii)
          END IF
        END DO
      END IF

! this is the number of paths/mixed paths/layers that will be output each chunk
      iStore=0
      DO ii=1,iNatm
        iStore=iStore+iaNumLayersOut(ii)
      END DO                                 !num of radiances to output

      iWhichStore=iNumPathsOut+iNumMixPathsOut !num of paths/MP to output
      iStore=iStore+iWhichStore                !total num to output

      print *,'For each k-comp file used in the atmos.x run '
      print *,'Total # of paths to output       = ',iNumPathsOut
      print *,'Total # of mixed paths to output = ',iNumMixPathsOut
      print *,'Total # of radiances to output   = ',iStore-iWhichStore
      print *,'GRAND total to be output         = ',iStore

! the following info is for each k-comp file

! cccccccc this is the summary info that HAS to be read in!!!!!!!!
      read(iIOUN) iTotal,iOutNumber     !!! typically iTotal = 89 (chunks) and iOutNumber = 5 (TwoSlab rads)
      read(iIOUN) (iaNumOut(iI),iI=1,iOutNumber)

      iJ=0
      DO iI=1,iOutNumber
        iJ=iJ+iaNumOut(iI)
        END DO
      IF (iJ .NE. iStore) THEN
        print *,'hmm : number of things to read  aintcha jivin'
        STOP
      END IF

      print *,'kMaxEntire,kMaxPts = ',kMaxEntire,kMaxPts
      DO iJ = 1,iTotal
        CALL ReadData2(iIOUN,iJ,iOutNumber,raFreq,raaEntire)
      END DO
      CLOSE(iIOUN)

      call output2(caInName,raFreq,raaEntire,iTotal,iOutNumber,iBinORNC)

      END PROGRAM

!************************************************************************
!!! typically iTotal = 89 (chunks) and iOutNumber = 5 (TwoSlab rads)
      SUBROUTINE output2(caInName,raFreq,raaEntire,iTotal,iOutNumber,iBinOrNC)

      use netcdf


      IMPLICIT NONE

      include 'convolve.param'

      CHARACTER*80 caInName
      REAL raaEntire(kMaxEntire,100),raFreq(kMaxEntire)
      INTEGER iTotal,iOutNumber,iBinOrNC

! local variables
      CHARACTER(LEN=80) caOutName
      INTEGER iI,iJ,iK,iIOUN,iFreqPts

!! netcdf
      INTEGER :: ncid,ndims

      iIOUN=11
      iFreqPts = kMaxPts * iTotal

! get the name of output file
      CALL GetOutName2(caInName,caOutName,iBinOrNC)

      print *,'opening caOutName ',caOutName

      IF (iBinOrNC .EQ. 1) THEN            !binary format for rdairs.m
! put in the header info ... iI tells how many spectra to expect
! only need to save this header when we are outputting the VERY first spectra
        
        ii = iOutNumber
        OPEN(UNIT=iIOUN,FILE=caOutName,STATUS='UNKNOWN',FORM='UNFORMATTED')
        write(iIOUN) 1,iFreqPts,iI
        write(iIOUN) (raFreq(iI),iI=1,iFreqPts)
        CLOSE(iIOUN)
! now output data!!!
        OPEN(UNIT=iIOUN,FILE=caOutName,STATUS='OLD',FORM='UNFORMATTED',ACCESS='APPEND')
        DO iK = 1,iOutNumber 
          write(iIOUN) (raaEntire(iJ,iK),iJ=1,iFreqPts)
        END DO
        CLOSE(iIOUN)

      ELSEIF (iBinOrNC .EQ. 0) THEN     !column text format for rdairs.m
! put in the header info ... iI tells how many spectra to expect
! only need to save this header when we are outputting the VERY first spectra
        OPEN(UNIT=iIOUN,FILE=caOutName,STATUS='UNKNOWN',FORM='FORMATTED')
        DO iK = 1,iOutNumber
          write(iIOUN,5050) -9999.0,-9999.0
          DO iJ=1,iFreqPts
            write(iIOUN,5050) raFreq(iJ),raaEntire(iJ,iK)
          END DO
        END DO
        CLOSE(iIOUN)
 5050 FORMAT(f10.4,'  ',1pe12.5)

      ELSEIF (iBinOrNC .EQ. -1) THEN     !column text format for rdairs.m

        !! https://home.chpc.utah.edu/~thorne/computing/Examples_netCDF.pdf for f90)        
        !! call f90_netcdf 

        !! https://www.unidata.ucar.edu/software/netcdf/examples/programs/sfc_pres_temp_wr.f for f77
        call f77_netcdf(caOutName,raFreq,raaEntire,iTotal,iOutNumber)

      END IF

      RETURN
      END   

!************************************************************************
! this takes the input file name, finds the root (before .dat or .datJAC) 
! and appends .datCON
      SUBROUTINE GetOutName2(caInName,caOutName,iBinOrNC)

      CHARACTER*80 caInName,caOutName
      INTEGER iBinOrNC
      
      INTEGER iI

      DO iI=1,80
        caOutName(iI:iI)=' '
      END DO

! search for the .dat in caInName
      iI=80
 15   CONTINUE
      IF ((caInName(iI:iI) .NE. '.') .AND. (iI .GT. 1)) THEN
        iI=iI-1
        GO TO 15
      END IF
      IF (iI .EQ. 1) THEN
        print *,'wierd input name!!'
        STOP
      END IF

      caOutName(1:iI)=caInName(1:iI)
      IF (iBinOrNC .EQ. 0) THEN
        caOutName(iI+1:iI+6) = 'datTXT'
      ELSEIF (iBinOrNC .GT. 0) THEN
        caOutName(iI+1:iI+6) = 'datBIN'
      ELSEIF (iBinOrNC .LT. 0) THEN
        caOutName(iI+1:iI+2) = 'nc'
      END IF

      RETURN
      END

!************************************************************************
       SUBROUTINE UPCASE(BUF)

!      Convert a string to upper case

       CHARACTER*(*) BUF

       INTEGER I
       INTEGER IC
       INTEGER J

       DO I=1,LEN(BUF)
          IC=ICHAR( BUF(I:I) )
          IF (IC .GE. ICHAR('a') .AND. IC .LE. ICHAR('z')) THEN
             J=IC + (ICHAR('A') - ICHAR('a'))
             BUF(I:I)=CHAR(J)
          ENDIF
       ENDDO

       RETURN
       END

!************************************************************************
!     This is part of the netCDF package.
!     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
!     See COPYRIGHT file for conditions of use.

!     This example writes some surface pressure and temperatures. It is
!     intended to illustrate the use of the netCDF fortran 77 API. The
!     companion program sfc_pres_temp_rd.f shows how to read the netCDF
!     data file created by this program.

!     This program is part of the netCDF tutorial:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

!     Full documentation of the netCDF Fortran 77 API can be found at:
!     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77

!     $Id: sfc_pres_temp_wr.f,v 1.9 2007/01/24 19:31:54 russ Exp $

!!! typically iTotal = 89 (chunks) and iOutNumber = 5 (TwoSlab rads)

      SUBROUTINE f77_netcdf(caOutName,raFreq,raaEntire,iTotal,iOutNumber)
      implicit none
      include 'netcdf.inc'
      include 'convolve.param'

      REAL rTemp,raFreq(kMaxEntire)
      REAL raaEntire(kMaxEntire,100)
      INTEGER iTotal,iOutNumber

!     This is the name of the data file we will create.
      character*(*) caOutname

!     We are writing 1D data, raFreq         so need one netCDF dimensions.
!     We are writing 1D data, P(z),T(z),Q(z) so need one netCDF dimensions.
!     We are writing 2D data, raaEntire      so need two netCDF dimensions.

!     Always check the return code of every netCDF function call. In
!     this example program, any retval which is not equal to nf_noerr
!     (0) will call handle_err, which prints a netCDF error message, and
!     then exits with a non-zero return code.
     
!     We are writing 2D data
      integer NDIMS
      parameter (NDIMS=2)
      integer NFREQ, NRADLAY

!     When we create netCDF files, variables and dimensions, we get an ID for each one.
      integer ncid, varid, dimids(NDIMS)
      integer x_dimid, y_dimid

!     This is the data array we will write. 
!      REAL data_out(NRADLAY, NFREQ) which is TRANSPOSE of raaEntire

!     Loop indexes, and error handling.
      integer x, y, retval

!!! typically iTotal = 89 (chunks) and iOutNumber = 5 (TwoSlab rads)
      NFREQ   = iTotal * kMaxPts
      NRADLAY = iOutNumber

!     Create some pretend data. If this wasn't an example program, we
!     would have some real data to write, for example, model output.
!      data_out = raaEntire(1:iOutnumber,1:iTotal)

!     Create the netCDF file. The nf_clobber parameter tells netCDF to
!     overwrite this file, if it already exists.
      retval = nf_create(caOutName, NF_CLOBBER, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     Define the dimensions. NetCDF will hand back an ID for each.
      retval = nf_def_dim(ncid, "x", NFREQ, x_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, "y", NRADLAY, y_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     The dimids array is used to pass the IDs of the dimensions of
!     the variables. Note that in fortran arrays are stored in
!     column-major format.
      dimids(2) = x_dimid
      dimids(1) = y_dimid

!     Define the variable. The type of the variable in this case is
!     NF_INT  (4-byte integer).
!     NF_REAL (4-byte real).
      retval = nf_def_var(ncid, "data", NF_REAL, NDIMS, dimids, varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     End define mode. This tells netCDF we are done defining metadata.
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

!     Write the pretend data to the file. Although netCDF supports
!     reading and writing subsets of data, in this case we write all the
!     data in one operation.
      retval = nf_put_var_real(ncid, varid, transpose(raaEntire(1:NFREQ,1:iOutnumber)))
      if (retval .ne. nf_noerr) call handle_err(retval)

!     Close the file. This frees up any internal netCDF resources
!     associated with the file, and flushes any buffers.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

      print *,'*** SUCCESS writing example file ',caOutName,iTotal,iOutNumber
      end


!************************************************************************
!SUBROUTINE f90_netcdf !! https://home.chpc.utah.edu/~thorne/computing/Examples_netCDF.pdf for f90)

! USE netcdf
! IMPLICIT NONE
! REAL(KIND=4), DIMENSION(NX), INTENT(IN) :: xpos
! REAL(KIND=4), DIMENSION(NY), INTENT(IN) :: ypos
! REAL(KIND=4), DIMENSION(NX,NY), INTENT(IN) :: idata
! INTEGER(KIND=4) :: ncid, x_dimid, y_dimid
! INTEGER(KIND=4) :: x_varid, y_varid, varid
! INTEGER(KIND=4), DIMENSION(2) :: dimids
! INTEGER(KIND=4), INTENT(IN) :: NX, NY
! CHARACTER(LEN=50), INTENT(IN) :: outfile
! !Create the netCDF file.
! CALL check(nf90_create(outfile, NF90_CLOBBER, ncid))
! !Define the dimensions.
! CALL check(nf90_def_dim(ncid, "lon", NX, x_dimid))
! CALL check(nf90_def_dim(ncid, "lat", NY, y_dimid))
! !Define coordinate variables
! CALL check(nf90_def_var(ncid, "lon", NF90_REAL, x_dimid, x_varid))
! CALL check(nf90_def_var(ncid, "lat", NF90_REAL, y_dimid, y_varid))
! dimids = (/ x_dimid, y_dimid /)
! !Define variable
! CALL check(nf90_def_var(ncid, "Perturbations", NF90_FLOAT, dimids, varid))
! CALL check(nf90_enddef(ncid)) !End Definitions
! !Write Data
! CALL check(nf90_put_var(ncid, x_varid, xpos))
! CALL check(nf90_put_var(ncid, y_varid, ypos))
! CALL check(nf90_put_var(ncid, varid, idata))
! CALL check(nf90_close(ncid))

! END SUBROUTINE writegrid
!=====================================================
      subroutine handle_err(errcode)
      implicit none
      include 'netcdf.inc'
      integer errcode

      print *, 'Error: ', nf_strerror(errcode)
      stop 2
      end

!************************************************************************
! https://stdlib.fortran-lang.org/proc/check.html
!subroutine xcheck(condition, msg, code, warn)
    !! version: experimental
    !!
    !! Checks the value of a logical condition
    !! ([Specification](../page/specs/stdlib_error.html#description))
    !!
    !!##### Behavior
    !!
    !! If `condition == .false.` and:
    !!
    !!   * No other arguments are provided, it stops the program with the default
    !!     message and exit code `1`;
    !!   * `msg` is provided, it prints the value of `msg`;
    !!   * `code` is provided, it stops the program with the given exit code;
    !!   * `warn` is provided and `.true.`, it doesn't stop the program and prints
    !!     the message.
    !!
    !!##### Examples
    !!
    !!* If `a /= 5`, stops the program with exit code `1`
    !!  and prints `Check failed.`
    !!``` fortran
    !!  call check(a == 5)
    !!```
    !!
    !!* As above, but prints `a == 5 failed`.
    !!``` fortran
    !!  call check(a == 5, msg='a == 5 failed.')
    !!```
    !!
    !!* As above, but doesn't stop the program.
    !!``` fortran
    !!  call check(a == 5, msg='a == 5 failed.', warn=.true.)
    !!```
    !!
    !!* As example #2, but stops the program with exit code `77`
    !!``` fortran
    !!  call check(a == 5, msg='a == 5 failed.', code=77)
    !!```

    !
    ! Arguments
    ! ---------

!    logical, intent(in) :: condition
!    character(*), intent(in), optional :: msg
!    integer, intent(in), optional :: code
!    logical, intent(in), optional :: warn
!    character(*), parameter :: msg_default = 'Check failed.'

!    if (.not. condition) then
!        if (optval(warn, .false.)) then
!            write(stderr,*) optval(msg, msg_default)
!        else
!            call error_stop(optval(msg, msg_default), optval(code, 1))
!        end if
!    end if

! end subroutine xcheck
!************************************************************************
     
