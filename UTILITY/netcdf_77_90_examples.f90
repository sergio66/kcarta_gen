c https://www.unidata.ucar.edu/software/netcdf/examples/programs/
c https://www.unidata.ucar.edu/software/netcdf/examples/programs/
c https://www.unidata.ucar.edu/software/netcdf/examples/programs/
c https://www.unidata.ucar.edu/software/netcdf/examples/programs/
c https://www.unidata.ucar.edu/software/netcdf/examples/programs/
c https://www.unidata.ucar.edu/software/netcdf/examples/programs/
c https://www.unidata.ucar.edu/software/netcdf/examples/programs/
c https://www.unidata.ucar.edu/software/netcdf/examples/programs/
c https://www.unidata.ucar.edu/software/netcdf/examples/programs/
c https://www.unidata.ucar.edu/software/netcdf/examples/programs/

c************************************************************************

C https://www.unidata.ucar.edu/software/netcdf/examples/programs/sfc_pres_temp_wr.f
C     This is part of the netCDF package.
C     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
C     See COPYRIGHT file for conditions of use.

C     This example writes some surface pressure and temperatures. It is
C     intended to illustrate the use of the netCDF fortran 77 API. The
C     companion program sfc_pres_temp_rd.f shows how to read the netCDF
C     data file created by this program.

C     This program is part of the netCDF tutorial:
C     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

C     Full documentation of the netCDF Fortran 77 API can be found at:
C     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77

C     $Id: sfc_pres_temp_wr.f,v 1.9 2007/01/24 19:31:54 russ Exp $

      program sfc_pres_temp_wr
      implicit none
      include 'netcdf.inc'

C     This is the name of the data file we will create.
      character*(*) FILE_NAME
      parameter (FILE_NAME='sfc_pres_temp.nc')
      integer ncid

C     We are writing 2D data, a 6 x 12 lat-lon grid. We will need two
C     netCDF dimensions.
      integer NDIMS
      parameter (NDIMS=2)
      integer NLATS, NLONS
      parameter (NLATS = 6, NLONS = 12)
      character*(*) LAT_NAME, LON_NAME
      parameter (LAT_NAME='latitude', LON_NAME='longitude')
      integer lon_dimid, lat_dimid

C     In addition to the latitude and longitude dimensions, we will also
C     create latitude and longitude netCDF variables which will hold the
C     actual latitudes and longitudes. Since they hold data about the
C     coordinate system, the netCDF term for these is: "coordinate
C     variables."
      real lats(NLATS), lons(NLONS)
      integer lat_varid, lon_varid
      real START_LAT, START_LON
      parameter (START_LAT = 25.0, START_LON = -125.0)

C     We will write surface temperature and pressure fields. 
      character*(*) PRES_NAME, TEMP_NAME
      parameter (PRES_NAME='pressure')
      parameter (TEMP_NAME='temperature')
      integer pres_varid, temp_varid
      integer dimids(NDIMS)

C     It's good practice for each variable to carry a "units" attribute.
      character*(*) UNITS
      parameter (UNITS = 'units')
      character*(*) PRES_UNITS, TEMP_UNITS, LAT_UNITS, LON_UNITS
      parameter (PRES_UNITS = 'hPa', TEMP_UNITS = 'celsius')
      parameter (LAT_UNITS = 'degrees_north')
      parameter (LON_UNITS = 'degrees_east')

C     We will create some pressure and temperature data to write out.
      real pres_out(NLONS, NLATS), temp_out(NLONS, NLATS)
      real SAMPLE_PRESSURE
      parameter (SAMPLE_PRESSURE = 900.0)
      real SAMPLE_TEMP
      parameter (SAMPLE_TEMP = 9.0)

C     Loop indices.
      integer lat, lon

C     Error handling.
      integer retval

C     Create pretend data. If this wasn't an example program, we would
C     have some real data to write, for example, model output.
      do lat = 1, NLATS
         lats(lat) = START_LAT + (lat - 1) * 5.0
      end do
      do lon = 1, NLONS
         lons(lon) = START_LON + (lon - 1) * 5.0
      end do
      do lon = 1, NLONS
         do lat = 1, NLATS
            pres_out(lon, lat) = SAMPLE_PRESSURE + 
     +           (lon - 1) * NLATS + (lat - 1)
            temp_out(lon, lat) = SAMPLE_TEMP + 
     +           .25 * ((lon - 1) * NLATS + (lat - 1))
         end do
      end do

C     Create the file. 
      retval = nf_create(FILE_NAME, nf_clobber, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Define the dimensions.
      retval = nf_def_dim(ncid, LAT_NAME, NLATS, lat_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_dim(ncid, LON_NAME, NLONS, lon_dimid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Define the coordinate variables. They will hold the coordinate
C     information, that is, the latitudes and longitudes. A varid is
C     returned for each.
      retval = nf_def_var(ncid, LAT_NAME, NF_REAL, 1, lat_dimid, 
     +     lat_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, LON_NAME, NF_REAL, 1, lon_dimid, 
     +     lon_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Assign units attributes to coordinate var data. This attaches a
C     text attribute to each of the coordinate variables, containing the
C     units.
      retval = nf_put_att_text(ncid, lat_varid, UNITS, len(LAT_UNITS), 
     +     LAT_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, lon_varid, UNITS, len(LON_UNITS), 
     +     LON_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Define the netCDF variables. The dimids array is used to pass the
C     dimids of the dimensions of the netCDF variables.
      dimids(1) = lon_dimid
      dimids(2) = lat_dimid
      retval = nf_def_var(ncid, PRES_NAME, NF_REAL, NDIMS, dimids, 
     +     pres_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_def_var(ncid, TEMP_NAME, NF_REAL, NDIMS, dimids, 
     +     temp_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Assign units attributes to the pressure and temperature netCDF
C     variables.
      retval = nf_put_att_text(ncid, pres_varid, UNITS, len(PRES_UNITS), 
     +     PRES_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_att_text(ncid, temp_varid, UNITS, len(TEMP_UNITS), 
     +     TEMP_UNITS)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     End define mode.
      retval = nf_enddef(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Write the coordinate variable data. This will put the latitudes
C     and longitudes of our data grid into the netCDF file.
      retval = nf_put_var_real(ncid, lat_varid, lats)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, lon_varid, lons)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Write the pretend data. This will write our surface pressure and
C     surface temperature data. The arrays of data are the same size as
C     the netCDF variables we have defined.
      retval = nf_put_var_real(ncid, pres_varid, pres_out)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_put_var_real(ncid, temp_varid, temp_out)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Close the file.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)
   
C     If we got this far, everything worked as expected. Yipee!
      print *,'*** SUCCESS writing example file sfc_pres_temp.nc!'
      end

      subroutine handle_err(errcode)
      implicit none
      include 'netcdf.inc'
      integer errcode

      print *, 'Error: ', nf_strerror(errcode)
      stop 2
      end
c************************************************************************
c https://www.unidata.ucar.edu/software/netcdf/examples/programs/sfc_pres_temp_rd.f

C     This is part of the netCDF package.
C     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
C     See COPYRIGHT file for conditions of use.

C     This is an example which reads some surface pressure and
C     temperatures. The data file read by this program is produced
C     comapnion program sfc_pres_temp_wr.f. It is intended to illustrate
C     the use of the netCDF fortran 77 API.

C     This program is part of the netCDF tutorial:
C     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

C     Full documentation of the netCDF Fortran 77 API can be found at:
C     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77

C     $Id: sfc_pres_temp_rd.f,v 1.8 2007/01/24 19:45:09 russ Exp $

      program sfc_pres_temp_rd
      implicit none
      include 'netcdf.inc'

C     This is the name of the data file we will read.
      character*(*) FILE_NAME
      parameter (FILE_NAME='sfc_pres_temp.nc')
      integer ncid

C     We are reading 2D data, a 6 x 12 lat-lon grid.
      integer NDIMS
      parameter (NDIMS=2)
      integer NLATS, NLONS
      parameter (NLATS = 6, NLONS = 12)
      character*(*) LAT_NAME, LON_NAME
      parameter (LAT_NAME='latitude', LON_NAME='longitude')
      integer lat_dimid, lon_dimid

C     For the lat lon coordinate netCDF variables.
      real lats(NLATS), lons(NLONS)
      integer lat_varid, lon_varid

C     We will read surface temperature and pressure fields. 
      character*(*) PRES_NAME, TEMP_NAME
      parameter (PRES_NAME='pressure')
      parameter (TEMP_NAME='temperature')
      integer pres_varid, temp_varid
      integer dimids(NDIMS)

C     To check the units attributes.
      character*(*) UNITS
      parameter (UNITS = 'units')
      character*(*) PRES_UNITS, TEMP_UNITS, LAT_UNITS, LON_UNITS
      parameter (PRES_UNITS = 'hPa', TEMP_UNITS = 'celsius')
      parameter (LAT_UNITS = 'degrees_north')
      parameter (LON_UNITS = 'degrees_east')
      integer MAX_ATT_LEN
      parameter (MAX_ATT_LEN = 80)
      character*(MAX_ATT_LEN) pres_units_in, temp_units_in
      character*(MAX_ATT_LEN) lat_units_in, lon_units_in
      integer att_len

C     Read the data into these arrays.
      real pres_in(NLONS, NLATS), temp_in(NLONS, NLATS)

C     These are used to calculate the values we expect to find.
      real START_LAT, START_LON
      parameter (START_LAT = 25.0, START_LON = -125.0)
      real SAMPLE_PRESSURE
      parameter (SAMPLE_PRESSURE = 900.0)
      real SAMPLE_TEMP
      parameter (SAMPLE_TEMP = 9.0)

C     We will learn about the data file and store results in these
C     program variables.
      integer ndims_in, nvars_in, ngatts_in, unlimdimid_in

C     Loop indices
      integer lat, lon

C     Error handling
      integer retval

C     Open the file. 
      retval = nf_open(FILE_NAME, nf_nowrite, ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     There are a number of inquiry functions in netCDF which can be
C     used to learn about an unknown netCDF file. NF_INQ tells how many
C     netCDF variables, dimensions, and global attributes are in the
C     file; also the dimension id of the unlimited dimension, if there
C     is one.
      retval = nf_inq(ncid, ndims_in, nvars_in, ngatts_in, 
     +     unlimdimid_in)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     In this case we know that there are 2 netCDF dimensions, 4 netCDF
C     variables, no global attributes, and no unlimited dimension.
      if (ndims_in .ne. 2 .or. nvars_in .ne. 4 .or. ngatts_in .ne. 0 
     +     .or. unlimdimid_in .ne. -1) stop 2

C     Get the varids of the latitude and longitude coordinate variables.
      retval = nf_inq_varid(ncid, LAT_NAME, lat_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_inq_varid(ncid, LON_NAME, lon_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Read the latitude and longitude data.
      retval = nf_get_var_real(ncid, lat_varid, lats)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_var_real(ncid, lon_varid, lons)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Check to make sure we got what we expected.
      do lat = 1, NLATS
         if (lats(lat) .ne. START_LAT + (lat - 1) * 5.0) stop 2
      end do
      do lon = 1, NLONS
         if (lons(lon) .ne. START_LON + (lon - 1) * 5.0) stop 2
      end do

C     Get the varids of the pressure and temperature netCDF variables.
      retval = nf_inq_varid(ncid, PRES_NAME, pres_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_inq_varid(ncid, TEMP_NAME, temp_varid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Read the surface pressure and temperature data from the file.
C     Since we know the contents of the file we know that the data
C     arrays in this program are the correct size to hold all the data.
      retval = nf_get_var_real(ncid, pres_varid, pres_in)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_get_var_real(ncid, temp_varid, temp_in)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     Check the data. It should be the same as the data we wrote.
      do lon = 1, NLONS
         do lat = 1, NLATS
             if (pres_in(lon, lat) .ne. SAMPLE_PRESSURE +
     +           (lon - 1) * NLATS + (lat - 1)) stop 2
             if (temp_in(lon, lat) .ne. SAMPLE_TEMP +
     +           .25 * ((lon - 1) * NLATS + (lat - 1))) stop 2
         end do
      end do

C     Each of the netCDF variables has a "units" attribute. Let's read
C     them and check them.

      retval = nf_get_att_text(ncid, lat_varid, UNITS, lat_units_in)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_inq_attlen(ncid, lat_varid, UNITS, att_len)
      if (retval .ne. nf_noerr) call handle_err(retval)
      if (lat_units_in(1:att_len) .ne. LAT_UNITS) stop 2
 
      retval = nf_get_att_text(ncid, lon_varid, UNITS, lon_units_in)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_inq_attlen(ncid, lon_varid, UNITS, att_len)
      if (retval .ne. nf_noerr) call handle_err(retval)
      if (lon_units_in(1:att_len) .ne. LON_UNITS) stop 2

      retval = nf_get_att_text(ncid, pres_varid, UNITS, pres_units_in)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_inq_attlen(ncid, pres_varid, UNITS, att_len)
      if (retval .ne. nf_noerr) call handle_err(retval)
      if (pres_units_in(1:att_len) .ne. PRES_UNITS) stop 2

      retval = nf_get_att_text(ncid, temp_varid, UNITS, temp_units_in)
      if (retval .ne. nf_noerr) call handle_err(retval)
      retval = nf_inq_attlen(ncid, temp_varid, UNITS, att_len)
      if (retval .ne. nf_noerr) call handle_err(retval)
      if (temp_units_in(1:att_len) .ne. TEMP_UNITS) stop 2

C     Close the file. This frees up any internal netCDF resources
C     associated with the file.
      retval = nf_close(ncid)
      if (retval .ne. nf_noerr) call handle_err(retval)

C     If we got this far, everything worked as expected. Yipee!
      print *,'*** SUCCESS reading example file sfc_pres_temp.nc!'
      end

      subroutine handle_err(errcode)
      implicit none
      include 'netcdf.inc'
      integer errcode

      print *, 'Error: ', nf_strerror(errcode)
      stop 2
      end
c************************************************************************ 
c************************************************************************ 
c************************************************************************ 

! This is part of the netCDF package.
! Copyright 2006 University Corporation for Atmospheric Research/Unidata.
! See COPYRIGHT file for conditions of use.

! This example writes some surface pressure and temperatures. It is
! intended to illustrate the use of the netCDF fortran 90 API. The
! companion program sfc_pres_temp_rd.f shows how to read the netCDF
! data file created by this program.

! This program is part of the netCDF tutorial:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

! Full documentation of the netCDF Fortran 90 API can be found at:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

! $Id: sfc_pres_temp_wr.f90,v 1.9 2007/01/24 19:32:10 russ Exp $

program sfc_pres_temp_wr
  use netcdf
  implicit none

  ! This is the name of the data file we will create.
  character (len = *), parameter :: FILE_NAME = "sfc_pres_temp.nc"
  integer :: ncid

  ! We are writing 2D data, a 6 x 12 lat-lon grid. We will need two
  ! netCDF dimensions.
  integer, parameter :: NDIMS = 2
  integer, parameter :: NLATS = 6, NLONS = 12
  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_NAME = "longitude"
  integer :: lat_dimid, lon_dimid

  ! In addition to the latitude and longitude dimensions, we will also
  ! create latitude and longitude netCDF variables which will hold the
  ! actual latitudes and longitudes. Since they hold data about the
  ! coordinate system, the netCDF term for these is: "coordinate
  ! variables."
  real :: lats(NLATS), lons(NLONS)
  integer :: lat_varid, lon_varid
  real, parameter :: START_LAT = 25.0, START_LON = -125.0

  ! We will write surface temperature and pressure fields. 
  character (len = *), parameter :: PRES_NAME="pressure"
  character (len = *), parameter :: TEMP_NAME="temperature"
  integer :: pres_varid, temp_varid
  integer :: dimids(NDIMS)

  ! It's good practice for each variable to carry a "units" attribute.
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: PRES_UNITS = "hPa"
  character (len = *), parameter :: TEMP_UNITS = "celsius"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"

  ! We will create some pressure and temperature data to write out.
  real :: pres_out(NLONS, NLATS), temp_out(NLONS, NLATS)
  real, parameter :: SAMPLE_PRESSURE = 900.0
  real, parameter :: SAMPLE_TEMP = 9.0

  ! Loop indices
  integer :: lat, lon

  ! Create pretend data. If this wasn't an example program, we would
  ! have some real data to write, for example, model output.
  do lat = 1, NLATS
     lats(lat) = START_LAT + (lat - 1) * 5.0
  end do
  do lon = 1, NLONS
     lons(lon) = START_LON + (lon - 1) * 5.0
  end do
  do lon = 1, NLONS
     do lat = 1, NLATS
        pres_out(lon, lat) = SAMPLE_PRESSURE + (lon - 1) * NLATS + (lat - 1)
        temp_out(lon, lat) = SAMPLE_TEMP + .25 * ((lon - 1) * NLATS + (lat - 1))
     end do
  end do

  ! Create the file. 
  call check( nf90_create(FILE_NAME, nf90_clobber, ncid) )

  ! Define the dimensions.
  call check( nf90_def_dim(ncid, LAT_NAME, NLATS, lat_dimid) )
  call check( nf90_def_dim(ncid, LON_NAME, NLONS, lon_dimid) )

  ! Define the coordinate variables. They will hold the coordinate
  ! information, that is, the latitudes and longitudes. A varid is
  ! returned for each.
  call check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
  call check( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )

  ! Assign units attributes to coordinate var data. This attaches a
  ! text attribute to each of the coordinate variables, containing the
  ! units.
  call check( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
  call check( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )

  ! Define the netCDF variables. The dimids array is used to pass the
  ! dimids of the dimensions of the netCDF variables.
  dimids = (/ lon_dimid, lat_dimid /)
  call check( nf90_def_var(ncid, PRES_NAME, NF90_REAL, dimids, pres_varid) )
  call check( nf90_def_var(ncid, TEMP_NAME, NF90_REAL, dimids, temp_varid) )

  ! Assign units attributes to the pressure and temperature netCDF
  ! variables.
  call check( nf90_put_att(ncid, pres_varid, UNITS, PRES_UNITS) )
  call check( nf90_put_att(ncid, temp_varid, UNITS, TEMP_UNITS) )

  ! End define mode.
  call check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. This will put the latitudes
  ! and longitudes of our data grid into the netCDF file.
  call check( nf90_put_var(ncid, lat_varid, lats) )
  call check( nf90_put_var(ncid, lon_varid, lons) )

  ! Write the pretend data. This will write our surface pressure and
  ! surface temperature data. The arrays of data are the same size as
  ! the netCDF variables we have defined.
  call check( nf90_put_var(ncid, pres_varid, pres_out) )
  call check( nf90_put_var(ncid, temp_varid, temp_out) )

  ! Close the file.
  call check( nf90_close(ncid) )
   
  ! If we got this far, everything worked as expected. Yipee! 
  print *,"*** SUCCESS writing example file sfc_pres_temp.nc!"

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  
end program sfc_pres_temp_wr
c************************************************************************
! This is part of the netCDF package.
! Copyright 2006 University Corporation for Atmospheric Research/Unidata.
! See COPYRIGHT file for conditions of use.

! This example writes some surface pressure and temperatures. It is
! intended to illustrate the use of the netCDF fortran 90 API. The
! companion program sfc_pres_temp_rd.f shows how to read the netCDF
! data file created by this program.

! This program is part of the netCDF tutorial:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

! Full documentation of the netCDF Fortran 90 API can be found at:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

! $Id: sfc_pres_temp_wr.f90,v 1.9 2007/01/24 19:32:10 russ Exp $

program sfc_pres_temp_wr
  use netcdf
  implicit none

  ! This is the name of the data file we will create.
  character (len = *), parameter :: FILE_NAME = "sfc_pres_temp.nc"
  integer :: ncid

  ! We are writing 2D data, a 6 x 12 lat-lon grid. We will need two
  ! netCDF dimensions.
  integer, parameter :: NDIMS = 2
  integer, parameter :: NLATS = 6, NLONS = 12
  character (len = *), parameter :: LAT_NAME = "latitude"
  character (len = *), parameter :: LON_NAME = "longitude"
  integer :: lat_dimid, lon_dimid

  ! In addition to the latitude and longitude dimensions, we will also
  ! create latitude and longitude netCDF variables which will hold the
  ! actual latitudes and longitudes. Since they hold data about the
  ! coordinate system, the netCDF term for these is: "coordinate
  ! variables."
  real :: lats(NLATS), lons(NLONS)
  integer :: lat_varid, lon_varid
  real, parameter :: START_LAT = 25.0, START_LON = -125.0

  ! We will write surface temperature and pressure fields. 
  character (len = *), parameter :: PRES_NAME="pressure"
  character (len = *), parameter :: TEMP_NAME="temperature"
  integer :: pres_varid, temp_varid
  integer :: dimids(NDIMS)

  ! It's good practice for each variable to carry a "units" attribute.
  character (len = *), parameter :: UNITS = "units"
  character (len = *), parameter :: PRES_UNITS = "hPa"
  character (len = *), parameter :: TEMP_UNITS = "celsius"
  character (len = *), parameter :: LAT_UNITS = "degrees_north"
  character (len = *), parameter :: LON_UNITS = "degrees_east"

  ! We will create some pressure and temperature data to write out.
  real :: pres_out(NLONS, NLATS), temp_out(NLONS, NLATS)
  real, parameter :: SAMPLE_PRESSURE = 900.0
  real, parameter :: SAMPLE_TEMP = 9.0

  ! Loop indices
  integer :: lat, lon

  ! Create pretend data. If this wasn't an example program, we would
  ! have some real data to write, for example, model output.
  do lat = 1, NLATS
     lats(lat) = START_LAT + (lat - 1) * 5.0
  end do
  do lon = 1, NLONS
     lons(lon) = START_LON + (lon - 1) * 5.0
  end do
  do lon = 1, NLONS
     do lat = 1, NLATS
        pres_out(lon, lat) = SAMPLE_PRESSURE + (lon - 1) * NLATS + (lat - 1)
        temp_out(lon, lat) = SAMPLE_TEMP + .25 * ((lon - 1) * NLATS + (lat - 1))
     end do
  end do

  ! Create the file. 
  call check( nf90_create(FILE_NAME, nf90_clobber, ncid) )

  ! Define the dimensions.
  call check( nf90_def_dim(ncid, LAT_NAME, NLATS, lat_dimid) )
  call check( nf90_def_dim(ncid, LON_NAME, NLONS, lon_dimid) )

  ! Define the coordinate variables. They will hold the coordinate
  ! information, that is, the latitudes and longitudes. A varid is
  ! returned for each.
  call check( nf90_def_var(ncid, LAT_NAME, NF90_REAL, lat_dimid, lat_varid) )
  call check( nf90_def_var(ncid, LON_NAME, NF90_REAL, lon_dimid, lon_varid) )

  ! Assign units attributes to coordinate var data. This attaches a
  ! text attribute to each of the coordinate variables, containing the
  ! units.
  call check( nf90_put_att(ncid, lat_varid, UNITS, LAT_UNITS) )
  call check( nf90_put_att(ncid, lon_varid, UNITS, LON_UNITS) )

  ! Define the netCDF variables. The dimids array is used to pass the
  ! dimids of the dimensions of the netCDF variables.
  dimids = (/ lon_dimid, lat_dimid /)
  call check( nf90_def_var(ncid, PRES_NAME, NF90_REAL, dimids, pres_varid) )
  call check( nf90_def_var(ncid, TEMP_NAME, NF90_REAL, dimids, temp_varid) )

  ! Assign units attributes to the pressure and temperature netCDF
  ! variables.
  call check( nf90_put_att(ncid, pres_varid, UNITS, PRES_UNITS) )
  call check( nf90_put_att(ncid, temp_varid, UNITS, TEMP_UNITS) )

  ! End define mode.
  call check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. This will put the latitudes
  ! and longitudes of our data grid into the netCDF file.
  call check( nf90_put_var(ncid, lat_varid, lats) )
  call check( nf90_put_var(ncid, lon_varid, lons) )

  ! Write the pretend data. This will write our surface pressure and
  ! surface temperature data. The arrays of data are the same size as
  ! the netCDF variables we have defined.
  call check( nf90_put_var(ncid, pres_varid, pres_out) )
  call check( nf90_put_var(ncid, temp_varid, temp_out) )

  ! Close the file.
  call check( nf90_close(ncid) )
   
  ! If we got this far, everything worked as expected. Yipee! 
  print *,"*** SUCCESS writing example file sfc_pres_temp.nc!"

contains
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check  
end program sfc_pres_temp_wr

************************************************************************

https://web.lmd.jussieu.fr/~lmdz/Distrib/Chgt_install/LMDZ20210529.trunk/netcdf-4.0.1/man4/netcdf-tutorial.html#simple_005fxy_005fwr_002ef
https://web.lmd.jussieu.fr/~lmdz/Distrib/Chgt_install/LMDZ20210529.trunk/netcdf-4.0.1/man4/netcdf-tutorial.html#simple_005fxy_005fwr_002ef
https://web.lmd.jussieu.fr/~lmdz/Distrib/Chgt_install/LMDZ20210529.trunk/netcdf-4.0.1/man4/netcdf-tutorial.html#simple_005fxy_005fwr_002ef

2.1.2.1 simple_xy_wr.f
     C     This is part of the netCDF package.
     C     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
     C     See COPYRIGHT file for conditions of use.
     
     C     This is a very simple example which writes a 2D array of
     C     sample data. To handle this in netCDF we create two shared
     C     dimensions, "x" and "y", and a netCDF variable, called "data".
     
     C     This example demonstrates the netCDF Fortran 77 API. This is part
     C     of the netCDF tutorial, which can be found at:
     C     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
     
     C     Full documentation of the netCDF Fortran 77 API can be found at:
     C     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77
     
     C     $Id: simple_xy_wr.f,v 1.11 2008/08/20 22:29:56 russ Exp $
     
           program simple_xy_wr
           implicit none
           include 'netcdf.inc'
     
     C     This is the name of the data file we will create.
           character*(*) FILE_NAME
           parameter (FILE_NAME='simple_xy.nc')
     
     C     We are writing 2D data, a 12 x 6 grid.
           integer NDIMS
           parameter (NDIMS=2)
           integer NX, NY
           parameter (NX = 6, NY = 12)
     
     C     When we create netCDF files, variables and dimensions, we get back
     C     an ID for each one.
           integer ncid, varid, dimids(NDIMS)
           integer x_dimid, y_dimid
     
     C     This is the data array we will write. It will just be filled with
     C     a progression of integers for this example.
           integer data_out(NY, NX)
     
     C     Loop indexes, and error handling.
           integer x, y, retval
     
     C     Create some pretend data. If this wasn't an example program, we
     C     would have some real data to write, for example, model output.
           do x = 1, NX
              do y = 1, NY
                 data_out(y, x) = (x - 1) * NY + (y - 1)
              end do
           end do
     
     C     Always check the return code of every netCDF function call. In
     C     this example program, any retval which is not equal to nf_noerr
     C     (0) will call handle_err, which prints a netCDF error message, and
     C     then exits with a non-zero return code.
     
     C     Create the netCDF file. The nf_clobber parameter tells netCDF to
     C     overwrite this file, if it already exists.
           retval = nf_create(FILE_NAME, NF_CLOBBER, ncid)
           if (retval .ne. nf_noerr) call handle_err(retval)
     
     C     Define the dimensions. NetCDF will hand back an ID for each.
           retval = nf_def_dim(ncid, "x", NX, x_dimid)
           if (retval .ne. nf_noerr) call handle_err(retval)
           retval = nf_def_dim(ncid, "y", NY, y_dimid)
           if (retval .ne. nf_noerr) call handle_err(retval)
     
     C     The dimids array is used to pass the IDs of the dimensions of
     C     the variables. Note that in fortran arrays are stored in
     C     column-major format.
           dimids(2) = x_dimid
           dimids(1) = y_dimid
     
     C     Define the variable. The type of the variable in this case is
     C     NF_INT (4-byte integer).
           retval = nf_def_var(ncid, "data", NF_INT, NDIMS, dimids, varid)
           if (retval .ne. nf_noerr) call handle_err(retval)
     
     C     End define mode. This tells netCDF we are done defining metadata.
           retval = nf_enddef(ncid)
           if (retval .ne. nf_noerr) call handle_err(retval)
     
     C     Write the pretend data to the file. Although netCDF supports
     C     reading and writing subsets of data, in this case we write all the
     C     data in one operation.
           retval = nf_put_var_int(ncid, varid, data_out)
           if (retval .ne. nf_noerr) call handle_err(retval)
     
     C     Close the file. This frees up any internal netCDF resources
     C     associated with the file, and flushes any buffers.
           retval = nf_close(ncid)
           if (retval .ne. nf_noerr) call handle_err(retval)
     
           print *,'*** SUCCESS writing example file simple_xy.nc!'
           end
     
           subroutine handle_err(errcode)
           implicit none
           include 'netcdf.inc'
           integer errcode
     
           print *, 'Error: ', nf_strerror(errcode)
           stop 2
           end
c************************************************************************
https://web.lmd.jussieu.fr/~lmdz/Distrib/Chgt_install/LMDZ20210529.trunk/netcdf-4.0.1/man4/netcdf-tutorial.html#simple_005fxy_005fwr_002ef
https://web.lmd.jussieu.fr/~lmdz/Distrib/Chgt_install/LMDZ20210529.trunk/netcdf-4.0.1/man4/netcdf-tutorial.html#simple_005fxy_005fwr_002ef
https://web.lmd.jussieu.fr/~lmdz/Distrib/Chgt_install/LMDZ20210529.trunk/netcdf-4.0.1/man4/netcdf-tutorial.html#simple_005fxy_005fwr_002ef


2.1.2.2 simple_xy_rd.f
     C     This is part of the netCDF package.
     C     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
     C     See COPYRIGHT file for conditions of use.
     
     C     This is a simple example which reads a small dummy array, from a
     C     netCDF data file created by the companion program simple_xy_wr.f.
     
     C     This is intended to illustrate the use of the netCDF fortran 77
     C     API. This example program is part of the netCDF tutorial, which can
     C     be found at:
     C     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
     
     C     Full documentation of the netCDF Fortran 77 API can be found at:
     C     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f77
     
     C     $Id: simple_xy_rd.f,v 1.8 2007/02/14 20:59:20 ed Exp $
     
           program simple_xy_rd
           implicit none
           include 'netcdf.inc'
     
     C     This is the name of the data file we will read.
           character*(*) FILE_NAME
           parameter (FILE_NAME='simple_xy.nc')
     
     C     We are reading 2D data, a 12 x 6 grid.
           integer NX, NY
           parameter (NX = 6, NY = 12)
           integer data_in(NY, NX)
     
     C     This will be the netCDF ID for the file and data variable.
           integer ncid, varid
     
     C     Loop indexes, and error handling.
           integer x, y, retval
     
     C     Open the file. NF_NOWRITE tells netCDF we want read-only access to
     C     the file.
           retval = nf_open(FILE_NAME, NF_NOWRITE, ncid)
           if (retval .ne. nf_noerr) call handle_err(retval)
     
     C     Get the varid of the data variable, based on its name.
           retval = nf_inq_varid(ncid, 'data', varid)
           if (retval .ne. nf_noerr) call handle_err(retval)
     
     C     Read the data.
           retval = nf_get_var_int(ncid, varid, data_in)
           if (retval .ne. nf_noerr) call handle_err(retval)
     
     C     Check the data.
           do x = 1, NX
              do y = 1, NY
                 if (data_in(y, x) .ne. (x - 1) * NY + (y - 1)) then
                    print *, 'data_in(', y, ', ', x, ') = ', data_in(y, x)
                    stop 2
                 end if
              end do
           end do
     
     C     Close the file, freeing all resources.
           retval = nf_close(ncid)
           if (retval .ne. nf_noerr) call handle_err(retval)
     
           print *,'*** SUCCESS reading example file ', FILE_NAME, '!'
           end
     
           subroutine handle_err(errcode)
           implicit none
           include 'netcdf.inc'
           integer errcode
     
           print *, 'Error: ', nf_strerror(errcode)
           stop 2
           end

c************************************************************************
Previous: simple_xy_wr.f90, Up: simple_xy in F90
2.1.3.2 simple_xy_rd.f90
     ! This is part of the netCDF package.
     ! Copyright 2006 University Corporation for Atmospheric Research/Unidata.
     ! See COPYRIGHT file for conditions of use.
     
     ! This is a simple example which reads a small dummy array, from a
     ! netCDF data file created by the companion program simple_xy_wr.f90.
     
     ! This is intended to illustrate the use of the netCDF fortran 90
     ! API. This example program is part of the netCDF tutorial, which can
     ! be found at:
     ! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
     
     ! Full documentation of the netCDF Fortran 90 API can be found at:
     ! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90
     
     ! $Id: simple_xy_rd.f90,v 1.10 2009/02/25 21:44:07 ed Exp $
     
     program simple_xy_rd
       use netcdf
       implicit none
     
       ! This is the name of the data file we will read.
       character (len = *), parameter :: FILE_NAME = "simple_xy.nc"
     
       ! We are reading 2D data, a 12 x 6 grid.
       integer, parameter :: NX = 6, NY = 12
       integer :: data_in(NY, NX)
     
       ! This will be the netCDF ID for the file and data variable.
       integer :: ncid, varid
     
       ! Loop indexes, and error handling.
       integer :: x, y
     
       ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
       ! the file.
       call check( nf90_open(FILE_NAME, NF90_NOWRITE, ncid) )
     
       ! Get the varid of the data variable, based on its name.
       call check( nf90_inq_varid(ncid, "data", varid) )
     
       ! Read the data.
       call check( nf90_get_var(ncid, varid, data_in) )
     
       ! Check the data.
       do x = 1, NX
          do y = 1, NY
             if (data_in(y, x) /= (x - 1) * NY + (y - 1)) then
                print *, "data_in(", y, ", ", x, ") = ", data_in(y, x)
                stop "Stopped"
             end if
          end do
       end do
     
       ! Close the file, freeing all resources.
       call check( nf90_close(ncid) )
     
       print *,"*** SUCCESS reading example file ", FILE_NAME, "! "
     
     contains
       subroutine check(status)
         integer, intent ( in) :: status
     
         if(status /= nf90_noerr) then
           print *, trim(nf90_strerror(status))
           stop 2
         end if
       end subroutine check
     end program simple_xy_rd

c************************************************************************
Next: simple_xy_rd.f90, Previous: simple_xy in F90, Up: simple_xy in F90
2.1.3.1 simple_xy_wr.f90
     !     This is part of the netCDF package.
     !     Copyright 2006 University Corporation for Atmospheric Research/Unidata.
     !     See COPYRIGHT file for conditions of use.
     
     !     This is a very simple example which writes a 2D array of
     !     sample data. To handle this in netCDF we create two shared
     !     dimensions, "x" and "y", and a netCDF variable, called "data".
     
     !     This example demonstrates the netCDF Fortran 90 API. This is part
     !     of the netCDF tutorial, which can be found at:
     !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial
     
     !     Full documentation of the netCDF Fortran 90 API can be found at:
     !     http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90
     
     !     $Id: simple_xy_wr.f90,v 1.10 2008/08/20 22:43:53 russ Exp $
     
     program simple_xy_wr
       use netcdf
       implicit none
     
       ! This is the name of the data file we will create.
       character (len = *), parameter :: FILE_NAME = "simple_xy.nc"
     
       ! We are writing 2D data, a 12 x 6 grid.
       integer, parameter :: NDIMS = 2
       integer, parameter :: NX = 6, NY = 12
     
       ! When we create netCDF files, variables and dimensions, we get back
       ! an ID for each one.
       integer :: ncid, varid, dimids(NDIMS)
       integer :: x_dimid, y_dimid
     
       ! This is the data array we will write. It will just be filled with
       ! a progression of integers for this example.
       integer :: data_out(NY, NX)
     
       ! Loop indexes, and error handling.
       integer :: x, y
     
       ! Create some pretend data. If this wasn't an example program, we
       ! would have some real data to write, for example, model output.
       do x = 1, NX
          do y = 1, NY
             data_out(y, x) = (x - 1) * NY + (y - 1)
          end do
       end do
     
       ! Always check the return code of every netCDF function call. In
       ! this example program, wrapping netCDF calls with "call check()"
       ! makes sure that any return which is not equal to nf90_noerr (0)
       ! will print a netCDF error message and exit.
     
       ! Create the netCDF file. The nf90_clobber parameter tells netCDF to
       ! overwrite this file, if it already exists.
       call check( nf90_create(FILE_NAME, NF90_CLOBBER, ncid) )
     
       ! Define the dimensions. NetCDF will hand back an ID for each.
       call check( nf90_def_dim(ncid, "x", NX, x_dimid) )
       call check( nf90_def_dim(ncid, "y", NY, y_dimid) )
     
       ! The dimids array is used to pass the IDs of the dimensions of
       ! the variables. Note that in fortran arrays are stored in
       ! column-major format.
       dimids =  (/ y_dimid, x_dimid /)
     
       ! Define the variable. The type of the variable in this case is
       ! NF90_INT (4-byte integer).
       call check( nf90_def_var(ncid, "data", NF90_INT, dimids, varid) )
     
       ! End define mode. This tells netCDF we are done defining metadata.
       call check( nf90_enddef(ncid) )
     
       ! Write the pretend data to the file. Although netCDF supports
       ! reading and writing subsets of data, in this case we write all the
       ! data in one operation.
       call check( nf90_put_var(ncid, varid, data_out) )
     
       ! Close the file. This frees up any internal netCDF resources
       ! associated with the file, and flushes any buffers.
       call check( nf90_close(ncid) )
     
       print *, "*** SUCCESS writing example file simple_xy.nc! "
     
     contains
       subroutine check(status)
         integer, intent ( in) :: status
     
         if(status /= nf90_noerr) then
           print *, trim(nf90_strerror(status))
           stop 2
         end if
       end subroutine check
     end program simple_xy_wr
