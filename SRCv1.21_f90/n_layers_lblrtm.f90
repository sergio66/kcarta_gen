! Copyright 2014
! University of Maryland Baltimore County
! All Rights Reserved

MODULE n_layers_lblrtm

USE basic_common
USE spline_and_sort_and_common
USE s_misc
USE freqfile
USE n_layers
USE n_mr_common

IMPLICIT NONE

CONTAINS

! this file does the LBLRTM TAPE5/TAPE6

! see eg http://shadow.eas.gatech.edu/~vvt/lblrtm/lblrtm_inst.html
!       http://www.ssec.wisc.edu/~paulv/Fortran90/AtmProfile/Modules.html
!       https://svn.ssec.wisc.edu/repos/uwphysret/trunk/mfiles/compute_F.m
!       JCHAR = 1-6           - default to value for specified model atmosphere
!              = " ",A         - volume mixing ratio (ppmv)
!              = B             - number density (cm-3)
!              = C             - mass mixing ratio (gm/kg)
!              = D             - mass density (gm m-3)
!              = E             - partial pressure (mb)
!              = F             - dew point temp (K) *H2O only*
!              = G             - dew point temp (C) *H2O only*
!              = H             - relative humidity (percent) *H2O only*
!              = I             - available for user definition

!      LBLRTM     iaGasUnits(iG) = 12   !!! assume hardcoded VMR
!      LBLRTM     pressures are in mb, which is what klayers wants

!      rP(min/max)KCarta is in N/m2 which is x100 what it would be in mb

! ESw.d ENw.d ESw.dEe ENw.dEe output format see http://www.enautica.pt/publico/professores/vfranco/formats.pdf


!************************************************************************
END MODULE n_layers_lblrtm
