PROGRAM get_varname_from_table

  ! Replaces in a given netcdf file all variables with varname "var<code>", 
  ! where<code> is the GRIB code number, by the "correct" varname. Besides, the
  ! attributes "long_name" and "units are created. The values of varname and 
  ! "longname", and "units" are looked up in the table <table>, where <table>
  ! is the table number of the variable.
  !
  ! This program should be called after a call of the ECHAM5 afterburner program
  ! *after* for variables with a table number not equal to 128 (= Local code 
  ! table for Standard ECHAM; "correct" varname and attribute long_varname is 
  ! determinedby *after* in that case.)
  !
  ! Input: 
  !
  ! $1 : Input file name (netcdf file to be transformed)
  !
  ! M. Schraner, ETH Zurich, 10.3.2009

  USE netcdf

  IMPLICIT NONE

  INTEGER :: i, status, ncid, nvartot, table, code
  CHARACTER(31) :: ncname, varname, varname_new, gridtype
  CHARACTER(80) :: longname, units

  ! Executable statements:

  CALL GETARG(1, ncname)

  ! Open netcdf file:   
  status=NF90_OPEN(TRIM(ncname), NF90_WRITE, ncid)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! Inquire number of variables:
  status=NF90_INQUIRE(ncid, nvariables=nvartot)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! Enter define mode:
  status = NF90_REDEF(ncid)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  DO i=1, nvartot

     ! Get variable name:
     status=NF90_INQUIRE_VARIABLE(ncid, i, varname)
     IF (status /= NF90_NOERR) CALL handle_err(status)

     ! Replace variable name in case that it begins with "var":
     IF (varname(1:3) .EQ. 'var') THEN

        ! Get code number from variable name:
        IF (LEN_TRIM(varname) .EQ. 5) THEN
           READ(varname(4:5),'(i2.2)') code
        ELSE
           READ(varname(4:6),'(i3.3)') code
        ENDIF

        ! Get table number and grid type:
        status=NF90_GET_ATT(ncid, i, 'table', table)
        IF (status /= NF90_NOERR) CALL handle_err(status)
        status=NF90_GET_ATT(ncid, i, 'grid_type', gridtype)
        IF (status /= NF90_NOERR) CALL handle_err(status)

        ! Get varname, long_varname, and units:
        CALL get_varname_long_varname_units

        ! Rename variable name:
        status = NF90_RENAME_VAR(ncid, i, TRIM(varname_new))
        IF (status /= NF90_NOERR) CALL handle_err(status)

        ! Delete attributes "table" and "grid type" (re-set below):
        status = NF90_DEL_ATT(ncid, i, 'table')
        IF (status /= NF90_NOERR) CALL handle_err(status)
        status = NF90_DEL_ATT(ncid, i, 'grid_type')
        IF (status /= NF90_NOERR) CALL handle_err(status)

        ! Set attributes "long_name" and "units":
        status=NF90_PUT_ATT(ncid, i, 'long_name', TRIM(longname))
        IF (status /= NF90_NOERR) CALL handle_err(status)
        status=NF90_PUT_ATT(ncid, i, 'units', TRIM(units))
        IF (status /= NF90_NOERR) CALL handle_err(status)

        ! Set attributes "code" and "table":
        status=NF90_PUT_ATT(ncid, i, 'code', code)
        IF (status /= NF90_NOERR) CALL handle_err(status)
        status=NF90_PUT_ATT(ncid, i, 'table', table)
        IF (status /= NF90_NOERR) CALL handle_err(status)
        status=NF90_PUT_ATT(ncid, i, 'grid_type', TRIM(gridtype))
        IF (status /= NF90_NOERR) CALL handle_err(status)

     ENDIF

  END DO

  ! End define mode:
  status=NF90_ENDDEF(ncid)
  IF (status /= NF90_NOERR) CALL handle_err(status)

  ! Close file:
  status = NF90_CLOSE(ncid)
  IF (status /= NF90_NOERR) CALL HANDLE_ERR(status)

CONTAINS

  SUBROUTINE get_varname_long_varname_units
    
    SELECT CASE (table)

    ! Codes in table 131 and 199 (chemistry table):
    CASE(131, 199)
       
       SELECT CASE (code)

       CASE(11)
          varname_new='O3'
          longname='O3 volume mixing ratio'
          units='mol/mol'
       CASE(12)
          varname_new='O'
          longname='O volume mixing ratio'
          units='mol/mol'
       CASE(13)
          varname_new='NO'
          longname='NO volume mixing ratio'
          units='mol/mol'
       CASE(14)
          varname_new='NO2'
          longname='NO2 volume mixing ratio'
          units='mol/mol'
       CASE(15)
          varname_new='HNO3'
          longname='HNO3 volume mixing ratio'
          units='mol/mol'
          
       CASE(16)
          varname_new='NO3'
          longname='NO3 volume mixing ratio'
          units='mol/mol'
       CASE(17)
          varname_new='N2O5'
          longname='N2O5 volume mixing ratio'
          units='mol/mol'
       CASE(18)
          varname_new='HNO4'
          longname='HNO4 volume mixing ratio'
          units='mol/mol'
       CASE(19)
          varname_new='HO2'
          longname='HO2 volume mixing ratio'
          units='mol/mol'
       CASE(20)
          varname_new='ClNO3'
          longname='ClNO3 volume mixing ratio'
          units='mol/mol'
 
       CASE(21)
          varname_new='ClO'
          longname='ClO volume mixing ratio'
          units='mol/mol'
       CASE(22)
          varname_new='N'
          longname='N volume mixing ratio'
          units='mol/mol'
       CASE(23)
          varname_new='OH'
          longname='OH volume mixing ratio'
          units='mol/mol'
       CASE(24)
          varname_new='H'
          longname='H volume mixing ratio'
          units='mol/mol'
       CASE(25)
          varname_new='Cl'
          longname='Cl volume mixing ratio'
          units='mol/mol'

       CASE(26)
          varname_new='HOCl'
          longname='HOCl volume mixing ratio'
          units='mol/mol'
       CASE(27)
          varname_new='ODSCLS'
          longname='ODSCLS volume mixing ratio'
          units='mol/mol'
       CASE(28)
          varname_new='ODSCLL'
          longname='ODSCLL volume mixing ratio'
          units='mol/mol'
       CASE(29)
          varname_new='N2O'
          longname='N2O volume mixing ratio'
          units='mol/mol'
       CASE(30)
          varname_new='CH4'
          longname='CH4 volume mixing ratio'
          units='mol/mol'
 
       CASE(31)
          varname_new='CO'
          longname='CO volume mixing ratio'
          units='mol/mol'
       CASE(32)
          varname_new='HCl'
          longname='HCl volume mixing ratio'
          units='mol/mol'
       CASE(33)
          varname_new='H2'
          longname='H2 volume mixing ratio'
          units='mol/mol'
       CASE(34)
          varname_new='H2O2'
          longname='H2O2 volume mixing ratio'
          units='mol/mol'
       CASE(35)
          varname_new='H2O'
          longname='H2O volume mixing ratio'
          units='mol/mol'

       CASE(36)
          varname_new='Cl2'
          longname='Cl2 volume mixing ratio'
          units='mol/mol'
       CASE(37)
          varname_new='Cl2O2'
          longname='Cl2O2 volume mixing ratio'
          units='mol/mol'
       CASE(38)
          varname_new='PSC1'
          longname='PSC1 volume mixing ratio'
          units='mol/mol'
       CASE(39)
          varname_new='PSC2'
          longname='PSC2 volume mixing ratio'
          units='mol/mol'
       CASE(40)
          varname_new='CH3O2H'
          longname='CH3O2H volume mixing ratio'
          units='mol/mol'
 
       CASE(41)
          varname_new='CH2O'
          longname='CH2O volume mixing ratio'
          units='mol/mol'
       CASE(42)
          varname_new='BrO'
          longname='BrO volume mixing ratio'
          units='mol/mol'
       CASE(43)
          varname_new='BrNO3'
          longname='BrNO3 volume mixing ratio'
          units='mol/mol'
       CASE(44)
          varname_new='BrCl'
          longname='BrCl volume mixing ratio'
          units='mol/mol'
       CASE(45)
          varname_new='HBr'
          longname='HBr volume mixing ratio'
          units='mol/mol'

       CASE(46)
          varname_new='HOBr'
          longname='HOBr volume mixing ratio'
          units='mol/mol'
       CASE(47)
          varname_new='ODSBR'
          longname='ODSBR volume mixing ratio'
          units='mol/mol'
       CASE(48)
          varname_new='CH3'
          longname='CH3 volume mixing ratio'
          units='mol/mol'
       CASE(49)
          varname_new='CH3O2'
          longname='CH3O2 volume mixing ratio'
          units='mol/mol'
       CASE(50)
          varname_new='CH3O'
          longname='CH3O volume mixing ratio'
          units='mol/mol'
 
       CASE(51)
          varname_new='HCO'
          longname='HCO volume mixing ratio'
          units='mol/mol'
       CASE(52)
          varname_new='Br'
          longname='Br volume mixing ratio'
          units='mol/mol'
       CASE(53)
          varname_new='OD'
          longname='OD volume mixing ratio'
          units='mol/mol'
       CASE(54)
          varname_new='TOTOZ'
          longname='Ozone column'
          units='DU'
       CASE(55)
          varname_new='TOTNO2'
          longname='NO2 Column'
          units='molec/cm2'

       CASE(56)
          varname_new='ClOx'
          longname='ClOx volume mixing ratio'
          units='mol/mol'
       CASE(57)
          varname_new='Cly'
          longname='Cly volume mixing ratio'
          units='mol/mol'
       CASE(58)
          varname_new='CCly'
          longname='CCly volume mixing ratio'
          units='mol/mol'
       CASE(59)
          varname_new='Bry'
          longname='Bry volume mixing ratio'
          units='mol/mol'
       CASE(60)
          varname_new='CBry'
          longname='CBry volume mixing ratio'
          units='mol/mol'
 
       CASE(61)
          varname_new='NOy'
          longname='NOy volume mixing ratio'
          units='mol/mol'
       CASE(62)
          varname_new='NOx'
          longname='NOx volume mixing ratio'
          units='mol/mol'
 
       CASE(95)
          varname_new='FAMFIXCL'
          longname='Cly-family fixer'
          units='-'
       CASE(96)
          varname_new='FAMFIXBR'
          longname='Bry-family fixer'
          units='-'
       CASE(97)
          varname_new='FAMFIXN'
          longname='NOy-family fixer'
          units='-'
 
       CASE(98)
          varname_new='SEDPSC1'
          longname='PSC1 sedimentation: HNO3 loss/gain'
          units='mol/cm3/s'
       CASE(99)
          varname_new='SEDPSC2'
          longname='PSC2 sedimentation: H2O loss/gain'
          units='mol/cm3/s'

       CASE(211)
          varname_new='RXN_1'
          longname='Ozone-depleting NOx cycles'
          units='mol/cm3/s'
       CASE(212)
          varname_new='RXN_2'
          longname='Ozone-depleting HOx cycles'
          units='mol/cm3/s'
       CASE(213)
          varname_new='RXN_3'
          longname='Ozone-depleting Cl cycles'
          units='mol/cm3/s'
       CASE(214)
          varname_new='RXN_4'
          longname='Ozone-depleting Br cycles'
          units='mol/cm3/s'
       CASE(215)
          varname_new='RXN_5'
          longname='Net chemical ozone production'
          units='mol/cm3/s'

       CASE(216)
          varname_new='RXN_6'
          longname='Sum of the HO2+NO and RO2+NO reactions'
          units='mol/cm3/s'
       CASE(217)
          varname_new='RXN_7'
          longname='Sum of OD+H2O, O3+HO2, O3+OH and C5H8+O3'
          units='mol/cm3/s'
       CASE(218)
          varname_new='RXN_8'
          longname='NO2+hv->NO+O'
          units='mol/cm3/s'
       CASE(219)
          varname_new='RXN_9'
          longname='OD+H2O->OH+OH'
          units='mol/cm3/s'
       CASE(220)
          varname_new='RXN_10'
          longname='Total OH loss'
          units='mol/cm3/s'

       CASE(240)
          varname_new='RXN_11'
          longname='CO+OH+M->H+CO2+M'
          units='mol/cm3/s'
       CASE(241)
          varname_new='RXN_12'
          longname='CH4+OH->CH3+H2O'
          units='mol/cm3/s'
       CASE(242)
          varname_new='RXN_13'
          longname='CH4+CL->CH3+HCL'
          units='mol/cm3/s'
       CASE(243)
          varname_new='RXN_14'
          longname='H2O2 production'
          units='mol/cm3/s'
       CASE(244)
          varname_new='RXN_15'
          longname='Production of all hydrogen peroxides'
          units='mol/cm3/s'

       CASE(245)
          varname_new='RXN_16'
          longname='HNO3 production'
          units='mol/cm3/s'
       CASE(246)
          varname_new='RXN_17'
          longname='Sum of RO2+NO reactions'
          units='mol/cm3/s'
       CASE(247)
          varname_new='RXN_18'
          longname='Sum of RO2+HO2 reactions'
          units='mol/cm3/s'
       CASE(248)
          varname_new='RXN_19'
          longname='Sum of RO2+RO2 reactions'
          units='mol/cm3/s'
       CASE(249)
          varname_new='RXN_20'
          longname='CH3CO3+NO2+M->PAN+M'
          units='mol/cm3/s'

      CASE(156)
          varname_new='RXN_21'
          longname='O2 photolysis'
          units='mol/cm3/s'
      CASE(157)
          varname_new='RXN_22'
          longname='Cl2O2 photolysis'
          units='mol/cm3/s'
      CASE(158)
          varname_new='RXN_23'
          longname='O3+hv->O(1D)'
          units='mol/cm3/s'
      CASE(159)
          varname_new='RXN_24'
          longname='Total O(1D) production'
          units='mol/cm3/s'
      CASE(160)
          varname_new='RXN_25'
          longname='Total ozone production'
          units='mol/cm3/s'

      CASE(161)
          varname_new='RXN_26'
          longname='Total ozone loss'
          units='mol/cm3/s'
      CASE(162)
          varname_new='RXN_27'
          longname='O3+OH->HO2+O2'
          units='mol/cm3/s'
      CASE(163)
          varname_new='RXN_28'
          longname='O3+HO2->OH+O2'
          units='mol/cm3/s'
      CASE(165)
          varname_new='RXN_29'
          longname='HO2+NO->NO2+OH'
          units='mol/cm3/s'
      CASE(166)
          varname_new='RXN_30'
          longname='CH3O2+NO'
          units='mol/cm3/s'

      CASE(170)
          varname_new='RXN_31'
          longname='O3+isoprene'
          units='mol/cm3/s'
      CASE(171)
          varname_new='RXN_32'
          longname='Total CH4 loss'
          units='mol/cm3/s'
      CASE(172)
          varname_new='RXN_33'
          longname='Total CO loss'
          units='mol/cm3/s'
      CASE(173)
          varname_new='RXN_34'
          longname='Total OH production'
          units='mol/cm3/s'
!      CASE(174)
!          varname_new='OZDEPOS'
!          longname='Ozone deposition rate'
!          units='m/s'

       CASE(221)
          varname_new='O3Losst'
          longname='Ozone loss from tracer'
          units='mol/mol/s'
       CASE(222)
          varname_new='O3Lossd'
          longname='Ozone loss diagnosed'
          units='mol/mol/s'
       CASE(224)
          varname_new='O3Prod'
          longname='Ozone production'
          units='mol/mol/s'
       CASE(225)
          varname_new='O3M1'
          longname='Pre-Chem ozone'
          units='mol/mol'
       CASE(226)
          varname_new='O3'
          longname='Post-chem ozone'
          units='mol/mol'
       CASE(227)
          varname_new='cO3'
          longname='Chem ozone tendency'
          units='mol/mol/s'
       CASE(228)
          varname_new='Reg'
          longname='Regions'
          units='none'
       END SELECT

    END SELECT

  END SUBROUTINE get_varname_long_varname_units

  SUBROUTINE HANDLE_ERR(status)
      
    INTEGER, INTENT(IN) :: status
      
    IF (status /= NF90_NOERR) THEN
       PRINT *, TRIM(NF90_STRERROR(status))
       STOP
    END IF
    
  END SUBROUTINE HANDLE_ERR  

END PROGRAM get_varname_from_table
