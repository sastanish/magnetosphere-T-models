module parser

  implicit none

  public : parse_input_file

contains

  subroutine parse_input_file(inFile,model,dipole,cells,bounds,dataFilePath&
                              outputFilePath, outputFlags)

    character(len=*),  intent(in)  :: file_path
    integer, intent(out) :: model, dipole
    logical, intent(out), dimension(8) :: outputFlags
    character(len=99), intent(out) :: outputFilePath
    character(len=99), intent(out) :: dataFilePath
    integer, intent(out) :: cells(3), bounds(6)

    integer:: fileUnit, ioStat


    ! outputFlags look like:
    

    ! Namelist definition===============================
    namelist /control/ &
        model, & ! 0-none, 1-TS04, 2-TA16
        dipole, & ! 0-none, 1-dipole08, 2-IGRF
        cells, & ! [nx, ny, nz]
        bounds, & ! [ x+,x-,y+,y-,z+,z- ]
        dataFilePath, & ! Path to data.txt
        outputFilePath, & ! output file. Each timestep will be outputFilePath_{lineNumber}.nc
        outputFlags ! [bvec,bmag,c2_t1_vec,c2_t1_mag,c2_t2_vec,c2_t2_mag,c2_t3_vec,c2_t3_mag]
      model = 1
      dipole = 1
      cells = [10, 10, 10]
      bounds = [ -5, 5, -5, 5, -5, 5 ]
      dataFilePath = "./data.txt"
      outputFilePath = "./ouput"
      outputFlags = [ .False., .False., .False., .False., .False., .False., .False., .False. ]
    ! Namelist definition===============================

    call open_inputfile(inFile, fileUnit, ioStat)
    if (iostat /= 0) then
      print *, "input file './control.nml' cannot be opened"
      return
    end if

    read (nml=ORDER, iostat=ioStat, unit=fileUnit)
    call close_inputfile(inFile, fileUnit, ioStat)
    if (iostat /= 0) then
      print *, "input file './control.nml' cannot be read"
      return
    end if
  end subroutine parse_input_file

  subroutine open_inputfile(file_path, file_unit, iostat)
    !! Check whether file exists, with consitent error message
    !! return the file unit
    character(len=*),  intent(in)  :: file_path
    integer,  intent(out) :: file_unit, iostat

    inquire (file=file_path, iostat=iostat)
    if (iostat /= 0) then
      write (stderr, '(3a)') 'Error: file "', trim(file_path), '" not found!'
    end if
    open (action='read', file=file_path, iostat=iostat, newunit=file_unit)
  end subroutine open_inputfile

  subroutine close_inputfile(file_path, file_unit, iostat)
    !! Check the reading was OK
    !! return error line IF not
    !! close the unit
    character(len=*),  intent(in)  :: file_path
    character(len=1000) :: line
    integer,  intent(in) :: file_unit, iostat

    if (iostat /= 0) then
      write (stderr, '(2a)') 'Error reading file :"', trim(file_path)
      write (stderr, '(a, i0)') 'iostat was:"', iostat
      backspace(file_unit)
      read(file_unit,fmt='(A)') line
      write(stderr,'(A)') &
          'Invalid line : '//trim(line)
    end if
    close (file_unit)
  end subroutine close_inputfile

end module parser
