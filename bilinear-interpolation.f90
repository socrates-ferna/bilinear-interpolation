PROGRAM bilinear_interpolation
    use iso_fortran_env, only: real64
    use subrutinas
    IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! INTERPOLACIÓN BILINEAL !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    REAL(real64), ALLOCATABLE, DIMENSION(:,:) :: input_array, bitmap
    REAL(real64), ALLOCATABLE, DIMENSION(:) :: exact
    INTEGER, ALLOCATABLE, DIMENSION(:) :: npxintx, npxinty, cumnpxintx,cumnpxinty
    REAL(real64), DIMENSION(3) :: Q11, Q22
    REAL(real64) :: Q21, Q12, L2, LINF
    REAL(real64) :: domainwidth, xleft, xright, dx, domainheight, ybot, ytop, dy, &
                         resto_intanterior, liy, lix
    INTEGER :: status, dim, idim, jdim, npxhigh, npxwide, npixels, i, j, k, arrelem, & 
                        columna, fila, iniciofilabit, finalfilabit
    CHARACTER(200) :: msg, filename
    CHARACTER(4) :: input_mode
    CHARACTER(16) :: format_str, func
    LOGICAL :: ex, save=.TRUE.
    PROCEDURE(SCHWEFEL2D), POINTER :: function_pointer
    NAMELIST / input_conf / xleft,xright,ytop,ybot,idim,jdim,npxwide,npxhigh

    function_pointer => SCHWEFEL2D !! Función utilizada por defecto, se podrían incluir otras y procedimientos para elegirlas
    input_mode = 'func' !! choose 'func' for self generating the exact and 'file' for reading an input file
    func = 'schwefel'
!----------------------------!
! BLOQUE I: ENTRADA DE DATOS !
!----------------------------!
    !explciar que está la opción de leer la función porque lo normal en programas de posproc
    !es leer archivos de resultados y dibujar las cosas desde ahí

    OPEN(1,FILE='input_conf.nml',ACTION='READ',IOSTAT=status,IOMSG=msg)
    READ(1,NML=input_conf)
    dim = idim * jdim
    ALLOCATE(input_array(dim,3))
    !IF (npxhigh < 10*idim) npxhigh = 10*idim  !! Por si el número de puntos interpolados solicitado es bajo
    !IF(npxwide < 10*jdim) npxwide = 10*jdim
    npixels = npxwide*npxhigh
    ALLOCATE(bitmap(npixels,3))
    domainwidth = xright - xleft
    dx = domainwidth / (jdim - 1)
    domainheight = ytop - ybot
    dy = domainheight / (idim - 1)
    bitmap(:,:) = 0.0   !!! LOS PIXEL VALUES SON INTEGER NECESITAS OTRO ARRAY O QUIZÁS UNA ESTRUCTURA ESPECIAL ¿O CONVERTIRLO A INT?

    IF (input_mode == 'read') THEN !! IMPORTANTE, LA INFO EN input_conf.nml DEBE SER CONCORDANTE CON EL ARRAY QUE SE LEA EN input_array
        OPEN(2,FILE='input_array',ACTION='read',IOSTAT=status,IOMSG=msg)
        DO i = 1,dim
            READ(2,'(3(ES13.6,TR2))') input_array(i,1), input_array(i,2), input_array(i,3)
            WRITE(*,'(A7,3(F8.2))') 'I read', input_array(i,1), input_array(i,2), input_array(i,3)
        END DO

        CLOSE(2)

    ELSE IF(input_mode == 'func') THEN
        DO i=1,idim
            DO j=1,jdim
                arrelem = jdim * (i-1) + j
                input_array(arrelem,1) = xleft + (j-1) * dx
                input_array(arrelem,2) = ytop - (i-1) * dy
                WRITE(*,'(2(A2,F8.2),2(A2,I0))') 'x=',input_array(arrelem,1),'y=',input_array(arrelem,2),'i',i,'j',j
            END DO
        END DO
        print*, 'grid done'
        CALL GENEXACT(exact,input_array(:,1),input_array(:,2), function_pointer)
        print*,'exact done'
        input_array(:,3) = exact(:)
        print*,input_array(:,1)
    END IF
!    GOTO 10

 !   2 CONTINUE
!    IF(input_mode == 'func') THEN
!        CALL GENEXACT(exact,input_array(:,1),input_array(:,2), function_pointer)
!        input_array(:,3) = exact (:) !Sobreescribir los valores de input leídos, nos ahorramos temporalmente generar la malla de puntos
!    END IF


  !  GOTO 10

!----------------------------------!
! BLOQUE II: BITMAP INITIALISATION !
!----------------------------------!
    !Asumimos que la malla de entrada está ordenada por filas
    !empezando por la esquina superior izquierda y terminando en la inferior derecha
    !Si no fuera así tendríamos que incluir aquí un algoritmo de ordenación de los elementos
    !o alguna herramienta que nos permitiera acceder a las cuaternas de puntos

    dx = domainwidth / npxwide
    dy = domainheight / npxhigh
    DO i=1,npxhigh
        DO j=1,npxwide
            arrelem = npxwide * (i-1) + j
            bitmap(arrelem,1) = xleft + (j-0.5D0)*dx
            bitmap(arrelem,2) = ytop - (i-0.5D0)*dy
        END DO
    END DO !YA TENEMOS LAS COORDENADAS DE LOS CENTROS DE CADA "PIXEL" INTERPOLADO, SON SIMBÓLICAS, NO SON NECESARIAS PARA CONSTRUIR EL .BMP
    ! SI QUEREMOS SACAR LOS VALORES PARA CUBRIR HASTA LOS LÍMITES EXACTOS DEL DOMINIO HAREMOS OTRO ARRAY CON LOS PUNTOS ADECUADAMENTE DISTRIBUIDOS
    !Contadores para caracterizar las subcuadrículas de la matriz
    ALLOCATE(npxintx(1:jdim-1))
    ALLOCATE(npxinty(1:idim-1))
    ALLOCATE(cumnpxintx(1:jdim-1))
    ALLOCATE(cumnpxinty(1:idim-1))
    print*,'allocations done'
    resto_intanterior = 0.0
    cumnpxintx(1) = 1

    DO j=1, jdim-1
        print*, resto_intanterior
        print*, input_array(j:j+1,1)
        lix = input_array(j+1,1) - input_array(j,1) + resto_intanterior
        print*, 'lix',lix
        npxintx(j) = NINT(lix / dx) !Entero más cercano de "píxeles" que caben en el subintervalo en x, me vale de índice
        resto_intanterior = lix/dx - npxintx(j)
        IF(j .gt. 1) THEN
            cumnpxintx(j) = cumnpxintx(j-1) + npxintx(j)
        END IF
    END DO
    resto_intanterior = 0.0
    print*,'x dir ok'
    cumnpxinty(1) = 0
    j=1
    DO i=1, dim - 2*jdim+1, jdim
        print*,'i=',i
        print*, resto_intanterior
        liy = input_array(i,2) - input_array(i+jdim,2) + resto_intanterior
        npxinty(j) = NINT(liy / dy)
        resto_intanterior = liy/dy - npxinty(j)
        IF (j .gt. 1) THEN
            cumnpxinty(j) = cumnpxinty(j-1) + npxinty(j)
        END IF
        j=j+1
    END DO
    IF(SUM(npxinty) > npxhigh) npxinty(idim-1) = npxinty(idim-1) - 1
    IF(SUM(npxintx) > npxwide) npxintx(jdim-1) = npxintx(jdim-1) - 1
    print*,'dx',dx,'dy',dy,'lix',lix,'liy',liy
    WRITE(*,'(3(I8))') cumnpxintx(:)
    WRITE(*,'(3(I8))') cumnpxinty(:)
    WRITE(*,'(3(I8))') npxintx(:)
    WRITE(*,'(3(I8))') npxinty(:)

    !!!! YA TIENES LOS VECTORES ACUMULATIVOS DE PÍXELES, TIENES QUE USAR ESO DENTRO DEL LOOP QUE RECORRE LA CUADRÍCULA PARA
    !!!! TENER EL COMIENZO Y EL FINAL DE SUB-ARRAY QUE TIENE QUE PASARLE A LA SUBRUTINA
    
    !ncuadriculas = (idim-1) * (jdim-1)  !! Definimos esta variable por deporte prácticamente
    

    !!mejor hazlo por filas y columnas
    columna = 1
    fila = 1
    !k = 0
    print*, 'llego al bucle'
    DO i=1,dim-jdim-1 !cada fila son los datos de una cuadrícula
        print*,'i',i

        IF (MOD(i,jdim) == 0) THEN
            fila = fila + 1
            columna = 1
            CYCLE
        END IF
        
        Q12 = input_array(i,3)
        Q11 = input_array(i+jdim,:)
        Q21 = input_array(i+jdim+1,3)
        Q22 = input_array(i+1,:)
        

        DO j=1,npxinty(fila) !!!REVISA ESTO PUEDE QUE HAYAS CAMBIADO FILAS Y COLUMNAS AQUÍ
            iniciofilabit= cumnpxinty(fila)*npxwide + (j-1)*npxwide+cumnpxintx(columna)   !esto parece que está bien para cada fila
            finalfilabit=iniciofilabit+npxintx(columna) - 1 ! el subindice 2 depende de la cuadrícula que estés haciendo
            !le pasamos la fila a la interpolación
            WRITE(*,*) 'iniciofilabit', iniciofilabit, 'finalfilabit', finalfilabit
            CALL BILINEAR(bitmap(iniciofilabit:finalfilabit,1:3),iniciofilabit,finalfilabit, &
                          Q11,Q22,Q12,Q21)
            format_str="(10(F8.3))" !!esto lo hago para probar que se pueden enchufar variables en la definición del formato, que no sabía si se podía
            WRITE(*,format_str) (bitmap(k,3), k=iniciofilabit,finalfilabit)
            WRITE(*,'(10(F8.3))') (bitmap(k,1), k=iniciofilabit,finalfilabit)
            WRITE(*,'(10(F8.3))') (bitmap(k,2), k=iniciofilabit,finalfilabit)
        END DO
        !k = k + npxinty(fila) * npxwide
        columna = columna + 1
        print*, columna
    END DO

    CALL ERRORES(bitmap,L2,LINF,function_pointer,save)
    PRINT*, 'L2=',L2
    PRINT*, 'LINF=',LINF

    INQUIRE(FILE='errors.csv',EXIST=ex)
    IF (ex) THEN
        OPEN(UNIT=100,FILE='errors.csv',STATUS='OLD', ACCESS='APPEND', ACTION='WRITE',IOSTAT=status,IOMSG=msg)

    ELSE
        OPEN(UNIT=100,FILE='errors.csv',STATUS='NEW', ACTION='WRITE',IOSTAT=status,IOMSG=msg)        
        WRITE(100,*) 'func,npxwide,npxhigh,idim,jdim,domheight,domwidth,L2,LINF'
    END IF

    005 FORMAT(A8,',',4(I0,','),4(F7.2,:,','))
    WRITE(100,005) TRIM(ADJUSTL(func)),npxwide,npxhigh,idim,jdim,domainheight,domainwidth,L2,LINF
    CLOSE(100)

    !!do another check as in errors.csv
    !!add pxwidth and height to filename ... plus bounding box?
    filename = TRIM(ADJUSTL(func)) // '.dat'
    OPEN(UNIT=200,FILE=filename,STATUS='NEW',ACTION='WRITE',IOSTAT=status,IOMSG=msg)
    010 FORMAT(6(ES10.3,:,','))

    DO i=1,npixels
        WRITE(200,010) bitmap(i,1),bitmap(i,2),bitmap(i,3) !!!whole bitmap output
    END DO
    CLOSE(200)


    !WRITE(idim_str,'(I2)') npxwide
    !format_str = "'" // idim_str // "(F8.3)"
    !WRITE(*,'(10(F8.3))') (bitmap(i,3), i=1,npixels)



END PROGRAM bilinear_interpolation 
