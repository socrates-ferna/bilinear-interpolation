PROGRAM bilinear_interpolation
    use iso_fortran_env, only: real64
    use subrutinas, only: BILINEAR
    IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! INTERPOLACIÓN BILINEAL !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    REAL(real64), ALLOCATABLE, DIMENSION(:,:) :: input_array, bitmap
    INTEGER, ALLOCATABLE, DIMENSION(:) :: npxintx, npxinty, cumnpxintx,cumnpxinty
    REAL(real64), DIMENSION(3) :: Q11, Q22
    REAL(real64) :: Q21, Q12
    REAL(real64) :: domainwidth, xleft, xright, dx, domainheight, ybot, ytop, dy, &
                         resto_intanterior, liy, lix
    INTEGER :: status, dim, idim, jdim, npxhigh, npxwide, npixels, i, j, k, arrelem, & 
                        columna, fila, iniciofilabit, finalfilabit
    CHARACTER(200) :: msg
    CHARACTER(1) :: idim_str
    CHARACTER(9) :: format_str

!----------------------------!
! BLOQUE I: ENTRADA DE DATOS !
!----------------------------!

    OPEN(1,FILE='input_array',ACTION='read',IOSTAT=status,IOMSG=msg)
    READ(1,001) idim
    READ(1,001) jdim
    dim = idim*jdim
    001 FORMAT(5X,I3)
    ALLOCATE(input_array(dim,3))
    DO i = 1,dim
        READ(1,'(3(ES13.6,TR2))') input_array(i,1), input_array(i,2), input_array(i,3)
        WRITE(*,'(A7,3(F8.2))') 'I read', input_array(i,1), input_array(i,2), input_array(i,3)
    END DO
    READ(1,'(8X,I4)') npxwide
    READ(1,'(8X,I4)') npxhigh
    WRITE(*,*) 'I read npxwide', npxwide, 'npxhigh', npxhigh
    CLOSE(1)

!----------------------------------!
! BLOQUE II: BITMAP INITIALISATION !
!----------------------------------!
    !Asumimos que la malla de entrada está ordenada por filas
    !empezando por la esquina superior izquierda y terminando en la inferior derecha
    !Si no fuera así tendríamos que incluir aquí un algoritmo de ordenación de los elementos
    !o alguna herramienta que nos permitiera acceder a las cuaternas de puntos
    npixels = npxwide*npxhigh
    ALLOCATE(bitmap(npixels,6))
    xright = MAXVAL(input_array(:,1))
    xleft = MINVAL(input_array(:,1))
    domainwidth = xright - xleft
    dx = domainwidth / npxwide
    ybot = MINVAL(input_array(:,2))
    ytop = MAXVAL(input_array(:,2))
    domainheight = ytop - ybot
    dy = domainheight / npxhigh
    bitmap(:,:) = 0.0D0   !!! LOS PIXEL VALUES SON INTEGER, O LOS CONVIERTES AL ESCRIBIR O TE HACES OTRO ARARY
    
    DO i=1,npxhigh
        DO j=1,npxwide
            arrelem = npxwide * (i-1) + j
            bitmap(arrelem,1) = xleft + (j-0.5D0)*dx
            bitmap(arrelem,2) = ytop - (i-0.5D0)*dy
        END DO
    END DO !YA TENEMOS LAS COORDENADAS DE LOS CENTROS DE CADA PIXEL, SON SIMBÓLICAS, NO SON NECESARIAS PARA CONSTRUIR EL BMP
    ! SI QUEREMOS SACAR LOS VALORES PARA CUBRIR HASTA LOS LÍMITES EXACTOS DEL DOMINIO HAREMOS OTRO ARRAY
    !Contadores para caracterizar las subcuadrículas de la matriz
    ALLOCATE(npxintx(1:idim-1))
    ALLOCATE(npxinty(1:jdim-1))
    ALLOCATE(cumnpxintx(1:idim-1))
    ALLOCATE(cumnpxinty(1:idim-1))

    resto_intanterior = 0.0
    cumnpxintx(1) = 1
    DO i=1, idim-1
        lix = input_array(i+1,1) - input_array(i,1) + resto_intanterior
        npxintx(i) = NINT(lix / dx) !Entero más cercano de "píxeles" que caben en el subintervalo en x, me vale de índice
        resto_intanterior = lix/dx - npxintx(i)
        IF(i .gt. 1) THEN
            cumnpxintx(i) = cumnpxintx(i-1) + npxintx(i)
        END IF
    END DO
    resto_intanterior = 0.0

    cumnpxinty(1) = 1
    DO j=1, dim - 2*idim+1, idim
        liy = input_array(j,2) - input_array(j+idim,2) + resto_intanterior
        npxinty(j) = NINT(liy / dy)
        resto_intanterior = liy/dy - npxinty(j)
        IF (j .gt. 1) THEN
            cumnpxinty(j) = cumnpxinty(j-1) + npxinty(j)
        END IF
    END DO
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
    DO i=1,dim-idim-1 !cada fila son los datos de una cuadrícula
        print*,'i',i

        IF (MOD(i,idim) == 0) THEN
            fila = fila + 1
            columna = 1
            CYCLE
        END IF
        
        Q12 = input_array(i,3)
        Q11 = input_array(i+idim,:)
        Q21 = input_array(i+idim+1,3)
        Q22 = input_array(i+1,:)
        

        DO j=1,npxinty(fila) 
            iniciofilabit=(j-1)*npxwide+cumnpxintx(columna)   !esto parece que está bien para cada fila
            finalfilabit=iniciofilabit+npxintx(columna) - 1 ! el subindice 2 depende de la cuadrícula que estés haciendo
            !le pasamos la fila a la interpolación
            WRITE(*,*) 'iniciofilabit', iniciofilabit, 'finalfilabit', finalfilabit
            CALL BILINEAR(bitmap(iniciofilabit:finalfilabit,1:3),iniciofilabit,finalfilabit, &
                          Q11,Q22,Q12,Q21)
            WRITE(*,'(10(F8.3))') (bitmap(k,3), k=iniciofilabit,finalfilabit)
            WRITE(*,'(10(F8.3))') (bitmap(k,1), k=iniciofilabit,finalfilabit)
            WRITE(*,'(10(F8.3))') (bitmap(k,2), k=iniciofilabit,finalfilabit)
        END DO
        columna = columna + 1
        print*, columna
    END DO

    !WRITE(idim_str,'(I2)') npxwide
    !format_str = "'" // idim_str // "(F8.3)"
    !WRITE(*,'(10(F8.3))') (bitmap(i,3), i=1,npixels)



END PROGRAM bilinear_interpolation 
