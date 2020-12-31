PROGRAM bilinear_interpolation
    use iso_fortran_env, only: real64
    use subrutinas
    IMPLICIT NONE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! INTERPOLACIÓN BILINEAL !!!!!!!!!!!!!!!!!!!!!!!!!
!
! IMPORTANTE: COMENTARIOS A LA IMPLEMENTACION
!
! -La idea inicial era intercomunicar este programa con un script copiado de stackoverflow de C que escribiera un .bmp
! - Por el alcance del trabajo se descartó la idea y con ella la generación del array con
!      los valores RGB que correspondieran a la escala de color elegida entre min(f(x,y)) y max(f(x,y))
!         implementar esta escala no tiene mucha dificultad así que es sólo una cuestión de tiempo
! -Como el .bmp se escribe como una sola línea de bits recorriendo la imagen por filas
!      desde la esquina superior izquierda hasta la inf derecha pensé en guardar el array de la misma manera, es decir,
!         3 columnas en el array de entrada y en el bitmap: x,y,f(x,y)... 
!   ... ¡¡¡Con lo fácil que hubiera sido usar un array de 3 dimensiones y luego RESHAPE!!! Y ni siquiera eso hacía falta
! - Me di cuenta tarde y ya tuve que tirar hasta el final con este programa tan contraintuitivo
!
! Resumen del programa:
! 1. Lectura de las dimensiones del dominio, de la malla de entrada y del bitmap de salida con la NAMELIST
! 2. Lectura del array en input_array o generación con la función que se elija, sólo está disponible la función Schwefel, pero añadir más es trivial
! 3. Como se pueden agrupar los elementos del bitmap según los cuatro puntos originales que les rodean, se crean varios arrays auxiliares:
!       3.1 npxintx y npxinty dicen cuántos píxeles entran en x y en y en cada una de las cuadrículas encerradas por 4 puntos de la malla original
!       3.2 cumnpxintx y cumnpxinty sirven de soporte para acceder a los elementos adecuados del bitmap
! 4. El bucle exterior del método recorre la malla original hasta llegar al punto superior izquierdo de la última subcuadrícula, toma la
!       información necesaria en cada vuelta y le pasa en el bucle interior a la subrutina BILINEAR
!           todas las "subfilas" de píxeles a interpolar
!
! -Hemos dejado WRITE y print comentados por todo el código por si fuera necesario echar mano de ellos
! -En general hemos intentado hacer un programa eficiente que ahorrara accesos a memoria e hiciera operaciones sobre trozos del
!     bitmap los más grandes posible, no parece que lo hayamos conseguido.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    REAL(real64), ALLOCATABLE, DIMENSION(:,:) :: input_array, bitmap
    REAL(real64), ALLOCATABLE, DIMENSION(:) :: exact
    INTEGER, ALLOCATABLE, DIMENSION(:) :: npxintx, npxinty, cumnpxintx,cumnpxinty
    REAL(real64), DIMENSION(3) :: Q11, Q22
    REAL(real64) :: Q21, Q12, L2, LINF
    REAL(real64) :: domainwidth, xleft, xright, dx, domainheight, ybot, ytop, dy, &
                         resto_intanterior, liy, lix
    INTEGER :: status, dim, idim, jdim, npxhigh, npxwide, npixels, i, j, k, arrelem, & 
                        columna, fila, iniciofilabit, finalfilabit, index
    CHARACTER(200) :: msg, filename
    CHARACTER(4) :: input_mode
    CHARACTER(16) :: format_str, func,i_str,j_str,pxw_str,pxh_str,aux_str
    LOGICAL :: ex, save=.TRUE.
    PROCEDURE(SCHWEFEL2D), POINTER :: function_pointer
    NAMELIST / input_conf / xleft,xright,ytop,ybot,idim,jdim,npxwide,npxhigh

    function_pointer => SCHWEFEL2D !! Función utilizada por defecto, se podrían incluir otras y procedimientos para elegirlas
    input_mode = 'func' !! choose 'func' for self generating the exact and 'file' for reading an input file
    func = 'schwefel'
!----------------------------!
! BLOQUE I: ENTRADA DE DATOS !
!----------------------------!

    OPEN(1,FILE='input_conf.nml',ACTION='READ',IOSTAT=status,IOMSG=msg)
    READ(1,NML=input_conf)
    dim = idim * jdim !j:columnas, i:filas
    ALLOCATE(input_array(dim,3))
    IF (npxhigh < 10*idim) npxhigh = 10*idim  !! Por si el número de puntos interpolados solicitado es bajo
    IF(npxwide < 10*jdim) npxwide = 10*jdim
    npixels = npxwide*npxhigh
    ALLOCATE(bitmap(npixels,3))
    domainwidth = xright - xleft
    dx = domainwidth / (jdim - 1)
    domainheight = ytop - ybot
    dy = domainheight / (idim - 1)
    bitmap(:,:) = 200.0 !!!Me resulta más fácil encontrar bugs si no inicializo los arrays con el valor 0, que suelo pasar por alto

    IF (input_mode == 'read') THEN !! IMPORTANTE, LA INFO EN input_conf.nml DEBE SER CONCORDANTE CON EL ARRAY QUE SE LEA EN input_array
        OPEN(2,FILE='input_array',ACTION='read',IOSTAT=status,IOMSG=msg)
        DO i = 1,dim
            READ(2,'(3(ES13.6,TR2))') input_array(i,1), input_array(i,2), input_array(i,3)
            !WRITE(*,'(A7,3(F8.2))') 'I read', input_array(i,1), input_array(i,2), input_array(i,3)
        END DO

        CLOSE(2)

    ELSE IF(input_mode == 'func') THEN
        DO i=1,idim
            DO j=1,jdim
                arrelem = jdim * (i-1) + j
                input_array(arrelem,1) = xleft + (j-1) * dx
                input_array(arrelem,2) = ytop - (i-1) * dy
                !WRITE(*,'(2(A2,F8.2),2(A2,I0))') 'x=',input_array(arrelem,1),'y=',input_array(arrelem,2),'i',i,'j',j
            END DO
        END DO
        !print*, 'grid done'
        filename = 'analytical'
        CALL INTTOSTRING(i_str,idim,index,aux_str)
        CALL APPENDSTRING(filename,i_str,'_')
        CALL INTTOSTRING(j_str,jdim,index,aux_str)
        CALL APPENDSTRING(filename,j_str,'_')
        CALL APPENDSTRING(filename,'.dat','')
        CALL GENEXACT(exact,input_array(:,1),input_array(:,2),function_pointer,save,filename)
        !print*,'exact done'
        input_array(:,3) = exact(:)
        !print*,input_array(:,1)
    END IF

!----------------------------------!
! BLOQUE II: BITMAP INITIALISATION !
!----------------------------------!
    !el caso de malla estructurada equiespaciada es el sencillo, si la malla fuera desestructurada sería necesaria una estrategia
    ! de algoritmos de búsqueda de vecinos para conocer los puntos más cercanos

    dx = domainwidth / npxwide  !sobrescribimos las variables dx y dy porque los valores anteriores no son necesarios ya
    dy = domainheight / npxhigh
    DO i=1,npxhigh
        DO j=1,npxwide
            arrelem = npxwide * (i-1) + j
            bitmap(arrelem,1) = xleft + (j-0.5D0)*dx
            bitmap(arrelem,2) = ytop - (i-0.5D0)*dy
        END DO
    END DO 
    !YA TENEMOS LAS COORDENADAS DE LOS CENTROS DE CADA "PIXEL" INTERPOLADO
    !SI QUISIÉRAMOA SACAR LOS VALORES PARA CUBRIR HASTA LOS LÍMITES EXACTOS DEL DOMINIO HARÍAMOS OTRO ARRAY CON LOS PUNTOS ADECUADAMENTE DISTRIBUIDOS
    
    !Contadores para caracterizar las subcuadrículas de la matriz
    ALLOCATE(npxintx(1:jdim-1))
    ALLOCATE(npxinty(1:idim-1))
    ALLOCATE(cumnpxintx(1:jdim-1))
    ALLOCATE(cumnpxinty(1:idim-1))
    !print*,'allocations done'
    resto_intanterior = 0.0
    cumnpxintx(1) = 1

    DO j=1, jdim-1
        !print*, resto_intanterior
        !print*, input_array(j:j+1,1)
        lix = input_array(j+1,1) - input_array(j,1) + resto_intanterior
        !print*, 'lix',lix
        npxintx(j) = NINT(lix / dx) !Entero más cercano de "píxeles" que caben en el subintervalo en x
        resto_intanterior = (lix/dx - npxintx(j)) * dx !Me guardo el resto para modificar el intervalo de la siguiente cuadrícula
        IF(j .gt. 1) THEN
            cumnpxintx(j) = cumnpxintx(j-1) + npxintx(j-1)
        END IF
    END DO

    resto_intanterior = 0.0
    !print*,'x dir ok'
    cumnpxinty(1) = 0
    j=1
    DO i=1, dim - 2*jdim+1, jdim
        !print*,'i=',i
        !print*, resto_intanterior
        liy = input_array(i,2) - input_array(i+jdim,2) + resto_intanterior
        npxinty(j) = NINT(liy / dy)
        resto_intanterior = (liy/dy - npxinty(j)) * dy
        IF (j .gt. 1) THEN
            cumnpxinty(j) = cumnpxinty(j-1) + npxinty(j-1)
        END IF
        j=j+1
    END DO

    IF(SUM(npxinty) > npxhigh) npxinty(idim-1) = npxinty(idim-1) - 1
    IF(SUM(npxintx) > npxwide) npxintx(jdim-1) = npxintx(jdim-1) - 1
    !print*,'dx',dx,'dy',dy,'lix',lix,'liy',liy
    !WRITE(*,'(10(I8))') cumnpxintx(:)
    !WRITE(*,'(10(I8))') cumnpxinty(:)
    !WRITE(*,'(10(I8))') npxintx(:)
    !WRITE(*,'(10(I8))') npxinty(:)

    format_str="(10(F8.3))"
    columna = 1
    fila = 1 !CORRESPONDEN A LA FILA Y COLUMNA DE CADA CUADRÍCULA DE PUNTOS DE LA MALLA ORIGINAL
    print*, 'llego al bucle'

!----------------------------------!
! BLOQUE III: BUCLE DEL MÉTODO     !
!----------------------------------!
    DO i=1,dim-jdim-1
        !print*,'i',i

        IF (MOD(i,jdim) == 0) THEN
            fila = fila + 1
            columna = 1
            CYCLE ! el último elemento de cada fila en la malla original no es el superior izquierdo de ninguna cuadrícula
        END IF
        
        Q12 = input_array(i,3)
        Q11 = input_array(i+jdim,:)
        Q21 = input_array(i+jdim+1,3)
        Q22 = input_array(i+1,:)
        

        DO j=1,npxinty(fila) 
            iniciofilabit= cumnpxinty(fila)*npxwide + (j-1)*npxwide+cumnpxintx(columna)
            finalfilabit=iniciofilabit+npxintx(columna) - 1
            !WRITE(*,*) 'iniciofilabit', iniciofilabit, 'finalfilabit', finalfilabit
            CALL BILINEAR(bitmap(iniciofilabit:finalfilabit,1:3),iniciofilabit,finalfilabit, &
                          Q11,Q22,Q12,Q21)
            
            !WRITE(*,format_str) (bitmap(k,3), k=iniciofilabit,finalfilabit)
            !WRITE(*,format_str) (bitmap(k,1), k=iniciofilabit,finalfilabit)
            !WRITE(*,format_str) (bitmap(k,2), k=iniciofilabit,finalfilabit)
        END DO

        columna = columna + 1
        !print*, columna
    END DO

!----------------------------------!
! BLOQUE IV: ARCHIVOS DE SALIDA    !
!----------------------------------!

    filename = 'analytical'
    CALL INTTOSTRING(pxw_str,npxwide,index,aux_str)
    CALL APPENDSTRING(filename,pxw_str,'_')


    CALL INTTOSTRING(pxh_str,npxhigh,index,aux_str)
    CALL APPENDSTRING(filename,pxh_str,'_')
    CALL APPENDSTRING(filename,'.dat','')
    
    CALL ERRORES(bitmap,L2,LINF,function_pointer,save,filename)



    PRINT*, 'L2=',L2
    PRINT*, 'LINF=',LINF

    INQUIRE(FILE='errors.csv',EXIST=ex)
    IF (ex) THEN
        OPEN(UNIT=100,FILE='errors.csv',STATUS='OLD', ACCESS='APPEND', ACTION='WRITE',IOSTAT=status,IOMSG=msg)

    ELSE
        OPEN(UNIT=100,FILE='errors.csv',STATUS='NEW', ACTION='WRITE',IOSTAT=status,IOMSG=msg)        
        WRITE(100,*) 'func,npxwide,npxhigh,idim,jdim,domheight,domwidth,L2,LINF'
    END IF

    005 FORMAT(A8,',',4(I0,','),4(F8.2,:,','))
    WRITE(100,005) TRIM(ADJUSTL(func)),npxwide,npxhigh,idim,jdim,domainheight,domainwidth,L2,LINF
    CLOSE(100)
    
    filename = TRIM(ADJUSTL(func))
    CALL INTTOSTRING(i_str,idim,index,aux_str)
    CALL APPENDSTRING(filename,i_str,'_')

    CALL INTTOSTRING(j_str,jdim,index,aux_str)
    CALL APPENDSTRING(filename,j_str,'_')
    
    CALL INTTOSTRING(pxw_str,npxwide,index,aux_str)
    CALL APPENDSTRING(filename,pxw_str,'_')

    CALL INTTOSTRING(pxh_str,npxhigh,index,aux_str)
    CALL APPENDSTRING(filename,pxh_str,'_')

    CALL APPENDSTRING(filename,'.dat','')

    CALL WRITEARRAY(bitmap(:,3),bitmap(:,1),bitmap(:,2),filename,200)

END PROGRAM bilinear_interpolation 
