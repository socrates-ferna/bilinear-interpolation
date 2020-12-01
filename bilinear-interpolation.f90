PROGRAM bilinear_interpolation
    use iso_fortran_env, only: real64

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! INTERPOLACIÓN BILINEAL !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    REAL(real64), ALLOCATABLE, DIMENSION(:,:) :: input_array, bitmap
    REAL(real64) :: domainwidth, xleft, xright, dx, domainheight, ybot, ytop, dy
    INTEGER :: status, dim, npxhigh, npxwide, i, j, arrelem
    CHARACTER(200) :: msg

!----------------------------!
! BLOQUE I: ENTRADA DE DATOS !
!----------------------------!

    OPEN(1,FILE='input_array',ACTION='read',IOSTAT=status,IOMSG=msg)
    READ(1,001) dim
    001 FORMAT(10X,I3)
    ALLOCATE(input_array(dim,3))
    DO i = 1,dim
        READ(1,'(3(ES13.6,TR2))') input_array(i,1), input_array(i,2), input_array(i,3)
        WRITE(*,*) 'I read', input_array(i,1), input_array(i,2), input_array(i,3)
    END DO
    READ(1,'(12X,I4)') npxwide
    READ(1,'(13X,I4)') npxhigh
    CLOSE(1)

!----------------------------------!
! BLOQUE II: BITMAP INITIALISATION !
!----------------------------------!
    !Asumimos que la malla de entrada está ordenada por filas
    !empezando por la esquina superior izquierda y terminando en la inferior derecha
    !Si no fuera así tendríamos que incluir aquí un algoritmo de ordenación de los elementos
    !o alguna herramienta que nos permitiera acceder a las cuaternas de puntos
    npixels = npxwide*npxhigh
    ALLOCATE(bitmap(npixels,5))
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
    
    






END PROGRAM bilinear_interpolation 
