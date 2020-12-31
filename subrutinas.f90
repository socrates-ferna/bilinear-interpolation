MODULE SUBRUTINAS
use iso_fortran_env, only: real64

IMPLICIT NONE

CONTAINS
    SUBROUTINE BILINEAR(fila, ifila, ffila, q11, q22, fq12, fq21)
        INTEGER, INTENT(IN) :: ifila, ffila
        REAL(real64), INTENT(INOUT), DIMENSION(ifila:ffila,3) :: fila
        REAL(real64), INTENT(IN) :: fq12, fq21
        REAL(real64) :: a0,a1,a2,a3,x1,x2,y2,y1,difx,dify, fq11, fq22
        REAL(real64), INTENT(IN), DIMENSION(3) :: q11, q22
        ! Calculamos los coeficientes que nos resuelven la interpolación
        x1 = q11(1)
        x2 = q22(1)
        y1 = q11(2) !! Estas asignaciones las hacemos para que
        y2 = q22(2) !! el código sea un poquito más legible
        difx = x1-x2
        dify = y1-y2
        fq11 = q11(3)
        fq22 = q22(3)

        a0 = fq11*x2*y2/((difx*dify)) + fq12*x2*y1/(difx*(-dify)) + &
            fq21*x1*y2/(difx*(-dify)) + fq22*x1*y1/(difx*dify)
        
        a1 = fq11*y2/(difx*(-dify)) + fq12*y1/(difx*dify) + fq21*y2/(difx*dify) +&
            fq22*y1/(difx*(-dify))

        a2 = fq11*x2/(difx*(-dify)) + fq12*x2/(difx*dify) + fq21*x1/(difx*dify) +&
            fq22*x1/(difx*(-dify))
        
        a3 = fq11/(difx*dify) + fq12/(difx*(-dify)) + fq21/(difx*(-dify)) + &
             fq22/(difx*dify)

        fila(ifila:ffila,3) = a0 + a1*fila(ifila:ffila,1) + a2*fila(ifila:ffila,2) + &
                              a3*fila(ifila:ffila,1)*fila(ifila:ffila,2)

        END SUBROUTINE BILINEAR

        REAL(real64) FUNCTION SCHWEFEL2D(x,y)
            REAL(real64) :: x,y
            SCHWEFEL2D=418.9829*2-x*SIN(SQRT(ABS(x)))-y*SIN(SQRT(ABS(y)))
        END FUNCTION SCHWEFEL2D

        SUBROUTINE GENEXACT(exact,x,y,fn,save,filename)
            REAL(real64), DIMENSION(:), ALLOCATABLE :: exact
            REAL(real64), DIMENSION(:) :: x,y
            INTEGER :: i
            CHARACTER(*) :: filename
            LOGICAL :: save
            PROCEDURE(SCHWEFEL2D),POINTER :: fn
            IF(.not. ALLOCATED(exact)) THEN
                ALLOCATE(exact(SIZE(x)))
            ELSE IF (ALLOCATED(exact)) THEN
                DEALLOCATE(exact)
                ALLOCATE(exact(SIZE(x)))
            END IF

            DO i=1,SIZE(exact)
                exact(i) = fn(x(i),y(i))
            END DO

            IF(save) THEN
                CALL WRITEARRAY(exact,x,y,filename,300)
            END IF
        END SUBROUTINE GENEXACT
        
        SUBROUTINE ERRORES(array, L2, LINF, fn,save,filename)
            CHARACTER(*) :: filename
            REAL(real64), DIMENSION(:,:), INTENT(IN) :: array
            REAL(real64), INTENT(OUT) :: L2, LINF
            REAL(real64), DIMENSION(:), ALLOCATABLE :: exact
            PROCEDURE(SCHWEFEL2D), POINTER :: fn
            LOGICAL :: save
            !INTEGER :: status, i
            !CHARACTER(200) :: msg

            CALL GENEXACT(exact,array(:,1),array(:,2),fn,save,filename)
            L2 = NORM2(array(:,3) - exact)
            LINF = MAXVAL(ABS(array(:,3) - exact))
        END SUBROUTINE ERRORES

        SUBROUTINE INTTOSTRING(str,int,i,aux_str)
            CHARACTER(*), INTENT(OUT) :: str
            INTEGER, INTENT(IN) :: int
            INTEGER, INTENT(OUT) :: i
            CHARACTER(10) :: i_str
            CHARACTER(*), INTENT(OUT) :: aux_str

            WRITE(str,'(I0)') int
            i = INDEX(ADJUSTL(str),' ')
            WRITE(i_str,'(I0)') i
            aux_str = "A" // TRIM(ADJUSTL(i_str))

        END SUBROUTINE INTTOSTRING

        SUBROUTINE APPENDSTRING(str,append,delim)
            CHARACTER(*), INTENT(INOUT) :: str
            CHARACTER(*), INTENT(IN) :: append,delim
            str = TRIM(ADJUSTL(str)) // TRIM(ADJUSTL(delim)) // TRIM(ADJUSTL(append))
        END SUBROUTINE

        SUBROUTINE WRITEARRAY(f,x,y,filename,unit)
            REAL(real64), INTENT(IN), DIMENSION(:) :: f,x,y
            CHARACTER(*), INTENT(IN) :: filename
            INTEGER :: status, i, unit
            CHARACTER(200) :: msg

            OPEN(UNIT=unit,FILE=filename,STATUS='NEW',ACTION='WRITE',IOSTAT=status,IOMSG=msg)
            010 FORMAT(6(ES10.3,:,','))
            DO i=1,SIZE(f)
                WRITE(unit,010) x(i),y(i),f(i)
            END DO
            CLOSE(unit)
        END SUBROUTINE

END MODULE SUBRUTINAS
