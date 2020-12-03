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
        !Hacemos esto porque es más barato que hacer la interpolación en dos pasos, primero la x y luego la y
        !a0, a1, a2, a3; f(x,y)~= a0+a1*x+a2*y+a3*xy
        x1 = q11(1)
        x2 = q22(1)
        y1 = q11(2) !! Estas asignaciones las hacemos para que
        y2 = q22(2) !! el código sea un poquito más legible, y para no acceder a los arrays tantas veces
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

        fila(:,3) = a0 + a1*fila(:,1) + a2*fila(:,2) + a3*fila(:,1)*fila(:,2)

        END SUBROUTINE BILINEAR
END MODULE SUBRUTINAS
