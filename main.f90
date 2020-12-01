
module Interpolación_bilineal
!funcion que calcula f(x,y) en funcion de 4 puntos conocidos 

use iso_fortran_env, only: real64
implicit none
contains

function interpolacion_bilineal_rectangulos(pxyin, fin, pxyout) result(fout)
!hace la interpolacion bilineal por rectangulos.
!hay que meterle las coordenadas de las puntas del rectangulo, los valores de la funcion en las esquinas
!y las coordenadas del punto a calcular

real(real64) , intent(in) :: pxyin(2,2) !la primera es elsubindice y la segunda si es x o y
real(re64) , intent(in) :: fin(2,2), pxyout(2)
real(rela64) ,allocatable :: fout(1)

!declaracion de variables

fout(1) = 1/((pxyin(2,1)-pxyin(1,1))*(pxyin(2,2)-pxyin(1,2)))*( &
(fin(1,1)* (pxyin(2,1)-pxyout(1))*(pxyin(2,2)-pxyout(2)) ) + &
(fin(2,1)*(-pxyin(1,1)+pxyout(1))*(pxyin(2,2)-pxyout(2)) ) + &
(fin(1,2)*(pxyin(2,1)-pxyout(1))*(-pxyin(1,2)+pxyout(2)) ) + &
(fin(2,2)*(-pxyin(1,1)+pxyout(1))*(-pxyin(1,2)+pxyout(2)) ) )

!la funcion per se

end function
end module