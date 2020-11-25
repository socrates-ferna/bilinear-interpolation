
program Interpolación_bilineal
!funcion que calcula f(x,y) en funcion de 4 puntos conocidos

use iso_Fortran_env, only: real64

real(real64) :: pxyin(2,2) !la primera es elsubindice y la segunda si es x o y
real(rela64) :: fin(2,2), pxyout(2), fout(1)


fout(1) = 1/((pxyin(2,1)-pxyin(1,1))*(pxyin(2,2)-pxyin(1,2)))*( &
(fin(1,1)* (pxyin(2,1)-pxyout(1))*(pxyin(2,2)-pxyout(2)) ) + &
(fin(2,1)*(-pxyin(1,1)+pxyout(1))*(pxyin(2,2)-pxyout(2)) ) + &
(fin(1,2)*(pxyin(2,1)-pxyout(1))*(-pxyin(1,2)+pxyout(2)) ) + &
(fin(2,2)*(-pxyin(1,1)+pxyout(1))*(-pxyin(1,2)+pxyout(2)) ) )
