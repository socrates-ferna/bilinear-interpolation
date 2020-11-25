
program Interpolación_bilineal
!funcion que calcula f(x,y) en funcion de 4 puntos conocidos

use iso_Fortran_env, only: real64

real(real64) :: pxyin(2,2) !la primera es elsubindice y la segunda si es x o y
real(rela64) :: fin(2,2), pxyout(2), fout(1)

real(real64) :: dist11, dist12, dist21, dist22

fout(1)  = (fin(1)*(pxyin(2,1)-pxyout(1))*(pxyin(2,2)-pxyout(2))/((pxyin(2,1)-pxyin(1,1))*(pxyin(2,2)-pxyin(1,2)))