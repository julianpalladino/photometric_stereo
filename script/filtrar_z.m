function [z_filtrado] = filtrar(ancho, alto, z)
  z_filtrado = zeros(ancho * alto, 1);
  ancho = ancho *2;
  k = 1;
  for i = 0:alto-1
    z_filtrado(k) = z(i*ancho+1); k++;
    z_filtrado(k) = z(i*ancho+2); k++;
    saltear = 1;
    for j = i*ancho+3:i*ancho+3+ancho-4
      if saltear
        saltear = !saltear;
        continue
      endif
      z_filtrado(k) = z(j); k++;
      saltear = !saltear;
    endfor
  endfor
end
