
set terminal png
set output "result.png"
set key right outside
f="stat.log"
p f u 1:2 t "potential E"   ,\
  f u 1:3 t "kinetic E" ,\
  f u 1:4 t "all E"       ,\
  f u 1:5 t "temperature" ,\
  f u 1:6 t "pressure"    

