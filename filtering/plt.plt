
set multiplot layout 2,2 title "Фильтрация изображений" scale 1,1
set palette gray
#set terminal qt
set yrange [*:*] reverse

set title "Исходное изображение"
#set size 1,1
set cbrange [0:255]
plot "image.txt" matrix with image

set title "Восстановленное изображение"
# set cbrange [0:255]
plot "image4.txt" matrix with image

set title "Спектр изображения"
set autoscale cb
#set size 1,1
plot "image2.txt" matrix with image

set title "Отфильтрованный спектр"
plot "image3.txt" matrix with image



unset multiplot

pause mouse close
