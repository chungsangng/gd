n=`expr $2`
  while expr $3 + 1 - $n
  do
    cp gd_3-psi.o$1$n.txt gd_3-psi.txt
    cp gd_3-aphi.o$1$n.txt gd_3-aphi.txt
    ./gd_3_B.x < gd_3.o$1$n.txt >  gd_3_B.o$1$n.txt
    rm gd_3-psi.txt
    rm gd_3-aphi.txt
    mv gd_3-BL.txt gd_3-BL.o$1$n.txt
    mv gd_3-BS.txt gd_3-BS.o$1$n.txt
    n=`expr $n + 1`
  done
