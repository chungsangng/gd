n=`expr $2`
  while expr $3 + 1 - $n
  do
    cp gd_3-psi.o$1$n.txt gd_3_chk-psi0.txt
    cp gd_3-aphi.o$1$n.txt gd_3_chk-aphi0.txt
    ./gd_3_chk.x < gd_3.o$1$n.txt >  gd_3_chk.o$1$n.txt
    rm gd_3_chk-psi0.txt
    rm gd_3_chk-aphi0.txt
    mv gd_3-psi.txt gd_3_chk-psi.o$1$n.txt
    mv gd_3-aphi.txt gd_3_chk-aphi.o$1$n.txt
    n=`expr $n + 1`
  done
