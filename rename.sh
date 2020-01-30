

fun_T(){

	#double T = 50;
	sed "s/double T = $1;/double T = $2;/g" Muonium0.h>Muonium_$2.h
	mv Muonium_$2.h Muonium0.h

}


fun_T $1 $2
