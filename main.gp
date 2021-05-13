dsa_pub = [Mod(16, 2359457956956300567038105751718832373634513504534942514002305419604815592073181005834416173), 589864489239075141759526437929708093408628376133735628500576354901203898018295251458604043, 2028727269671031475103905404250865899391487240939480351378663127451217489613162734122924934];

check(s,dsa_pub) = {
  my(h,r,g,q,X);
  [h,r,s] = s;
  [g,q,X] = dsa_pub;
  X=Mod(X,g.mod);
  lift( (g^h*X^r)^(1/s % q) ) % q == r;
}

/************************************************************************************/
/*************** tentative de résoudre le challenge en s'inspirant du td 7 **********/

\\ Un parcours de la liste à la recherche d'une collision sur la valeur de r s'étant révélé infructueux (aucun message n'a été signé avec le même k), on va chercher un autre type de collision, faisant intervenir un a fixé.


\\ On cherche a tel que deux messages m et m' ont été signés avec k et k+a.
\\ C'est à dire, r(m')=g^(k+a)=g^a * g^k = g^a *(r(m))
find_couple(l)={
	my(i,j,a,g);
	g=dsa_pub[1];
	borne=2^20;
	gpow=g;
	a=18172;
	while( a<borne,
	i=1;
	gpow=g^a;
	\\print("Dans le while, gpow = ",gpow," ...\n");
	while(i< #l,
		 j=i+1;
		 while(j<= #l,
		 	    if(l[j][2]==(l[i][2]*gpow),/*print("trouvé !\n")*/;return([i,j,a]));
			    j=j+1);
			    i=i+1);
			    
			    a=a+1;
	);
	print("Aucun couple ne correspond...\n");
	return ([-1,-1,-1]);

}



\\ Une fois que l'on a la collision, on retrouve k (voir les calculs plus bas)
find_k(i,j,a,list)={
	my(g,q,X);
	[g,q,X] = dsa_pub;	
	s_k2=list[j][3];
	s_k1=list[i][3];
	r_k2=list[j][2];
	r_k1=list[i][2];
	h_k2=list[j][1];
	h_k1=list[i][1];
	print("g : ", g);
	print("q : ", q);
	my(k);
	k=(g^a*s_k1 - s_k2)^(-1) * (h_k1*g^a - h_k2 + a*s_k2);
	\\print("k : ", lift(k));
	\\lift(k);
	return (lift(k));
}


/*
Calcul de k à partir de r,r' et s,s' obtenus avec deux messages m et m' signés respectivement avec k et k+a.
ks * g^a  - (k+a)*s'   =  g^a*h(m) + xr * g^a  - h(m') - xr'
ks * g^a  - (k+a)*s'   =  g^a*h(m) + xr * g^a  - h(m') - x (g^a * r)
ks * g^a  - (k+a)*s'   =  g^a*h(m) - h(m')
k*s *g^a  - ks' - a*s' =  g^a*h(m) - h(m')
k (s*g^a - s') - a*s' =  g^a*h(m) - h(m')
==> k (s*g^a - s') =  g^a*h(m) - h(m') + a*s'
Un calcul qui semble correct car seul k est inconnu, et les autres entités sont inversibles.
*/



\\ A partir de k et du message m signé avec k, on retrouve x
find_x(i,k,list)={
	my(q,x,s,r,h);
	q = dsa_pub[2];
	s=list[i][3];
	r=list[i][2];
	h=list[i][1];
	x=lift(Mod(lift((s*k-h)*r^(-1)),q));
}



\\ Cependant, cette méthode, discutée avec Marina, ne fonctionne pas. Il doit y avoir une erreur de calcul quelque part.
\\ Pour l'instant on laisse cette idée de côté, et on va tenter un brute force sur k, qui n'est pas dans un intervalle trop grand.


list=readvec("input.txt");
\\[i,j,a]=find_couple(list);
\\print("i : "i, "\nj : ", j,"\na : ",a);
\\ka=find_k(i,j,a,list);
\\print(ka);
\\print(find_x(i,ka,list));


/***********************************************************************************/
/***************************** Méthode moins délicate ******************************/


\\ Comme il n'existe pas deux messages ayant été signés avec le même k, on passe en mode brute force.
\\ On va tirer aléatoirement des k dans l'espace indiqué, calculer g^k et chercher dans la liste un r qui correspond.
\\ On aura alors la valeur de k, et on pourra donc retrouver x.


\\ première version avec une liste, trop longue
find_k_brute_force(l)={
	my(i,g,k);
	g=dsa_pub[1];
	gpow=g;
	trouve=0;
	while( trouve==0,
	       k==random(10^10-1)+1;
	       gpow=g^k;
	       for(i=1, i< #l,
		 	    if(l[i][2]==gpow,trouve=1;return([i,k]));
		);
	);
}


\\ On essaye d'utiliser une map pour gagner en rapidité
find_k_brute_force_bis(l)={
	map =Map();
	for(i=1,#l,mapput(map,l[i][2],l[i]));
	trouve=0;
	[g,q,X]=dsa_pub;
	while(trouve==0,
		k=random(10^10-1)+1;
			r_k=lift(Mod(lift(g^k),q));
			if(mapisdefined(map,r_k),trouve=1;[h,r,s]=mapget(map,r_k)));
	\\Maintenant qu'on a la valeur de k pour un message donné, on retrouve x :
	x=lift(Mod(lift((s*k-h)*r^(-1)),q));
	print(x);
}


\\print("ici !");
find_k_brute_force_bis(list);
\\print(x);

