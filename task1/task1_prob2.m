x0= [0 0 0 0 0 0 0];
options = optimset('LargeScale','off','Display','iter');
[x,fval,exitflag,output] = fmincon(@objfun2,x0,[],[],[],[],[],[],@confun2,options);