x0 = [0 0 0 0];
options = optimset('LargeScale','off','Display','iter');
[x,fval,exitflag,output] = fmincon(@objfun1,x0,[],[],[],[],[],[],@confun1,options);