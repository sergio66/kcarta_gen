beta = 0.00001 : 0.1 : 1;

mu = 1; xyz = mu*beta - log(mu*beta + 1); xyx = xyz./beta./beta; y1 = xyz;
mu = 0; xyz = mu*beta - log(mu*beta + 1); xyx = xyz./beta./beta; y0 = xyz;

plot(beta,y1-y0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

beta = 0.00001 : 0.1 : 1;

mu = 1; xyz = beta - log(beta + 1); xyx = xyz./beta./beta; y1 = xyz;
mu = 0; xyz = beta - log(1);        xyx = xyz./beta./beta; y0 = xyz;

plot(beta,y1-y0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xyz = beta - log(beta + 1);

xyz = xyz./beta./beta;

plot(beta, xyz)
