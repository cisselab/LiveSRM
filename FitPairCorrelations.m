function [rho_null,sigma_null,rho_exp,sigma_exp,A,xi]=FitPairCorrelations(r,g)

    fnull = @(parameters)null_model_vector_error(r,g,parameters(1),parameters(2));
    fexp = @(parameters)exponential_model_vector_error(r,g,parameters(1),parameters(2),parameters(3),parameters(4));

    null_parameters_initial_guess=[1,1];
    null_parameters_lower_bounds = [0,0];
    null_parameters_upper_bounds = [Inf,Inf];
    
    exp_parameters_initial_guess = [1,1,1,1];
    exp_parameters_lower_bounds = [0,0,0,0];
    exp_parameters_upper_bounds = [Inf,Inf,Inf,Inf];
    
    
    parametersnull=lsqnonlin(fnull,null_parameters_initial_guess,null_parameters_lower_bounds,null_parameters_upper_bounds);
    rho_null=parametersnull(1);
    sigma_null=parametersnull(2);
    
    parametersexp=lsqnonlin(fexp,exp_parameters_initial_guess,exp_parameters_lower_bounds,exp_parameters_upper_bounds);
    rho_exp=parametersexp(1);
    sigma_exp= parametersexp(2);
    A = parametersexp(3);
    xi = parametersexp(4);

%     figure
%     plot(r,g,'.k')
%     hold on
%     plot(r,1+1/4/pi/sigma_null^2/rho_null*exp(-r.^2/4/sigma_null^2),'g')
%     plot(r,1+1/4/pi/sigma_exp^2/rho_exp*exp(-r.^2/4/sigma_exp^2)+A*exp(-r/xi).*(erf(r/2/sigma_exp-sigma_exp/xi)+erf(sigma_exp/xi)),'r')
%     

end

function vector_error=null_model_vector_error(r,g,rho,sigma)
    vector_error=g-(1+1/4/pi/sigma^2/rho*exp(-r.^2/4/sigma^2));
end

function vector_error = exponential_model_vector_error(r,g,rho,sigma,A,xi)

    vector_error = g-(1+1/4/pi/sigma^2/rho*exp(-r.^2/4/sigma^2)+A*exp(-r/xi).*(erf(r/2/sigma-sigma/xi)+erf(sigma/xi)));

end