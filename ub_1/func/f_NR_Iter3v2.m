function [xn]=f_NR_Iter(func,x0,B,FR)

maxiter=500;
nx=size(FR,1);
rpms=(B.omega/(2*pi/60));
if rpms < 300
    nloads=1; 
    disp('lower than 300 rpm')
elseif rpms >= 300 && rpms < 5000
    disp('now 300-5000 rpm')
    nloads=2;
elseif rpms >= 5000 && rpms < 70000
    disp('over than 5000-70,000 rpm')
    nloads=5;
end

  for jj=1:nloads,
    fact=jj/nloads;
       fact2=fact^2;
        Rscale=diag(ones(nx,1)*fact2);       % scaling matrix
        FRloads=Rscale*FR;    
       disp(['****** Loading step: ' num2str(fact2),'/1  ******'])
   
    for ii=1:maxiter
    % Calculation of the bearing stiffness matrix------------------------------- 
    % bearing force vector at given operation point 
    [EqS]=feval(func,x0,B,FRloads);
    %h=eps^(1/2);   
    h=1e-9;
    % loop over columns of the Jacobian
    for jcol=1:length(x0)
        xh = x0; % Set all elements to x0 values
        % Perturb jcol element of Disph
         xh(jcol) = x0(jcol)+h;
        % Calculate variated Force Fh
        %[Eqh]=feval(func,xh,B,FR);
        [EqhS]=feval(func,xh,B,FRloads);
        % Form column of Jacobian by finite differences
        Jac(:,jcol) = (EqhS-EqS)/h;
    end
    % Solve new displacements
    xn=x0-Jac\EqS; % using temporary variable ee
    x0=xn;
    if ii == maxiter
        disp('Maximum iteration rounds found. Did not iterate to solution')
        break
    end
    % Set convergence criteria
    if ii==1; Re=0.001; end; %ANSYS uses factor 0.005 %was 0.001
    %disp(['Iteration number: ' num2str(ii)])
    %disp(['Convergence Norm = ' num2str(norm(EqS),'%0.5e') '    Criterion =' num2str(Re,'%0.5e')] )
    if norm(EqS)<Re
       %disp(['Solution CONVERGED at iteration number ' num2str(ii)])
        %disp(['Convergence Norm = ' num2str(norm(Eq),'%0.5e') '    Criterion = ' num2str(Re,'%0.5e')] )
        %disp(Disp') % display coordinates of the end tip
        break
    end
   end
  end
           
    
end
