function D = Cost_function(stepAt, theta)
global optTol oBP BP rBP PosScale NegScale Domain
oBP = stepAt;
BP  = log10(stepAt);%  Breakpoint of the function
rBP = floor((BP-eps)*100)/100; % Rounded Breakpoint of the function
SymmetricDomain = all(stepAt == 10.^(mean(Domain,2))); % Matlba's multidimenstional integrators work very fast on symmetric domains
if isnan(PosScale)
    [PosScale, NegScale]  = Pos_Class_Ratio();
end
switch size(Domain,1)
    case 1
        %% Domain has one dimension integrate over it
        
        dx = @(x)Cost_1D(theta, x, []);
        D  = quadgk(dx,Domain(1,1), Domain(1,2),'Waypoints', BP, 'AbsTol', 0, 'RelTol', 10^(optTol+1));
        
    case 2
        %% Domain has two dimention integrate over it
        dx = @(x,y)GeneralCost(theta, x, y, []);
        D  = Pej_quad2d(dx,Domain(1,1), Domain(1,2),Domain(2,1), Domain(2,2),10^(optTol+1), BP );
        
    case 3
        %% Domain has three dimention integrate over it
        dz = @(z)Cost_3D(theta, Domain(1:2,:), z, []);
        D  = Pej_Quad(dz,Domain(3,1), Domain(3,2), 10^(optTol+1));
       
    otherwise
        %% Domain has more than two dimension integrate over it the last two
        
end

end

function dz = Cost_3D(theta, Domain, Z, Constants)
%% Total cost for a give parameter set
global optTol BP
RelativeTol = 10.^(optTol-length(Constants)+1); % Calculate inner integrals with lower persision
dz = nan(size(Z));
for i = 1:size(Z,2)
    dy = @(x,y)GeneralCost(theta, x, y, Z(i)*ones(size(x)));
    dz(i)  = Pej_quad2d(dy,Domain(1,1), Domain(1,2),Domain(2,1), Domain(2,2),RelativeTol, BP);
end
end

function dx = Cost_1D(theta, Xl, Constants)
global Functioninstance FunctionArray FunctionWeight oBP PosScale NegScale
Px = BackProb_log([Xl; repmat(Constants,1,size(Xl,2))]);
X       = 10.^Xl; clear Xl
theta   = 10.^theta;

if ~isempty(Constants)
    % There are Other fixed variables provided.
    Constants = 10.^Constants;
    X = [X; repmat(Constants,1,size(X,2))];
end
dx = zeros(1, size(X,2));
for f = 1:size(FunctionArray,1)
    Functioninstance = zeros(size(FunctionArray,2), size(FunctionArray,3));
    Functioninstance(:,:) = FunctionArray(f, :,:); % Note of stupidity: Functioninstance is defined as global in both IdealFun(X) and Evaluate_FunctionCnt_Minimal(X, theta)
    Ix = Evaluate_Function_Minimal(X>oBP);    In = ~Ix;
    CAn = log10(Evaluate_FunctionCnt_Minimal(X, theta));
    CAn(Ix) =  -PosScale(f) * FunctionWeight(f) * CAn(Ix);
    CAn(In) =   NegScale(f) * FunctionWeight(f) * CAn(In);
    dx = dx + CAn;
end

dx = dx .* Px;
end


function dx = GeneralCost(theta, Xl, Yl, Zl)
global Functioninstance FunctionArray FunctionWeight oBP PosScale NegScale
myXl = [Xl(:)'; Yl(:)'; Zl(:)'];
Px = BackProb_log(myXl);
X       = 10.^myXl; clear myXl
theta   = 10.^theta;

dx = zeros(1,size(X,2));
for f = 1:size(FunctionArray,1)
    Functioninstance = zeros(size(FunctionArray,2), size(FunctionArray,3));
    Functioninstance(:,:) = FunctionArray(f, :,:); % Note of stupidity: Functioninstance is defined as global in both IdealFun(X) and Evaluate_FunctionCnt_Minimal(X, theta)
    Ix = Evaluate_Function_Minimal(X>oBP);    In = ~Ix;
    CAn = log10(Evaluate_FunctionCnt_Minimal(X, theta));
    CAn(Ix) =  -PosScale(f) * FunctionWeight(f) * CAn(Ix);
    CAn(In) =   NegScale(f) * FunctionWeight(f) * CAn(In);
    dx = dx + CAn;
end
 
dx = dx .* Px;
dx = reshape(dx, size(Xl));
end


function X = Pej_Quad(dx, DomainL, DomainH, RelativeTol)
global rBP
X = quadl(dx,DomainL, rBP,RelativeTol ) + quadl(dx,rBP+.001, DomainH,RelativeTol );
end

% function dy = Cost_2D(theta, Domain,  Y, Constants)
% %% Total cost for a give parameter set
% global optTol BP
% RelativeTol = 10.^(optTol-length(Constants)+1); % Calculate inner integrals with lower persision
% dy = nan(size(Y));
% 
% for i = 1:size(Y,2)
%     dx = @(x)Cost_1D(theta, x, [Y(i); Constants]);
%     dy(i) = quadgk(dx,Domain(1,1), Domain(1,2),'Waypoints', BP, 'AbsTol', 0, 'RelTol', RelativeTol);
% end
% end