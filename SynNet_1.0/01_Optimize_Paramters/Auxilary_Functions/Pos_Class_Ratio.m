function [PosScale, NegScale] = Pos_Class_Ratio()
global Domain Functioninstance FunctionArray BP optTol_Mass Class_Mass_Balance

if ~Class_Mass_Balance
    disp('Class mass balancing deactive.')
    PosScale = ones(size(FunctionArray,1),1);
    NegScale = ones(size(FunctionArray,1),1);
    return
end
fprintf('Calculating Positive/Negative class ratio')
optTol_Mass = 1E-4;
m1 = ones(size(FunctionArray,1),1);
m2 = ones(size(FunctionArray,1),1);

for f = 1:size(FunctionArray,1)
    fprintf(' F%d', f);
    Functioninstance = zeros(size(FunctionArray,2), size(FunctionArray,3));
    Functioninstance(:,:) = FunctionArray(f, :,:); % Note of stupidity: Functioninstance is defined as global in both IdealFun(X) and Evaluate_FunctionCnt_Minimal(X, theta)
    
    
    switch max(abs(Functioninstance(:)))
        case 1
            %% Domain has one dimension integrate over it
            m1(f) = quadgk(@(x)Ix_1D(x, [])       ,Domain(1,1), Domain(1,2),'Waypoints', BP, 'AbsTol', 0, 'RelTol', optTol_Mass);
            if m1(f)>.49;
                m2(f) = quadgk(@(x)In_1D(x, [])       ,Domain(1,1), Domain(1,2),'Waypoints', BP, 'AbsTol', 0, 'RelTol', optTol_Mass);
            end
        case 2
            %% Domain has two dimention integrate over it
            dm1 = @(x,y)Ix_General(x, y, []);
            dm2 = @(x,y)In_General(x, y, []);
            m1(f)  = Pej_quad2d(dm1,Domain(1,1), Domain(1,2),Domain(2,1), Domain(2,2),optTol_Mass, BP );
            if m1(f)>.49;
                m2(f)  = Pej_quad2d(dm2,Domain(1,1), Domain(1,2),Domain(2,1), Domain(2,2),optTol_Mass, BP );
            end
        case 3
            %% Domain has three dimention integrate over it
            dm1    = @(z)Ix_3D(Domain(1:2,:), z, []);
            dm2    = @(z)In_3D(Domain(1:2,:), z, []);
            m1(f)  = Pej_Quad(dm1,Domain(3,1), Domain(3,2), optTol_Mass);
            if m1(f)>.49;
                m2(f)  = Pej_Quad(dm2,Domain(3,1), Domain(3,2), optTol_Mass);
            end
            
        otherwise
            %% Domain has more than two dimension integrate over it the last two
            
    end
    m1(f) = max(m1(f), 0); % Just in case integrator produces tiny negative values
    m2(f) = max(m2(f), 0); % Just in case integrator produces tiny negative values
    if m1(f)>m2(f)
        % edit the large one to complement the small one, cuz small one is more
        % accurate.
        m1(f) = 1-m2(f);
    else
        m2(f) = 1-m1(f);
    end
end

PosScale = 1 ./ (m1) ;
NegScale = 1 ./ (m2);
fprintf('. done!\n');
end


function dx = Ix_General(Xl, Yl, Zl)
global BP

myXl = [Xl(:)'; Yl(:)'; Zl(:)'];
dx = Evaluate_Function_Minimal(myXl>BP).*BackProb_log(myXl);
dx = reshape(dx, size(Xl));
end

function dx = In_General(Xl, Yl, Zl)
global BP

myXl = [Xl(:)'; Yl(:)'; Zl(:)'];
dx = (1 - Evaluate_Function_Minimal(myXl>BP)).*BackProb_log(myXl);
dx = reshape(dx, size(Xl));
end


function dx = In_1D(Xl, Constants)
global BP

if ~isempty(Constants)
    % There are Other fixed variables provided.
    Xl = [Xl; repmat(Constants,1,size(Xl,2))];
end
dx = (1 - Evaluate_Function_Minimal(Xl>BP)).*BackProb_log(Xl);
end

function dz = In_3D(xyDomain, zl, Constants)
%% Total cost for a give parameter set
global optTol_Mass BP
dz = nan(size(zl));

for i = 1:size(zl,2)
    dy = @(x,y)In_General(x, y, zl(i)*ones(size(x)));
    dz(i) = Pej_quad2d(dy, xyDomain(1,1), xyDomain(1,2),xyDomain(2,1), xyDomain(2,2),optTol_Mass, BP);
end
end

function dx = Ix_1D(Xl, Constants)
global BP

if ~isempty(Constants)
    % There are Other fixed variables provided.
    Xl = [Xl; repmat(Constants,1,size(Xl,2))];
end
dx = Evaluate_Function_Minimal(Xl>BP).*BackProb_log(Xl);
end

function dz = Ix_3D(xyDomain, zl, Constants)
%% Total cost for a give parameter set
global optTol_Mass BP
dz = nan(size(zl));

for i = 1:size(zl,2)
    dy = @(x,y)Ix_General(x, y, zl(i)*ones(size(x)));
    dz(i) = Pej_quad2d(dy, xyDomain(1,1), xyDomain(1,2),xyDomain(2,1), xyDomain(2,2),optTol_Mass, BP);
end
end


function X = Pej_Quad(dx, DomainL, DomainH, RelativeTol)
global rBP BP
X = quadl(dx,DomainL, rBP,RelativeTol ) + quadl(dx,rBP+.001, DomainH,RelativeTol );
%X = quadgk(dx, DomainL, DomainH,'Waypoints', BP, 'AbsTol', 0, 'RelTol', RelativeTol);
end


% function dy = In_2D(xDomain, Yl, Constants)
% %% Total cost for a give parameter set
% global optTol_Mass
% dy = nan(size(Yl));
%
% for i = 1:size(Yl,2)
%     dy(i) = Pej_Quad(@(x)In_1D(x, [Yl(i); Constants]),xDomain(1), xDomain(2), optTol_Mass);
% end
% end

% function dy = Ix_2D(xDomain, Yl, Constants)
% %% Total cost for a give parameter set
% global optTol_Mass
% dy = nan(size(Yl));
%
% for i = 1:size(Yl,2)
%     dy(i) = Pej_Quad(@(x)Ix_1D(x, [Yl(i); Constants]),xDomain(1), xDomain(2), optTol_Mass);
% end
% end