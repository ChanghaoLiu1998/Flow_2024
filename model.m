% -------------------------------------------------------------------------
% Set parameters 
% -------------------------------------------------------------------------

maxx1 = max(vcmplx(1).bc(:,1)); 
minx1 = min(vcmplx(1).bc(:,1));
maxx2 = max(vcmplx(1).bc(:,2)); 
minx2 = min(vcmplx(1).bc(:,2));
maxx3 = max(vcmplx(1).bc(:,3)); 
minx3 = min(vcmplx(1).bc(:,3));

smallTol = 1e-4;

for i = 1:vcmplx(4).num(4).val
    counting(i) = length(nonzeros(vcmplx(4).bndop(3).indx(i,:)));
end
dimension = max(counting);

for i = 1:vcmplx(3).num(3).val
    counting(i) = length(nonzeros(vcmplx(3).bndop(2).indx(i,:)));
end
dimensionedge = max(counting);

%--------------------------------------------------------------------------
% -------------------Assignment of voids-----------------------------------
%--------------------------------------------------------------------------
poresdata_void = xlsread("voids.xlsx");
Volume_void = poresdata_void(:, 1);

poresdata_VE = xlsread("VE.xlsx");
Volume_VE(:, 1) = poresdata_VE(:, 1);
Volume_VE(:, 2) = poresdata_VE(:, 2);

Length_2cell = vcmplx(2).vvol;
[a, idx] = sort(vcmplx(3).vvol, 'ascend');    
c = zeros(vcmplx(3).num(3).val, 1);
d = c;
grain = vcmplx(4).vvol;
grain_conductivity = 1 * ones(sum(vcmplx(4).num(4).val), dimension);

closest_indices = zeros(size(Volume_void));
used_indices = false(size(grain));

for i = 1:length(Volume_void)
    differences = abs(grain - Volume_void(i));
    differences(used_indices) = Inf;
    
    [~, closest_index] = min(differences);
    closest_indices(i) = closest_index;

    used_indices(closest_index) = true;
end

vf_void = (sum(grain(closest_indices))) / (maxx1 * maxx2 * maxx3);

thick = zeros(vcmplx(3).num(3).val, 1);

closest_indices = sortrows(closest_indices, 'ascend');
surfacearea = vcmplx(3).vvol;
for i = 1:length(Volume_void)
    area = surfacearea(nonzeros(vcmplx(4).bndop(3).indx(closest_indices(i), :)));
    number = nonzeros(vcmplx(4).bndop(3).indx(closest_indices(i), :));

    for j = 1:length(area)
        axislength(i, j) = norm([vcmplx(4).bc(closest_indices(i), 1), vcmplx(4).bc(closest_indices(i), 2), vcmplx(4).bc(closest_indices(i), 3)] ...
            - [vcmplx(3).bc(number(j), 1), vcmplx(3).bc(number(j), 2), vcmplx(3).bc(number(j), 3)]);
        
        grain_conductivity(closest_indices(i), j) = ((area(j) / pi)) / (8 * 1 * (10^(-3)));
        thick(number(j)) = axislength(i, j);
    end
end

thickness_edge = zeros(vcmplx(2).num(2).val, 1);
for i = 1:length(closest_indices)
    tmp = nonzeros(vcmplx(4).bndop(2).indx(closest_indices(i), :));
    for j = 1:length(tmp)
        thickness_edge(tmp(j)) = norm([vcmplx(4).bc(closest_indices(i), 1), vcmplx(4).bc(closest_indices(i), 2), vcmplx(4).bc(closest_indices(i), 3)] ...
            - [vcmplx(2).bc(tmp(j), 1), vcmplx(2).bc(tmp(j), 2), vcmplx(2).bc(tmp(j), 3)]);
    end
end

%--------------------------------------------------------------------------
% --------Assignment of VE--------------------------------------------------
%--------------------------------------------------------------------------

place = find(thick);

b = vcmplx(3).vvol;
for i = 1:length(closest_indices)
    indexxxx = nonzeros(vcmplx(4).bndop(3).indx(closest_indices(i), :));
    b(indexxxx) = 0;
end

selected_volumes = [];
accumulated_volume = 0;
vf_VE = 0;

random_values = [];  

while vf_VE < vf - vf_void
    random_value = rand;
    random_values = [random_values; random_value];  
    selected_volume = interp1(Volume_VE(:, 2), Volume_VE(:, 1), random_value, 'linear', 'extrap');
    
    accumulated_volume = accumulated_volume + selected_volume;
    selected_volumes = [selected_volumes; selected_volume];
    
    vf_VE = accumulated_volume / (maxx1 * maxx2 * maxx3);
end

porosity = vf_VE + vf_void;

for i = 1:length(selected_volumes)
    lower_bound = (pi / 4) * (5 * selected_volumes(i) / pi)^(2 / 3); 
    upper_bound = (pi / 4) * (20 * selected_volumes(i) / pi)^(2 / 3); 
    
    valid_b_values = b(b > lower_bound); 
    
    if ~isempty(valid_b_values)
        selected_2cell(i, 1) = min(valid_b_values);
        
        selected_idx = find(b == selected_2cell(i), 1, 'first'); 
        
        selected_index(i, 1) = selected_idx;
        
        b(selected_idx) = -Inf;
    else
        selected_2cell(i, 1) = NaN;
        selected_index(i, 1) = NaN;
    end
end

thickness = zeros(vcmplx(3).num(3).val, 1);
for i = 1:length(selected_index)
    thick(selected_index(i), 1) = selected_volumes(i) / selected_2cell(i);
    thickness(selected_index(i), 1) = selected_volumes(i) / selected_2cell(i);
end

for i = 1:vcmplx(3).num(3).val
    tmp = nonzeros(vcmplx(3).bndop(2).indx(i, :));
    if thick(i) ~= 0
        thickness_edge(tmp) = max(thickness_edge(tmp), thickness(i));
    end
end

% -------------------------------------------------------------------------
% STEP 3 (large pore allocation in pristine material)
% -------------------------------------------------------------------------

for i = 1:vcmplx(3).num(3).val
    c(i) = ((thick(i))^2) / (12 * 1 * (10^(-3)));
end

for i = 1:length(c)
    if c(i) == 0
        c(i) = 1;
    end
end

% -----------------------------------------------------------------------
% Build diagonals of the diffusion coefficient matrix
% -----------------------------------------------------------------------

k1mat = ones(2 * vcmplx(2).num(2).val, 1); 
for i = 1:length(thickness_edge)
    if thickness_edge(i) ~= 0
        k1mat(2 * i - 1, 1) = ((thickness_edge(i, 1) / 2)^2) / (8 * (1 * (10^(-3))));  
        k1mat(2 * i, 1) = ((thickness_edge(i, 1) / 2)^2) / (8 * (1 * (10^(-3))));
    else
        k1mat(2 * i - 1, 1) = 1; 
        k1mat(2 * i, 1) = 1; 
    end
end

k2mat = ones(sum(vcmplx(3).num(2).val), 1);  
ct = 0;
for i = 1:vcmplx(3).num(3).val
    for j = 1:size(vcmplx(3).bndop(2).indx, 2)
        if (vcmplx(3).bndop(2).indx(i, j) ~= 0)
            ct = ct + 1;
            k2mat(ct, 1) = c(i);
        end  
    end
end

k3mat = ones(sum(vcmplx(4).num(3).val), 1);   
ct = 0;
for i = 1:vcmplx(4).num(4).val
    num_nonzeros = length(nonzeros(vcmplx(4).bndop(3).indx(i, :)));
    for j = 1:num_nonzeros
        if (vcmplx(4).bndop(3).indx(i, j) ~= 0)
            ct = ct + 1;
            k3mat(ct, 1) = grain_conductivity(i, j);
        end  
    end
end

k = spdiags([k1mat; k2mat; k3mat], 0, fcmplx(2).num(2).val, fcmplx(2).num(2).val);


