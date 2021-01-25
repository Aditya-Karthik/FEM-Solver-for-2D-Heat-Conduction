% 2D STEADY STATE HEAT CONDUCTION %

% DOMAIN (Rectangular)
width = 10;
height = 10;
domain_coord = [0 0;width 0;width height;0 height]; % Clockwise from Left Bottom about origin

% _______________________________________________________________________________________________________________________________________________%
% MESHING

no_of_nodes_per_element = 3;
no_of_elements_x = 10;                                                           % Number of elements in the X-Direction
no_of_elements_y = 10;                                                           % Number of elements in the Y-Direction
no_of_elements = 2*no_of_elements_x*no_of_elements_y;
no_of_nodes_x = no_of_elements_x + 1;
no_of_nodes_y = no_of_elements_y + 1;

spacing_x = distance(domain_coord(1,:),domain_coord(2,:))/no_of_elements_x;
spacing_y = distance(domain_coord(1,:),domain_coord(4,:))/no_of_elements_y;

no_of_nodes = (no_of_nodes_x)*(no_of_nodes_y);
node_coord = zeros(no_of_nodes,2);

% Assigning co-ordinates to nodes
for i=1:no_of_nodes_x                           
    for j=1:no_of_nodes_y
        node_coord(i+((j-1)*no_of_nodes_x),1) = (domain_coord(1,1)+(i-1)*spacing_x);
        node_coord(i+((j-1)*no_of_nodes_x),2) = (domain_coord(1,2)+(j-1)*spacing_y);
    end
end

% Element connectivity
element_connectivity = zeros((no_of_elements),3);
counter = 1;
for i=1:no_of_elements_y
    for j=1:no_of_elements_x
        element_connectivity(counter,1) = j + (i-1)*no_of_nodes_x;
        element_connectivity(counter,2) = 1 + j + (i-1)*no_of_nodes_x;        
        element_connectivity(counter,3) = 1 + j + i*no_of_nodes_x;
        counter = counter + 1;
    end
    
    for k = 1:no_of_elements_x
        element_connectivity(counter,1) = k + (i-1)*no_of_nodes_x;
        element_connectivity(counter,2) = 1 + k + i*no_of_nodes_x;        
        element_connectivity(counter,3) = k + i*no_of_nodes_x;
        counter = counter + 1;
    end
end


left_nodes = 1:no_of_nodes_x:1+(no_of_elements_y*no_of_nodes_x);
right_nodes = no_of_nodes_x:no_of_nodes_x:no_of_nodes_x*(no_of_elements_y+1);
bottom_nodes = 1:1:no_of_nodes_x;
top_nodes = (no_of_nodes_x*no_of_elements_y)+1:1:no_of_nodes_x*(no_of_elements_y+1);

left_elements = (no_of_nodes_x):(2*no_of_elements_x):no_of_elements-no_of_elements_x+1;
right_elements = no_of_elements_x:(2*no_of_elements_x):no_of_elements-no_of_elements_x;
bottom_elements = 1:1:no_of_elements_x;
top_elements = no_of_elements-no_of_elements_x+1:1:no_of_elements;
elements = {bottom_elements;right_elements;top_elements;left_elements};
nodes = {bottom_nodes;right_nodes;top_nodes;left_nodes};
A = sparse(no_of_elements);                                                              % Area of elements

% _______________________________________________________________________________________________________________________________________________%
% MATERIAL PROPERTIES
k_e = zeros(no_of_elements);     % Thermal Conductivity
k_val = 1;                       % W/mK
k_e(:) = k_val;
k_nodes = zeros(no_of_nodes,1);
k_nodes(:) = k_val;
% _______________________________________________________________________________________________________________________________________________%
 
% BOUNDARY CONDITIONS (Temperature, Heat Flux or Convection BC)
    BC = zeros(4);              
    % BC(i) : i = 1 (Bottom)
    % BC(i) : i = 2 (Right)
    % BC(i) : i = 3 (Top)
    % BC(i) : i = 4 (Left)
    % BC = 1 ---> Temperature BC
    % BC = 2 ---> Heat Flux BC
    % BC = 3 ---> Convection BC
      BC(1) = 2;
      BC(2) = 1;
      BC(3) = 2;
      BC(4) = 1;
      
% Bottom   
    T_bottom = rand;     % Temperature  
    q_bottom = 0;       % Heat Flux
    h_bottom = rand;     % Convection
    T_inf_bottom = rand; % Convection 
   
   
% Right   
    T_right = 100;       % Temperature  
    q_right = rand;      % Heat Flux
    h_right = rand;      % Convection
    T_inf_right = rand;  % Convection    
    

% Top  
    T_top = rand;        % Temperature  
    q_top = 0;           % Heat Flux
    h_top = rand;        % Convection
    T_inf_top = rand;    % Convection 
    
 
% Left   
    T_left = 30;        % Temperature  
    q_left = rand;       % Heat Flux
    h_left = rand;       % Convection
    T_inf_left = rand;   % Convection 
    
   
T_wall = [T_bottom T_right T_top T_left];
q = [q_bottom q_right q_top q_left];
h = [h_bottom h_right h_top h_left];
T_inf = [T_inf_bottom T_inf_right T_inf_top T_inf_left];

% Heat Generation Elements   
heat_gen_elements = [];    % Central Elements have been assumed to be heat sources
Q = 0;
   
% _______________________________________________________________________________________________________________________________________________%
  
% SOLUTION

T = sparse(no_of_nodes,1);

Kg = zeros(no_of_nodes,no_of_nodes);

for(k=1:no_of_elements)
    
    for(i=1:no_of_nodes_per_element)
        for(j=1:no_of_nodes_per_element)
            x(i,j) = node_coord(element_connectivity(k,i),1) - node_coord(element_connectivity(k,j),1);
            y(i,j) = node_coord(element_connectivity(k,i),2) - node_coord(element_connectivity(k,j),2);
        end                
    end
    
    J = [x(1,3) y(1,3);x(2,3) y(2,3)]; 
    A_e = 0.5*det(J);
    A(k) = A_e;
    
    B = [y(2,3) y(3,1) y(1,2) ;
         x(3,2) x(1,3) x(2,1) ;];
     
    B = (1/det(J))*B;
    B_T = transpose(B);
    
    K_e = k_e(k)*A_e*(B_T*B);
    
    temp = [(element_connectivity(k,1));(element_connectivity(k,2));(element_connectivity(k,3))];    
    Kg(temp,temp) = Kg(temp,temp) + K_e;    
end

RHS = zeros(no_of_nodes,1);

% Heat Generation
for(i=heat_gen_elements)
       
    r_Q = (Q*A(i)/3)*[1;1;1];
    temp = [(element_connectivity(i,1));(element_connectivity(i,2));(element_connectivity(i,3))];
    RHS(temp) = RHS(temp) + r_Q;
end

% Heat Flux BCs
edges = [1 2;2 3;3 2;1 3];
M{1} = [1;1;0];
M{2} = [0;1;1];
M{3} = [0;1;1];
M{4} = [1;0;1];
for(j=1:1:4)
    if(BC(j) ~= 2)        
        continue 
    end       
    for(i=elements{j})        
        l_edge = distance(node_coord(element_connectivity(i,edges(j,1)),:),node_coord(element_connectivity(i,edges(j,2)),:));
        r_q = ((q(j)*l_edge)/2)*M{j};
        temp = [(element_connectivity(i,1));(element_connectivity(i,2));(element_connectivity(i,3))];
        RHS(temp) = RHS(temp) - r_q;
    end    
end

% Convection BCs

edges = [1 2;2 3;3 2;1 3];
M{1} = [1;1;0];
M{2} = [0;1;1];
M{3} = [0;1;1];
M{4} = [1;0;1];
for(j=1:1:4)
    if(BC(j) ~= 3)        
        continue 
    end       
    for(i=elements{j})   
        l_edge = distance(node_coord(element_connectivity(i,edges(j,1)),:),node_coord(element_connectivity(i,edges(j,2)),:));
        r_inf = ((h(j)*T_inf(j)*l_edge)/2)*M{j};
        temp = [(element_connectivity(i,1));(element_connectivity(i,2));(element_connectivity(i,3))];
        RHS(temp) = RHS(temp) + r_inf;
    end    
end

edges = [1 2;2 3;3 2;1 3];
M{1} = [1 2 0;2 1 0;0 0 0];
M{2} = [0 0 0;0 2 1;0 1 2];
M{3} = [0 0 0;0 2 1;0 1 2];
M{4} = [1 0 2;0 0 0;2 0 1];
for(j=1:1:4)
    if(BC(j) ~= 3)         
        continue 
    end       
    for(i=elements{j})    
        l_edge = distance(node_coord(element_connectivity(i,edges(j,1)),:),node_coord(element_connectivity(i,edges(j,2)),:));
        h_t = ((h(j)*l_edge)/6)*M{j};
        temp = [(element_connectivity(i,1));(element_connectivity(i,2));(element_connectivity(i,3))];
        Kg(temp,temp) = Kg(temp,temp) + h_t;
    end    
end


% Temperature BCs
known_nodes = [];
for(i=1:4)
   if(BC(i) ~= 1)               
       continue 
   end   
   T(nodes{i}) = T_wall(i);
   known_nodes = cat(2,known_nodes,nodes{i});
   temp_edge = zeros(no_of_nodes,1);
   temp_edge(nodes{i}) = T_wall(i);
   RHS = RHS - (Kg*temp_edge);
end


unknown_nodes = setdiff((1:no_of_nodes),known_nodes);
RHS = RHS(unknown_nodes);
T(unknown_nodes) = Kg(unknown_nodes,unknown_nodes)\RHS;
% _______________________________________________________________________________________________________________________________________________%
% POST PROCESSING

% Mesh
figure('Name','Mesh','NumberTitle','off');
trimesh(element_connectivity,node_coord(:,1),node_coord(:,2));
hold on;

% Heat Flux
boundary_nodes = union(bottom_nodes,right_nodes);
boundary_nodes = union(boundary_nodes,top_nodes);
boundary_nodes = union(boundary_nodes,left_nodes);
heat_flux_nodes = setdiff(1:no_of_nodes,boundary_nodes);
no_of_heat_flux_nodes = (4 + no_of_nodes - 2*(no_of_nodes_x+no_of_nodes_y));
q_x = zeros(no_of_nodes,1);
q_x(heat_flux_nodes,1) = -k_nodes(heat_flux_nodes).*(T(heat_flux_nodes+1) - T(heat_flux_nodes))/(spacing_x);
q_y = zeros(no_of_nodes,1);
q_y(heat_flux_nodes,1) = -k_nodes(heat_flux_nodes).*(T(heat_flux_nodes+no_of_nodes_x) - T(heat_flux_nodes))/(spacing_y);
figure('Name','Heat Flux','NumberTitle','off');
axis([0 width 0 height]);
hold on;
quiver(node_coord(heat_flux_nodes,1),node_coord(heat_flux_nodes,2),q_x(heat_flux_nodes),q_y(heat_flux_nodes));

% Temperature Distribution
figure('Name','Temperature Distribution','NumberTitle','off');
x = node_coord(:,1);
y = node_coord(:,2);
z = T(1:1:no_of_nodes);
scatter3(x,y,z,'*');
hold on;











