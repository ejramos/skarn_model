function plot_pluton(pluton_dof,rock_dof,Grid)

% author: Evan J. Ramos
% date:   16 May 2016

% Description: This function outlines the pluton in the head/streamline
% plot figure.

% Input: pluton_dof --> degrees of freedom associated with the pluton
%                       boundary (cell centers)
%        rock_dof   --> degrees of freedom associated with the locations
%                       around the pluton 

outline             = zeros(Grid.Ny,Grid.Nx);
[Xc, Yc]            = meshgrid(Grid.xc,Grid.yc);

%Locating the minimum x, maximum x, and maximum y values 
outline(pluton_dof) = 1;
x_dists             = zeros(length(pluton_dof),1);
y_dists             = zeros(length(pluton_dof),1);
count               = 1;

for i = 1:Grid.Ny
    for j = 1:Grid.Nx
        
        if outline(i,j) == 1
            x_dists(count) = Grid.xc(j);
            y_dists(count) = Grid.Ly - Grid.yc(i);
            count          = count + 1;
        end
        
    end
end

xmin                = min(x_dists);
xmax                = max(x_dists);
ymin                = min(y_dists);
ymax                = max(y_dists);

%top of pluton
line([xmin xmax],[ymax-ymin ymax-ymin],...
     'Color',[0 0 0],'LineWidth',2)
%left side of pluton
if xmin ~= Grid.xc(1)
    line([xmin xmin],ymax - [ymin ymax],...
         'Color',[0 0 0],'LineWidth',2)
end
%right side of pluton
if xmax ~= Grid.xc(end)
    line([xmax xmax],ymax - [ymin ymax],...
         'Color',[0 0 0],'LineWidth',2)
end

%points around pluton
if ~isempty(rock_dof)
    for j = 1:length(rock_dof)
        htex = text(Xc(rock_dof(j)),Yc(rock_dof(j)),char(64+j));
        set(htex,'FontWeight','bold','Color','r','FontSize',20)
    end
end

end