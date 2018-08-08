function plot_flownet(N,H,PSI,head,stream,Grid)

% author: Evan Ramos
% date:   28 October 2015

% Description: Plots a flownet with an equal number of equally spaced head
% contours and streamlines.

% Input: N = numer of contours & streamlines
% h = Grid.N by 1 column vector of heads
% PSI = Ny+1 by Nx+1 matrix containing the stream function
% head = string specifying the linestyle for the head contours
% stream = string specifying the linestyle for the streamlines
% Grid = structure containing all information about the grid.

[Xc,Yc] = meshgrid(Grid.xc,Grid.yc);     % 2D coords of cell centers
[Xp,Yp] = meshgrid(Grid.xf,Grid.yf);     % 2D coords of cell corners


% heads
contour(Xc,Yc,H,N,'LineColor','Blue','LineStyle',head)
hold on

% stream lines
contour(Xp,Yp,PSI,N,'LineColor','Red','LineStyle',stream,'Linewidth',1)

end