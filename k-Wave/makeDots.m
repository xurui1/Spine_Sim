function Dots = makeDots(Nx,Ny,Ndots,coords)
%MAKEDOTS

% force integer values
Nx = round(Nx);
Ny = round(Ny);
Ndots = round(Ndots);
coords = round(coords);


% define literals
MAGNITUDE = 1;

% create empty matrix
Dots = zeros(Nx, Ny);

for i = 1:Ndots
    x = coords(i,1);
    y = coords(i,2);
    Dots(x,y) = MAGNITUDE;
end


