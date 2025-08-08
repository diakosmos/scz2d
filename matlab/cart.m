
iterations = 1000;

nu = 1;     % ???/

imax = 15; % zones in x direction
jmax = 15;  % zones in y dithrection

dt = 0.001;

dx = 1/(jmax-1);
dy = 1 / (imax-1);

x1d = linspace(0,1,imax);
y1d = linspace(0,1,jmax);
[x,y]=meshgrid(x1d,y1d);

func = sin(x) .* cos(2 .* y) + y.^2;   % just some random initialization
for i=1:imax
    for j=1:jmax
        func2(i,j) = func(i,jmax + 1 - j);
    end
end
func = func2;

while iterations > 0
    iterations = iterations - 1;

            d_x_func = diff(func,1,1)./dx;
         
    for j=1:jmax
        d_x_func(1,j) = 0; 
        d_x_func(imax,j) = 0;    
    end
    
    for i=1:imax
        for j=(1+1):jmax
            d_y_func(i,j) = (func(i,j) - func(i,j-1))./dy;
        end
    end
    for i=1:imax
        d_y_func(i,1) = 0;
        d_y_func(i,jmax) = 0;
    end
    
    for i=1:imax
        for j=2:(jmax-1)
            d_x_x_func(i,j) = (func(i,j+1) - 2 * func(i,j) + func(i,j-1)) ./ (dx * dx);
        end
    end
    for i=1:imax
        d_x_x_func(i,1) = (func(i,2) - func(i,1)) ./ dx ./ (dx / 2);
        d_x_x_func(i,jmax) = (func(i,jmax) - func(i,jmax-1)) ./ dx ./ (dx /2);
    end
    for i=2:(imax-1)
        for j=1:jmax            
            d_y_y_func(i,j) = (func(i+1,j) - 2 * func(i,j) + func(i-1,j)) ./ (dy * dy);
        end
    end
    for j=1:jmax
        d_y_y_func(1,j) = (func(2,j) - func(1,j)) ./ dy ./ (dy/2);
        d_y_y_func(imax,j) = (func(imax,j) - func(imax-1,j)) ./ dy ./ (dy / 2);
    end

    dtfunc = nu .* (d_x_x_func + d_y_y_func);
    
    func = func + dtfunc * dt;
    
    mesh(x,y,func);
    pause(0.1);
    
end