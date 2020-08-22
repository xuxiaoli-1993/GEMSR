% read and plot gems residual
clearvars
close all

fileID = fopen('../tests/Riemann/3_np1_largedomain/gems.res.dat','r');
ndim = 2;

if(ndim == 2)
    formatSpec = '%i %e %e %e %e %e';
    sizeA = [6 Inf];
elseif(ndim == 1)
    formatSpec = '%i %e %e %e %e';
    sizeA = [5 Inf];
end
    
A = fscanf(fileID,formatSpec,sizeA);

n = A(1,:);

if(ndim == 2)
    dp = A(3,:);
    du = A(4,:);
    dv = A(5,:);
    dT = A(6,:);
elseif(ndim == 1)
    dp = A(3,:);
    du = A(4,:);
    dT = A(5,:);
end

fig1 = figure(1);
set(fig1, 'Position',[100,100,1500,400]);

subplot(1,3,1);
plot(n, dp, '-o', 'LineWidth', 2);
xlabel('nadv'); ylabel('dp/p');

subplot(1,3,2);
hold all
if(ndim == 2)
    plot(n, du, 'r--o', 'LineWidth', 2);
    plot(n, dv, 'b-o', 'LineWidth', 2);
    xlabel('nadv'); ylabel('dV/V');
    legend('du/u', 'dv/v');
elseif(ndim == 1)
    plot(n, du, 'r', 'LineWidth', 2);
    xlabel('nadv'); ylabel('dV/V');
end

subplot(1,3,3);
plot(n, dT, '-o', 'LineWidth', 2);
xlabel('nadv'); ylabel('dT/T');


        
        
        
        
        
        
        
