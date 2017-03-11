% Function to convert all fig files in current directory to eps
function fig2eps_gray()

% Assumptions and modifications
% - all fig files are in current directory
% 

% All fig files in current directory
figs = dir('*.fig');
nfigs = length(figs);
if nfigs == 0
    disp('No figure files to convert');
else
    for i = 1:nfigs
        % Open figure
        open(figs(i).name);
        % Save as colour eps, remove .fig
        name = strtok(figs(i).name, '.');
        % Export to gray colormap
        style = hgexport('factorystyle');
        style.Color = 'gray';
        hgexport(gcf, [name '.eps'], style);
        %saveas(gcf, name, 'epsc');
    end
    close all
end