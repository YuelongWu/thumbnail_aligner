function visualize_mesh(Ms, gry)
    if nargin < 2
        gry = false;
    end
    clf;
    hold on;
    for k = 1:numel(Ms)
        if gry
            cc = hsv2rgb([rand(1), 0, 0.6 + 0.35*rand(1)]);
        else
            cc = hsv2rgb([rand(1), 0.6 + 0.35*rand(1,2)]);
        end
        triplot(Ms{k}.TR.ConnectivityList, double(Ms{k}.TR.Points(:,1)), ...
            double(Ms{k}.TR.Points(:,2)),'Color',cc);
        text(mean(double(Ms{k}.TR.Points(:,1))), mean(double(Ms{k}.TR.Points(:,2))), num2str(k))
    end
end