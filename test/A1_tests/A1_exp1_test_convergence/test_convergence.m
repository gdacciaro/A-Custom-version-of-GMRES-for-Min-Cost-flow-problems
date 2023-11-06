path_to_root = "../../../";
experiment_title = "exp_1";
addpath(path_to_root)
format long;
seed = 42;
filenames       = ["graphs/net8_8_3.dmx","graphs/net10_8_3.dmx", "graphs/net12_8_3.dmx"];
reorth_flags    = [false, true];
threshold       = 1e-10;
debug           = false;
colors          = ["#0072BD","#D95319"];

file_path = experiment_title+"_results.csv";
fileID = fopen(file_path, 'w');
fprintf(fileID, "file_name;cond;det;reorth;relative residual;number of iterations;time\n");

for i = 1:length(filenames)

    [E, D, b] = utility_read_matrix(path_to_root+filenames(i), seed, debug);
    [c,d] = calculate_det_and_cond(D,E);
    starting_point = b;
    
    res_history = {};

    string_list = split(path_to_root+filenames(i), "/");
    name = string_list(end);
    tmp = split(name, '.');
    name = tmp(1);


    for j = 1:length(reorth_flags)
        tic;
        [x, r_rel, residuals, break_flag, k] = our_gmres(D, E, NaN, b, starting_point, threshold, reorth_flags(j), debug);
        execution_time = toc;

        res_history{end+1} = residuals;

        fprintf(fileID,"%s;%e;%e;%d;%e;%d;%f\n", name, c,d,reorth_flags(j), r_rel, k,execution_time);
    end
    
    plot_file_name = experiment_title+"_"+name+".png";
    plot_res(res_history,plot_file_name, threshold);
end

fclose(fileID);

function plot_res(residuals, filename, threshold)
    colors = ["#D95319", "#4DBEEE"];
    figure;

    p = semilogy(cell2mat(residuals(1)), 'LineWidth',2);
    p.Color = colors(1);

    hold on;
    for i =2:numel(residuals)
       p = semilogy(cell2mat(residuals(i)), 'LineWidth',2);
       p.Color = colors(i);
    end
    yline(threshold,'--','Threshold','LineWidth',3);
    legend(["Without reorth.","With reorth."]);
    xlabel('iteration');
    ylabel('residual');
    hold off;
    if ~isempty(filename)
        saveas(gcf, filename);
    end
end


function [c,d] = calculate_det_and_cond(D,E)
    dim = size(D, 1) + size(E, 1);

    A = zeros(dim, dim);
    A(1:size(D, 1), 1:size(D, 1)) = diag(D);
    A(size(D, 1)+1:end, 1:size(E, 2)) = E;
    A(1:size(D, 1), size(E, 2)+1:end) = E';
    
    c = cond(A);
    d = det(A);
end
