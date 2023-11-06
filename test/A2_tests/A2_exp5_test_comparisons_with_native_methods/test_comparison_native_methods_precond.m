path_to_root = "../../../";
experiment_title = "exp_5";
addpath(path_to_root)
format long;
seed = 42;
filenames = ["graphs/net10_8_3.dmx", "graphs/net10_8_3.dmx", "graphs/net12_8_3.dmx"];
threshold = 1e-10;
debug = false;
trials = 5;

file_path = experiment_title+"_results.csv";
fileID = fopen(file_path, 'w');
fprintf(fileID, "file_name;cond;det;our relative residual;our number of iterations;our time;GMRES relative residual;GMRES number of iterations;MINRES time;GMRES relative residual;MINRES number of iterations;MINRES time;  creation time of S;\n");

for i = 1:length(filenames)
    filename = filenames(i);
    [E, ~, b] = utility_read_matrix(path_to_root+filename, seed, debug);
    
    D = ones(size(E,2), 1);
    
    [S, P, total_time_S] = create_preconditioner(D,E); 
    S = sparse(S);
    starting_point  = b;
    
    [A,c,d] = calculate_det_and_cond(D,E);
    A = sparse(A);
    residuals = {};

    total_time = 0;
    for trial=1:trials
        tic;
        [~, our_r_rel, our_res_vec, break_flag, our_k] = our_gmres(D, E, S, b, starting_point, threshold, true, debug);
        trial_time = toc;
        total_time = total_time + trial_time;
    end
    our_time = total_time/trials;

    residuals{1} = our_res_vec;

    dim = size(D, 1) + size(E, 1);
    
    total_time = 0;
    for trial=1:trials
        tic;
        [~, ~, gmres_r_rel, gmres_n_iter, res_vec] = gmres(A, b, [], threshold, dim, P',P,starting_point);
        trial_time = toc;
        total_time = total_time + trial_time;
    end
    gmres_time = total_time/trials;
    gmres_k = gmres_n_iter(2);

    residuals{2} = res_vec/norm(b);

    total_time = 0;
    for trial=1:trials
        tic;
        [~, ~, minres_r_rel, minres_n_iter, res_vec] = minres(A, b, threshold, dim, P',P,starting_point);
        trial_time = toc;
        total_time = total_time + trial_time;
    end
    minres_time = total_time/trials;

    residuals{3} = res_vec/norm(b);
    

    string_list = split(filename, "/");
    name = string_list(end);
    tmp = split(name, '.');
    name = tmp(1);
    plot_file_name = experiment_title+"_"+name+"_"+"comparison_with_native_methods.png";
    plot_res(residuals, plot_file_name, threshold);
   
    fprintf(fileID,"%s;%e;%e;%e;%d;%f;%e;%d;%f;%e;%d;%f;%f;\n", name,c,d, our_r_rel, our_k, our_time, ...
                    gmres_r_rel, gmres_k, gmres_time, ...
                    minres_r_rel, minres_n_iter, minres_time, total_time_S);
end
fclose(fileID);

function plot_res(residuals, filename, threshold)
    colors = ["#0072BD","#D95319","#EDB120"];
    figure;
    
    p = semilogy(cell2mat(residuals(1)), 'LineWidth',2);
    p.Color = colors(1);

    hold on;
    for i =2:numel(residuals)
       p = semilogy(cell2mat(residuals(i)), 'LineWidth',2);
       p.Color = colors(i);
    end
    yline(threshold,'--','Threshold','LineWidth',3);
    legend(["Our method","MATLAB GMRES","MATLAB MINRES"]);
    xlabel('iteration');
    ylabel('residual');
    hold off;
    if ~isempty(filename)
        saveas(gcf, filename);
    end
end

function [A, c,d] = calculate_det_and_cond(D,E)
    dim = size(D, 1) + size(E, 1);

    A = zeros(dim, dim);
    A(1:size(D, 1), 1:size(D, 1)) = diag(D);
    A(size(D, 1)+1:end, 1:size(E, 2)) = E;
    A(1:size(D, 1), size(E, 2)+1:end) = E';
    
    c = cond(A);
    d = det(A);
end

