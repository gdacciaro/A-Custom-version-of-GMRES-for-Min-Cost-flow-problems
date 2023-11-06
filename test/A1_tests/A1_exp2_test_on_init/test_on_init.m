path_to_root = "../../../";
experiment_title = "exp_2";
addpath(path_to_root)
format long;
seed = 42;
filenames = ["graphs/net10_8_3.dmx"];
reorth_flag = true;
init_mode = ["random", "random_between_interval","identity","2_dist_eig","5_dist_eig","all_diff"];
colors = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#4DBEEE","#A2142F"];
threshold = 1e-10;
debug = false;

file_path = experiment_title+"_results.csv";
fileID = fopen(file_path, 'w');
fprintf(fileID, "mode;cond;det;relative residual;number of iterations;time\n");

for j = 1:length(filenames)
    filename = filenames(j);
    for i = 1:length(init_mode)
        
        [E, ~, b] = utility_read_matrix(path_to_root+filename, seed, debug);
    
        D = init_D(size(E,2), init_mode(i));

        [c,d] = calculate_det_and_cond(D,E);
    
        starting_point  = b;
        
        tic;
        [x, r_rel, residuals, break_flag, k] = our_gmres(D, E, NaN, b, starting_point, threshold, reorth_flag, debug);
        execution_time = toc;
        
        fprintf(fileID,"%s;%e;%e;%e;%d;%e\n",init_mode(i),c,d,r_rel, k,execution_time);

        string_list = split(filename, "/");
        name = string_list(end);
        tmp = split(name, '.');
        name = tmp(1);
        plot_file_name = experiment_title+"_"+name+"_"+init_mode(i)+".png";
        plot_res(residuals,colors(i), plot_file_name,threshold);
    end

end

fclose(fileID);

function plot_res(residuals, color, filename,threshold)
    figure;
    p = semilogy(residuals, 'LineWidth',2);
    p.Color=color;
    xlabel('iteration');
    ylabel('residual');
    yline(threshold,'--','Threshold','LineWidth',3);
    if ~isempty(filename)
        saveas(gcf, filename);
    end
end 


function [D] = init_D(dim, mode)
    switch mode
        case "random"
            D = rand(dim, 1);
        case "random_between_interval"
            D = rand(dim, 1);  
            D = D*4+1;
        case "identity"
            D = ones(dim, 1);
        case "2_dist_eig"
            D = ones(dim, 1);    %In this case the first half are all ones,
            D(round((dim/2)+1):end) = 2; % while the second half are all twos.
        case "5_dist_eig"
            distinct_values = 5;
            D = ones(dim, 1);
            for i = 2:distinct_values
                start_index = round(((dim/distinct_values)*(i-1))+1);
                end_index = round((dim/distinct_values)*(i));
               
                D(start_index:end_index) = i;
            end
        case "all_diff"
            D = 1:1.5:(1 + (dim-1)*1.5);
            D = D';
        otherwise
            disp("Init mode not valid, D initialized with uniform distribution")
            D = rand(dim, 1);
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
