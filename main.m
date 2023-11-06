seed = 42;

disp("=======================================");
disp("Select the graph:");
disp("1 : 256 nodes,    2048 edges ");
disp("2 : 1024 nodes,   8192 edges ");
disp("3 : 4096 nodes,   32768 edges ");
disp("=======================================");
x = input("Select the probem [1,2,3]: ","s");

if x == "1"
    filename = "net8_8_3.dmx";
elseif  x == "2"
    filename = "net10_8_3.dmx";
elseif  x == "3"
    filename = "net12_8_3.dmx"; 
else
    disp("Error: please, select a valid choice: [1,2,3].")
    return
end

disp("=======================================");
disp("Select the initialization approach of D:");
disp("1 : all ones");
disp("2 : two distinct values");
disp("3 : five distinct values");
disp("4 : all distinct values");
disp("5 : random in [0,1]");
disp("6 : random in [4,5]");
disp("=======================================");
x = input("Select the problem [1,2,3,4,5,6]: ","s");

[E, ~, b] = utility_read_matrix("graphs/"+filename, seed, false);

D_size = size(E, 2);
if x == "1"
    D = ones(D_size, 1);
elseif  x == "2"
    D = ones(D_size, 1);    %In this case the first half are all ones,
    D(round((D_size/2)+1):end) = 2; % while the second half are all twos.
elseif  x == "3"
    distinct_values = 5;
    D = ones(D_size, 1);
    for i = 2:distinct_values
        start_index = round(((D_size/distinct_values)*(i-1))+1);
        end_index = round((D_size/distinct_values)*(i));
       
        D(start_index:end_index) = i;
    end
elseif  x == "4"
    D = (1:1.5:(1 + (D_size-1)*1.5))';
elseif  x == "5" 
    D = rand(D_size, 1);
elseif  x == "6"
    D = rand(D_size, 1);  
    D = D*4+1;
else
    disp("Error: please, select a valid choice: [1,2,3,4,5,6].")
    return
end

starting_point = b;
threshold = 1e-10;
reorth_flag = true;
debug = true;

disp("=======================================");
x = input("Preconditioner? [y/n]: ","s");

if x == "y"
    [S, ~, creation_time_S] = create_preconditioner(D,E, false);
    S = sparse(S);
elseif x == "n"
    S = NaN;
else
    disp("Error: please, select a valid choice: (y,n)")
    return
end

tic;
[x, r_rel, residuals, break_flag, k] = our_gmres(D, E, S, b, starting_point, threshold, reorth_flag, debug);
exec_time =  toc;
fprintf("Res. Rel: %e | Iter: %d | Trial time : %f\n",  r_rel, k, exec_time);
