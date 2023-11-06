function run_everything()
    folder = 'A1_tests';
    recursive_read_folder(folder)
    
    folder = 'A2_tests';
    recursive_read_folder(folder)
    
    disp("End.")
end

function recursive_read_folder(folder)
    fileList = dir(folder);
    for i = 1:length(fileList)
        fileName = fileList(i).name;
        
        if strcmp(fileName, '.') || strcmp(fileName, '..')
            continue;
        end

        pathFile = fullfile(folder, fileName);
        
        if isfile(pathFile)
            str = strtrim(fileName); 
            result = endsWith(str, ".m");
            if result
                launch(pathFile,fileName)
            end

        elseif isfolder(pathFile)
            recursive_read_folder(pathFile)
        end
    end
end

function launch(pathFile,fileName)
    disp(">>> Running experiment:"+fileName)
    run(pathFile);
end
