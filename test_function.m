clear
addpath ..
alpha = -2;
gamma = 3;
simulation_name = '8_18_21_alpha_-2_gamma_3';

delete *.out

num = [60,90,50,60,35,30,20,60,35,50,30,40,25,30,20,40,25];
total_samples = [1904,712,712,388,388,155,155,376,376,280,280,237,237,137,137,227,227];
names = {'metabric','brca','brca','hnsc','hnsc','laml','laml','luad','luad','lusc','lusc','ov','ov','paad','paad','sarc','sarc'};
Nc = [2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2];

for i = 1:length(num)
    intervals = 1:((total_samples(i)-1)/num(i)):total_samples(i);
    intervals = round(intervals);
    intervals = [intervals(1:(end-1))' intervals(2:end)'];
    intervals(1:end-1,2)=intervals(1:end-1,2)-1;
    for j = 1:num(i)
        counter = 1;
        for k = intervals(j,1):intervals(j,2)
            if isfile(sprintf('%s%i_%i_%.1f_%i.mat',names{1,i},Nc(i),alpha,gamma,k))
                test_existence(counter) = 1;
                counter = counter + 1;
            else
                test_existence(counter) = 0;
                counter = counter + 1;
            end
        end
        if all(test_existence)
           rerun(i,j) = 0;
        else
           rerun(i,j) = 1;
        end
        clear test_existence
    end
end
fileID = fopen(sprintf('script.sh'),'w');
for i = 1:length(num)
    intervals = 1:((total_samples(i)-1)/num(i)):total_samples(i);
    intervals = round(intervals);
    intervals = [intervals(1:(end-1))' intervals(2:end)'];
    intervals(1:end-1,2)=intervals(1:end-1,2)-1;
    for j = 1:num(i)
        if rerun(i,j) == 1
            fprintf(fileID, sprintf('sbatch -J "%s_%i" --partition=large --nodes=1 --exclusive --output=%s_%i.out --time=6:00:00 --mem=230gb --wrap "matlab -nodisplay -nosplash -nodesktop -r ''cd %s;vec_curv_fun(%i,%i,%i,%i,%i);exit;''"\n',names{1,i},j,names{1,i},j,simulation_name,intervals(j,1),intervals(j,2),alpha,gamma,i));
            %fprintf(fileID, sprintf('sleep 1\n'));
        end
    end
end
fclose(fileID);

unix2dos('script.sh',logical(1));

