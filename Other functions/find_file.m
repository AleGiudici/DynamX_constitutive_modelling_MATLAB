function file_name = find_file(extension, protocol, test_name, directory)

A=dir(directory);

if(protocol == 0) % static test
    locator = find(test_name=='_');
    if(test_name(1:2)=='PD')
        if(size(test_name,2)==6)
            target = [test_name(1:locator(1)+1),'.',test_name(locator(1)+2:end)];
        else
            target = [test_name(1:locator(1)),'0.',test_name(locator(1)+1:end)];
        end
    elseif(test_name(1:2)=='FL')
        target = test_name;
    end
elseif(protocol==1)
    if((test_name(end-9)=='5')&(test_name(end-10)~='_'))
        target = [test_name(1:end-3-7),'.',test_name(end-2-7)];
    else
        target = test_name(1:end-3-7);
    end
end

check = 0;
i = 0;
while(strcmp(A(i+1).name(1),'.'))
    i = i+1;
end
% i=3;

while((i<size(A,1))&(check == 0))
    i = i+1;
    if(strcmp(A(i).name(end+1-size(extension,2):end),extension))
        name_file = A(i).name;
        j = 1;
        
        while((check==0)&(j<=size(name_file,2)-size(target,2)+1))
            part_name = name_file(j:j+size(target,2)-1);
            check = strcmp(part_name,target);
            j = j+1;
        end
    end
end

file_name = A(i).name;