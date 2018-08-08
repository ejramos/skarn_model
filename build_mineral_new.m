function Mineral = build_mineral_new(filename,Grid)

%author: Evan J. Ramos
%date:   12 Jan 2017

%Description: This function creates the Mineral structure, which contains
%properties of which affect the stable isotope transport equation. It
%further reads in a permeability matrix and from the permeability
%prescribed to that node, it assigns are certain assemblage of minerals
%associated with a rock type that has a similar permeability.

%%%%%% FILE FORMAT %%%%%%
%1    xxxxx    mineral1      mineral2      mineral3       mineral4   ...
%2      D       xxxxx         xxxxx          xxxxx         xxxxx
%3      E       xxxxx         xxxxx          xxxxx         xxxxx
%4      F       xxxxx         xxxxx          xxxxx         xxxxx
%5      Ea      xxxxx         xxxxx          xxxxx         xxxxx
%6      A0      xxxxx         xxxxx          xxxxx         xxxxx
%7     A_bar    xxxxx         xxxxx          xxxxx         xxxxx
%8      nu      xxxxx         xxxxx          xxxxx         xxxxx
%9      xm      xxxxx         xxxxx          xxxxx         xxxxx
%10    d18O     xxxxx         xxxxx          xxxxx         xxxxx

%%%%%
%The values for each of the minerals come from these sources:
%Northrop and Clayton (1966)
%O'Neil and Taylor (1967)
%O'Neil (1969)
%Cole and Ohmoto (1986)
%Clayton et al. (1989)
%Person et al. (2007) --> effective surface area for reactions (A_bar)


%% Create Mineral structure
fid   = fopen(filename);
count = 1;

if fid == -1
    disp('Error opening file')
else
    while feof(fid) ~= 1
        
        %read in individual lines of the file as strings
        aline = fgetl(fid);
        
        if count == 1
        %Special treatment for the first line of the file. Creating a cell
        %array of all of the mineral names which will ultimately become the
        %primary fields to the Mineral structure.
        
            [~,rest]  = strtok(aline);
            min_names = strtrim(rest);
            
            while ~isempty(strfind(min_names,'  '))
                min_names = strrep(min_names,'  ',' ');
            end
            
            blank_loc = strfind(min_names,' ');
            min_total = length(blank_loc)+1;
            mins      = cell(min_total,1);
            
            for i = 1:min_total
                if i == 1
                    mins{i} = strtrim(min_names(1:blank_loc(1)));
                elseif i == min_total
                    mins{i} = strtrim(min_names(blank_loc(end):end));
                else
                    mins{i} = strtrim(min_names(blank_loc(i-1):blank_loc(i)));
                end
            end
            
        else
        %Create the subfields for each of minerals
        [subfield, vals] = strtok(aline);
        vals = str2num(vals);
        for j = 1:length(vals)
            Mineral.(mins{j}).(subfield) = vals(j);
            Mineral.(mins{j}).phi        = -ones(Grid.N,1);
            %preallocated as -1s as an indicator that the value of that
            %element in the vector has not been reassigned (see the
            %function make_medium.m to make sense of it)
        end
            
        end
        
        count = count + 1;
    end
end

fclose(fid);
end