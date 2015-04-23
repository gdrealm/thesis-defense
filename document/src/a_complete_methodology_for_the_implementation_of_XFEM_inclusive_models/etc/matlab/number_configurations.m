function[] = number_configurations()
% This function computes all possible combinations of level-set function
% values at the corner nodes to obtain the number of 3D triangulation
% combinations possible.

% x,y,z coordinates of element.
% This will form a 3D cubic element cube.
ex = [0 1 1 0 0 1 1 0];
ey = [0 0 1 1 0 0 1 1];
ez = [0 0 0 0 1 1 1 1];

%% Sweep over all possible level-set configurations.

% Twelve edges. Maximum number of configurations:
maxc = 2^12;
icase = zeros(maxc,1);
hexsect = cell(127,8);

% Set initial configuration, k
% This loop calculates the total number of possible configurations of
% edge intersections
% Routine uses simple negative/positive level set values
k = 0;
for i1 = -1:2:1
    for i2 = -1:2:1
        for i3 = -1:2:1
            for i4 = -1:2:1
                for i5 = -1:2:1
                    for i6 = -1:2:1
                        for i7 = -1:2:1
                            for i8 = -1:2:2
                                levs = [i1 i2 i3 i4 i5 i6 i7 i8];
                                [nsct,isct,xsct,ysct,zsct,levs] = xfem8isct(ex,ey,ez,levs);
                                % isct is a 12 x 1 matrix representing each edge of the 3D cube element
                                % Since each edge can have one or zero intersections, the isct vector becomes a binary display
                                cbin = sprintf('%d%d%d%d%d%d%d%d%d%d%d%d',isct(1),isct(2),isct(3),isct(4),isct(5),isct(6),isct(7),isct(8),isct(9),isct(10),isct(11),isct(12));
                                % cbin is a binary number displayed as a string, and converted into the decimal cdec.
                                cdec = bin2dec(cbin);
                                % If icase is equal to zero, it means it is a new unique configuration and therefore,
                                % number of total configurations k should
                                % be increased by one
                                if icase(cdec+1) == 0 && nsct > 0
                                    k = k+1;
                                    [xtet,ytet,ztet,ctet,ptet,pnd,plist,Tetp1,Tetp2,tet1G,tet1L,tet2G,tet2L] = xfem8tet(isct,xsct,ysct,zsct,ex,ey,ez,levs);
                                    hexsect{k,1} = cdec;
                                    hexsect{k,2} = xtet;
                                    hexsect{k,3} = ytet;
                                    hexsect{k,4} = ztet;
                                    hexsect{k,5} = ctet;
                                    hexsect{k,6} = ptet;
                                    hexsect{k,7} = pnd;
                                    hexsect{k,8} = plist;
                                    [max(plist) min(plist)];
                                    icase(cdec+1)=1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

fprintf('Number of configurations = %d\n',k);

end